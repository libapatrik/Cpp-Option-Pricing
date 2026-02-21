#include "cppfm/market/SviCalibration.h"
#include "cppfm/calibration/LinearAlgebra.h"
#include "cppfm/calibration/Optimizer.h"
#include <cmath>
#include <numeric>

SviParams::SviParams(double a, double b, double rho, double m, double sigma)
	: a(a), b(b), rho(rho), m(m), sigma(sigma)
{
}

double SviParams::dw(double k) const
{
	double dkm = k - m;
	return b * (rho + dkm / std::sqrt(dkm * dkm + sigma * sigma));
}

double SviParams::d2w(double k) const
{
	double dkm = k - m;
	double s2 = dkm * dkm + sigma * sigma;
	return b * sigma * sigma / (s2 * std::sqrt(s2));
}

double SmileParams::gFunction(double k) const
{
	// when g(k) >= 0 for all k, the implied density is non-negative - no butterfly arbitrage
	// We check this in the builder
	// if k=0, as long as w(0) > 0, the division is fine (4th constraint)
	double w = totalVariance(k);
	double wp = dw(k);
	double wpp = d2w(k);

	// Gatheral (2004) Eq. for risk-neutral density
	double term1 = 1.0 - k * wp / (2.0 * w);
	double term2 = wp * wp / 4.0 * (1.0 / w + 0.25);
	return term1 * term1 - term2 + wpp / 2.0;
}

// ============================================================================
// SVI
// ============================================================================

SviParams::SviParams(const std::vector<double> &v) { fromVector(v); }

double SviParams::totalVariance(double k) const
{
	double dkm = k - m;
	return a + b * (rho * dkm + std::sqrt(dkm * dkm + sigma * sigma));
}

std::vector<double> SviParams::toVector() const
{
	return {a, b, rho, m, sigma};
}

void SviParams::fromVector(const std::vector<double> &v)
{
	a = v[0];
	b = v[1];
	rho = v[2];
	m = v[3];
	sigma = v[4];
}

std::vector<double> SviParams::lowerBounds() const
{
	return {-1.0, 1e-8, -0.999, -2.0, 1e-6};
}

std::vector<double> SviParams::upperBounds() const
{
	return {1.0, 5.0, 0.999, 2.0, 5.0};
}

int SviParams::nParams() const { return 5; }

std::unique_ptr<SmileParams> SviParams::clone() const
{
	return std::make_unique<SviParams>(*this);
}

SmileFitResult calibrateSmile(const SmileParams &guess,
							  const std::vector<double> &strikes,
							  const std::vector<double> &marketVols,
							  double forward, double T, const LMOptions &opts)
{
	int n = static_cast<int>(strikes.size());

	// convert to SVI space
	std::vector<double> k(n), wMkt(n);
	for (int i = 0; i < n; ++i)
	{
		k[i] = std::log(strikes[i] / forward);
		wMkt[i] = marketVols[i] * marketVols[i] * T;
	}

	// clone guess so the lambda can mutate it
	auto prototype = guess.clone();
	auto residuals = [&](const std::vector<double> &p) -> std::vector<double>
	{
		prototype->fromVector(p);
		std::vector<double> r(n);
		for (int i = 0; i < n; ++i)
			r[i] = prototype->totalVariance(k[i]) - wMkt[i];
		return r;
	};

	LevenbergMarquardt lm;
	LMResult result = lm.solve(residuals,
							   guess.toVector(),
							   guess.lowerBounds(),
							   guess.upperBounds(),
							   nullptr,
							   opts);

	auto finalR = residuals(result.params);
	double sumSq = std::inner_product(finalR.begin(), finalR.end(), finalR.begin(), 0.0);

	SmileFitResult fit;
	fit.params = guess.clone();
	fit.params->fromVector(result.params);
	fit.rmse = std::sqrt(sumSq / n);
	fit.converged = result.converged;
	return fit;
}

bool SviParams::satisfiesConstraints() const
{
	if (b < 0.0)
		return false;
	if (std::abs(rho) >= 1.0)
		return false; // |rho| < 1
	if (sigma <= 0.0)
		return false;
	// non-negative total variance: a + b*sigma*sqrt(1-rho^2) >= 0
	if (a + b * sigma * std::sqrt(1.0 - rho * rho) < 0.0)
		return false;
	// R Lee wing slopes: b*(1+|rho|) <= 2
	if (b * (1.0 + std::abs(rho)) > 2.0)
		return false;
	return true;
}

// ============================================================================
// SSVI
// ============================================================================

SsviParams::SsviParams(double rho, double eta, double gamma)
	: rho(rho), eta(eta), gamma(gamma) {}

SsviParams::SsviParams(const std::vector<double> &v)
{
	fromVector(v);
}

double SsviParams::phi() const
{
	return phi(theta);
}

double SsviParams::phi(double th) const
{
	// power-law parametrization
	return eta / (std::pow(th, gamma) * std::pow(1.0 + th, 1.0 - gamma));
}

double SsviParams::totalVariance(double k) const
{
	double p = phi();
	double pk = p * k;
	return (theta / 2.0) * (1.0 + rho * pk + std::sqrt((pk + rho) * (pk + rho) + (1.0 - rho * rho)));
}

std::vector<double> SsviParams::toVector() const
{
	return {rho, eta, gamma};
}
void SsviParams::fromVector(const std::vector<double> &v)
{
	rho = v[0];
	eta = v[1];
	gamma = v[2];
}

std::vector<double> SsviParams::lowerBounds() const
{
	return {-0.999, 0.001, 0.0};
}
std::vector<double> SsviParams::upperBounds() const
{
	return {0.999, 5.0, 1.0};
}

int SsviParams::nParams() const
{
	return 3;
}

double SsviParams::dw(double k) const
{
	double p = phi();
	double pk = p * k;
	double disc = std::sqrt((pk + rho) * (pk + rho) + (1.0 - rho * rho));
	return (theta / 2.0) * p * (rho + (pk + rho) / disc);
}

double SsviParams::d2w(double k) const
{
	double p = phi();
	double pk = p * k;
	double inner = (pk + rho) * (pk + rho) + (1.0 - rho * rho);
	double disc = std::sqrt(inner);
	return (theta / 2.0) * p * p * (1.0 - rho * rho) / (disc * disc * disc);
}

bool SsviParams::satisfiesConstraints() const
{
	if (std::abs(rho) >= 1.0)
		return false;
	if (eta <= 0.0)
		return false;
	if (gamma < 0.0 || gamma > 1.0)
		return false;
	if (eta * (1.0 + std::abs(rho)) > 2.0)
		return false; // large-theta wing bound
	return true;
}

std::unique_ptr<SmileParams> SsviParams::clone() const
{
	return std::make_unique<SsviParams>(*this);
}

bool SsviParams::satisfiesButterflyConstraint(double th) const
{
	// Thm 4.3
	return th * phi(th) * (1.0 + std::abs(rho)) < 4.0;
}

std::unique_ptr<VolatilitySurface> SsviCalibrationResult::buildSurface(
	const std::vector<double> &strikeGrid,
	const std::vector<double> &forwards,
	const DiscountCurve &discountCurve) const
{
	VolatilitySurfaceBuilder builder;
	builder.setDiscountCurve(discountCurve);

	SsviParams local = params;

	for (int i = 0; i < static_cast<int>(maturities.size()); ++i)
	{
		local.theta = thetas[i];
		double T = maturities[i];
		for (double K : strikeGrid)
		{
			double k = std::log(K / forwards[i]);
			double iv = std::sqrt(local.totalVariance(k) / T);
			builder.setVolatility(K, T, iv);
		}
	}
	return builder.build();
}

SsviCalibrationResult calibrateSsvi(const std::vector<std::vector<double>> &strikesPerMaturity,
									const std::vector<std::vector<double>> &volsPerMaturity,
									const std::vector<double> &forwards,
									const std::vector<double> &maturities,
									const LMOptions &opts)
{
	return {}; // * work in progress, on paper
}
