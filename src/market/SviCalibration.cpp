#include "cppfm/market/SviCalibration.h"
#include "cppfm/calibration/LinearAlgebra.h"
#include "cppfm/calibration/Optimizer.h"
#include "cppfm/utils/InterpolationSchemes.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

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

double SsviParams::dwdtheta(double k) const
{
	// chain rule: dw/dtheta = w/theta + k * dw/dk * (dphi/dtheta)/phi
	// (dphi/dtheta)/phi = h = -gamma/theta - (1-gamma)/(1+theta)
	if (theta < 1e-10)
		return 0.0;
	double h = -gamma / theta - (1.0 - gamma) / (1.0 + theta);
	return totalVariance(k) / theta + k * dw(k) * h;
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

// enforce monotonicity - want to ensure calendar-abr-free
static std::vector<double> enforceMonotonicity(const std::vector<double> &values)
{
	int n = static_cast<int>(values.size());
	if (n <= 1)
		return values;

	// each block (sum, count)
	std::vector<std::pair<double, int>> blocks;

	for (int i = 0; i < n; ++i)
	{
		blocks.push_back({values[i], 1});
		// merge backward while we have violation
		while (blocks.size() >= 2)
		{
			auto &curr = blocks[blocks.size() - 1];
			auto &prev = blocks[blocks.size() - 2];
			if (curr.first / curr.second >= prev.first / prev.second)
				break;
			// merge into prev
			prev.first += curr.first;
			prev.second += curr.second;
			blocks.pop_back();
		}
	}

	// expand back
	std::vector<double> result;
	result.reserve(n);
	for (auto &[sum, count] : blocks)
	{
		double avg = sum / count;
		for (int j = 0; j < count; ++j)

		{
			result.push_back(avg);
		}
	}
	return result;
}

// linear interpolation to get ATM vol when SVI fit fails
static double interpolateAtmVol(const std::vector<double> &strikes,
								const std::vector<double> &vols, double forward)
{
	// strikes may be provided as unsorted
	assert(std::is_sorted(strikes.begin(), strikes.end()));
	for (int i = 0; i < static_cast<int>(strikes.size()) - 1; ++i)
	{
		if (strikes[i] <= forward && strikes[i + 1] >= forward)
		{
			double t = (forward - strikes[i]) / (strikes[i + 1] - strikes[i]);
			return vols[i] + t * (vols[i + 1] - vols[i]);
		}
	}
	// forward outside strike range
	if (forward <= strikes.front())
		return vols.front();
	return vols.back();
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
									const SsviCalibrationOptions &opts)
{
	int nMat = static_cast<int>(maturities.size());

	// stage 1: per-slice SVI -> thetas + warm-start rho
	std::vector<double> thetas(nMat);
	std::vector<std::vector<double>> logStrikes(nMat);
	std::vector<std::vector<double>> wMkt(nMat);
	double sumRho = 0.0;
	int nConverged = 0;

	for (int i = 0; i < nMat; ++i)
	{
		SviParams sviGuess;
		auto fit = calibrateSmile(sviGuess, strikesPerMaturity[i], volsPerMaturity[i],
								  forwards[i], maturities[i]);

		if (fit.converged)
		{
			thetas[i] = fit.params->totalVariance(0.0);
			auto *svi = dynamic_cast<SviParams *>(fit.params.get());
			if (svi)
			{
				sumRho += svi->rho;
				++nConverged;
			}
		}
		else
		{
			// fallback: interpolate ATM vol from market quotes
			double atmVol = interpolateAtmVol(strikesPerMaturity[i],
											  volsPerMaturity[i], forwards[i]);
			thetas[i] = atmVol * atmVol * maturities[i];
		}

		int nK = static_cast<int>(strikesPerMaturity[i].size());
		logStrikes[i].resize(nK);
		wMkt[i].resize(nK);
		for (int j = 0; j < nK; ++j)
		{
			logStrikes[i][j] = std::log(strikesPerMaturity[i][j] / forwards[i]);
			wMkt[i][j] = volsPerMaturity[i][j] * volsPerMaturity[i][j] * maturities[i];
		}
	}

	// enforce monotonicity
	if (opts.enforceMonotonicity)
		thetas = enforceMonotonicity(thetas);

	// nudge the flat regions so dthetadT > 0 for local vol
	for (int i = 1; i < nMat; ++i)
	{
		if (thetas[i] <= thetas[i - 1])
			thetas[i] = thetas[i - 1] + 1e-8;
	}

	// stage 2: global SSVI fit
	SsviParams local;

	// normalized residuals -- each slice weighted equally in relative terms
	auto residuals = [&](const std::vector<double> &p) -> std::vector<double>
	{
		local.fromVector(p);
		std::vector<double> r;
		for (int i = 0; i < nMat; ++i)
		{
			local.theta = thetas[i];
			double scale = 1.0 / thetas[i];
			for (int j = 0; j < static_cast<int>(logStrikes[i].size()); ++j)
				r.push_back((local.totalVariance(logStrikes[i][j]) - wMkt[i][j]) * scale);
		}
		return r;
	};

	// analytical Jacobian -- dw/d{rho, eta, gamma}
	JacobianFunc jacobian = [&](const std::vector<double> &p) -> Matrix
	{
		local.fromVector(p);
		Matrix J;
		for (int i = 0; i < nMat; ++i)
		{
			local.theta = thetas[i];
			double th = thetas[i];
			double ph = local.phi();
			double rho_val = local.rho;
			double eta_val = local.eta;
			double logRatio = std::log((1.0 + th) / th);
			double scale = 1.0 / th;

			for (int j = 0; j < static_cast<int>(logStrikes[i].size()); ++j)
			{
				double k = logStrikes[i][j];
				double pk = ph * k;
				double D = std::sqrt((pk + rho_val) * (pk + rho_val) + (1.0 - rho_val * rho_val));

				// dw/dk reused by eta and gamma partials
				double dw_dk = (th / 2.0) * ph * (rho_val + (pk + rho_val) / D);

				double dw_drho = (th / 2.0) * pk * (1.0 + 1.0 / D);
				double dw_deta = (k / eta_val) * dw_dk;
				double dw_dgamma = k * logRatio * dw_dk;

				J.push_back({dw_drho * scale, dw_deta * scale, dw_dgamma * scale});
			}
		}
		return J;
	};

	// grid search - pick best starting point
	double warmRho = (nConverged > 0) ? sumRho / nConverged : -0.3; // warm starter on the rhos
	std::vector<SsviParams> candidates;

	if (opts.useGridSearch)
	{
		candidates.push_back(SsviParams(warmRho, 0.5, 0.5));
		candidates.push_back(SsviParams(-0.3, 0.3, 0.3));
		candidates.push_back(SsviParams(-0.7, 0.8, 0.7));
	}
	else
	{
		candidates.push_back(SsviParams());
	}

	auto evalNorm = [&](const SsviParams &c) -> double
	{
		auto r = residuals(c.toVector());
		return std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
	};

	SsviParams bestGuess = candidates[0];
	double bestNorm = evalNorm(bestGuess);
	for (size_t c = 1; c < candidates.size(); ++c)
	{
		double norm = evalNorm(candidates[c]);
		if (norm < bestNorm)
		{
			bestNorm = norm;
			bestGuess = candidates[c];
		}
	}

	// LM solve
	LevenbergMarquardt lm;
	LMResult lmResult = lm.solve(residuals, bestGuess.toVector(),
								 bestGuess.lowerBounds(), bestGuess.upperBounds(),
								 jacobian, opts.lmOpts);

	// RMSE in original total-variance units (undo normalization)
	local.fromVector(lmResult.params);
	double sumSq = 0.0;
	int totalPts = 0;
	for (int i = 0; i < nMat; ++i)
	{
		local.theta = thetas[i];
		for (int j = 0; j < static_cast<int>(logStrikes[i].size()); ++j)
		{
			double r = local.totalVariance(logStrikes[i][j]) - wMkt[i][j];
			sumSq += r * r;
			++totalPts;
		}
	}

	SsviCalibrationResult result;
	result.params.fromVector(lmResult.params);
	result.thetas = thetas;
	result.maturities = maturities;
	result.rmse = std::sqrt(sumSq / totalPts);
	result.converged = lmResult.converged;

	// check Gatheral & Jacquier arb-free conditions
	bool arbFree = result.params.satisfiesConstraints();
	for (int i = 0; i < nMat && arbFree; ++i)
		arbFree = result.params.satisfiesButterflyConstraint(thetas[i]);
	result.arbitrageFree = arbFree;

	return result;
}

// ============================================================================
// SsviSurface
// ============================================================================

SsviSurface::SsviSurface(const SsviCalibrationResult &result,
						 const std::vector<double> &forwards)
	: _params(result.params),
	  _maturities(result.maturities),
	  _thetas(result.thetas),
	  _forwards(forwards)
{
	_thetaInterp = std::make_unique<CubicSplineInterpolation>(
		_maturities, _thetas,
		CubicSplineInterpolation::BoundaryType::Natural,
		ExtrapolationType::Linear);

	_forwardInterp = std::make_unique<CubicSplineInterpolation>(
		_maturities, _forwards,
		CubicSplineInterpolation::BoundaryType::Natural,
		ExtrapolationType::Linear);
}

double SsviSurface::thetaAt(double T) const
{
	return std::max((*_thetaInterp)(T), 1e-10);
}

double SsviSurface::dthetadT(double T) const
{
	// floor at 0 -- cubic spline can overshoot, negative would mean calendar arb
	return std::max(_thetaInterp->derivative(T), 0.0);
}

double SsviSurface::forwardAt(double T) const
{
	return (*_forwardInterp)(T);
}

double SsviSurface::localVariance(double k, double T) const
{
	// sigma^2_LV = (dw/dtheta * dtheta/dT) / g(k)
	SsviParams local = _params;
	local.theta = thetaAt(T);

	double dwdt = local.dwdtheta(k);
	double dtdt = dthetadT(T);
	double g = std::max(local.gFunction(k), 1e-8);

	return std::max(dwdt * dtdt / g, 1e-10);
}

double SsviSurface::localVolatility(double spot, double T) const
{
	double F = forwardAt(T);
	double k = std::log(spot / F);
	return std::sqrt(localVariance(k, T));
}

double SsviSurface::impliedVolatility(double strike, double T) const
{
	double F = forwardAt(T);
	double k = std::log(strike / F);

	SsviParams local = _params;
	local.theta = thetaAt(T);
	return std::sqrt(local.totalVariance(k) / T);
}

SsviSurface SsviCalibrationResult::buildAnalyticalSurface(const std::vector<double> &forwards) const
{
	return SsviSurface(*this, forwards);
}
