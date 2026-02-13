#include "cppfm/market/SviCalibration.h"
#include "cppfm/calibration/Optimizer.h"
#include <cmath>
#include <numeric>

SviParams::SviParams(double a, double b, double rho, double m, double sigma)
	: a(a), b(b), rho(rho), m(m), sigma(sigma)
{
}

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
