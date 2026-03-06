/**
	Heston Model Calibration
	- Fits (v0, kappa, vbar, sigma_v, rho) to market implied vols
	- Uses HestonCF + COSPricer for fast pricing, LM for optimization
	- Vega-weighted price residuals: (C_model - C_mkt) / vega ≈ IV diff
*/

#ifndef CPPFM_HESTON_CALIBRATION_H
#define CPPFM_HESTON_CALIBRATION_H

#include <cppfm/calibration/Optimizer.h>
#include <string>
#include <vector>

struct HestonParams
{
	double v0 = 0.04;		// initial variance
	double kappa = 1.5;		// mean reversion speed
	double vbar = 0.04;		// long-term variance
	double sigma_v = 0.3;	// vol of vol
	double rho = -0.7;		// spot-vol correlation

	HestonParams() = default;
	HestonParams(double v0, double kappa, double vbar, double sigma_v, double rho);
	HestonParams(const std::vector<double> &v);

	std::vector<double> toVector() const;
	void fromVector(const std::vector<double> &v);

	// box constraints for LM
	static std::vector<double> lowerBounds();
	static std::vector<double> upperBounds();

	// 2*kappa*vbar > sigma_v^2
	bool satisfiesFellerCondition() const;

	static constexpr int nParams() { return 5; }
};

// per-maturity market data
struct HestonSliceData
{
	double T;
	std::vector<double> strikes;
	std::vector<double> marketIVs;
};

struct HestonCalibrationResult
{
	HestonParams params;
	double rmse = 0.0;
	std::vector<double> sliceRmse; // per-maturity RMSE
	int iterations = 0;
	bool converged = false;
	std::string message;
};

// main entry point: calibrate across all maturities simultaneously
HestonCalibrationResult calibrateHeston(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	const HestonParams &guess = {},
	const LMOptions &opts = {});

// coarse grid search for a decent starting point
HestonParams hestonGridSearch(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	int gridPoints = 3);

// parallel version of grid search
HestonParams hestonGridSearchParallel(
	const std::vector<HestonSliceData> &slices,
	double S0, double r,
	int gridPoints = 3);

#endif
