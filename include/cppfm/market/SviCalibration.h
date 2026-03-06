#ifndef CPPFM_SVICALIBRATION_H
#define CPPFM_SVICALIBRATION_H

#include "cppfm/calibration/Optimizer.h"
#include "cppfm/market/DiscountCurve.h"
#include "cppfm/market/VolatilitySurface.h"
#include "cppfm/utils/InterpolationSchemes.h"
#include <memory>
#include <vector>

// forward declarations
class VolatilitySurface;
class DiscountCurve;
class SsviSurface;

// base class, abstract

class SmileParams
{
public:
	virtual ~SmileParams() = default;				  // hold unique_ptr<SmileParams>
	virtual double totalVariance(double k) const = 0; // each derived class defines its own formula
	double impliedVol(double k, double T) const
	{
		return std::sqrt(totalVariance(k) / T);
	}

	// LM needs vector<double>
	virtual std::vector<double> toVector() const = 0;
	virtual void fromVector(const std::vector<double> &v) = 0;

	// toVector() gives LM initial guesses
	// fromVector() updates the object after LM finds a solution

	// bounds for the LM's inputs
	virtual std::vector<double> lowerBounds() const = 0;
	virtual std::vector<double> upperBounds() const = 0;

	virtual int nParams() const = 0;
	// Dims of param vector. SVI=5, SSVI=3,4; idea: could freeze some params

	virtual std::unique_ptr<SmileParams> clone() const = 0;
	// for calibrateSmile() it takes const SmileParams& but needs mutable copy
	// for the residual function lambda

	// analytical derivatives
	virtual double dw(double k) const = 0;
	virtual double d2w(double k) const = 0;
	// Gatheral density must be >= 0 for no butterfly arbitrage
	virtual double gFunction(double k) const; // model agnostic, uses dw/d2w
};

// ============================================================================
// SVI
// ============================================================================
class SviParams : public SmileParams
{
public:
	// SVI - 5 params - serve as initial guesses
	double a = 0.04;
	double b = 0.01;
	double rho = -0.3;
	double m = 0.0; // shift
	double sigma = 0.2;

	SviParams() = default;
	SviParams(double a, double b, double rho, double m, double sigma);
	SviParams(const std::vector<double> &v);

	// w(k) = a + b * (rho*(k-m) + sqrt((k-m)^2 + sigma^2))
	double totalVariance(double k) const override;

	std::vector<double> toVector() const override;
	void fromVector(const std::vector<double> &v) override;

	std::vector<double> lowerBounds() const override;
	std::vector<double> upperBounds() const override;
	int nParams() const override;

	std::unique_ptr<SmileParams> clone() const override;
	// default copy constructor copies those 5 doubles

	// analytical derivatives
	double dw(double k) const override;
	double d2w(double k) const override;

	// post-calibration checks (not checked during LM)
	bool satisfiesConstraints() const;
};

// The results from the calibrateSmile()
struct SmileFitResult
{
	std::unique_ptr<SmileParams> params; // pointer to the fitted smile
	double rmse = 0.0;					 // error metric
	bool converged = false;				 // if LM hits the tolerance or max iters
};

// Calibration
SmileFitResult calibrateSmile(const SmileParams &guess,
							  const std::vector<double> &strikes,
							  const std::vector<double> &marketVols,
							  double forward, double T,
							  const LMOptions &opts = {});

// ============================================================================
// SSVI
// ============================================================================

class SsviParams : public SmileParams
{
public:
	// 3 params
	double rho = -0.5;
	double eta = 0.1;
	double gamma = 0.5;

	// per slice
	double theta = 0.04;

	SsviParams() = default;
	SsviParams(double rho, double eta, double gamma);
	SsviParams(const std::vector<double> &v);

	// phi(theta) = eta / (theta^gamma * (1 + theta)^(1-gamma))
	double phi() const;
	double phi(double th) const;

	// SSVI total variance
	double totalVariance(double k) const override;

	std::vector<double> toVector() const override; // { rho, eta, gamma}
	void fromVector(const std::vector<double> &v) override;

	std::vector<double> lowerBounds() const override;
	std::vector<double> upperBounds() const override;
	int nParams() const override;

	std::unique_ptr<SmileParams> clone() const override;

	// analytical derivatives in k
	double dw(double k) const override;
	double d2w(double k) const override;
	double dwdtheta(double k) const; // dw/dtheta for Dupire local vol
	// gFunction is inherited

	// Gatheral constraints
	// 1. theta * phi(theta) * (1 + |rho|) < 4    (butterfly, Thm 4.3)
	// 2. eta * (1 + |rho|) <= 2                   (large-theta limit)

	bool satisfiesConstraints() const;
	bool satisfiesButterflyConstraint(double th) const;
};

// calibration
struct SsviCalibrationResult
{
	SsviParams params;
	std::vector<double> thetas;
	std::vector<double> maturities;
	double rmse = 0.0;
	bool converged = false;
	bool arbitrageFree = false; // Gatheral & Jacquier constraints satisfied

	// build VolatilitySurface from the fitted SSVI
	std::unique_ptr<VolatilitySurface> buildSurface(
		const std::vector<double> &strikeGrid,
		const std::vector<double> &forwards,
		const DiscountCurve &discountCurve) const;

	// analytical Dupire local vol -- no grid, no finite differences
	SsviSurface buildAnalyticalSurface(const std::vector<double> &forwards) const;
};

struct SsviCalibrationOptions
{
	LMOptions lmOpts = {};
	bool enforceMonotonicity = true; // PAVA on extracted thetas
	int nGridPoints = 3;	  // number of starting points for grid search (noting 0 = skip)
};

// this gives: mkt data -> per-slice SVI -> extract theta -> global SSVI -> result
SsviCalibrationResult calibrateSsvi(
	const std::vector<std::vector<double>> &strikesPerMaturity,
	const std::vector<std::vector<double>> &volsPerMaturity,
	const std::vector<double> &forwards,
	const std::vector<double> &maturities,
	const SsviCalibrationOptions &opts = {});
// ============================================================================
// SsviSurface — analytical Dupire local vol from SSVI parametric form
// ============================================================================

class SsviSurface
{
public:
	SsviSurface(const SsviCalibrationResult &result,
				const std::vector<double> &forwards);

	SsviSurface(SsviSurface&&) = default;
	SsviSurface& operator=(SsviSurface&&) = default;

	double localVolatility(double spot, double time) const;
	double impliedVolatility(double strike, double time) const;
	double localVariance(double k, double time) const; // k-space, for diagnostics

private:
	SsviParams _params;
	std::vector<double> _maturities, _thetas, _forwards;
	std::unique_ptr<InterpolationScheme> _thetaInterp;
	std::unique_ptr<InterpolationScheme> _forwardInterp;

	double thetaAt(double T) const;
	double dthetadT(double T) const;
	double forwardAt(double T) const;
};

#endif
