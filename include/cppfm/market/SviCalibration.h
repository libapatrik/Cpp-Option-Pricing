#ifndef CPPFM_SVICALIBRATION_H
#define CPPFM_SVICALIBRATION_H

#include "cppfm/calibration/Optimizer.h"
#include "cppfm/market/DiscountCurve.h"
#include "cppfm/market/VolatilitySurface.h"
#include <memory>
#include <vector>

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

	// build VolatilitySurface from the fitted SSVI
	std::unique_ptr<VolatilitySurface> buildSurface(
		const std::vector<double> &strikeGrid,
		const std::vector<double> &forwards,
		const DiscountCurve &discountCurve) const;
};

SsviCalibrationResult calibrateSsvi(
	const std::vector<std::vector<double>> &strikesPerMaturity,
	const std::vector<std::vector<double>> &volsPerMaturity,
	const std::vector<double> &forwards,
	const std::vector<double> &maturities,
	const LMOptions &opts = {});

#endif
