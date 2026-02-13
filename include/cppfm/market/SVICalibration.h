#ifndef CPPFM_SVICALIBRATION_H
#define CPPFM_SVICALIBRATION_H

#include "cppfm/calibration/Optimizer.h"
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
	double m = 0.0; // not shifted
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

#endif
