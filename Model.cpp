#include "Model.h"
#include "Utils.h"
#include "VolatilitySurface.h"
#include <cmath>

// Model
Model::Model(double initValue)
	: _initValue(initValue)
{
}

Model::Model(const Model& model)
	: _initValue(model._initValue)
{
}

Model& Model::operator=(const Model& model)
{
	if (this != &model)
	{
		_initValue = model._initValue;
	}
	return *this;
}

// ============================================================================
// BlackScholesModel Implementation
// ============================================================================
BlackScholesModel::BlackScholesModel(double spot, const DiscountCurve& discountCurve, double volatility)
	: Model(spot), _drift(discountCurve.rate()), _volatility(volatility)
	// Actually it calls the copy constructor of each data member
{
}

BlackScholesModel::BlackScholesModel(const BlackScholesModel& model)
	: Model(model), _drift(model._drift), _volatility(model._volatility)
{
}

BlackScholesModel * BlackScholesModel::clone() const
{
	return new BlackScholesModel(*this);
	// create new pointer copy construct of current object -> returns a copy of this method
}


BlackScholesModel& BlackScholesModel::operator=(const BlackScholesModel & model)
{
	if (this != &model)
	{
		Model::operator=(model);
		_drift = model._drift;
		_volatility = model._volatility;
	}

	return *this;
}


double BlackScholesModel::drift(double time, double assetPrice) const
{
	return _drift * assetPrice;
}

double BlackScholesModel::diffusion(double time, double assetPrice) const
{
	return _volatility * assetPrice;
}



// ============================================================================
// DupireModel Implementation
// ============================================================================
DupireModel::DupireModel(double spot, const VolatilitySurface& volSurface)
	: Model(spot), _volSurfacePtr(volSurface.clone().release())
{
/* DupireModel volSurface.clone() makes independent copy of the VolatilitySurface
 * DupireModel has its own VolatilitySurface copy, independent of the original VolatilitySurface
 * Risk-free rate is extracted from the VolatilitySurface's DiscountCurve
 * NOTE: clone() returns unique_ptr<VolatilitySurface>, we release() to get raw pointer
 * and unique_ptr<const VolatilitySurface> constructor takes ownership
 */
}


DupireModel::DupireModel(const DupireModel& model)
	: Model(model), _volSurfacePtr(model._volSurfacePtr->clone().release())
{
}

DupireModel* DupireModel::clone() const
{
	return new DupireModel(*this);
}

DupireModel& DupireModel::operator=(const DupireModel& model)
{
	if (this != &model)
	{
		Model::operator=(model);
		// Exception-safe assignment: clone() returns unique_ptr, release() to get raw pointer
		_volSurfacePtr.reset(model._volSurfacePtr->clone().release()); // handles delete and assignment
	}
	return *this;
}

bool DupireModel::operator==(const Model& model) const
{
	const auto* dupireModel = dynamic_cast<const DupireModel*>(&model);
	if (!dupireModel)
		return false;
	
	return (_initValue == dupireModel->_initValue) && 
		   (*_volSurfacePtr == *dupireModel->_volSurfacePtr);
		   // dereferencing the unique_ptr to get the VolatilitySurface object
}

double DupireModel::drift(double time, double assetPrice) const
{
	// LECTURE NOTES CORRECTION: Use DiscountCurve properly for time-dependent rates
	// For time-dependent rates, we need the instantaneous rate r(t) at time t
	// This is computed as: r(t) = -d/dt[log(B(t))] where B(t) is the discount factor
	// For numerical approximation: r(t) ≈ -[log(B(t+ε)) - log(B(t))]/ε
	// const double eps = 1e-6;
	const DiscountCurve& discountCurve = _volSurfacePtr->discountCurve();
	// double discount_t = discountCurve.discount(time);
	// double discount_t_plus = discountCurve.discount(time + eps);
	// double riskFreeRate = -std::log(discount_t_plus / discount_t) / eps;
	double riskFreeRate = discountCurve.instantaneousRate(time);
	
	return riskFreeRate * assetPrice; // under risk-neutral measure -> dS/S = r(t)dt + σ(S,t)dW
}

double DupireModel::diffusion(double time, double assetPrice) const
{
	// Local volatility from Dupire formula
	double localVol = localVolatility(assetPrice, time);
	return localVol * assetPrice;
}



double DupireModel::localVolatility(double spot, double time) const
{
	// Delegate the VolatilitySurface's implementation of the Dupire formula.
	return _volSurfacePtr->localVolatility(spot, time);
/* 
 * Market Data: VolatilitySurface contains market implied volatilities
 * Local Volatility: Dupire formula converts implied to local volatility
 * Pricing: Local volatility used in SDE for path simulation
 */
}
