#include "Model.h"
#include "VolatilitySurface.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

// ============================================================================
// ModelBase Implementation
// ============================================================================
ModelBase::ModelBase(double initValue)
	: _initValue(initValue)
{
}

ModelBase::ModelBase(const ModelBase& model)
	: _initValue(model._initValue)
{
}

ModelBase& ModelBase::operator=(const ModelBase& model)
{
	if (this != &model)
	{
		_initValue = model._initValue;
	}
	return *this;
}

// ============================================================================
// Model1D Implementation
// ============================================================================
Model1D::Model1D(double initValue)
	: ModelBase(initValue)
{
}

// ============================================================================
// Model2D Implementation
// ============================================================================
Model2D::Model2D(double initValue)
	: ModelBase(initValue)
{
}

// ============================================================================
// BlackScholesModel Implementation
// ============================================================================
BlackScholesModel::BlackScholesModel(double spot, const DiscountCurve& discountCurve, double volatility)
	: Model1D(spot), _drift(discountCurve.rate()), _volatility(volatility)
	// Actually it calls the copy constructor of each data member
{
}

BlackScholesModel::BlackScholesModel(const BlackScholesModel& model)
	: Model1D(model), _drift(model._drift), _volatility(model._volatility)
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
		Model1D::operator=(model);
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
	: Model1D(spot), _volSurfacePtr(volSurface.clone().release())
{
/* DupireModel volSurface.clone() makes independent copy of the VolatilitySurface
 * DupireModel has its own VolatilitySurface copy, independent of the original VolatilitySurface
 * Risk-free rate is extracted from the VolatilitySurface's DiscountCurve
 * NOTE: clone() returns unique_ptr<VolatilitySurface>, we release() to get raw pointer
 * and unique_ptr<const VolatilitySurface> constructor takes ownership
 */
}


DupireModel::DupireModel(const DupireModel& model)
	: Model1D(model), _volSurfacePtr(model._volSurfacePtr->clone().release())
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
		Model1D::operator=(model);
		// Exception-safe assignment: clone() returns unique_ptr, release() to get raw pointer
		_volSurfacePtr.reset(model._volSurfacePtr->clone().release()); // handles delete and assignment
	}
	return *this;
}

bool DupireModel::operator==(const Model1D& model) const
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

// ============================================================================
// HestonModel Implementation
// ============================================================================

HestonModel::HestonModel(double spot, const DiscountCurve& discountCurve, double v0, double kappa, double vbar, double sigma_v, double rho)
		: Model2D(spot), _riskFreeRate(discountCurve.rate()), _v0(v0), _kappa(kappa), _vbar(vbar), _sigma_v(sigma_v), _rho(rho)
{
	validateParameters();
}

HestonModel::HestonModel(const HestonModel& model)
	: Model2D(model), _riskFreeRate(model._riskFreeRate), _v0(model._v0), _kappa(model._kappa), _vbar(model._vbar), _sigma_v(model._sigma_v), _rho(model._rho)
{
}

HestonModel* HestonModel::clone() const
{
	return new HestonModel(*this);
}

HestonModel& HestonModel::operator=(const HestonModel& model)
{
	if (this != &model) {
		Model2D::operator=(model);
		_riskFreeRate = model._riskFreeRate;
		_v0 = model._v0;
		_kappa = model._kappa;
		_vbar = model._vbar;
		_sigma_v = model._sigma_v;
		_rho = model._rho;
	}
	return *this;
}

bool HestonModel::operator==(const Model2D& model) const 
{
	const auto* hestonModel = dynamic_cast<const HestonModel*>(&model);
	if (!hestonModel) {
		return false;
	}
	// Compare all parameters with tolerance for floating-point comparison
	const double tol = 1e-15;
	return (std::abs(_initValue - hestonModel->_initValue) < tol) && 
	       (std::abs(_riskFreeRate - hestonModel->_riskFreeRate) < tol) &&
	       (std::abs(_v0 - hestonModel->_v0) < tol) &&
	       (std::abs(_kappa - hestonModel->_kappa) < tol) &&
	       (std::abs(_vbar - hestonModel->_vbar) < tol) &&
	       (std::abs(_sigma_v - hestonModel->_sigma_v) < tol) &&
	       (std::abs(_rho - hestonModel->_rho) < tol);
}

// ============================================================================
// Heston 2D Interface Implementation
// ============================================================================

std::pair<double, double> HestonModel::drift2D(double time, double assetPrice, double variance) const
{
	// Asset price drift: dS/dt = r*S  
	double drift_S = _riskFreeRate * assetPrice;
	
	// Variance drift (CIR): dV/dt = κ(v̄ - V)
	double drift_V = _kappa * (_vbar - variance);
	
	return {drift_S, drift_V};
}

std::pair<double, double> HestonModel::diffusion2D(double time, double assetPrice, double variance) const
{
	// Protect against negative variance (numerical safety)
	double V_safe = std::max(0.0, variance);
	double sqrt_V = std::sqrt(V_safe);
	
	// Asset price diffusion: σ_S = √V * S
	double diffusion_S = sqrt_V * assetPrice;
	
	// Variance diffusion (CIR): σ_V = σ_v * √V  
	double diffusion_V = _sigma_v * sqrt_V;
	
	return {diffusion_S, diffusion_V};
}

bool HestonModel::satisfiesFellerCondition() const
{
	// Feller condition: 2κv̄ ≥ σ_v²
	// If satisfied, the CIR process for variance stays strictly positive
	return (2.0 * _kappa * _vbar) >= (_sigma_v * _sigma_v);
}


void HestonModel::validateParameters() const
{
	// ========================================================================
	// PARAMETER VALIDATION - Basic Constraints (THROW on failure)
	// ========================================================================
	
	// Validate initial variance v0 > 0
	if (_v0 <= 0.0) {
		throw std::invalid_argument("HestonModel: Initial variance v0 must be strictly positive. Got v0 = " + std::to_string(_v0));
	}
	
	// Validate mean reversion speed κ > 0
	if (_kappa <= 0.0) {
		throw std::invalid_argument("HestonModel: Mean reversion speed kappa (κ) must be strictly positive. Got κ = " + std::to_string(_kappa));
	}
	
	// Validate long-term variance v̄ > 0
	if (_vbar <= 0.0) {
		throw std::invalid_argument("HestonModel: Long-term variance vbar (v̄) must be strictly positive. Got v̄ = " + std::to_string(_vbar));
	}
	
	// Validate vol-of-vol σ_v > 0
	if (_sigma_v <= 0.0) {
		throw std::invalid_argument("HestonModel: Vol-of-vol sigma_v (σ_v) must be strictly positive. Got σ_v = " + std::to_string(_sigma_v));
	}
	
	// Validate correlation ρ ∈ [-1, 1]
	if (_rho < -1.0 || _rho > 1.0) {
		throw std::invalid_argument("HestonModel: Correlation rho (ρ) must be in [-1, 1]. Got ρ = " + std::to_string(_rho));
	}
	
	// ========================================================================
	// FELLER CONDITION - Warn if violated but ALLOW construction
	// ========================================================================
	// Feller condition: 2κv̄ ≥ σ_v²
	// If violated, variance process V(t) can reach zero with positive probability
	// Some discretization schemes (e.g., BK/QE, BK/TG) can handle this gracefully
	
	if (!satisfiesFellerCondition()) {
		std::cerr << "WARNING [HestonModel]: Feller condition NOT satisfied!\n"
		          << "  Required: 2*κ*v̄ ≥ σ²_v\n"
		          << "  Current:  2*" << _kappa << "*" << _vbar << " = " << (2.0 * _kappa * _vbar)
		          << " < σ²_v = " << (_sigma_v * _sigma_v) << "\n"
		          << "  Implication: Variance process can reach zero.\n"
		          << "  Use appropriate discretization schemes (BK/QE, BK/TG) that handle V->0.\n"
		          << std::endl;
	}
}