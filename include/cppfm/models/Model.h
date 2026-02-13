// #pragma once

#ifndef MODEL_H
#define MODEL_H

#include <memory>

// Forward declaration
class DiscountCurve;
class VolatilitySurface;

/**
 * TODO:
 *  Allows for extension to multidimensional models in the future
 */

/** NOTES:
 * Getters - if just access the data member, use i.e. rate()
 * 		   - if need to compute something, use i.e.
 * computeDupireLocalVolatility()
 *
 */

/**
 * ============================================================================
 * MODEL HIERARCHY - Pure Design
 * ============================================================================
 *
 * ModelBase:
 *   ├── Model1D: Single state variable models (S)
 *   │     ├── BlackScholesModel
 *   │     └── DupireModel
 *   │
 *   └── Model2D: Two state variable models (S, V)
 *         ── HestonModel
 *
 * ModelBase - Minimal common base for all financial models
 * Provides only the essential shared interface
 *
 * Abstract class = class that have at least 1 virtual pure method
 * "Incomplete" class -> cannot instantiate this class
 */

class ModelBase
{
public:
	ModelBase(double initValue);
	ModelBase(const ModelBase &model);
	ModelBase &operator=(const ModelBase &model);
	virtual ~ModelBase() = default; // ALWAYS DECLARE BASE CLASS DESTRUCTOR
									// VIRTUAL -> Avoid memory leak

	// Getters here
	inline double initValue() const
	{
		return _initValue;
	} // do not modify any data member

	// private:
protected:
	double _initValue; // S_0 - initial asset price
};

// ============================================================================
// 1D MODELS - Single state variable (asset price S)
// ============================================================================
// The Financial Models we consider are defined by Ito processes only
// Drift term and diffusion terms are functions of 2 variables - deterministic

/**
 * Model1D - Base class for models with single state variable
 * SDE: dS = μ(t, S)dt + σ(t, S)dW
 */
class Model1D : public ModelBase
{
public:
	Model1D(double initValue);
	virtual ~Model1D() = default;

	virtual double drift(double time,
						 double assetPrice) const = 0; // Virtual Pure method
	virtual double
	diffusion(double time, double assetPrice) const = 0; // Virtual Pure method

	// Note: Risk-free rates are accessed through DiscountCurve objects

	// RULE: If a pointer of this class is used as a data member of another
	// class [here PathSimulator] then we need to delegate the copy construction
	// and pointing to another method This method is: clone()
	virtual Model1D *clone() const = 0; // Virtual Pure method

	virtual bool operator==(const Model1D &model) const = 0; // == operator
};

// For backward compatibility - keep old name as alias
using Model = Model1D;

// ============================================================================
// 2D MODELS - Two state variables (asset price S, auxiliary variable V)
// ============================================================================
// 2D Models are defined by coupled Ito processes
// Both state variables evolve according to correlated SDEs

/**
 * Model2D - Base class for models with two state variables
 * SDE system:
 *   dS = μ_S(t, S, V)dt + σ_S(t, S, V)dW^S
 *   dV = μ_V(t, S, V)dt + σ_V(t, S, V)dW^V
 *   with correlation d⟨W^S, W^V⟩ = ρdt
 */
class Model2D : public ModelBase
{
public:
	Model2D(double initValue);
	virtual ~Model2D() = default;

	/**
	 * Drift coefficients for 2D system
	 * @return pair{μ_S, μ_V}
	 */
	virtual std::pair<double, double>
	drift2D(double time, double assetPrice,
			double auxVariable) const = 0; // Virtual Pure method

	/**
	 * Diffusion coefficients for 2D system
	 * @return pair{σ_S, σ_V}
	 */
	virtual std::pair<double, double>
	diffusion2D(double time, double assetPrice,
				double auxVariable) const = 0; // Virtual Pure method

	// RULE: Same clone() pattern as Model1D for polymorphic copying
	virtual Model2D *clone() const = 0; // Virtual Pure method
	virtual bool operator==(const Model2D &model) const = 0; // == operator

	/**
	 * Correlation coefficient between the two Brownian motions
	 * @return ρ in [-1, 1]
	 */
	virtual double correlation() const = 0; // Virtual Pure method
};

// ============================================================================
// 1D MODEL IMPLEMENTATIONS
// ============================================================================
/**
 *  Black-Scholes Model
 *  SDE: dS = rS dt + σS dW - constant volatility model
 */
class BlackScholesModel : public Model1D // "Is a" relationship
{
public:
	// Default constructor
	BlackScholesModel() = delete;
	// Constructor with parameters
	BlackScholesModel(double spot, const DiscountCurve &discountCurve,
					  double volatility);
	// Copy constructor
	BlackScholesModel(const BlackScholesModel &model);
	// Clone method
	BlackScholesModel *clone() const override;

	// Copy Assignment Operator
	BlackScholesModel &operator=(const BlackScholesModel &model);

	bool
	operator==(const Model1D &model) const override // required by base class
	{ // dynamic_cast - Derived1 vs. Derived2; safe downcasting in polymorphic
	  // hierarchies.
		// Check at runtime if a Model reference is actually a BlackScholesModel
		// object. CRTP - useful for static polymorphism, does not support
		// runtime polymorphism (i.e. lose ability to use base class
		// pointers/references for derived objects)

		// check if the Model reference passed as an argument is actually a
		// BlackScholesModel
		const auto *bsModel = dynamic_cast<const BlackScholesModel *>(
			&model);  // try to convert Model to BlackScholesModel
		if (!bsModel) // if the cast fails; ie not of type BlackScholesModel
			return false;
		// Compare private data members	of the two BlackScholesModel objects.
		// -> to access the members of the bsModel pointer
		// Check if _initvalue, _drift, _volatility are equal between current
		// object (this) and bsModel object
		return (_initValue == bsModel->_initValue) &&
			   (_drift == bsModel->_drift) &&
			   (_volatility ==
				bsModel->_volatility); // access members of the pointers
		// return is true - if the two objects are of the same type
		// BlackScholesModel and their data members have identical values
	} // need to be smarter

	// Destructor
	~BlackScholesModel() override = default;

	// Virtual -> You CAN override base class virtual methods
	// Virtual pure -> You HAVE TO override base class pure virtual methods
	double drift(double time, double assetPrice) const override;
	double diffusion(double time, double assetPrice) const override;
	// Getters for model parameters - accessors
	double drift() const { return _drift; }
	double volatility() const { return _volatility; }

private: // default constructor will call the default constructor for each
		 // data member
	// double _initValue;
	double _drift;
	double _volatility;
};

/*
 * Dupire Local Volatility Model
 * SDE: dS = rS dt + sigma(S, t) S dW
 */
class DupireModel : public Model
{
public:
	// ~DupireModel() override = default; // Destructor
	/**
	 * Constructor with volatility surface
	 * @param spot Initial spot price
	 * @param volSurface Volatility surface for local volatility
	 */
	DupireModel(double spot, const VolatilitySurface &volSurface);

	DupireModel(const DupireModel &model); // Copy constructor
	DupireModel *clone() const override;   // Clone method
	DupireModel &
	operator=(const DupireModel &model); // Copy assignment operator
	bool operator==(const Model1D &model) const override; // == operator

	// Virtual methods from base class
	double drift(double time, double assetPrice) const override;
	double diffusion(double time, double assetPrice) const override;
	// double drift() const;
	// double volatility() const; // Volatility is accessed through the
	// VolatilitySurface Getters
	const VolatilitySurface &volatilitySurface() const
	{
		return *_volSurfacePtr; // pointer to the VolatilitySurface
	}

	/**
	 * Get local volatility at given spot and time
	 * @param spot Current asset price
	 * @param time Current time
	 * @return Local volatility
	 */
	double localVolatility(double spot, double time) const;

private:
	// const VolatilitySurface* _volSurfacePtr; // remove the raw pointer
	std::unique_ptr<const VolatilitySurface> _volSurfacePtr;

	// internal Local Vol Solver as a data member
};

// ============================================================================
// 2D MODEL IMPLEMENTATIONS
// ============================================================================

/**
 * Heston Stochastic Volatility Model
 *
 * 2D SDE system (S: asset price, V: variance):
 *   dS = r*S*dt + √V*S*dW^S
 *   dV = κ(v̄ - V)*dt + σ_v*√V*dW^V
 *   d⟨W^S, W^V⟩ = ρ*dt
 */
class HestonModel : public Model2D
{
public:
	/**
	 * Constructor for Heston model
	 * @param spot S_0: Initial asset price
	 * @param discountCurve: Provides risk-free rate r(t)
	 * @param v0: Initial variance v(0) > 0
	 * @param kappa: κ > 0, mean reversion speed
	 * @param vbar: θ > 0, long-term variance mean
	 * @param sigma_v: σ_v > 0, volatility of volatility (vol-of-vol)
	 * @param rho: ρ in [-1,1], correlation between W_X and W_V
	 */
	HestonModel(double spot, const DiscountCurve &discountCurve, double v0,
				double kappa, double vbar, double sigma_v, double rho);

	HestonModel(const HestonModel &model); // copy constructor
	HestonModel *clone() const override;   // clone method
	HestonModel &
	operator=(const HestonModel &model); // copy assignment operator
	bool operator==(const Model2D &model)
		const override; // equality == operator; for polymorphic comparison
	~HestonModel() override = default; // destructor

	// Implement Model2D interface
	std::pair<double, double> drift2D(double time, double assetPrice,
									  double variance) const override;
	std::pair<double, double> diffusion2D(double time, double assetPrice,
										  double variance) const override;
	double correlation() const override { return _rho; }

	// Accessors for parameters
	inline double v0() const { return _v0; }
	inline double kappa() const { return _kappa; }
	inline double vbar() const { return _vbar; }
	inline double sigma_v() const { return _sigma_v; }
	inline double rho() const { return _rho; }
	inline double riskFreeRate() const { return _riskFreeRate; }

	/**
	 * Check if Feller condition is satisfied: 2κv̄ ≥ σ_v²
	 * If satisfied, the variance process stays strictly positive almost surely
	 * @return true if Feller condition is satisfied
	 */
	bool satisfiesFellerCondition() const;

private:
	void validateParameters() const; // validate parameters for Heston model

	double _riskFreeRate;
	double _v0;
	double _kappa;
	double _vbar;
	double _sigma_v;
	double _rho;
};

#endif // MODEL_H
