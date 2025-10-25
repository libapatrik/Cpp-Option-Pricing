// #pragma once

#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include <vector>
#include <memory>

// Forward declaration
class DiscountCurve;

/** TODO:
		Can delegate the computation to the DiscountCurve; just reuse the ε - eps
		Add Black-Scholes choice for put call.
		Create h/cpp for finite difference methods and finite element methods
		
		Don't do the numerical approximation in the model; delegate to the DiscountCurve
		DupireModel
		VolatilitySurface -> Greeks -> Calibration [+interpolation methods] = Study of Volatility
		Pricing under DupireModel

		VolSurface = implied from market data/prices
		Give me prices/data I will reconstruct the VolSurf which is being priced in by the traders

		yield-curve and VolSurface = backbone of the market data
		IVolSurf is viewed as initial vol for the SVM
		DiscountCurve is viewed as initial rate for StochasticIntRates

		We assume we got data and so we can calibrate - we just design the class/blueprints
		Most important: How we will represent the data; list of tenors, strikes and ImpVols matrix [raw data]
		Model needs to have continuous -> ex/interpolations schemes! Determines the robustness of the model

		HOMEWORK:
			1. See the notes: interpolation cubic spline - along strikes for 1 tenor, then do Thomas algo per
			VolSurface - give me a value I need to interpolate, then go to InterpScheme and get it done.
			InterpolationSchemes.hpp/cpp with multiple methods
		Careful: Cubic spline only works when we interpolate the strikes, interpolating along the tenor we will need information set of curves and interpolate in between the curves.

		Market data is for market data team
		Pricing lib must have proper yield-curve, IV-surface, then interpolate and get the model, then becomes the input for my calibration.
		What is market data, modelling assumption + interpolations!

		Stochastic Interest Rates

		HestonModel

		Generalise the BlackScholesModel for non-constant rates

	PDE: 
		Grid construction
		Boundary conditions (Dirichlet, Neumann, Robin)
		Numerical derivatives (1st, 2nd order, mixed)
		Time-stepping schemes (Explicit, Implicit, Crank-Nicolson)
		Tridiagonal solvers (Thomas algorithm)
		Stability analysis
	
	American option pricing:

*/

/** NOTES:
 * Getters - if just access the data member, use i.e. rate()
 * 		   - if need to compute something, use i.e. computeDupireLocalVolatility()
 *
 */


// Abstract class = class that have at least 1 virtual pure method
// "Incomplete" class -> cannot instantiate this class
class Model
{
public:
	Model(double initValue);
	Model(const Model& model);
	Model& operator=(const Model& model);
	virtual bool operator==(const Model& model) const = 0; // == operator
	virtual ~Model() = default; // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

	virtual double drift(double time, double assetPrice) const = 0; // Virtual Pure method
	virtual double diffusion(double time, double assetPrice) const = 0; // Virtual Pure method
	
	// Note: Risk-free rates are accessed through DiscountCurve objects
	// No virtual methods needed for risk-free rate or volatility so don't need this ?

	// RULE: If a pointer of this class is used as a data member of another class [here PathSimulator]
	// then we need to delegate the copy construction and pointing to another method
	// This method is: clone()

	virtual Model* clone() const = 0; // Virtual Pure method

	// Getters here
	inline double initValue() const // do not modify any data member
	{
		return _initValue;
	}

//private:
protected:
	double _initValue;
};


// The Financial Models we consider are defined by Ito processes only
// Drift term and diffusion terms are functions of 2 variables - deterministic

class BlackScholesModel : public Model  // "Is a" relationship 
{
public:
	// Default constructor
	BlackScholesModel() = delete;
	// Constructor with parameters
	BlackScholesModel(double spot, const DiscountCurve& discountCurve, double volatility);
	// Copy constructor
	BlackScholesModel(const BlackScholesModel& model);
	// Clone method
	BlackScholesModel* clone() const override;

	// Copy Assignment Operator
	BlackScholesModel& operator=(const BlackScholesModel& model);

	bool operator==(const Model& model) const override // required by base class
	{	// dynamic_cast - Derived1 vs. Derived2; safe downcasting in polymorphic hierarchies.
		// Check at runtime if a Model reference is actually a BlackScholesModel object.
		// CRTP - useful for static polymorphism, does not support runtime polymorphism
		// (i.e. lose ability to use base class pointers/references for derived objects)

		// check if the Model reference passed as an argument is actually a BlackScholesModel
		const auto* bsModel = dynamic_cast<const BlackScholesModel*>(&model); // try to convert Model to BlackScholesModel
		if (!bsModel) // if the cast fails; ie not of type BlackScholesModel
			return false;
		// Compare private data members	of the two BlackScholesModel objects.
		// -> to access the members of the bsModel pointer
		// Check if _initvalue, _drift, _volatility are equal between current object (this) and bsModel object
		return (_initValue == bsModel->_initValue) && (_drift == bsModel->_drift) && (_volatility == bsModel->_volatility); // access members of the pointers
		// return is true - if the two objects are of the same type BlackScholesModel and their data members have identical values
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


private: // default constructor will call the default constructor for each data member
	//double _initValue;
	double _drift;
	double _volatility;
};


// For Stochastic interest rate models: Hull-White model
class HullWhiteModel : public Model 
{   // for modelling with stochastic interest rate
	// simulate, short rate r(t) as Hull-White process
	// in MC discounting will be pathwise exp(int_0^T r(t) dt)
public:

protected:

};

// Forward declaration
class VolatilitySurface;
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
	DupireModel(double spot, const VolatilitySurface& volSurface);

	DupireModel(const DupireModel& model); // Copy constructor
	DupireModel* clone() const override; // Clone method
	DupireModel& operator=(const DupireModel& model); // Copy assignment operator
	bool operator==(const Model& model) const override; // == operator
	
	// Virtual methods from base class
	double drift(double time, double assetPrice) const override;
	double diffusion(double time, double assetPrice) const override;
	// double drift() const;
	// double volatility() const; // Volatility is accessed through the VolatilitySurface
	
	// Getters
	const VolatilitySurface& volatilitySurface() const 
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


};



#endif //MODEL_H
