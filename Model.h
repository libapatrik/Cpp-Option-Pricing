// #pragma once

#ifndef MODEL_HPP
#define MODEL_HPP

#include <cmath>
#include <vector>

/* Notes
TO DO:
1. override specifier in BlackScholesModel::operator==
The base class expects a reference to Model, while the derived class uses a reference to BlackScholesModel.
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
	BlackScholesModel(double spot, double mu, double sigma);
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

	// Getters 
	inline double riskFreeRate() const 
	{ 
		return _drift; 
	}
	
	inline double volatility() const 
	{ 
		return _volatility; 
	}

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


/* TODO:
	DupireModel


	HestonModel - time grid, RNG
		COS method - to invert the ch.f.
		Broadie-Kaya exact simulation
		CIR Integrated Variance process

	Stochastic Interest rate
*/

#endif
