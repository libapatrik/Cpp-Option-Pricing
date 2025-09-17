// #pragma once

#ifndef MODEL_HPP
#define MODEL_HPP

#include <cmath>
#include <vector>

/* NOTES:
Base Class Model
	Derived Class BlackScholesModel
	Derived Class DupireModel
	Derived Class HestonModel

Base Class Financial Instrument
	Derived BASE Class Option
		Derived Class EuropeanOption
		Derived Class AmericanOption
		Derived Class AsianOption

	Derive Class Other Instruments

Base Class PricingMethod
	Derived Class MonteCarloPricing
	Derived Class FiniteDifferencePricing
	Derived Class COSMethodPricing

Base Class PathSimulator
	Derived Class Euler-Maruyana PathSimulator
	Derived Class Milstein PathSimulator
	Derived Class Broadie-Kaya PathSimulator scheme


*/


// Abstraction
// What is a model?

// 1) Initial Value
// 2) Drift term and volatility term


// Abstract class = class that have at least 1 virtual pure method
// "Incomplete" class -> cannot instantiate this class
class Model
{
public:
	Model(double initValue);
	Model(const Model& model);
	Model& operator=(const Model& model);
	virtual ~Model() = default; // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

	virtual double drift(double time, double assetPrice) const = 0; // Virtual Pure method
	virtual double diffusion(double time, double assetPrice) const = 0; // Virtual Pure method

	// RULE: If a pointer of this class is used as a data member of another class [here PathSimulator]
	// then we need to delegate the copy construction and pointing to another method
	// This method is: clone()

	virtual Model* clone() const = 0; // Virtual Pure method

	// Getters here
	inline double initValue() const { // do not modify any data member
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

	// Destructor
	~BlackScholesModel() override = default;

	// Virtual -> You CAN override base class virtual methods
	// Virtual pure -> You HAVE TO override base class pure virtual methods
	double drift(double time, double assetPrice) const override;
	double diffusion(double time, double assetPrice) const override;

private: // default constructor will call the default constructor for each data member
	//double _initValue;
	double _drift;
	double _volatility;
};




#endif
