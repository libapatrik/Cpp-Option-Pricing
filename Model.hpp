// #pragma once

#ifndef MODEL_HPP
#define MODEL_HPP

#include <cmath>
#include <vector>


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

	virtual Model* clone() const = 0; // Virtual Pure method


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


class PathSimulator
{
public:
	// Constructor
	PathSimulator(const double& initValue,
					const std::vector<double>& timePoints,
					const Model& model);

	// Destructor
	virtual ~PathSimulator();

protected:
	double _initValue;
	std::vector<double> _timePoints;
	const Model* _model; // Create Pure Virtual method inside base class Model

};

// We need to delete the pointer to the _model to avoid memory leak
inline PathSimulator::~PathSimulator() {
	delete _model;

	// Calls Destructor methond of Model to destroy *_model
	// Model's destructor is virtual -> Calls the derived class destructor if needed
	// Make sure that the right Model Derived class destructor is called

}


#endif