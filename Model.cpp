#include "Model.hpp"
#include "Utils.h"
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

// BlackScholesModel
BlackScholesModel::BlackScholesModel(double spot, double mu, double sigma)
	: Model(spot), _drift(mu), _volatility(sigma)
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


/*double Model::drift(double time, double assetPrice) const
{
	return 0.0;
}

double Model::diffusion(double time, double assetPrice) const
{
	return 0.0;
}*/



// PathSimulator
PathSimulator::PathSimulator(const double &initValue, const std::vector<double> &timePoints, const Model &model)
	: _initValue(initValue), _timePoints(timePoints), _model(model.clone())
{
}


// BlackScholesFormula
// TODO: declare optType, assetPrice etc, as private members of a class Model? Or BlackScholesModel? Access them via getters
double BlackScholesFormula(double& drift, double& volatility, double& spot, double& time, double& assetPrice, double& strike, std::string& optType)
{
	// Implement the Black-Scholes formula here
	double tmp = volatility * sqrt(time);
	double d1 = (log(spot/strike) + (drift + 0.5 * pow(volatility, 2)) * time) / tmp;
	double d2 = d1 - tmp;
	return spot * Utils::stdNormCdf(d1) - strike * exp(-drift * time) * Utils::stdNormCdf(d2);
}

/*
Under each model - i.e. BSM, Dupire, Heston etc.
Have different option classes: Class EuropeanOption, Class AmericanOption, Class AsianOption, Class BarrierOption, etc.?
Puts and calls.
Then for these options we could have different pricing methods: BlackScholesFormula, MonteCarloSimulation, FiniteDifferenceMethod, etc.?
*/