#include "Model.hpp"


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