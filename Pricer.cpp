//
// Created by Patrik  Liba on 18/09/2025.
//

#include "Pricer.h"
#include "Utils.h"
#include "DiscountCurve.h"
#include "BlackScholesFormulas.h"
#include <cmath>



Pricer::Pricer(const Model &model, const DiscountCurve &discountCurve)
    :   _modelPtr(model.clone()), _discountCurvePtr(discountCurve.clone()) // add to memory
{
}

Pricer::~Pricer()
{   // delete the memory allocated with close
    if (_modelPtr) {  // Check if Null before deleting pointer
        delete _modelPtr;
        _modelPtr = nullptr;  // Prevent double deletion; ensuring that the pointer does not point to deallocated member location
    }
    if (_discountCurvePtr) {  // Check if Null before deleting pointer
        delete _discountCurvePtr;
        _discountCurvePtr = nullptr;  // Prevent double deletion
    }
}

BlackScholesPricer::BlackScholesPricer() = default;

BlackScholesPricer::BlackScholesPricer(const BlackScholesModel &model, const DiscountCurve &discountCurve)
    : Pricer(model, discountCurve)
{
    _bsModelPtr = dynamic_cast<const BlackScholesModel*>(_modelPtr); // set it once
    if (!_bsModelPtr)
        throw std::runtime_error("BlackScholesPricer requires BlackScholesModel");
}

BlackScholesPricer::~BlackScholesPricer() = default;
BlackScholesPricer * BlackScholesPricer::clone() const
{
    return new BlackScholesPricer(*this);
}

double BlackScholesPricer::price(const EuropeanOptionPayoff &option) const
{
    if (_bsModelPtr == nullptr)
        throw std::runtime_error("BlackScholesPricer requires BlackScholesModel");

    return BlackScholesFormulas::price(
        _bsModelPtr->initValue(), 
        option.strike(), 
        *_discountCurvePtr,
        _bsModelPtr->volatility(), 
        option.maturity(), 
        option.type()
    );
}


MonteCarloPricer::MonteCarloPricer(const Model& model, const DiscountCurve& discountCurve, const PathSimulator& simulator, size_t numPaths) // with params
    : Pricer(model, discountCurve),  _numPaths(numPaths), _simulatorPtr(&simulator) // explicitly passing the model to base class constructor
{  // NOTE: PathSimulator may need clone method for safety
}

MonteCarloPricer* MonteCarloPricer::clone() const
{
    return new MonteCarloPricer(*this); // destroyed by Pricer destructor!
}

MonteCarloPricer::~MonteCarloPricer() = default;

double MonteCarloPricer::price(const EuropeanOptionPayoff &option) const
{
    // Simulate paths, sum the payoffs, take average, discount
    double T = option.maturity();
    double r = -std::log(_discountCurvePtr->discount(T)) / T;
    // if rate is stochastic, then switch to pathwise exp(int_0^T r(t) dt); Model.h - HullWhiteModel

    size_t numPaths = _numPaths > 0 ? _numPaths : 10000;
    double sumPayoff = 0.0;
    for (size_t p = 0; p < numPaths; ++p)
    {
        std::vector<double> path = _simulatorPtr->path();
        double ST = path.back(); // last element of the path
        // Compute payoff at maturity
        double payoff = (option.type() == Option::Type::Call) // condition ? value_if_true : value_if_false
            ? std::max(ST - option.strike(), 0.0)
            : std::max(option.strike() - ST, 0.0);
        sumPayoff += payoff; // sum of payoffs
    }

    double avgPayoff = sumPayoff / static_cast<double>(numPaths); // average payoff
    return std::exp(-r * T) * avgPayoff; // discount 
}


// double FDPricer