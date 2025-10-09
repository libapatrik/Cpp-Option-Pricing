//
// Created by Patrik  Liba on 18/09/2025.
//

#include "Pricer.h"
#include "Utils.h"
#include "DiscountCurve.h"
#include <cmath>

using namespace std;

Pricer::Pricer(const Model &model, const DiscountCurve &discountCurve)
    :   _modelPtr(model.clone()), _discountCurvePtr(discountCurve.clone()) // add to memory
{
}

Pricer::~Pricer()
{   // delete the memory allocated with close
    if (_modelPtr) {  // Check if Null before deleting pointer
        delete _modelPtr;
        _modelPtr = nullptr;  // Prevent double deletion 
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
}

BlackScholesPricer::~BlackScholesPricer() = default;
BlackScholesPricer * BlackScholesPricer::clone() const
{
    return new BlackScholesPricer(*this);
}

static double normCdf(double x)
{
    return Utils::stdNormCdf(x);
}

double BlackScholesPricer::d1(double spot, double strike, double r, double sigma, double T)
{
    return (log(spot / strike) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
}

double BlackScholesPricer::d2(double d1Val, double sigma, double T)
{
    return d1Val - sigma * sqrt(T);
}

double BlackScholesPricer::price(const EuropeanOptionPayoff &option) const
{
    if (_bsModelPtr == nullptr)
        throw "BlackScholesPricer requires BlackScholesModel";

    double S0 = _bsModelPtr->initValue();
    double T = option.maturity();
    // discount = e^(-r*T)  =>  r = -ln(P(0,T))/T
    double r = -log(_discountCurvePtr->discount(T)) / T; // infer the rate from the discount curve
    double sigma = _bsModelPtr->getVolatility();
    double K = option.strike();

    if (T <= 0.0)
        return max( (option.getType() == Option::Type::Call ? S0 - K : K - S0), 0.0 ); // condition ? value_if_true : value_if_false

    double d1Val = d1(S0, K, r, sigma, T);
    double d2Val = d2(d1Val, sigma, T);

    if (option.getType() == Option::Type::Call)
    {
        return S0 * normCdf(d1Val) - K * exp(-r * T) * normCdf(d2Val);
    }
    else // Put via put-call parity or direct formula
    {
        return K * exp(-r * T) * normCdf(-d2Val) - S0 * normCdf(-d1Val);
    }
}


MonteCarloPricer::MonteCarloPricer(const Model& model, const DiscountCurve& discountCurve, const PathSimulator& simulator, size_t numPaths) // with params
    : Pricer(model, discountCurve),  _numPaths(numPaths), _simulatorPtr(&simulator) // explicitly passing the model to base class constructor
{  // NOTE: PathSimulator may need clone method for safety
}

MonteCarloPricer * MonteCarloPricer::clone() const
{
    return new MonteCarloPricer(*this); // destroyed by Pricer destructor!
}

MonteCarloPricer::~MonteCarloPricer() = default;

double MonteCarloPricer::price(const EuropeanOptionPayoff &option) const
{
    // Simulate paths, sum the payoffs, take average, discount
    double T = option.maturity();
    double r = -log(_discountCurvePtr->discount(T)) / T;
    // if rate is stochastic, then switch to pathwise exp(int_0^T r(t) dt); Model.h - HullWhiteModel

    size_t numPaths = _numPaths > 0 ? _numPaths : 10000;
    double sumPayoff = 0.0;
    for (size_t p = 0; p < numPaths; ++p)
    {
        vector<double> path = _simulatorPtr->path();
        double ST = path.back(); // last element of the path
        // Compute payoff at maturity
        double payoff = (option.getType() == Option::Type::Call) // condition ? value_if_true : value_if_false
            ? max(ST - option.strike(), 0.0)
            : max(option.strike() - ST, 0.0);
        sumPayoff += payoff; // sum of payoffs
    }

    double avgPayoff = sumPayoff / static_cast<double>(numPaths); // average payoff
    return exp(-r * T) * avgPayoff; // discount 
}
