//
// Created by Patrik  Liba on 18/09/2025.
//

#include "Pricer.h"

Pricer::Pricer(const Model &model)
    :   _modelPtr(model.clone()) // add to memory
{
}

Pricer::~Pricer()
{           // delete the memory allocated with close
    delete _modelPtr; // the 'new' called pointer in Pricer is deleted here - when the Pricer object is destroyed
}

double Pricer::price(const EuropeanOptionPayoff &option, double spot, double volatility, double r, double T,
                     int numPaths) const
{
}

MonteCarloPricer::MonteCarloPricer(const Model& model)
    : Pricer(model) // explicitly passing the model to base class constructor
{
}

MonteCarloPricer * MonteCarloPricer::clone() const
{
    return new MonteCarloPricer(*this); // destroyed by Pricer destructor!
}

MonteCarloPricer::~MonteCarloPricer() = default;

double MonteCarloPricer::price(const EuropeanOptionPayoff &option, double spot, double volatility, double r, double T, int numPaths) const
{
    // Has no protecteds (yet)
}
