//
// Created by Patrik  Liba on 18/09/2025.
#include <cppfm/market/FinancialInstrument.h>


Option::Option(Option::Type type)
    : _type(type)
{
}


EuropeanOptionPayoff::EuropeanOptionPayoff(Type type)
    : Option(type) // call base constructor
{
}

EuropeanOptionPayoff::EuropeanOptionPayoff(Type type, double strike, double maturity)
    : Option(type), _strike(strike), _maturity(maturity)
{
    if (strike <= 0.0)
        throw "Strike must be positive";
    if (maturity <= 0.0)
        throw "Maturity must be positive";
}

Option::Type EuropeanOptionPayoff::type() const
{
    return _type;
}

EuropeanOptionPayoff * EuropeanOptionPayoff::clone() const
{
    return new EuropeanOptionPayoff(*this); // copy construct of current object
}

