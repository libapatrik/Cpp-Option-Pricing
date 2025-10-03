//
// Created by Patrik  Liba on 18/09/2025.
#include "FinancialInstrument.h"


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
}

Option::Type EuropeanOptionPayoff::getType() const
{
    return _type;
}

EuropeanOptionPayoff * EuropeanOptionPayoff::clone() const
{
    return new EuropeanOptionPayoff(*this); // copy construct of current object
}

