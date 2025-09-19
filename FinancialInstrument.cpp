//
// Created by Patrik  Liba on 18/09/2025.
#include "FinancialInstrument.h"


EuropeanOptionPayoff::EuropeanOptionPayoff(Type type)
    : Option(type) // call base constructor
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

