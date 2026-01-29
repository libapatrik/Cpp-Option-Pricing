#include <cppfm/market/DiscountCurve.h>
#include <cmath>


FlatDiscountCurve::FlatDiscountCurve(double rate)
    : _rate(rate)
{
}

double FlatDiscountCurve::discount(double time) const
{
    return std::exp(-_rate * time);
}

double FlatDiscountCurve::rate() const
{
    return _rate;
}

FlatDiscountCurve* FlatDiscountCurve::clone() const
{
    return new FlatDiscountCurve(*this);
}


