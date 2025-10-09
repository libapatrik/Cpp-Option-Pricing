#ifndef CPPFM_DISCOUNTCURVE_H
#define CPPFM_DISCOUNTCURVE_H

#include <cmath>

class DiscountCurve // for modelling with deterministic rates
{
public:
    virtual double discount(double time) const = 0;
    virtual DiscountCurve* clone() const = 0;
    virtual double rate() const = 0;  // Virtual getter for rate information
    virtual ~DiscountCurve() = default;
};

class FlatDiscountCurve : public DiscountCurve
{
public:
    FlatDiscountCurve(double rate);

    double discount(double time) const override;

    double rate() const override;
    
    FlatDiscountCurve* clone() const override;

private:
    double _rate;
};

#endif //CPPFM_DISCOUNTCURVE_H
