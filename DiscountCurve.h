#ifndef CPPFM_DISCOUNTCURVE_H
#define CPPFM_DISCOUNTCURVE_H

#include <cmath>

class DiscountCurve // for modelling with deterministic rates
{
public:
    virtual double discount(double time) const = 0;
    virtual DiscountCurve* clone() const = 0;
    virtual double rate() const = 0;  // Virtual getter for rate information

    virtual double instantaneousRate(double time, double eps=1e-6) const {
        double discount_t = discount(time);
        double discount_t_plus = discount(time + eps);
        return -std::log(discount_t_plus / discount_t) / eps;
    }

    virtual ~DiscountCurve() = default;
};

class FlatDiscountCurve : public DiscountCurve
{
public:
    FlatDiscountCurve(double rate);
    double discount(double time) const override;
    double rate() const override;
    double instantaneousRate(double time, double eps=1e-6) const override {
        return _rate; 
    }
    FlatDiscountCurve* clone() const override;
    

private:
    double _rate;
};

#endif //CPPFM_DISCOUNTCURVE_H
