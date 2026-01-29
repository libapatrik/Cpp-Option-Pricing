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
        // Use central difference for second-order accuracy: r(t) = -d/dt[log(B(t))]
        // Central difference has O(eps^2) error vs forward difference O(eps)
        
        // Handle boundary case at t=0
        if (time < eps) {
            // Use forward difference at t=0 (no choice)
            double discount_t = discount(time);
            double discount_t_plus = discount(time + eps);
            return -std::log(discount_t_plus / discount_t) / eps;
        }
        
        // Central difference for interior points
        double discount_t_plus = discount(time + eps);
        double discount_t_minus = discount(time - eps);
        return -std::log(discount_t_plus / discount_t_minus) / (2.0 * eps);
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
