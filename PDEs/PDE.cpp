#include "PDE.h"
#include <stdexcept>

// ============================================================================
// ConstantCoefficientPDE Implementation
// ============================================================================

ConstantCoefficientPDE::ConstantCoefficientPDE(double a, double b, double c, double d, std::function<double(double)> initialCond)
    : _a(a), _b(b), _c(c), _d(d), _u0(initialCond)
{
    if (_a < 0.0) {
        throw std::invalid_argument("PDE: Diffusion coefficient must be non-negative");
    }
    if (!_u0) {
        throw std::invalid_argument("PDE: Initial condition function is required");
    }
}

double ConstantCoefficientPDE::diffusion(double x, double t) const
{
    return _a;
}

double ConstantCoefficientPDE::convection(double x, double t) const
{
    return _b;
}

double ConstantCoefficientPDE::reaction(double x, double t) const
{
    return _c;
}

double ConstantCoefficientPDE::source(double x, double t) const
{
    return _d;
}

double ConstantCoefficientPDE::initialCondition(double x) const
{
    return _u0(x);
}

ConstantCoefficientPDE* ConstantCoefficientPDE::clone() const
{
    return new ConstantCoefficientPDE(*this);
}

// ============================================================================
// VariableCoefficientPDE Implementation
// ============================================================================
// Constructor
VariableCoefficientPDE::VariableCoefficientPDE(
    std::function<double(double, double)> diffusion_fn,
    std::function<double(double, double)> convection_fn,
    std::function<double(double, double)> reaction_fn,
    std::function<double(double, double)> source_fn,
    std::function<double(double)> initialCond)
    : _a_fn(diffusion_fn), _b_fn(convection_fn), 
      _c_fn(reaction_fn), _d_fn(source_fn), _u0(initialCond)
{
    if (!_a_fn || !_b_fn || !_c_fn || !_d_fn) {
        throw std::invalid_argument("PDE: All coefficient functions are required");
    }
    if (!_u0) {
        throw std::invalid_argument("PDE: Initial condition function is required");
    }
}

double VariableCoefficientPDE::diffusion(double x, double t) const
{
    return _a_fn(x, t);
}

double VariableCoefficientPDE::convection(double x, double t) const
{
    return _b_fn(x, t);
}

double VariableCoefficientPDE::reaction(double x, double t) const
{
    return _c_fn(x, t);
}

double VariableCoefficientPDE::source(double x, double t) const
{
    return _d_fn(x, t);
}

double VariableCoefficientPDE::initialCondition(double x) const
{
    return _u0(x);
}

VariableCoefficientPDE* VariableCoefficientPDE::clone() const
{
    return new VariableCoefficientPDE(*this);
}

