#include "BoundaryConditions.h"
#include <stdexcept>

// ============================================================================
// DirichletBC Implementation
// ============================================================================

/**
 *  TODO: Add Cauchy boundary condition class
 */

DirichletBC::DirichletBC(std::function<double(double)> lowerFunc,
                         std::function<double(double)> upperFunc)
    : _lowerFunc(lowerFunc), _upperFunc(upperFunc)
{
    if (!_lowerFunc || !_upperFunc) {
        throw std::invalid_argument("DirichletBC: Boundary functions cannot be null");
    }
}

DirichletBC::DirichletBC(double lowerValue, double upperValue)
    : _lowerFunc([lowerValue](double t) { return lowerValue; }),
      _upperFunc([upperValue](double t) { return upperValue; })
{
}

double DirichletBC::lowerValue(double t) const
{
    return _lowerFunc(t);
}

double DirichletBC::upperValue(double t) const
{
    return _upperFunc(t);
}

DirichletBC* DirichletBC::clone() const
{
    return new DirichletBC(*this);
}
