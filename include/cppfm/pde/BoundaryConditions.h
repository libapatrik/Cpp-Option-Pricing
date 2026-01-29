#ifndef CPPFM_BOUNDARYCONDITIONS_H
#define CPPFM_BOUNDARYCONDITIONS_H

#include <functional>

/**
 * @class BoundaryCondition
 * @brief Abstract base class for PDE boundary conditions
 * 
 * Defines the behavior of the solution at the spatial boundaries:
 *   - Lower boundary: x = x_min
 *   - Upper boundary: x = x_max
 * 
 * Currently supported:
 *   - Dirichlet: u(x_boundary, t) = f(t)
 * 
 */
class BoundaryCondition
{
public:
    virtual ~BoundaryCondition() = default;
    
    virtual double lowerValue(double t) const = 0; // PVM
    virtual double upperValue(double t) const = 0; // PVM
    
    virtual BoundaryCondition* clone() const = 0;  // PVM
};

/**
 * @class DirichletBC
 * @brief Dirichlet boundary conditions: u(x_boundary, t) = f(t)
 * Specifies the value of u at the boundaries.
 */
class DirichletBC : public BoundaryCondition
{
public:
    /**
     * Constructor with time-dependent boundary functions
     * @param lowerFunc Value at x_min: u(x_min, t) = lowerFunc(t)
     * @param upperFunc Value at x_max: u(x_max, t) = upperFunc(t)
     */
    DirichletBC(std::function<double(double)> lowerFunc, std::function<double(double)> upperFunc);
    
    DirichletBC(double lowerValue, double upperValue);
    
    double lowerValue(double t) const override;
    double upperValue(double t) const override;
    
    DirichletBC* clone() const override;

private:
    std::function<double(double)> _lowerFunc;
    std::function<double(double)> _upperFunc;
};

#endif //CPPFM_BOUNDARYCONDITIONS_H
