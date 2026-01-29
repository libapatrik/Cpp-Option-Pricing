#ifndef CPPFM_PDE_H
#define CPPFM_PDE_H

#include <vector>
#include <functional>

/**
 * @class PDE
 * 
 * We want to solve the following PDE:
 *   ∂u/∂t = a(x,t)*∂^2u/∂x^2 + b(x,t)*∂u/∂x + c(x,t)*u + d(x,t)
 * 
 * With constant coefficients and variable coefficients.
 *   - a(x,t): diffusion coefficient (must be ≥ 0 for parabolic PDEs)
 *   - b(x,t): convection/drift coefficient
 *   - c(x,t): reaction coefficient
 *   - d(x,t): source term
 * 
 * The PDE is solved with:
 *   - Initial condition: u(x, 0) = u_0(x)
 *   - Boundary conditions at x_min and x_max (see BoundaryConditions.h)
 */
class PDE
{
public:
    virtual ~PDE() = default;
    
    // PDE Coefficients at point (x,t)
    virtual double diffusion(double x, double t) const = 0;    // a(x,t)
    virtual double convection(double x, double t) const = 0;   // b(x,t)
    virtual double reaction(double x, double t) const = 0;     // c(x,t)
    virtual double source(double x, double t) const = 0;       // d(x,t)
    
    // Initial condition: u(x, t=0)
    virtual double initialCondition(double x) const = 0;
    
    // Clone pattern for polymorphic copying
    virtual PDE* clone() const = 0; // this returns a raw pointer
    // virtual std::unique_ptr<PDE> clone() const = 0; 
};

/**
 * @class ConstantCoefficientPDE
 * @brief PDE with constant coefficients (independent of x and t)
 * 
 * Simplified case: ∂u/∂t = a*∂2u/∂x2 + b*∂u/∂x + c*u + d
 * 
 * Example: Heat equation with ∂u/∂t = α*∂2u/∂x2 has a=α, b=c=d=0
 */
class ConstantCoefficientPDE : public PDE
{
public:
    /**
     * Constructor
     * @param a Diffusion coefficient (must be ≥ 0)
     * @param b Convection coefficient
     * @param c Reaction coefficient
     * @param d Source term
     * @param initialCond Initial condition function u(x,0)
     */
    ConstantCoefficientPDE(double a, double b, double c, double d,
                          std::function<double(double)> initialCond);
    
    double diffusion(double x, double t) const override;
    double convection(double x, double t) const override;
    double reaction(double x, double t) const override;
    double source(double x, double t) const override;
    double initialCondition(double x) const override;
    
    ConstantCoefficientPDE* clone() const override;

private:
    double _a, _b, _c, _d;
    std::function<double(double)> _u0;
};

/**
 * @class VariableCoefficientPDE
 * @brief PDE with space/time-dependent coefficients
 * 
 * General case: Coefficients can depend on both x and t
 * 
 * Example: Black-Scholes PDE where a(S,t) = 1/2*σ^2*S^2, b(S,t) = r*S, c(S,t) = -r, d(S,t) = 0
 */
class VariableCoefficientPDE : public PDE
{
public:
    /**
     * Constructor with coefficient functions
     * @param diffusion_fn a(x,t)
     * @param convection_fn b(x,t)
     * @param reaction_fn c(x,t)
     * @param source_fn d(x,t)
     * @param initialCond u(x,0)
     */
    VariableCoefficientPDE(
        std::function<double(double, double)> diffusion_fn, // returns double 
        std::function<double(double, double)> convection_fn,
        std::function<double(double, double)> reaction_fn,
        std::function<double(double, double)> source_fn,
        std::function<double(double)> initialCond);
    
    double diffusion(double x, double t) const override;
    double convection(double x, double t) const override;
    double reaction(double x, double t) const override;
    double source(double x, double t) const override;
    double initialCondition(double x) const override;
    
    VariableCoefficientPDE* clone() const override;

private:
    std::function<double(double, double)> _a_fn;
    std::function<double(double, double)> _b_fn;
    std::function<double(double, double)> _c_fn;
    std::function<double(double, double)> _d_fn;
    std::function<double(double)> _u0;
};

#endif //CPPFM_PDE_H

