#ifndef CPPFM_SOLVER_H
#define CPPFM_SOLVER_H

#include <cppfm/pde/PDE.h>
#include <cppfm/pde/Grid.h>
#include <cppfm/pde/BoundaryConditions.h>
#include <cppfm/utils/InterpolationSchemes.h>
#include <vector>
#include <memory>

/**
 * =============================================================================
 * SOLVER: θ-METHOD IMPLEMENTATION
 * =============================================================================
 * The θ-method discretizes the PDE:
 *   ∂u/∂t = L(u)  where L is the spatial operator, u unknown function, t time (so time-dependent PDE)
 * 
 * Discretise time domain into steps of size Δt: t_n = n*Δt, u_n current soltion at time t_n, u_{n+1} next solution at time t_{n+1}
 *   (u^{n+1} - u^n)/Δt = θ·L(u^{n+1}) + (1-θ)·L(u^n)
 * 
 * Rearranging for unknown u^{n+1} LHS and knowns u^n RHS:
 *   (I - θΔt·L)u^{n+1} = (I + (1-θ)Δt·L)u^n
 *
 */

/**
 * @class Solver
 * @brief Base class for θ-method PDE solver with cubic spline interpolation
 * 
 */
class Solver
{
public:
    /**
     * Constructor
     * @param pde The PDE to solve (defines coefficients and initial condition)
     * @param grid Spatial and temporal discretization
     * @param bc Boundary conditions (Dirichlet)
     */
    Solver(const PDE& pde, const Grid& grid, const BoundaryCondition& bc);
    
    virtual ~Solver();
    
    // Prevent copying
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;
    
    /**
     * Solve the PDE from t=0 to t=T
     * @return Solution at final time T for all spatial points
     */
    virtual std::vector<double> solve() = 0; // pure virtual method 
    
    // Accessors
    /**
     *  Grid& grid() const { return *_grid; } - is bad!
     *  use const twice (1st) "returned grid cannot be modified"
     *  (2nd) "this function does not modify the object"
     */

    const Grid& grid() const { return *_grid; }
    const PDE& pde() const { return *_pde; }
    const BoundaryCondition& boundaryCondition() const { return *_bc; }
    const std::vector<double>& solution() const { return _solution; }
    
    // ========================================================================
    // Greeks and Value Queries (using Cubic Spline Interpolation)
    // ========================================================================
    
    double valueAt(double x) const; // value at spot price
    double derivativeAt(double x) const; // delta at spot price
    double secondDerivativeAt(double x) const; // gamma at spot price

protected:
    // Core components
    std::unique_ptr<PDE> _pde;              // PDE to solve
    std::unique_ptr<Grid> _grid;            // Spatial and temporal discretization
    std::unique_ptr<BoundaryCondition> _bc; // Boundary conditions
    std::vector<double> _solution;          // Solution at final time
    
    // Cubic spline for smooth interpolation ("lazy" construction = only built when needed)
    mutable std::unique_ptr<CubicSplineInterpolation> _interpolator;
    

    std::vector<double> computeInitialCondition() const;
    void applyBoundaryConditions(std::vector<double>& u, double t) const;
    
    // Interpolator 
    const CubicSplineInterpolation& getInterpolator() const;
    
    size_t findGridIndex(double x) const;
};

/**
 * @class ThetaMethodSolver
 * @brief Unified θ-method finite difference solver
 * 
 * NOTES:
 * =====================
 * This class implements the general θ-method discretization.
 * 
 * For θ = 0 (Explicit):
 *   - Directly evaluate u^{n+1} = u^n + Δt·L(u^n)
 *   - No linear system solve required
 *   - WARNING: Conditionally stable (CFL condition must be checked)
 * 
 * For θ > 0 (Implicit/Semi-implicit):
 *   - Solves (I - θΔt·L)u^{n+1} = RHS
 *   - Results in tridiagonal system (solved via Thomas algorithm)
 *   - RHS includes explicit part: (I + (1-θ)Δt·L)u^n
 * 
 * 3 cases:
 * ==========================
 * - For smooth solutions: θ = 0.5 (Crank-Nicolson)
 * - For maximum stability: θ = 1.0 (Implicit Euler)
 */
class ThetaMethodSolver : public Solver
{
public:
    /**
     * Enumeration for commonly used schemes
     */
    enum class Scheme {
        Explicit,         // θ = 0.0 (conditionally stable)
        Implicit,         // θ = 1.0 (unconditionally stable, 1st order)
        CrankNicolson     // θ = 0.5 (unconditionally stable, 2nd order)
    };
    
    /**
     * Constructor
     * @param pde The PDE to solve
     * @param grid Spatial and temporal discretization
     * @param bc Boundary conditions
     * @param scheme One of the three predefined schemes (default: CrankNicolson)
     */
    ThetaMethodSolver(const PDE& pde, const Grid& grid, 
                      const BoundaryCondition& bc, 
                      Scheme scheme = Scheme::CrankNicolson);
    
    std::vector<double> solve() override;     // solve the PDE using θ-method
    double theta() const { return _theta; }   // returns the θ param we used {0, 0.5, 1}
    Scheme scheme() const { return _scheme; } // returns the scheme type (for reference)

private:
    double _theta;      // θ parameter
    Scheme _scheme;     // Scheme type (for reference)
    
    void checkStability() const; // check CFL stability condition for explicit scheme
    
    /**
     * Take one θ-method time step (unified implementation)
     * 
     * ALGORITHM:
     * ==========
     * Solves the linear system: (I - θΔt·L)u^{n+1} = (I + (1-θ)Δt·L)u^n
     * 
     * This single method handles ALL schemes:
     * - θ = 0.0 (Explicit):      LHS becomes identity, full explicit RHS
     * - θ = 0.5 (Crank-Nicolson): Symmetric implicit/explicit weighting  
     * - θ = 1.0 (Implicit):       Full implicit LHS, simplified RHS
     * 
     * The method builds a tridiagonal system and solves it using the
     * Thomas algorithm (O(N) complexity). For θ=0, the system degenerates
     * to an identity matrix which is solved trivially.
     * 
     * @param u_old Solution at time t^n
     * @param t_old Time t^n
     * @param t_new Time t^{n+1} = t^n + Δt
     * @return Solution at time t^{n+1}
     */
    std::vector<double> timeStep(const std::vector<double>& u_old, 
                                  double t_old, double t_new);
};

#endif //CPPFM_SOLVER_H

