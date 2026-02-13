#include <cppfm/pde/Solver.h>
#include <cppfm/utils/Utils.h>
#include <cppfm/utils/InterpolationSchemes.h>
#include <algorithm>
#include <stdexcept>
#include <cmath>

/** SOLVE: ∂u/∂t = a(x,t)·∂²u/∂x² + b(x,t)·∂u/∂x + c(x,t)·u + d(x,t)
 * NOTES: See Lecture CF_Fang 10 in Dropbox for reference
 */


Solver::Solver(const PDE& pde, const Grid& grid, const BoundaryCondition& bc)
    : _pde(pde.clone()),
      _grid(grid.clone()),
      _bc(bc.clone())          // each clone returns a unique pointer
{
    // Pre-allocate solution vector
    size_t N_x = _grid->numSpotPoints();
    _solution.resize(N_x, 0.0); // initialise solution vector with zeros
}

Solver::~Solver() = default;

std::vector<double> Solver::computeInitialCondition() const
{
    size_t N_x = _grid->numSpotPoints();  // _grid is unique pointer to Grid, so we use -> to access numSpotPoints()
    std::vector<double> u0(N_x);          // initialise u0 vector with zeros
    
    // Evaluate initial condition at each grid point
    for (size_t i = 0; i < N_x; ++i) {
        double x = _grid->spot(i);      
        u0[i] = _pde->initialCondition(x);
    }
    
    return u0;
}

void Solver::applyBoundaryConditions(std::vector<double>& u, double t) const
{
    // Dirichlet boundary conditions: directly set boundary values
    u[0] = _bc->lowerValue(t);
    u[u.size() - 1] = _bc->upperValue(t); // u.size() - 1 could cause issues if u is empty (size_t is unsigned, so -1 is very large positive number)
                                          // but we are fine here because u is pre-allocated vector
}

size_t Solver::findGridIndex(double x) const
{
    const std::vector<double>& spots = _grid->spots();
    
    if (x < spots.front() || x > spots.back()) { // spots.back() == spots[spots.size() - 1]
        throw std::out_of_range("Solver::findGridIndex: x outside grid bounds");
    }
    
    auto it = std::upper_bound(spots.begin(), spots.end(), x);
    
    if (it == spots.end()) {
        return spots.size() - 2;
    }
    if (it == spots.begin()) {
        return 0;
    }
    
    return std::distance(spots.begin(), it) - 1;
}

const CubicSplineInterpolation& Solver::getInterpolator() const
{
    if (!_interpolator) {
        // Lazy construction: build cubic spline only when first needed
        // This gives us smooth C^2 continuous interpolation for accurate Greeks
        const std::vector<double>& spots = _grid->spots();
        _interpolator = std::make_unique<CubicSplineInterpolation>(           // creates a unique pointer to a new CubicSplineInterpolation object
            spots, _solution, CubicSplineInterpolation::BoundaryType::Natural // passes the spots, solution, and boundary type to the constructor
        ); // creates a new CubicSplineInterpolation object and assigns it to the _interpolator unique pointer
    }      // if _interpolator is not null
    return *_interpolator; // return reference to owned object
}

// ============================================================================
// Value and Greek (with Cubic Spline)
// ============================================================================

double Solver::valueAt(double x) const
{
    if (_solution.empty()) { // check if solution is empty
        throw std::runtime_error("Solver::valueAt: No solution available - call solve() first.");
    }
    
    // Use cubic spline for C^2 continuity
    return getInterpolator().interpolate(x); 
    // getInterpolator() returns const CubicSplineInterpolation&
    // call interpolate() method on that reference to return the interpolated value at x
}

double Solver::derivativeAt(double x) const
{
    if (_solution.empty()) {
        throw std::runtime_error("Solver::derivativeAt: No solution available.");
    }
    
    // Cubic spline provides analytical derivative -> smooth Delta
    return getInterpolator().derivative(x);
}

double Solver::secondDerivativeAt(double x) const
{
    if (_solution.empty()) {
        throw std::runtime_error("Solver::secondDerivativeAt: No solution available.");
    }
    
    // Cubic spline provides analytical second derivative -> smooth Gamma
    return getInterpolator().secondDerivative(x);
}

// ============================================================================
// ThetaMethodSolver Implementation
// ============================================================================

ThetaMethodSolver::ThetaMethodSolver(const PDE& pde, const Grid& grid, const BoundaryCondition& bc, Scheme scheme)
    : Solver(pde, grid, bc), _scheme(scheme)
{
    // Map scheme enum to θ value
    switch (scheme) { // scheme is an enum value (Explicit, Implicit, CrankNicolson)
        case Scheme::Explicit:
            _theta = 0.0;
            break;
        case Scheme::Implicit:
            _theta = 1.0;
            break;
        case Scheme::CrankNicolson:
            _theta = 0.5;
            break;
    }
}

void ThetaMethodSolver::checkStability() const
{
    // Only explicit EULER method (θ = 0) needs stability check
    if (_theta > 0.0) return; // if θ > 0, return (no stability check needed)
    
    // Find maximum diffusion coefficient
    const std::vector<double>& spots = _grid->spots();
    size_t N_x = _grid->numSpotPoints();
    double max_a = 0.0;
    
    for (size_t i = 0; i < N_x; ++i) {
        double a = _pde->diffusion(spots[i], 0.0);
        max_a = std::max(max_a, a);
    }
    
    // Find minimum grid spacing
    double min_dx = spots[1] - spots[0];
    for (size_t i = 1; i < N_x - 1; ++i) {
        double dx = spots[i+1] - spots[i];
        min_dx = std::min(min_dx, dx);
    }
    
    // Check CFL condition: μ = a*Δt/Δx² ≤ 0.5
    double dt = _grid->dt();
    double mu = max_a * dt / (min_dx * min_dx);
    
    if (mu > 0.5) {
        throw std::runtime_error(
            "CFL stability condition violated for Explicit scheme!\n"
            "  μ = " + std::to_string(mu) + " > 0.5\n"
            "  Solution: Increase time steps or use Implicit/CrankNicolson scheme."
        );
    }
}

std::vector<double> ThetaMethodSolver::solve()
{
    // STEP 0: Check stability for explicit scheme
    checkStability();
    
    // STEP 1: Initialise with initial condition u(x, 0)
    std::vector<double> u = computeInitialCondition();   // by value! So we copy the vector
    applyBoundaryConditions(u, 0.0);                     // by reference! So we modify the original vector in place
    
    // STEP 2: Time marching loop
    // Advance solution from t=0 to t=T in N_t steps
    size_t N_t = _grid->numTimeSteps();
    
    for (size_t n = 0; n < N_t; ++n) {
        double t_old = _grid->time(n);
        double t_new = _grid->time(n + 1);
        
        // Take one time step using θ-method
        u = timeStep(u, t_old, t_new);
    }
    
    // STEP 3: Store final solution and invalidate interpolator cache
    _solution = u;
    _interpolator.reset();  
    // Deletes the unique pointer to the CubicSplineInterpolation object and sets it to nullptr
    // This forces the next call to getInterpolator() to rebuild the spline with the new solution
    
    return _solution;  // Returns copy of the solution vector
}

std::vector<double> ThetaMethodSolver::timeStep(const std::vector<double>& u_old, double t_old, double t_new)
{
    /**
     * UNIFIED θ-METHOD IMPLEMENTATION
     * ================================
     * Formula: (I - θΔt·L)u^{n+1} = (I + (1-θ)Δt·L)u^n
     * 
     * This single implementation handles ALL schemes:
     * - θ = 0.0: Explicit (Forward Euler)     - LHS becomes identity
     * - θ = 0.5: Crank-Nicolson               - Symmetric weighting
     * - θ = 1.0: Implicit (Backward Euler)    - RHS simplifies
     * 
     * LINEAR SYSTEM: A·u^{n+1} = b
     * - LHS matrix A = I - θΔt·L  (tridiagonal, identity when θ=0)
     * - RHS vector b = (I + (1-θ)Δt·L)u^n + Δt·d
     * 
     * STRATEGY:
     * 1. Build tridiagonal system using finite difference stencil
     * 2. Compute RHS with (1-θ) weighting on explicit part
     * 3. Solve using Thomas algorithm (O(N) complexity)
     * 4. Boundary conditions embedded in system
     */
    
    size_t N_x = _grid->numSpotPoints();
    double dt = _grid->dt();
    
    // Allocate tridiagonal system: A·u_new = rhs
    std::vector<double> lower(N_x, 0.0);  // Lower-diagonal
    std::vector<double> diag(N_x, 1.0);   // Main diagonal (starts as identity)
    std::vector<double> upper(N_x, 0.0);  // Upper-diagonal
    std::vector<double> rhs(N_x, 0.0);    // Right-hand side
    
    /**
     * COEFFICIENT EVALUATION TIME:
     * - θ = 0.0 (explicit): t_old
     * - θ = 0.5 (CN): t_mid = (t_old + t_new)/2
     * - θ = 1.0 (implicit): t_new
     * - General: θ-weighted interpolation
     * 
     * Q: Are we interpolating time here? Or what is going on?
     * - No. We compute at which single time point to use to evalute the PDE
     */
    double t_eval = (1.0 - _theta) * t_old + _theta * t_new;
    
    // BUILD SYSTEM FOR INTERIOR POINTS: i = 1, 2, ..., N_x-2
    for (size_t i = 1; i < N_x - 1; ++i) {     // loop through interior grid points
        double x = _grid->spot(i);             // x is the grid point index
        
        // Evaluate PDE coefficients at θ-weighted time
        // x varies, t_eval is constant for each grid point
        double a = _pde->diffusion(x, t_eval);
        double b = _pde->convection(x, t_eval);
        double c = _pde->reaction(x, t_eval);
        double d = _pde->source(x, t_eval);
        
        // Grid spacing for non-uniform grids
        double h_plus = _grid->spot(i+1) - _grid->spot(i);    // forward step
        double h_minus = _grid->spot(i) - _grid->spot(i-1);   // backward step
        double h_sum = h_plus + h_minus;
        
        /**
         * FINITE DIFFERENCE STENCIL FOR UNIFORM AND LOG-SPACED GRIDS:
         * =================================================
         * The spatial operator L discretized gives a 3-point stencil:
         *   L(u_i) = α·u[i-1] + β·u[i] + γ·u[i+1] + δ
         * 
         * Where:
         *   α = a/Δx² - b/(2Δx)      (coefficient of u[i-1])
         *   β = -2a/Δx² + c          (coefficient of u[i])
         *   γ = a/Δx² + b/(2Δx)      (coefficient of u[i+1])
         *   δ = d                    (source term)
         *
         * MORE DERIVATIONS FOR NON-UNIFORM GRIDS:
         *   L(u) = a·∂²u/∂x² + b·∂u/∂x + c·u
                 = a·[coeff_{i-1}·u_{i-1} + coeff_i·u_i + coeff_{i+1}·u_{i+1}]
                   + b·[coeff'_{i-1}·u_{i-1} + coeff'_i·u_i + coeff'_{i+1}·u_{i+1}]
                   + c·u_i
            
         *  Second derivative coefficients:
                Coefficient of u_{i-1}: 2a/(h_minus * (h_plus + h_minus))
                Coefficient of u_i: -2a/(h_plus * h_minus)
                Coefficient of u_{i+1}: 2a/(h_plus * (h_plus + h_minus))
         * First derivative coefficients (central difference, non-uniform):
                Coefficient of u_{i-1}: -b*h_plus / (h_minus * (h_plus + h_minus))
                Coefficient of u_i: b*(h_plus - h_minus) / (h_plus * h_minus)
                Coefficient of u_{i+1}: b*h_minus / (h_plus * (h_plus + h_minus))
         *  u_{i-1} 
         *  alpha = 2a/(h_minus * h_sum) - b*h_plus/(h_minus * h_sum)
                  = (2a - b*h_plus) / (h_minus * h_sum)

         *  u_i
            beta = -2a/(h_plus * h_minus) + b*(h_plus - h_minus)/(h_plus * h_minus) + c
                 = (-2a + b*(h_plus - h_minus))/(h_plus * h_minus) + c

         *  u_{i+1}
            gamma = 2a/(h_plus * h_sum) + b*h_minus/(h_plus * h_sum)
                  = (2a + b*h_minus) / (h_plus * h_sum)
         */
        double alpha = (2.0 * a - b * h_plus) / (h_minus * h_sum);
        double beta = (-2.0 * a + b * (h_plus - h_minus)) / (h_plus * h_minus) + c;
        double gamma = (2.0 * a + b * h_minus) / (h_plus * h_sum);
        
        /**
         * LEFT-HAND SIDE: (I - θΔt·L)u^{n+1}
         * ===================================
         * Tridiagonal matrix A:
         *   A[i,i-1] = -θΔt·α
         *   A[i,i]   = 1 - θΔt·β
         *   A[i,i+1] = -θΔt·γ
         * 
         * Note: For θ=0, this becomes identity matrix 
         */
        lower[i] = -_theta * dt * alpha;
        diag[i] = 1.0 - _theta * dt * beta;
        upper[i] = -_theta * dt * gamma;
        
        // RHS: u^n + (1-θ)Δt·(α u[i-1] + β u[i] + γ u[i+1]) + Δt·d
        // source d gets full Δt weight, not (1-θ) — otherwise backward Euler loses it entirely
        double explicit_contribution = (1.0 - _theta) * dt * (
            alpha * u_old[i-1] + beta * u_old[i] + gamma * u_old[i+1]
        );
        rhs[i] = u_old[i] + explicit_contribution + dt * d;
    }
    
    /**
     * BOUNDARY CONDITIONS (Dirichlet)
     * ================================
     * Set first and last equations to enforce boundary values
     */
    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = _bc->lowerValue(t_new);
    
    lower[N_x-1] = 0.0;
    diag[N_x-1] = 1.0;
    rhs[N_x-1] = _bc->upperValue(t_new);
    
    /**
     * SOLVE TRIDIAGONAL SYSTEM
     * ========================
     * Thomas algorithm handles all θ cases efficiently (O(N) complexity)
     * For θ=0, solves identity system (very fast)
     */
    return ThomasAlgorithm::solve(lower, diag, upper, rhs);
} // ThetaMethodSolver::timeStep()
