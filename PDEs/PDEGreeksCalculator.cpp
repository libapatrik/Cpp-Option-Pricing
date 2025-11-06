#include "PDEGreeksCalculator.h"
#include "../Utils.h"  // For NumericalDerivatives
#include <stdexcept>
#include <cmath>

PDEGreeksCalculator::PDEGreeksCalculator(
    const PDE& pde,
    const Grid& grid,
    const BoundaryCondition& bc,
    double sigma_base,
    double r_base,
    ThetaMethodSolver::Scheme scheme)
        : _grid(grid.clone()),
        _bc(bc.clone()),
        _scheme(scheme),
        _sigma_base(sigma_base),
        _r_base(r_base),
        _spot(0.0),
        _solved(false)
{
    // Create and solve base case
    _base_solver = std::make_unique<ThetaMethodSolver>(pde, *_grid, *_bc, _scheme);
}

void PDEGreeksCalculator::solve(double spot)
{
    _spot = spot;
    _base_solver->solve();
    _solved = true;
}

double PDEGreeksCalculator::value() const
{
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator: Must call solve() first");
    }
    return _base_solver->valueAt(_spot);
}

// ============================================================================
// Direct Greeks (from spatial derivatives)
// ============================================================================
double PDEGreeksCalculator::delta() const
{
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::delta: Must call solve() first");
    }
    
    // Delta comes directly from first spatial derivative
    return _base_solver->derivativeAt(_spot);
}

double PDEGreeksCalculator::gamma() const
{
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::gamma: Must call solve() first");
    }
    
    // Gamma comes directly from second spatial derivative
    return _base_solver->secondDerivativeAt(_spot);
}

// ============================================================================
// Parameter Greeks
// ============================================================================
double PDEGreeksCalculator::vega(const PDEFactory& pdeFactory) const
{
    /**
     * VEGA VIA CENTRAL DIFFERENCE
     * ===========================
     * Vega = ∂V/∂σ = [V(σ+h) - V(σ-h)] / (2h)
     * Step size h is handled automatically by NumericalDerivatives utility
     * The spatial domain should scale as S_max = exp(σ√T) to capture the proper diffusion range.
     
     */
    
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::vega: Must call solve() first");
    }
    
    // Create function that solves PDE for given volatility
    auto valueFunction = [&](double sigma) {
        auto pde = pdeFactory(sigma, _r_base);
        return solvePDEAt(*pde, _spot, sigma); // use sigma-scaled grid
    };
    
    // Compute derivative using reusable utility (h = 1e-4 by default)
    return NumericalDerivatives::firstDerivative(valueFunction, _sigma_base);
}

double PDEGreeksCalculator::rho(const PDEFactory& pdeFactory) const
{
    /**
     * RHO VIA CENTRAL DIFFERENCE
     * ==========================
     * Rho = ∂V/∂r = [V(r+h) - V(r-h)] / (2h)
     * Step size h is handled automatically by NumericalDerivatives utility
     */
    
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::rho: Must call solve() first");
    }
    
    // Create function that solves PDE for given interest rate
    auto valueFunction = [&](double rate) {
        auto pde = pdeFactory(_sigma_base, rate);
        return solvePDEAt(*pde, _spot);
    };
    
    // Compute derivative using reusable utility (h = 1e-4 by default)
    return NumericalDerivatives::firstDerivative(valueFunction, _r_base);
}

double PDEGreeksCalculator::theta() const
{
    /**
     * THETA APPROXIMATION
     * ===================
     * Theta = -∂V/∂T = -[V(T) - V(T-dt)] / dt
     * 
     * We approximate using the PDE relationship:
     * ∂V/∂t = 0.5*σ²*S²*∂²V/∂S² + r*S*∂V/∂S - r*V
     * 
     * For Black-Scholes PDE, theta can be computed from
     * the other Greeks using the PDE itself!
     */
    
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::theta: Must call solve() first");
    }
    
    // Use Black-Scholes PDE relationship
    double S = _spot;
    double Delta = delta();
    double Gamma = gamma();
    double V = value();

    // From PDE: ∂V/∂t = -0.5*σ²*S²*Γ - r*S*Δ + r*V
    // Theta = -∂V/∂t (negative by convention)
    double theta_value = 0.5 * _sigma_base * _sigma_base * S * S * Gamma 
                        + _r_base * S * Delta 
                        - _r_base * V;
    
    return -theta_value;  // Negative sign for theta convention
}

// ============================================================================
// Cross-Derivatives 
// ============================================================================
double PDEGreeksCalculator::vanna(const PDEFactory& pdeFactory) const
{
    /**
     * VANNA VIA CENTRAL DIFFERENCE OF DELTA
     * =====================================
     * Vanna = ∂^2V/∂S∂σ = ∂Δ/∂σ = vega(σ) = [Δ(σ+h) - Δ(σ-h)] / (2h)
     */
    
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::vanna: Must call solve() first");
    }
    // Create function that computes delta for given volatility
    auto deltaFunction = [&](double sigma) {
        auto pde = pdeFactory(sigma, _r_base);
        return solvePDEDelta(*pde, _spot, sigma);
    };
    
    // Compute derivative of delta wrt sigma using reusable utility
    return NumericalDerivatives::firstDerivative(deltaFunction, _sigma_base);
}

double PDEGreeksCalculator::volga(const PDEFactory& pdeFactory) const
{
    /**
     * VOLGA VIA SECOND-ORDER FINITE DIFFERENCE
     * ========================================
     * Volga = ∂^2V/∂σ^2 = ∂ν/∂σ = [V(σ+h) - 2*V(σ) + V(σ-h)] / h²   - solve 3 PDEs
     * Volga = [vega(σ+h) - vega(σ-h)] / (2h) - nested function worse
     * Bump size h is handled automatically by NumericalDerivatives in Utils.cpp
     * TODO: We need to update the grid
     *   - grid is setup for _sigma_base = 0.2, but then we solve PDEs at σ = 0.2 +/- h 
     * so, we need new S_max = K + exp(n * σ * sqrt(T))
     */

    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::volga: Must call solve() first");
    }
    
    // Create function that solved PDE for a given sigma with scaled grid
    auto valueFunction = [&](double sigma) {
        auto pde = pdeFactory(sigma, _r_base);
        return solvePDEAt(*pde, _spot, sigma);
    };
    
    // Compute second derivative using reusable utility (h = 1e-4 by default)
    /// TODO: Volga blows up, why? and how to fix this. [see test_pde_greeks.cpp]
    // return NumericalDerivatives::secondDerivative(valueFunction, _sigma_base);
    return NumericalDerivatives::secondDerivative(valueFunction, _sigma_base);
}

// ============================================================================
// Compute All Greeks at Once
// ============================================================================

GreekResults PDEGreeksCalculator::computeAllGreeks(const PDEFactory& pdeFactory) const
{
    if (!_solved) {
        throw std::runtime_error("PDEGreeksCalculator::computeAllGreeks: Must call solve() first");
    }
    
    GreekResults greeks;
    
    // Value and direct spatial derivatives (free!)
    greeks.value = value();
    greeks.delta = delta();
    greeks.gamma = gamma();
    // Theta from PDE relationship (efficient!)
    greeks.theta = theta();
    // Parameter sensitivities (require additional solves)
    greeks.vega = vega(pdeFactory);
    greeks.rho = rho(pdeFactory);
    // Cross-derivatives (expensive!)
    greeks.vanna = vanna(pdeFactory);
    greeks.volga = volga(pdeFactory);
    
    return greeks;
}

// ============================================================================
// Helper Methods
// ============================================================================
// 2-parameter version for rho
double PDEGreeksCalculator::solvePDEAt(const PDE& pde, double spot) const
{
    // Solve PDE and return value at spot
    ThetaMethodSolver solver(pde, *_grid, *_bc, _scheme);
    solver.solve();
    return solver.valueAt(spot);
}
// 3-parameter version for vega/vanna/volga
double PDEGreeksCalculator::solvePDEAt(const PDE& pde, double spot, double sigma) const
{
    // Solve PDE and return value at spot
    auto scaled_grid = createScaledGrid(sigma);
    ThetaMethodSolver solver(pde, *scaled_grid, *_bc, _scheme);
    solver.solve();
    return solver.valueAt(spot);
}

// 2-parameter version: uses base grid
double PDEGreeksCalculator::solvePDEDelta(const PDE& pde, double spot) const
{
    ThetaMethodSolver solver(pde, *_grid, *_bc, _scheme);
    solver.solve();
    return solver.derivativeAt(spot);
}

// 3-parameter version: uses scaled grid
double PDEGreeksCalculator::solvePDEDelta(const PDE& pde, double spot, double sigma) const
{
    // Solve PDE and return delta at spot
    auto scaled_grid = createScaledGrid(sigma);
    ThetaMethodSolver solver(pde, *scaled_grid, *_bc, _scheme);
    solver.solve();
    return solver.derivativeAt(spot);
}

std::unique_ptr<Grid> PDEGreeksCalculator::createScaledGrid(double sigma) const
{
    /**
     * But, scale the grid in S/K?
     * What about Neuman BCs? 
     */
    double sigma_ratio = sigma / _sigma_base;
    auto bounds = _grid->getBounds();
    double S_min_base = bounds.first.first;
    double S_max_base = bounds.first.second;
    double S_max_scaled = S_min_base + (S_max_base - S_min_base) * sigma_ratio;
    
    return std::make_unique<Grid>(
        S_min_base,
        S_max_scaled,
        _grid->numSpotPoints(),
        _grid->timeMax(),
        _grid->numTimeSteps(),
        _grid->spacingType()
    );
}