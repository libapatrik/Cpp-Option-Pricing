#ifndef CPPFM_PDEGREEKSCALCULATOR_H
#define CPPFM_PDEGREEKSCALCULATOR_H

#include "Solver.h"
#include "PDE.h"
#include "Grid.h"
#include "BoundaryConditions.h"
#include <functional>


/**
 * @struct GreekResults
 * @brief Container for all computed Greeks
 */
struct GreekResults
{
    double value;   // Option value
    double delta;   // ∂V/∂S
    double gamma;   // ∂²V/∂S²
    double vega;    // ∂V/∂σ
    double theta;   // -∂V/∂T (time decay)
    double rho;     // ∂V/∂r
    double vanna;   // ∂²V/∂S∂σ
    double volga;   // ∂²V/∂σ²
};

class PDEGreeksCalculator // base class for all PDEGreeks calculators
{
public:
    /**
     * Type for PDE factory function
     * Takes (sigma, r) and returns a PDE
     */
    using PDEFactory = std::function<std::unique_ptr<PDE>(double sigma, double r)>;
    
    /**
     * Constructor
     * @param pde Base PDE (with base parameters)
     * @param grid Spatial and temporal discretization
     * @param bc Boundary conditions
     * @param sigma_base Base volatility
     * @param r_base Base interest rate
     * @param scheme Numerical scheme
     */
    PDEGreeksCalculator(
        const PDE& pde,
        const Grid& grid,
        const BoundaryCondition& bc,
        double sigma_base,
        double r_base,
        ThetaMethodSolver::Scheme scheme = ThetaMethodSolver::Scheme::CrankNicolson
    );
    
    /**
     * Solve base PDE at given spot
     * @param spot Spot price for evaluation
     */
    void solve(double spot);
    
    /**
     * Compute all Greeks using bump-and-revalue
     * @param pdeFactory Function that creates PDE with (sigma, r) parameters
     * @return Struct containing all Greeks
     * 
     * Note: Bump sizes handled automatically by NumericalDerivatives utility
     */
    GreekResults computeAllGreeks(const PDEFactory& pdeFactory) const;
    
    // Individual Greek computations (for flexibility)
    double delta() const;
    double gamma() const;
    double vega(const PDEFactory& pdeFactory) const;
    double rho(const PDEFactory& pdeFactory) const;
    double vanna(const PDEFactory& pdeFactory) const;
    double volga(const PDEFactory& pdeFactory) const;
    double theta() const;
    
    // Accessors
    double value() const;
    bool isSolved() const { return _solved; } 

private:
    std::unique_ptr<Grid> _grid;
    std::unique_ptr<BoundaryCondition> _bc;
    ThetaMethodSolver::Scheme _scheme;
    
    double _sigma_base;
    double _r_base;
    double _spot;
    bool _solved;
    
    // Store base solver
    mutable std::unique_ptr<ThetaMethodSolver> _base_solver;
    
    // Helper: Solve PDE and return value at spot
    double solvePDEAt(const PDE& pde, double spot) const;
    
    // Helper: Solve PDE and return delta at spot
    double solvePDEDelta(const PDE& pde, double spot) const;

    // Helper: Sigma PDEs for vega/volga
    double solvePDEAt(const PDE& pde, double spot, double sigma) const;
    double solvePDEDelta(const PDE& pde, double spot, double sigma) const;
    std::unique_ptr<Grid> createScaledGrid(double sigma) const;
};

#endif //CPPFM_PDEGREEKSCALCULATOR_H

