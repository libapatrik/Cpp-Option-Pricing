#include <gtest/gtest.h>
#include "../PDEs/PDE.h"
#include "../PDEs/Grid.h"
#include "../PDEs/BoundaryConditions.h"
#include "../PDEs/Solver.h"
#include <cmath>
#include <iostream>
#include <iomanip>


/// TODO: FIX ALL OF THIS!

// ============================================================================
// Test 1: Heat Equation with Analytical Solution
// ============================================================================

/**
 * Heat equation: ∂u/∂t = α*∂2u/∂x2
 *
 * This maps to our general PDE with:
 *   a(x,t) = α (diffusion coefficient)
 *   b(x,t) = 0 (no convection)
 *   c(x,t) = 0 (no reaction)
 *   d(x,t) = 0 (no source)
 *
 * Analytical solution for u(x,0) = sin(πx) with u(0,t)=u(1,t)=0:
 *   u(x,t) = sin(πx)*exp(-α*π2*t)
 */
TEST(PDESolverTest, HeatEquation_AnalyticalSolution)
{
    std::cout << "\n=== Heat Equation Test ===" << std::endl;

    // Problem parameters
    const double alpha = 0.1;  // Thermal diffusivity
    const double L = 1.0;      // Domain length [0, L]
    const double T = 0.1;      // Final time

    // Initial condition: u(x,0) = sin(π*x)
    auto u0 = [](double x) { return std::sin(M_PI * x); };

    // Analytical solution: u(x,t) = sin(π*x)*exp(-α*π2*t)
    auto u_exact = [alpha](double x, double t) {
        return std::sin(M_PI * x) * std::exp(-alpha * M_PI * M_PI * t);
    };

    // Define PDE: ∂u/∂t = α*∂2u/∂x2
    ConstantCoefficientPDE pde(alpha, 0.0, 0.0, 0.0, u0);

    // Boundary conditions: u(0,t) = u(1,t) = 0
    DirichletBC bc(0.0, 0.0);

    // Grid: For CFL stability with explicit method, need μ = α*dt/dx² ≤ 0.5
    // With dx = 0.01, need dt ≤ 0.0005, so N_t ≥ 200
    Grid grid(0.0, L, 101, T, 250);  // Increased time steps for stability

    std::cout << "Grid: " << grid.numSpotPoints() << " spatial points, "
              << grid.numTimeSteps() << " time steps" << std::endl;
    std::cout << "dt = " << grid.dt() << ", dx ≈ " << grid.spot(1) - grid.spot(0) << std::endl;

    // Test Explicit Solver
    std::cout << "\n--- Explicit Solver ---" << std::endl;
    ThetaMethodSolver explicitSolver(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);
    auto solution_explicit = explicitSolver.solve();

    // Compare with analytical solution at t=T
    double max_error_explicit = 0.0;
    for (size_t i = 0; i < grid.numSpotPoints(); ++i) {
        double x = grid.spot(i);
        double u_numerical = solution_explicit[i];
        double u_analytical = u_exact(x, T);
        double error = std::abs(u_numerical - u_analytical);
        max_error_explicit = std::max(max_error_explicit, error);
    }

    std::cout << "Max error (Explicit): " << max_error_explicit << std::endl;

    // Test Implicit Solver
    std::cout << "\n--- Implicit Solver ---" << std::endl;
    ThetaMethodSolver implicitSolver(pde, grid, bc, ThetaMethodSolver::Scheme::Implicit);
    auto solution_implicit = implicitSolver.solve();

    double max_error_implicit = 0.0;
    for (size_t i = 0; i < grid.numSpotPoints(); ++i) {
        double x = grid.spot(i);
        double u_numerical = solution_implicit[i];
        double u_analytical = u_exact(x, T);
        double error = std::abs(u_numerical - u_analytical);
        max_error_implicit = std::max(max_error_implicit, error);
    }

    std::cout << "Max error (Implicit): " << max_error_implicit << std::endl;

    // Test Crank-Nicolson Solver
    std::cout << "\n--- Crank-Nicolson Solver ---" << std::endl;
    ThetaMethodSolver cnSolver(pde, grid, bc, ThetaMethodSolver::Scheme::CrankNicolson);
    auto solution_cn = cnSolver.solve();

    double max_error_cn = 0.0;
    for (size_t i = 0; i < grid.numSpotPoints(); ++i) {
        double x = grid.spot(i);
        double u_numerical = solution_cn[i];
        double u_analytical = u_exact(x, T);
        double error = std::abs(u_numerical - u_analytical);
        max_error_cn = std::max(max_error_cn, error);
    }

    std::cout << "Max error (Crank-Nicolson): " << max_error_cn << std::endl;

    // Print sample values
    std::cout << "\nSample comparison at x = 0.5:" << std::endl;
    std::cout << "Analytical:      " << u_exact(0.5, T) << std::endl;
    std::cout << "Explicit:        " << explicitSolver.valueAt(0.5) << std::endl;
    std::cout << "Implicit:        " << implicitSolver.valueAt(0.5) << std::endl;
    std::cout << "Crank-Nicolson:  " << cnSolver.valueAt(0.5) << std::endl;

    // Assertions
    EXPECT_LT(max_error_explicit, 0.01) << "Explicit solver error too large";
    EXPECT_LT(max_error_implicit, 0.01) << "Implicit solver error too large";
    EXPECT_LT(max_error_cn, 0.001) << "Crank-Nicolson solver error too large";

    // Crank-Nicolson should be most accurate (second order in time)
    EXPECT_LT(max_error_cn, max_error_explicit);
    EXPECT_LT(max_error_cn, max_error_implicit);
}

// ============================================================================
// Test 2: Black-Scholes European Call Option
// ============================================================================

/**
 * Black-Scholes PDE: ∂V/∂t + 1/2*σ^2 * S2*∂2V/∂S2 + rS*∂V/∂S - rV = 0
 *
 * Rewrite as: ∂V/∂t = 1/2*σ^2S2*∂2V/∂S2 + rS*∂V/∂S - rV
 *
 * Maps to general PDE:
 *   a(S,t) = 1/2*σ^2*S^2
 *   b(S,t) = r*S
 *   c(S,t) = -r
 *   d(S,t) = 0
 *
 * Initial condition (payoff at t=0): V(S,0) = max(S - K, 0)
 * Boundary conditions:
 *   V(0,t) = 0
 *   V(S_max,t) = S_max - K*exp(-r*t)
 */
TEST(PDESolverTest, BlackScholes_EuropeanCall)
{
    std::cout << "\n=== Black-Scholes European Call Test ===" << std::endl;

    // Market parameters
    const double S0 = 100.0;    // Current spot price
    const double K = 100.0;     // Strike
    const double r = 0.05;      // Risk-free rate
    const double sigma = 0.2;   // Volatility
    const double T = 1.0;       // Maturity

    const double S_max = 3.0 * K;  // Upper boundary for S

    std::cout << "Parameters: S0=" << S0 << ", K=" << K << ", r=" << r
              << ", σ=" << sigma << ", T=" << T << std::endl;

    // Define Black-Scholes PDE coefficients
    auto a_fn = [sigma](double S, double t) { return 0.5 * sigma * sigma * S * S; };
    auto b_fn = [r](double S, double t) { return r * S; };
    auto c_fn = [r](double S, double t) { return -r; };
    auto d_fn = [](double S, double t) { return 0.0; };

    // Initial condition: Call payoff at t=0
    auto payoff = [K](double S) { return std::max(S - K, 0.0); };

    VariableCoefficientPDE pde(a_fn, b_fn, c_fn, d_fn, payoff);

    // Boundary conditions for call option
    auto lower_bc = [](double t) { return 0.0; };  // V(0,t) = 0
    auto upper_bc = [K, r, S_max, T](double t) {
        // V(S_max,t) = S_max - K*exp(-r*(T-t))
        double tau = T - t;  // Time to maturity
        return S_max - K * std::exp(-r * tau);
    };

    DirichletBC bc(lower_bc, upper_bc);

    // Grid: For CFL stability with explicit method:
    // max_a = 0.5*σ²*S_max² = 0.5*0.04*90000 = 1800
    // dx = S_max/200 = 1.5
    // Need dt ≤ 0.5*dx²/max_a = 0.5*2.25/1800 ≈ 0.000625
    // So N_t ≥ T/0.000625 = 1600
    Grid grid(0.0, S_max, 201, T, 2000);  // Increased time steps for stability

    std::cout << "Grid: " << grid.numSpotPoints() << " spatial points, "
              << grid.numTimeSteps() << " time steps" << std::endl;

    // Solve with ExplicitSolver
    std::cout << "\n--- Solving with Explicit ---" << std::endl;
    ThetaMethodSolver solverE(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);

    // Solve with ImplicitSolver
    std::cout << "\n--- Solving with ImplicitSolver ---" << std::endl;
    ThetaMethodSolver solverI(pde, grid, bc);

    // Solve with Crank-Nicolson (most accurate)
    std::cout << "\n--- Solving with Crank-Nicolson ---" << std::endl;
    ThetaMethodSolver solverC(pde, grid, bc);
    auto solutionE = solverE.solve();
    auto solutionI = solverI.solve();
    auto solutionC = solverC.solve();


    // Get price at S0
    double pde_priceE = solverE.valueAt(S0);
    double pde_priceI = solverI.valueAt(S0);
    double pde_priceC = solverC.valueAt(S0);

    // Black-Scholes analytical formula for comparison
    auto N = [](double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    };

    double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double bs_price = S0 * N(d1) - K * std::exp(-r * T) * N(d2);

    std::cout << "\nResults at S = " << S0 << ":" << std::endl;
    std::cout << "PDE Price Explicit:            " << std::setprecision(6) << pde_priceE << std::endl;
    std::cout << "PDE Price Implicit:            " << std::setprecision(6) << pde_priceI << std::endl;
    std::cout << "PDE Price Crank-Nicolson:      " << std::setprecision(6) << pde_priceC << std::endl;
    std::cout << "Black-Scholes Price:           " << std::setprecision(6) << bs_price << std::endl;

    std::cout << "\nPrice differences:" << std::endl;
    std::cout << "Difference Explicit:                 " << std::abs(pde_priceE - bs_price) << std::endl;
    std::cout << "Relative Error Explicit:             " << 100.0 * std::abs(pde_priceE - bs_price) / bs_price << "%" << std::endl;
    std::cout << "Difference Implicit:                 " << std::abs(pde_priceI - bs_price) << std::endl;
    std::cout << "Relative Error Implicit:             " << 100.0 * std::abs(pde_priceI - bs_price) / bs_price << "%" << std::endl;
    std::cout << "Difference Crank-Nicolson:           " << std::abs(pde_priceC - bs_price) << std::endl;
    std::cout << "Relative Error Crank-Nicolson:       " << 100.0 * std::abs(pde_priceC - bs_price) / bs_price << "%" << std::endl;

    // Delta and Gamma from PDE solution
    double pde_deltaE = solverE.derivativeAt(S0);
    double pde_gammaE = solverE.secondDerivativeAt(S0);
    double pde_deltaI = solverI.derivativeAt(S0);
    double pde_gammaI = solverI.secondDerivativeAt(S0);
    double pde_deltaC = solverC.derivativeAt(S0);
    double pde_gammaC = solverC.secondDerivativeAt(S0);

    // Analytical Greeks
    double bs_delta = N(d1);
    double bs_gamma = std::exp(-0.5 * d1 * d1) / (S0 * sigma * std::sqrt(2.0 * M_PI * T));

    std::cout << "\nGreeks at S = " << S0 << ":" << std::endl;
    std::cout << "Delta (PDE) Explicit:        " << pde_deltaE << std::endl;
    std::cout << "Delta (PDE) Implicit:        " << pde_deltaI << std::endl;
    std::cout << "Delta (PDE) Crank-Nicolson:  " << pde_deltaC << std::endl;
    std::cout << "Delta (BS):                  " << bs_delta << std::endl;
    std::cout << "Gamma (PDE) Explicit:        " << pde_gammaE << std::endl;
    std::cout << "Gamma (PDE) Implicit:        " << pde_gammaI << std::endl;
    std::cout << "Gamma (PDE) Crank-Nicolson:  " << pde_gammaC << std::endl;
    std::cout << "Gamma (BS):                  " << bs_gamma << std::endl;

    // Assertions
    double rel_errorE = std::abs(pde_priceE - bs_price) / bs_price;
    EXPECT_LT(rel_errorE, 0.005) << "PDE Explicit price differs from Black-Scholes by more than 0.5%";

    double delta_errorE = std::abs(pde_deltaE  - bs_delta);
    EXPECT_LT(delta_errorE, 0.01) << "Delta Explicit error too large";

    double rel_errorI = std::abs(pde_priceI - bs_price) / bs_price;
    EXPECT_LT(rel_errorI, 0.005) << "PDE Explicit price differs from Black-Scholes by more than 0.5%";

    double delta_errorI = std::abs(pde_deltaI  - bs_delta);
    EXPECT_LT(delta_errorI, 0.01) << "Delta Explicit error too large";

    double rel_errorC = std::abs(pde_priceC - bs_price) / bs_price;
    EXPECT_LT(rel_errorC, 0.005) << "PDE Explicit price differs from Black-Scholes by more than 0.5%";

    double delta_errorC = std::abs(pde_deltaC  - bs_delta);
    EXPECT_LT(delta_errorC, 0.01) << "Delta Explicit error too large";
}

// ============================================================================
// Test 3: Convergence Test
// ============================================================================
TEST(PDESolverTest, Convergence_HeatEquation)
{
    std::cout << "\n=== Convergence Test ===" << std::endl;

    const double alpha = 0.1;
    const double L = 1.0;
    const double T = 0.1;

    auto u0 = [](double x) { return std::sin(M_PI * x); };
    auto u_exact = [alpha, T](double x) {
        return std::sin(M_PI * x) * std::exp(-alpha * M_PI * M_PI * T);
    };

    ConstantCoefficientPDE pde(alpha, 0.0, 0.0, 0.0, u0);
    DirichletBC bc(0.0, 0.0);

    std::cout << "\nRefining grid to test convergence:" << std::endl;
    std::cout << std::setw(10) << "N_x"
              << std::setw(15) << "Error (CN)"
              << std::setw(15) << "Ratio" << std::endl;

    double prev_error = 0.0;

    for (int refinement = 0; refinement < 4; ++refinement) {
        size_t N_x = 50 * (1 << refinement) + 1;  // 51, 101, 201, 401
        size_t N_t = 50 * (1 << refinement);      // 50, 100, 200, 400

        Grid grid(0.0, L, N_x, T, N_t);
        ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::CrankNicolson);
        solver.solve();

        // Compute L2 error
        double error_sum = 0.0;
        for (size_t i = 0; i < N_x; ++i) {
            double x = grid.spot(i);
            double u_num = solver.solution()[i];
            double u_ana = u_exact(x);
            error_sum += (u_num - u_ana) * (u_num - u_ana);
        }
        double error = std::sqrt(error_sum / N_x);

        double ratio = (prev_error > 0) ? prev_error / error : 0.0;

        std::cout << std::setw(10) << N_x
                  << std::setw(15) << std::scientific << error
                  << std::setw(15) << std::fixed << ratio << std::endl;

        prev_error = error;

        if (refinement > 0) {
            // Crank-Nicolson is second order, so error should decrease by ~4 when doubling resolution
            EXPECT_GT(ratio, 3.0) << "Convergence rate too slow";
        }
    }
}

// ============================================================================
// Test 4: Diffusion with Source Term
// ============================================================================

TEST(PDESolverTest, DiffusionWithSource)
{
    std::cout << "\n=== Diffusion with Source Term ===" << std::endl;

    // ∂u/∂t = α*∂2u/∂x2 + f(x,t)
    // With constant source: f(x,t) = 1

    const double alpha = 0.1;
    const double L = 1.0;
    const double T = 0.5;

    auto u0 = [](double x) { return 0.0; };  // Start with u = 0

    // PDE: ∂u/∂t = α*∂2u/∂x2 + 1
    ConstantCoefficientPDE pde(alpha, 0.0, 0.0, 1.0, u0);
    DirichletBC bc(0.0, 0.0);  // Keep boundaries at 0

    Grid grid(0.0, L, 101, T, 100);

    std::cout << "Solving diffusion equation with constant source term" << std::endl;

    ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::CrankNicolson);
    auto solution = solver.solve();

    // With source term, solution should be positive in interior
    double u_mid = solver.valueAt(0.5);

    std::cout << "Solution at x=0.5, t=" << T << ": " << u_mid << std::endl;

    EXPECT_GT(u_mid, 0.0) << "Solution should be positive with positive source";

    // Print profile
    std::cout << "\nSolution profile at t=" << T << ":" << std::endl;
    for (double x = 0.0; x <= L; x += 0.2) {
        std::cout << "  x=" << x << ": u=" << solver.valueAt(x) << std::endl;
    }
}

// ============================================================================
// Test 5: Grid Functionality
// ============================================================================

TEST(PDESolverTest, Grid_Uniform)
{
    std::cout << "\n=== Grid Test (Uniform) ===" << std::endl;

    Grid grid(0.0, 10.0, 11, 1.0, 10);

    EXPECT_EQ(grid.numSpotPoints(), 11);
    EXPECT_EQ(grid.numTimeSteps(), 10);
    EXPECT_DOUBLE_EQ(grid.spotMin(), 0.0);
    EXPECT_DOUBLE_EQ(grid.spotMax(), 10.0);
    EXPECT_DOUBLE_EQ(grid.timeMax(), 1.0);
    EXPECT_DOUBLE_EQ(grid.dt(), 0.1);

    // Check uniform spacing
    for (size_t i = 0; i < grid.numSpotPoints() - 1; ++i) {
        double ds = grid.spot(i+1) - grid.spot(i);
        EXPECT_NEAR(ds, 1.0, 1e-10);
    }

    std::cout << "Uniform grid test passed" << std::endl;
}

TEST(PDESolverTest, Grid_LogSpaced)
{
    std::cout << "\n=== Grid Test (Log-Spaced) ===" << std::endl;

    Grid grid(1.0, 100.0, 21, 1.0, 10, Grid::SpacingType::LogSpaced);

    EXPECT_DOUBLE_EQ(grid.spot(0), 1.0);
    EXPECT_DOUBLE_EQ(grid.spot(grid.numSpotPoints()-1), 100.0);

    // Log spacing should have smaller steps near S_min
    double ds_low = grid.spot(1) - grid.spot(0);
    double ds_high = grid.spot(20) - grid.spot(19);

    std::cout << "Spacing near S_min: " << ds_low << std::endl;
    std::cout << "Spacing near S_max: " << ds_high << std::endl;

    EXPECT_LT(ds_low, ds_high) << "Log spacing should be tighter near S_min";

    std::cout << "Log-spaced grid test passed" << std::endl;
}

// int main(int argc, char** argv)
// {
//     std::cout << std::fixed << std::setprecision(6);
//     testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }

