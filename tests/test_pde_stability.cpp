#include <gtest/gtest.h>
#include <cppfm/pde/PDE.h>
#include <cppfm/pde/Grid.h>
#include <cppfm/pde/BoundaryConditions.h>
#include <cppfm/pde/Solver.h>
#include <cmath>
#include <iostream>
#include <iomanip>

/**
 * Stability Analysis for Explicit Euler Method
 *
 * For the heat equation ∂u/∂t = α·∂²u/∂x², the explicit method is stable if:
 *
 *   μ = α·Δt/Δx² ≤ 0.5
 *
 * This is called the CFL (Courant-Friedrichs-Lewy) condition.
 */
TEST(PDEStabilityTest, ExplicitEuler_StabilityCondition)
{
    std::cout << "\n=== Explicit Euler Stability Analysis ===" << std::endl;
    std::cout << "\nStability condition: mu = α·Δt/Δx² ≤ 0.5\n" << std::endl;

    const double alpha = 0.1;
    const double L = 1.0;
    const double T = 0.1;

    auto u0 = [](double x) { return std::sin(M_PI * x); };
    auto u_exact = [alpha, T](double x) {
        return std::sin(M_PI * x) * std::exp(-alpha * M_PI * M_PI * T);
    };

    ConstantCoefficientPDE pde(alpha, 0.0, 0.0, 0.0, u0);
    DirichletBC bc(0.0, 0.0);

    std::cout << std::setw(10) << "N_x"
              << std::setw(10) << "N_t"
              << std::setw(12) << "Δx"
              << std::setw(12) << "Δt"
              << std::setw(12) << "μ"
              << std::setw(15) << "Max Error"
              << std::setw(12) << "Status" << std::endl;
    std::cout << std::string(83, '-') << std::endl;

    // Test 1: Unstable (μ = 1.0, too large!) - should throw exception
    {
        size_t N_x = 101;
        double dx = L / (N_x - 1);
        double dt_unstable = 1.0 * dx * dx / alpha;  // μ = 1.0
        size_t N_t = static_cast<size_t>(T / dt_unstable);
        double mu = alpha * dt_unstable / (dx * dx);

        Grid grid(0.0, L, N_x, T, N_t);
        ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);
        
        // Solver should throw exception for unstable CFL condition
        bool exception_thrown = false;
        try {
            solver.solve();
        } catch (const std::runtime_error& e) {
            exception_thrown = true;
            std::string msg(e.what());
            EXPECT_TRUE(msg.find("CFL") != std::string::npos) << "Exception should mention CFL condition";
        }

        std::cout << std::setw(10) << N_x
                  << std::setw(10) << N_t
                  << std::setw(12) << std::fixed << std::setprecision(6) << dx
                  << std::setw(12) << dt_unstable
                  << std::setw(12) << mu
                  << std::setw(15) << "N/A"
                  << std::setw(12) << "PREVENTED" << std::endl;

        EXPECT_TRUE(exception_thrown) << "Should throw exception for μ=1.0 > 0.5";
    }

    // Test 2: Marginally stable (μ = 0.49, just below limit)
    {
        size_t N_x = 101;
        double dx = L / (N_x - 1);
        double dt_marginal = 0.49 * dx * dx / alpha;  // μ = 0.49 (avoid floating point issues at boundary)
        size_t N_t = static_cast<size_t>(T / dt_marginal) + 1;  // +1 to ensure dt is not too large
        double mu = alpha * dt_marginal / (dx * dx);

        Grid grid(0.0, L, N_x, T, N_t);
        ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);
        solver.solve();

        double max_error = 0.0;
        for (size_t i = 0; i < N_x; ++i) {
            double x = grid.spot(i);
            double error = std::abs(solver.solution()[i] - u_exact(x));
            max_error = std::max(max_error, error);
        }

        std::cout << std::setw(10) << N_x
                  << std::setw(10) << N_t
                  << std::setw(12) << std::fixed << dx
                  << std::setw(12) << dt_marginal
                  << std::setw(12) << mu
                  << std::setw(15) << std::scientific << max_error
                  << std::setw(12) << "MARGINAL" << std::endl;

        EXPECT_LT(max_error, 0.01) << "Should be stable with μ=0.49";
    }

    // Test 3: Stable (μ = 0.25)
    {
        size_t N_x = 101;
        double dx = L / (N_x - 1);
        double dt_stable = 0.25 * dx * dx / alpha;  // μ = 0.25
        size_t N_t = static_cast<size_t>(T / dt_stable);
        double mu = alpha * dt_stable / (dx * dx);

        Grid grid(0.0, L, N_x, T, N_t);
        ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);
        solver.solve();

        double max_error = 0.0;
        for (size_t i = 0; i < N_x; ++i) {
            double x = grid.spot(i);
            double error = std::abs(solver.solution()[i] - u_exact(x));
            max_error = std::max(max_error, error);
        }

        std::cout << std::setw(10) << N_x
                  << std::setw(10) << N_t
                  << std::setw(12) << std::fixed << dx
                  << std::setw(12) << dt_stable
                  << std::setw(12) << mu
                  << std::setw(15) << std::scientific << max_error
                  << std::setw(12) << "STABLE" << std::endl;

        EXPECT_LT(max_error, 0.01) << "Should be stable with μ=0.25";
    }

    std::cout << "\n✓ Explicit method requires μ ≤ 0.5 for stability" << std::endl;
    std::cout << "✓ Implicit and Crank-Nicolson are unconditionally stable" << std::endl;
}

/**
 * Demonstration: Why Implicit methods are preferred despite being more expensive
 */
TEST(PDEStabilityTest, ImplicitVsExplicit_Efficiency)
{
    std::cout << "\n=== Implicit vs Explicit: Efficiency Comparison ===" << std::endl;

    const double alpha = 0.1;
    const double L = 1.0;
    const double T = 1.0;  // Longer time horizon

    auto u0 = [](double x) { return std::sin(M_PI * x); };

    ConstantCoefficientPDE pde(alpha, 0.0, 0.0, 0.0, u0);
    DirichletBC bc(0.0, 0.0);

    size_t N_x = 101;
    double dx = L / (N_x - 1);

    std::cout << "\nFixed spatial resolution: N_x = " << N_x << ", Δx = " << dx << std::endl;
    std::cout << "Time horizon: T = " << T << std::endl;

    // Explicit: Must use small time step
    double dt_explicit = 0.4 * dx * dx / alpha;  // Safe explicit time step
    size_t N_t_explicit = static_cast<size_t>(T / dt_explicit);

    // Implicit: Can use much larger time step
    double dt_implicit = 10.0 * dt_explicit;  // 10x larger!
    size_t N_t_implicit = static_cast<size_t>(T / dt_implicit);

    std::cout << "\nExplicit method:" << std::endl;
    std::cout << "  Δt = " << dt_explicit << " (limited by stability)" << std::endl;
    std::cout << "  N_t = " << N_t_explicit << " time steps required" << std::endl;

    std::cout << "\nImplicit method:" << std::endl;
    std::cout << "  Δt = " << dt_implicit << " (no stability limit)" << std::endl;
    std::cout << "  N_t = " << N_t_implicit << " time steps required" << std::endl;

    std::cout << "\nTime step ratio: " << static_cast<double>(N_t_explicit) / N_t_implicit << "x" << std::endl;
    std::cout << "\n✓ Implicit methods can use much larger time steps!" << std::endl;
    std::cout << "✓ This often makes them faster overall despite the tridiagonal solve" << std::endl;
}

/**
 * Black-Scholes Stability: Variable coefficients
 */
TEST(PDEStabilityTest, BlackScholes_VariableCoefficients)
{
    std::cout << "\n=== Black-Scholes with Variable Coefficients ===" << std::endl;

    const double S0 = 100.0;
    const double K = 100.0;
    const double r = 0.05;
    const double sigma = 0.2;
    const double T = 1.0;
    const double S_max = 3.0 * K;

    // Black-Scholes coefficients: a(S) = 0.5*σ²*S²
    // The stability condition depends on S: worst case at S=S_max
    auto a_fn = [sigma](double S, double t) { return 0.5 * sigma * sigma * S * S; };
    auto b_fn = [r](double S, double t) { return r * S; };
    auto c_fn = [r](double S, double t) { return -r; };
    auto d_fn = [](double S, double t) { return 0.0; };
    auto payoff = [K](double S) { return std::max(S - K, 0.0); };

    VariableCoefficientPDE pde(a_fn, b_fn, c_fn, d_fn, payoff);

    auto lower_bc = [](double t) { return 0.0; };
    auto upper_bc = [K, r, S_max, T](double t) {
        return S_max - K * std::exp(-r * (T - t));
    };
    DirichletBC bc(lower_bc, upper_bc);

    std::cout << "\nFor Black-Scholes, diffusion coefficient varies with S:" << std::endl;
    std::cout << "  a(S) = 0.5·σ²·S²" << std::endl;
    std::cout << "  Maximum at S_max = " << S_max << std::endl;
    std::cout << "  a(S_max) = " << a_fn(S_max, 0.0) << std::endl;

    size_t N_x = 201;
    double dS = S_max / (N_x - 1);

    double a_max = a_fn(S_max, 0.0);
    double dt_stable = 0.4 * dS * dS / a_max;

    std::cout << "\nFor stability with explicit method:" << std::endl;
    std::cout << "  ΔS = " << dS << std::endl;
    std::cout << "  Required: Δt ≤ " << 0.5 * dS * dS / a_max << std::endl;
    std::cout << "  Safe choice: Δt ≤ " << dt_stable << std::endl;

    size_t N_t_stable = static_cast<size_t>(T / dt_stable) + 1;
    std::cout << "  This requires N_t ≥ " << N_t_stable << " time steps!" << std::endl;

    // Test with safe time step
    Grid grid(0.0, S_max, N_x, T, N_t_stable);
    std::cout << "\nTesting explicit solver with safe time step..." << std::endl;

    ThetaMethodSolver solver(pde, grid, bc, ThetaMethodSolver::Scheme::Explicit);
    solver.solve();

    double price = solver.valueAt(S0);
    std::cout << "Option price at S=" << S0 << ": " << price << std::endl;

    // Should be close to Black-Scholes
    auto N = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
    double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double bs_price = S0 * N(d1) - K * std::exp(-r * T) * N(d2);

    double rel_error = std::abs(price - bs_price) / bs_price;
    std::cout << "Black-Scholes price: " << bs_price << std::endl;
    std::cout << "Relative error: " << 100.0 * rel_error << "%" << std::endl;

    EXPECT_LT(rel_error, 0.05) << "Should match Black-Scholes within 5%";

    std::cout << "\n✓ Explicit method works when stability condition is satisfied" << std::endl;
}

// int main(int argc, char** argv)
// {
//     std::cout << std::fixed << std::setprecision(6);
//     testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }

