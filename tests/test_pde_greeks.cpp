#include <gtest/gtest.h>
#include "../PDEs/PDEGreeksCalculator.h"
#include "../PDEs/PDE.h"
#include "../PDEs/Grid.h"
#include "../PDEs/BoundaryConditions.h"
#include "../BlackScholesFormulas.h"
#include <cmath>

// ============================================================================
// Test Configuration
// ============================================================================

// TODO: Define test tolerances
// Define as constexpr - the value is known at compile time, and
// can be used in lambda functions without capture, and 
// they are guaranteed to exist before any object is constructed.

/// NOTE: OPTION_TOL is chosen such that tests fail, and we can see the differences.

constexpr double VALUE_TOL = 1e-3;   // For option value
constexpr double OPTION_TOL = 1e-3;  // For option value

constexpr double S0 = 100.0;     
constexpr double K = 100.0;
constexpr double r = 0.05;        
constexpr double sigma = 0.2;     
constexpr double T = 1.0;         

// ============================================================================
// Test Fixture
// ============================================================================
/// NOTE: For pdfFactory see: PDEs/PDEGreeksCalculator.h
class PDEGreeksTest : public ::testing::Test
{
protected:
    // TODO: Add setup code that runs before each test
    void SetUp() override
    {  
        // TODO: Create PDE factory for call options
        // This function should create a Black-Scholes PDE with given (sigma, r)
        // Formula: ∂V/∂t = 0.5*σ^2*S^2* ∂^2V/∂S^2 + r*S*∂V/∂S - r*V
        // Call Payoff: V = max(S - K, 0)
        pdeFactoryCall = [](double sig, double rate)->std::unique_ptr<PDE> {   // Black-Scholes pricing PDE
            return std::make_unique<VariableCoefficientPDE>(
                [sig](double S, double t) { return 0.5 * sig * sig * S * S; }, // diffusion function
                [rate](double S, double t) { return rate * S; },               // convection function
                [rate](double S, double t) { return -rate; },                  // reaction function
                [](double S, double t) { return 0.0; },                        // source function
                [](double S) { return std::max(S - K, 0.0); }                  // call payoff function
            );
        };

        pdeFactoryPut = [](double sig, double rate)->std::unique_ptr<PDE> {
            return std::make_unique<VariableCoefficientPDE>(
                [sig](double S, double t) { return 0.5 * sig * sig * S * S; }, 
                [rate](double S, double t) { return rate * S; },
                [rate](double S, double t) { return -rate; },
                [](double S, double t) { return 0.0; },
                [](double S) { return std::max(K - S, 0.0); }
            );
        };
    
        /// TODO: Create grid
        /**
         * ISSUE: When perturbing σ from 0.20 to 0.21, the diffusion range expands
         * With narrow domain [K/2, 2K], the upper boundary at S=2K contaminates
         * the solution at higher σ, creating asymmetric errors in the centered
         * finite difference for volga.
         *
         * TRADEOFF: Wider domain [K/5, 5K] removes boundary effects BUT increases
         * step size with fixed N_S = 4001:
         *   - Old: h = (2K - K/2)/4001 ≈ 0.0375  =>  h² = 0.0014 (finer mesh)
         *   - New: h = (5K - K/5)/4001 ≈ 0.12    =>  h² = 0.014  (10x coarser mesh)
         *
         * (Δσ)² = 10^-6 in finite differences.
         * So, coarser mesh on wider domain achieves correct volga.
         */

        S_min = 0.2 * K;   // K/5 for wide domain
        S_max = 5.0 * K;   // 5K
        N_S = 4001;        // Very fine grid for accurate second derivatives
        N_t = 4000;
        
        // Universal grid for both call and put
        grid = std::make_unique<Grid>(S_min, S_max, N_S, T, N_t, 
                                      Grid::SpacingType::Uniform); // Uniform or LogSpaced
                                      /// TODO: LogSpace is broken for Greeks, need to fix it.
        /**
        TODO: Create boundary conditions for European call
         *  - Lower (S=0): V = 0 (worthless if spot is zero)
         *  - Upper (S=S_max): V = S_max - K*e^(-r*(T-t)) (deep ITM)
         
         *  S_max = this->S_max is init-capture or generalised lambda capture
         *  Creates new variable S_max inside lambda, initialised with value of this->S_max
         *  This is required, because S_max is not static member variable, so we need to explicitly capture it
         *  this-> is accessing the member variable S_max of the current object
         */
        // Create boundary conditions for European call
        bcCall = std::make_unique<DirichletBC>(
           [](double t) { return 0.0; },
           [S_max = this->S_max](double t) { // S_max captured here
                return S_max - K * std::exp(-r * (T -t)); 
            }
        );
        // Create boundary conditions for European put
        bcPut = std::make_unique<DirichletBC>(
            [S_min = this->S_min](double t) {      // S_max captured here
                return K * std::exp(-r * (T - t)) - S_min;  // lower = Ke^(-r*(T-t)) - S_min
            },
            [](double t) { return 0.0; } // upper = 0
        );
        // Create base PDE
        basePDECall = pdeFactoryCall(sigma, r);
        // Create calculator for call
        calculatorCall = std::make_unique<PDEGreeksCalculator>(*basePDECall, *grid, *bcCall, sigma, r);
        calculatorCall->solve(S0); // Solve at spot S0

        // Create base PDE for put
        basePDEPut = pdeFactoryPut(sigma, r);
        // Create calculator for put
        calculatorPut = std::make_unique<PDEGreeksCalculator>(*basePDEPut, *grid, *bcPut, sigma, r);
        calculatorPut->solve(S0); // Solve at spot S0
    }
    
    // Member variables for reuse Call
    PDEGreeksCalculator::PDEFactory pdeFactoryCall;
    std::unique_ptr<Grid> grid;
    std::unique_ptr<BoundaryCondition> bcCall;
    std::unique_ptr<PDE> basePDECall;
    std::unique_ptr<PDEGreeksCalculator> calculatorCall;
    

    // Member variables for reuse Put
    PDEGreeksCalculator::PDEFactory pdeFactoryPut;      
    std::unique_ptr<BoundaryCondition> bcPut;
    std::unique_ptr<PDE> basePDEPut;
    std::unique_ptr<PDEGreeksCalculator> calculatorPut;


    double S_min, S_max;
    size_t N_S, N_t;
};

// ============================================================================
// TEST: All Greeks at Once
// ============================================================================
TEST_F(PDEGreeksTest, AllGreeksMatchAnalytical)
{
    std::cout << "TEST 1: PDE Greeks for European Call:" << std::endl;
    auto greeksCall = calculatorCall->computeAllGreeks(pdeFactoryCall);
    // Compare PDE Greeks to analytical Black-Scholes formulas for European Call
    EXPECT_NEAR(greeksCall.value, BlackScholesFormulas::callPrice(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksCall.delta, BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Call), OPTION_TOL);
    EXPECT_NEAR(greeksCall.gamma, BlackScholesFormulas::gamma(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksCall.vega, BlackScholesFormulas::vega(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksCall.theta, BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Call), OPTION_TOL);
    EXPECT_NEAR(greeksCall.rho, BlackScholesFormulas::rho(S0, K, r, sigma, T, Option::Type::Call), OPTION_TOL);
    EXPECT_NEAR(greeksCall.vanna, BlackScholesFormulas::vanna(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksCall.volga, BlackScholesFormulas::volga(S0, K, r, sigma, T), OPTION_TOL);

    
    std::cout << "\nTEST 2: PDE Greeks for European Put:" << std::endl;
    auto greeksPut = calculatorPut->computeAllGreeks(pdeFactoryPut);
    // Compare PDE Greeks to analytical Black-Scholes formulas for European Put
    EXPECT_NEAR(greeksPut.value, BlackScholesFormulas::putPrice(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksPut.delta, BlackScholesFormulas::delta(S0, K, r, sigma, T, Option::Type::Put), OPTION_TOL);
    EXPECT_NEAR(greeksPut.gamma, BlackScholesFormulas::gamma(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksPut.vega, BlackScholesFormulas::vega(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksPut.theta, BlackScholesFormulas::theta(S0, K, r, sigma, T, Option::Type::Put), OPTION_TOL);
    EXPECT_NEAR(greeksPut.rho, BlackScholesFormulas::rho(S0, K, r, sigma, T, Option::Type::Put), OPTION_TOL);
    EXPECT_NEAR(greeksPut.vanna, BlackScholesFormulas::vanna(S0, K, r, sigma, T), OPTION_TOL);
    EXPECT_NEAR(greeksPut.volga, BlackScholesFormulas::volga(S0, K, r, sigma, T), OPTION_TOL);
    
}
/** Why is volga bad?
 *  NOTE: Vanna seems fine.
    Volga is bad:
        - Solves 3 PDEs, each solution carries it numerical error
        - Division by Δσ^2 makes denominator small, leads to large errors
        - Requires much finer grid relative to other Greeks
 */



// ============================================================================
// TEST: Delta is Positive for Call
// ============================================================================
TEST_F(PDEGreeksTest, DeltaIsPositiveForCall)
{
    // Correct if they do not print anything
    double delta = calculatorCall->delta();
    EXPECT_GT(delta, 0.0);
    EXPECT_LT(delta, 1.0);
    EXPECT_NEAR(delta, 0.6, 0.15);
}

// ============================================================================
// TEST: Gamma is Positive
// ============================================================================
TEST_F(PDEGreeksTest, GammaIsPositive)
{
    double gamma = calculatorCall->gamma();
    EXPECT_GT(gamma, 0.0);
}

// ============================================================================
// TEST 10: Vega is Positive
// ============================================================================
TEST_F(PDEGreeksTest, VegaIsPositive)
{
    double vega = calculatorCall->vega(pdeFactoryCall);
    EXPECT_GT(vega, 0.0);
}

// Does delta converge with finer grid?
TEST(PDEGreeksConvergence, DeltaConvergesWithFinerGrid)
{
    // TODO: Test that delta gets more accurate with finer grids
    // Check that error vs analytical decreases
}

