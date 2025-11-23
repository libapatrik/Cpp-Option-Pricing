#include "Utils.h"
#include <cmath>
#include <stdexcept>

// Standard normal CDF
double Utils::stdNormCdf(double x)
{
    return 0.5 * erfc(-x / sqrt(2.0));
}

// ============================================================================
// Tridiagonal System Solver Implementation
// ============================================================================
// Used in: Cubic Spline Interpolation, PDEs (implicit, crank-nicolson)


std::vector<double> ThomasAlgorithm::solve(
    const std::vector<double>& lower,
    const std::vector<double>& diag,
    const std::vector<double>& upper,
    const std::vector<double>& rhs)
{
    // Validate input sizes
    size_t n = rhs.size();
    
    if (lower.size() != n || diag.size() != n || upper.size() != n) {
        throw std::invalid_argument("TridiagonalSolver::solve: All input vectors must have the same size");
    }
    
    if (n == 0) {
        throw std::invalid_argument("TridiagonalSolver::solve: Cannot solve empty system");
    }
    
    if (n == 1) {
        // Trivial case: b[0]*x[0] = d[0]
        if (std::abs(diag[0]) < 1e-14) {
            throw std::runtime_error(
                "TridiagonalSolver::solve: Singular matrix (diagonal element is zero)");
        }
        return {rhs[0] / diag[0]};
    }

    // Diagonal dominance, M-matrices etc.
    
    // ========================================================================
    // THOMAS ALGORITHM (General Implementation)
    // ========================================================================
    
    // Temporary storage for modified coefficients
    std::vector<double> c_prime(n, 0.0);  // Modified upper diagonal
    std::vector<double> d_prime(n, 0.0);  // Modified RHS
    std::vector<double> x(n, 0.0);        // Solution vector
    
    // ------------------------------------------------------------------------
    /// PROCESS: Mx = r -> LUx = r -> (1) Ly = r -> (2) Ux = y
    // Step 1: Forward Elimination
    // ------------------------------------------------------------------------
    // Transform the system to upper triangular form
    // Original:  b[i]*x[i] + c[i]*x[i+1] = d[i] - a[i]*x[i-1]
    // Modified:  x[i] + c'[i]*x[i+1] = d'[i]
    
    // First row (i = 0): b[0]*x[0] + c[0]*x[1] = d[0]
    if (std::abs(diag[0]) < 1e-14) {
        throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row 0");
    }
    c_prime[0] = upper[0] / diag[0];         // γ1 = c1/b1
    d_prime[0] = rhs[0] / diag[0];           // ρ1 = r1/b1
    
    // Interior rows (i = 1, 2, ..., n-1)
    for (size_t i = 1; i < n; ++i) {
        // Denominator after eliminating x[i-1]
        double m = diag[i] - lower[i] * c_prime[i-1];
        //          b[i] -  a[i] * γ[i-1]
        if (std::abs(m) < 1e-14) {
            throw std::runtime_error("TridiagonalSolver::solve: Singular matrix at row " + std::to_string(i));
        }
        
        c_prime[i] = upper[i] / m;
        // c[i] = c[i+1] / (b[i] - a[i] * γ[i-1])
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i-1]) / m;
        // d'[i] = (d[i] - a[i] * ρ[i-1]) / (b[i] - a[i] * γ[i-1])
    }
    
    // ------------------------------------------------------------------------
    // Step 2: Back Substitution
    // ------------------------------------------------------------------------
    // Solve the upper triangular system from bottom to top
    // x[i] = d'[i] - c'[i]*x[i+1]
    
    // Last equation (i = n-1): x[n-1] = d'[n-1]
    x[n-1] = d_prime[n-1];
    
    // Work backwards (i = n-2, n-3, ..., 0)
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
        //     ρ[i+1]     - γ[i+1]     * x[i+2]
        //  (Code i=k corresponds to PDF row k+1)
    }
    
    return x;
}       // end of ThomasAlgorithm::solve


// ============================================================================
// Numerical Derivatives Implementation
// ============================================================================

double NumericalDerivatives::firstDerivative(std::function<double(double)> f, double x, double h)
{
    // Scale-invariant step size: h_scaled = h * max(|x|, 1.0)
    // This ensures relative perturbation is consistent across different magnitudes
    double h_scaled = h * std::max(std::abs(x), 1.0);
    
    // Central difference: f'(x) = [f(x+h) - f(x-h)] / (2h)
    return (f(x + h_scaled) - f(x - h_scaled)) / (2.0 * h_scaled);
}

double NumericalDerivatives::secondDerivative(std::function<double(double)> f, double x, double h)
{
    // Scale-invariant step size
    double h_scaled = h * std::max(std::abs(x), 1.0);
    
    // Central second difference: f''(x) = [f(x+h) - 2f(x) + f(x-h)] / h^2
    return (f(x + h_scaled) - 2.0 * f(x) + f(x - h_scaled)) / (h_scaled * h_scaled);
}


