#ifndef CPPFM_UTILS_H
#define CPPFM_UTILS_H

#include <vector>
#include <functional>

class Utils
{
public: // expose methods which we want to expose
    // static method = class method and not a object/instance method!
    static double stdNormCdf(double x);
// private - methods those intermediate methods
private:

};

// ============================================================================
// Tridiagonal System Solver - Thomas Algorithm
// ============================================================================

/**
 * @namespace TridiagonalSolver
 * @brief Efficient solver for tridiagonal linear systems using Thomas Algorithm
 * 
 * Solves the tridiagonal system:
 *   a[i]·x[i-1] + b[i]·x[i] + c[i]·x[i+1] = d[i]
 * 
 * for i = 0, 1, ..., n-1
 * 
 * Where:
 *   - a[i]: lower diagonal (subdiagonal)
 *   - b[i]: main diagonal
 *   - c[i]: upper diagonal (superdiagonal)
 *   - d[i]: right-hand side
 * 
 * Boundary conditions are implicit:
 *   - a[0] is not used (no x[-1])
 *   - c[n-1] is not used (no x[n])
 * 
 * Complexity: O(n) - much more efficient than general Gaussian elimination O(n³)
 * 
 * Used in:
 *   - Cubic spline interpolation (computing second derivatives)
 *   - Implicit PDE solvers (solving for u^{n+1})
 *   - Crank-Nicolson schemes
 */
 
class ThomasAlgorithm 
{
public:
    static std::vector<double> solve(      // purely utility function without any state = static
        const std::vector<double>& lower,
        const std::vector<double>& diag,
        const std::vector<double>& upper,
        const std::vector<double>& rhs
    );

private:
    // This is just utility class, so we don't need to instantiate it
    /**
     *   Someone could do this mistake:
     *   ThomasAlgorithm solver;  // This compiles but is POINTLESS
     *   auto result = ThomasAlgorithm::solve(...);  // The object 'solver' is never used
     */
    ThomasAlgorithm() = delete; // prevent instantiation
};

class NumericalDerivatives
{
public:
    // Scale-invariant numerical derivatives: step size adapts to magnitude of x
    static double firstDerivative(std::function<double(double)> f, double x, double h = 1e-4);
    static double secondDerivative(std::function<double(double)> f, double x, double h = 1e-4);

private:
    NumericalDerivatives() = delete; // Same thing to prevent instantiation as in ThomasAlgorithm
};




#endif //CPPFM_UTILS_H