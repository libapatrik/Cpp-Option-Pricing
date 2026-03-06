#ifndef CPPFM_UTILS_H
#define CPPFM_UTILS_H

#include <vector>
#include <functional>
#include <complex>
#include <cppfm/third_party/pcg/pcg_random.hpp>


class Utils
{
public: // expose methods which we want to expose
    // static method = class method and not an object/instance method!
    static double stdNormCdf(double x);

    static double stdNormPdf(double x);

    /**
     * Standard Normal Characteristic Function
     * ChF of N(0,1): φ(ω) = exp(-ω²/2)
     *
     * @param omega Frequency parameter
     * @return Complex characteristic function value
     */
    static std::complex<double> stdNormChF(double omega);

    /** Moro's Inverse Normal CDF Algorithm
     * Inverse of the standard normal CDF: Φ^(-1)(u)
     *
     * @param u Uniform random variable in (0, 1)
     * @return Standard normal random variable Z ~ N(0,1)
     */
    static double inverseNormalCDF(double u);

    // private - methods those intermediate methods
private:
    Utils() = delete; // delete constructor; everything is static
};




// ============================================================================
// Tridiagonal System Solver - Thomas Algorithm
// ============================================================================

/**
 * Thomas Algorithm: Solves the tridiagonal system:
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
 *   - Crank-Nicholson schemes
 */

class ThomasAlgorithm
{
public:
    static std::vector<double> solve(   // purely utility function without any state = static
        const std::vector<double>& lower,
        const std::vector<double>& diag,
        const std::vector<double>& upper,
        const std::vector<double>& rhs
    );

private:
    ThomasAlgorithm() = delete; // delete constructor; everything is static
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

/**
 * CIR (Cox-Ingersoll-Ross) Process Sampler
 *
 * Exact sampling from CIR process using noncentral chi-squared distribution.
 * Used in Broadie-Kaya scheme for exact variance simulation.
 *
 * CIR process: dV = κ(θ - V)dt + σ√V dW
 */
class CIRSampler
{
public:
    /**
     * Sample V(t) given V(s) for s < t using exact distribution
     *
     * @param kappa Mean reversion speed κ
     * @param vbar Long-term mean θ (v-bar)
     * @param sigma Volatility of volatility σ
     * @param v_s Initial variance V(s)
     * @param s Initial time
     * @param t Final time
     * @param randomEngine Random engine for sampling
     * @return Sample of V(t)
     *
     * Uses noncentral chi-squared distribution with:
     *   δ = 4κθ/σ²  (degrees of freedom)
     *   λ = 2c·exp(-κτ)·V(s)  (noncentrality parameter)
     *   where c = 2κ/(σ²(1-exp(-κτ))) and τ = t-s
     */

    /**
     * CIR sampling with pre-generated uniform (for antithetic sampling support)
     * @param u Uniform random variable in (0,1) for CDF inversion
     */
    static double sampleCIR(double kappa, double vbar, double sigma, double v_s, double s, double t, double u);

private:
    CIRSampler() = delete;
};

#endif //CPPFM_UTILS_H
