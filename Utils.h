#ifndef CPPFM_UTILS_H
#define CPPFM_UTILS_H

#include <vector>
#include <functional>
#include <complex>


/**
 *  NOTES:
 *  (1) There are 2 implementations of Newton's method for CDF inversion: 
 *      (a) that I developed from scratch for Computational Finance
 *      (b) an optimized version that uses precomputed coefficients .
 *
 *  (2) Boost library used in: CIRSampler (sampling from non-chi-squared distribution), ChFIntegratedVarianceProcess (Bessel functions)
 *
 * 
 *  TODO:
 *  (1) Optimize std::exp/log/sqrt with Accelerate.h
 *  (2) Implement generic n-dim Cholesky decomposition to decorrelate the Brownian motions
 *  (3) Develop a linspace function as in Python
 *  (4) Develop a creating a table function?
 */


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
 * COS Method: Fourier-Cosine Expansion Method for Option Pricing
 * 
 * Recovers probability density and cumulative distribution functions
 * from characteristic functions using Fourier-cosine series expansion.
 * 
 * Reference: Fang & Oosterlee (2008), "A Novel Pricing Method for European Options"
 */
class Transforms
{
public:
    /**
     * Precomputed Fourier coefficients for fast CDF/PDF evaluation
     * Avoids recomputing F_k in each Newton iteration
     */
    struct Coefficients {std::vector<double> F_k;      // Fourier coefficients
                         std::vector<double> omega;    // Frequencies ω_k = kπ/(b-a)
                         double a;                     // Lower truncation bound
                         double b;                     // Upper truncation bound
                         size_t N;                     // Number of terms
    };
    
    /**
     * Precompute Fourier coefficients F_k from characteristic function
     * 
     * This is the expensive part - computing ChF at N frequencies.
     * By caching these, we can evaluate CDF/PDF very quickly.
     * 
     * @param a Lower truncation bound
     * @param b Upper truncation bound
     * @param N Number of terms in cosine expansion
     * @param chf Characteristic function
     * @return Precomputed coefficients
     */
    static Coefficients precomputeCoefficients(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf);
    
    /**
     * Fast CDF evaluation using precomputed coefficients
     * @param coeffs Precomputed Fourier coefficients
     * @param x Point at which to evaluate CDF
     * @return CDF value
     */
    static double evaluateCDF(const Coefficients& coeffs, double x);
    
    /**
     * Fast PDF evaluation using precomputed coefficients
     * @param coeffs Precomputed Fourier coefficients
     * @param x Point at which to evaluate PDF
     * @return PDF value
     */
    static double evaluatePDF(const Coefficients& coeffs, double x);
    
public:
    /**
     * Recover CDF from characteristic function
     * 
     * @param a Lower truncation bound
     * @param b Upper truncation bound
     * @param N Number of terms in cosine expansion
     * @param chf Characteristic function evaluated at frequencies omega
     * @param x Points at which to evaluate CDF
     * @return CDF values at x
     */
    static std::vector<double> recoverCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x);
    
    /**
     * Recover PDF from characteristic function
     * 
     * @param a Lower truncation bound
     * @param b Upper truncation bound
     * @param N Number of terms in cosine expansion
     * @param chf Characteristic function evaluated at frequencies omega
     * @param x Points at which to evaluate PDF
     * @return PDF values at x
     */
    static std::vector<double> recoverPDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x);
    
    /**
     * Invert CDF using Newton-Raphson method (Original Version)
     * 
     * Given a target probability p, find x such that F(x) = p
     * where F is the CDF recovered via COS method.
     * 
     * Uses Newton-Raphson: x_{n+1} = x_n - (F(x_n) - p) / f(x_n)
     * where f is the PDF (derivative of CDF).
     * 
     * @param a Lower truncation bound
     * @param b Upper truncation bound
     * @param N Number of terms in cosine expansion
     * @param chf Characteristic function
     * @param p Target probability in (0, 1)
     * @param max_iter Maximum number of Newton iterations
     * @param tol Convergence tolerance
     * @return x such that CDF(x) ≈ p
     * 
     * Algorithm:
     * 1. Evaluate CDF at multiple initial points to find good starting guess
     * 2. Use Newton-Raphson with COS-recovered CDF and PDF
     * 3. Apply safeguards for bounds and small derivatives
     * 
     * Note: This version recomputes F_k in each iteration (clear but not optimal).
     *       For better performance, use invertCDF_Optimized().
     */
    static std::pair<double, size_t> invertCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf,
                                               double p, size_t max_iter = 100, double tol = 1e-8);
    
    /**
     * Invert CDF using Optimized Newton-Raphson (Coefficient Caching)
     * 
     * Optimized version that precomputes F_k coefficients once,
     * then uses fast evaluation in Newton iterations.
     * 
     * Performance: 5-10x faster than invertCDF() due to coefficient caching.
     * 
     * @param a Lower truncation bound
     * @param b Upper truncation bound
     * @param N Number of terms in cosine expansion
     * @param chf Characteristic function
     * @param p Target probability in (0, 1)
     * @param max_iter Maximum number of Newton iterations
     * @param tol Convergence tolerance
     * @return x such that CDF(x) ≈ p
     * 
     * Algorithm:
     * 1. Precompute F_k coefficients ONCE (expensive)
     * 2. Fast initial guess using cached coefficients
     * 3. Newton-Raphson with fast CDF/PDF evaluation (cheap)
     * 4. Safeguards for bounds and small derivatives
     * 
     * Complexity:
     * - Precomputation: O(N) ChF evaluations (expensive)
     * - Per Newton iteration: O(N) arithmetic operations (cheap)
     * - Total: O(N) + O(iterations × N) ≈ O(N × iterations)
     *   vs original O(N × iterations × 2) since we avoid recomputing F_k
     */
     static std::pair<double, size_t> invertCDF_Optimized(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, 
                                                          double p, size_t max_iter = 100, double tol = 1e-8);

private:
    Transforms() = delete;
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

/**
 * Characteristic Function of Integrated Variance
 * 
 * Computes the characteristic function of ∫[s,t] V(u)du for the CIR process.
 * Used in Broadie-Kaya scheme with COS method to sample asset price.
 * 
 * Reference: Broadie & Kaya (2006), "Exact Simulation of Stochastic Volatility"
 */
class ChFIntegratedVariance
{
public:
    /**
     * Compute characteristic function of integrated variance
     * 
     * @param omega Frequency parameter
     * @param kappa Mean reversion speed κ
     * @param vbar Long-term mean θ
     * @param sigma Volatility of volatility σ (or γ in some papers)
     * @param v_s Variance at start V(s)
     * @param v_t Variance at end V(t)
     * @param tau Time interval τ = t - s
     * @return Complex characteristic function value
     * 
     * Uses modified Bessel functions of the first kind I_ν(z)
     */
    static std::complex<double> compute(double omega, double kappa, double vbar, double sigma, double v_s, double v_t, double tau);

private:
    /**
     * Modified Bessel Function I_ν(z) for complex argument
     * - we need to manually implement this as Boost does not support complex Bessel functions
     * Equivalent to scipy.special.iv(nu, z) in Python
     * Uses power series for small |z|, asymptotic expansion for large |z|
     * 
     * @param nu Order of Bessel function
     * @param z Complex argument
     * @return I_ν(z)
     */
    static std::complex<double> modifiedBesselI(double nu, std::complex<double> z);
    
    ChFIntegratedVariance() = delete;
};

#endif //CPPFM_UTILS_H