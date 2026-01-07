#ifndef CPPFM_UTILS_H
#define CPPFM_UTILS_H

#include <vector>
#include <functional>
#include <complex>
#include "fast_rng/pcg_random.hpp"


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
 *  (1.1) Make the PDE solvers more efficient 
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


/**
 * Heston Characteristic Function for Log-Returns
 *
 * φ(u) = E[exp(iu · ln(S_T/S_0))]
 *
 * Uses the "Little Heston Trap" formulation (Albrecher et al. 2007) for numerical stability.
 * Branch selection enforces Re(d) > 0 to avoid discontinuities from complex sqrt.
 */
class HestonCF
{
public:
    /**
     * @param kappa Mean reversion speed κ
     * @param vbar Long-term variance θ
     * @param sigma Vol-of-vol σ (xi in some notation)
     * @param rho Correlation ρ ∈ [-1, 1]
     * @param v0 Initial variance V(0)
     * @param r Risk-free rate
     * @param T Time to maturity
     */
    HestonCF(double kappa, double vbar, double sigma, double rho, double v0, double r, double T);

    /**
     * Evaluate φ(u) for log-return X = ln(S_T/S_0)
     */
    std::complex<double> operator()(double u) const;

    double maturity() const { return _T; }

private:
    double _kappa, _vbar, _sigma, _rho, _v0, _r, _T;
};


/**
 * COS Method for European Option Pricing
 *
 * Fang & Oosterlee (2008): "A Novel Pricing Method for European Options
 * Based on Fourier-Cosine Series Expansions"
 *
 * V = e^{-rT} · Σ' Re[φ(kπ/(b-a)) · e^{ikπ(x-a)/(b-a)}] · V_k
 *
 * where V_k are payoff coefficients computed from χ_k (for e^x) and ψ_k (for 1).
 */
class COSPricer
{
public:
    enum class OptionType { Call, Put };

    /**
     * Price European option using COS method
     *
     * @param S0 Spot price
     * @param K Strike
     * @param r Risk-free rate
     * @param T Maturity
     * @param chf Characteristic function of log-return ln(S_T/S_0)
     * @param type Call or Put
     * @param N Number of Fourier terms (128-256 typical)
     * @param L Truncation parameter: domain = [x0 - L·σ·√T, x0 + L·σ·√T]
     * @param sigma Volatility estimate for domain sizing (use √v0 or √vbar for Heston)
     */
    static double price(double S0, double K, double r, double T,
                        const std::function<std::complex<double>(double)>& chf,
                        OptionType type = OptionType::Call,
                        size_t N = 256, double L = 10.0, double sigma = 0.0);

    // Convenience wrappers
    static double callPrice(double S0, double K, double r, double T,
                            const std::function<std::complex<double>(double)>& chf,
                            size_t N = 256, double L = 10.0, double sigma = 0.0);

    static double putPrice(double S0, double K, double r, double T,
                           const std::function<std::complex<double>(double)>& chf,
                           size_t N = 256, double L = 10.0, double sigma = 0.0);

    // =========================================================================
    // Vectorized pricing - multiple strikes at same maturity
    // =========================================================================
    // Key optimization: Characteristic function φ(ω_k) depends only on T, not K.
    // By precomputing φ(ω_k) once and reusing for all strikes, we get massive
    // speedup when pricing many options at the same maturity.
    //
    // Speedup: O(N_strikes × N_terms) CF evaluations → O(N_terms) CF evaluations

    /**
     * Price multiple European options at the same maturity using vectorized COS
     *
     * @param S0 Spot price
     * @param strikes Vector of strike prices
     * @param r Risk-free rate
     * @param T Maturity (same for all strikes)
     * @param chf Characteristic function of log-return
     * @param type Call or Put
     * @param N Number of Fourier terms
     * @param L Truncation parameter
     * @param sigma Volatility estimate for domain sizing
     * @return Vector of option prices corresponding to each strike
     */
    static std::vector<double> prices(double S0, const std::vector<double>& strikes,
                                       double r, double T,
                                       const std::function<std::complex<double>(double)>& chf,
                                       OptionType type = OptionType::Call,
                                       size_t N = 256, double L = 10.0, double sigma = 0.0);

    static std::vector<double> callPrices(double S0, const std::vector<double>& strikes,
                                           double r, double T,
                                           const std::function<std::complex<double>(double)>& chf,
                                           size_t N = 256, double L = 10.0, double sigma = 0.0);

    static std::vector<double> putPrices(double S0, const std::vector<double>& strikes,
                                          double r, double T,
                                          const std::function<std::complex<double>(double)>& chf,
                                          size_t N = 256, double L = 10.0, double sigma = 0.0);

    /**
     * Implied volatility via Newton-Raphson on Black-Scholes
     */
    static double impliedVol(double price, double S0, double K, double r, double T,
                             bool isCall, double tol = 1e-8, size_t maxIter = 100);

private:
    /**
     * χ_k(c,d) = ∫[c,d] e^x · cos(kπ(x-a)/(b-a)) dx
     * Closed form from Fang & Oosterlee Eq. (22)
     */
    static double chi(size_t k, double a, double b, double c, double d);

    /**
     * ψ_k(c,d) = ∫[c,d] cos(kπ(x-a)/(b-a)) dx
     * Closed form from Fang & Oosterlee Eq. (23)
     */
    static double psi(size_t k, double a, double b, double c, double d);

    static double bsCall(double S, double K, double r, double T, double vol);
    static double bsVega(double S, double K, double r, double T, double vol);

    COSPricer() = delete;
};


/**
 * Analytical Local Volatility for Heston Model
 *
 * Computes Dupire local variance directly from Heston characteristic function
 * using COS method for pricing and numerical derivatives on prices.
 *
 * This approach matches the van der Stoep paper methodology:
 * - Direct finite differences on Fourier-generated prices
 * - NO intermediate implied volatility interpolation
 * - NO IV surface artifacts
 *
 * Dupire Formula (zero dividend):
 *   σ²_LV(K,T) = 2 · (∂C/∂T + rK·∂C/∂K) / (K² · ∂²C/∂K²)
 *
 * Reference: Gatheral (2006), "The Volatility Surface", Chapter 3
 */
class HestonLocalVol
{
public:
    /**
     * @param kappa Mean reversion speed κ
     * @param vbar Long-term variance θ
     * @param sigma Vol-of-vol σ
     * @param rho Correlation ρ ∈ [-1, 1]
     * @param v0 Initial variance V(0)
     * @param r Risk-free rate
     * @param S0 Spot price
     */
    HestonLocalVol(double kappa, double vbar, double sigma, double rho,
                   double v0, double r, double S0);

    /**
     * Compute local variance σ²_LV(K,T) at given strike and maturity
     * Delegates to localVolsAtTime for efficiency.
     */
    double localVariance(double K, double T) const;

    /**
     * Compute local volatility σ_LV(K,T) = √(σ²_LV)
     * Delegates to localVolsAtTime for efficiency.
     */
    double localVol(double K, double T) const;

    /**
     * Compute call price using COS method
     * Delegates to callPrices for efficiency.
     */
    double callPrice(double K, double T) const;

    /**
     * Set COS method parameters
     * @param N Number of Fourier terms (default 256)
     * @param L Truncation parameter (default 12)
     */
    void setCOSParams(size_t N, double L);

    // =========================================================================
    // Vectorized batch computation (primary implementation)
    // =========================================================================
    // All single-point methods delegate to these batch methods.
    // Uses vectorized COS pricing for optimal performance.

    /**
     * Compute local volatilities for all spots at same maturity
     * Primary implementation using vectorized COS.
     *
     * @param spots Vector of spot prices (typically the grid points)
     * @param T Time to maturity
     * @return Vector of local volatilities σ_LV for each spot
     */
    std::vector<double> localVolsAtTime(const std::vector<double>& spots, double T) const;

    /**
     * Compute call prices for multiple strikes at same maturity
     * Primary implementation using vectorized COS.
     */
    std::vector<double> callPrices(const std::vector<double>& strikes, double T) const;

private:
    double _kappa, _vbar, _sigma, _rho, _v0, _r, _S0;
    size_t _N = 256;    // COS terms
    double _L = 12.0;   // Truncation parameter
};

#endif //CPPFM_UTILS_H