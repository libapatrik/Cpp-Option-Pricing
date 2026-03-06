#ifndef CPPFM_COS_H
#define CPPFM_COS_H

#include <complex>
#include <functional>
#include <vector>

enum class NewtonMethod
{
    Original,
    Optimized
};

// COS Method: Fourier-Cosine Expansion for density/CDF recovery
// Fang & Oosterlee (2008), "A Novel Pricing Method for European Options"
class Transforms
{
public:
    // precomputed Fourier coefficients for fast CDF/PDF evaluation
    struct Coefficients {std::vector<double> F_k;      // Fourier coefficients
                         std::vector<double> omega;    // Frequencies omega_k = k*pi/(b-a)
                         double a;                     // Lower truncation bound
                         double b;                     // Upper truncation bound
                         size_t N;                     // Number of terms
    };

    // precompute F_k from characteristic function (the expensive part)
    static Coefficients precomputeCoefficients(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf);

    // fast CDF/PDF evaluation from precomputed coefficients
    static double evaluateCDF(const Coefficients& coeffs, double x);
    static double evaluatePDF(const Coefficients& coeffs, double x);

public:
    // recover CDF/PDF from characteristic function at points x
    static std::vector<double> recoverCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x);
    static std::vector<double> recoverPDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf, const std::vector<double>& x);

    // Newton-Raphson CDF inversion: find x s.t. F(x) = p
    // Original version recomputes F_k each iteration
    static std::pair<double, size_t> invertCDF(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf,
                                               double p, size_t max_iter = 100, double tol = 1e-8);

    // Optimized version: precomputes F_k once, 5-10x faster
     static std::pair<double, size_t> invertCDF_Optimized(double a, double b, size_t N, const std::function<std::complex<double>(double)>& chf,
                                                          double p, size_t max_iter = 100, double tol = 1e-8);

private:
    Transforms() = delete;
};

// Characteristic Function of Integrated Variance for CIR process
// Used in Broadie-Kaya scheme with COS method to sample asset price
// Broadie & Kaya (2006), "Exact Simulation of Stochastic Volatility"
class ChFIntegratedVariance
{
public:
    // phi(omega) = E[exp(i*omega * integral_s^t V(u)du)]
    // Uses modified Bessel functions I_nu(z)
    static std::complex<double> compute(double omega, double kappa, double vbar, double sigma, double v_s, double v_t, double tau);

private:
    // I_nu(z) for complex z - Boost doesn't support complex Bessel
    // power series for small |z|, asymptotic for large |z|
    static std::complex<double> modifiedBesselI(double nu, std::complex<double> z);

    ChFIntegratedVariance() = delete;
};


// Heston CF for log-returns: phi(u) = E[exp(iu * ln(S_T/S_0))]
// "Little Heston Trap" formulation (Albrecher et al. 2007) for numerical stability
class HestonCF
{
public:
    HestonCF(double kappa, double vbar, double sigma, double rho, double v0, double r, double T);

    // evaluate phi(u) for log-return X = ln(S_T/S_0)
    std::complex<double> operator()(double u) const;

    double maturity() const { return _T; }

private:
    double _kappa, _vbar, _sigma, _rho, _v0, _r, _T;
};

#endif // CPPFM_COS_H
