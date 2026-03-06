#ifndef CPPFM_COSPRICER_H
#define CPPFM_COSPRICER_H

#include <complex>
#include <functional>
#include <vector>

// COS Method for European Option Pricing
// Fang & Oosterlee (2008)
// V = e^{-rT} * sum' Re[phi(k*pi/(b-a)) * e^{ik*pi*(x-a)/(b-a)}] * V_k
class COSPricer
{
public:
    enum class OptionType { Call, Put };

    // price single European option
    static double price(double S0, double K, double r, double T,
                        const std::function<std::complex<double>(double)>& chf,
                        OptionType type = OptionType::Call,
                        size_t N = 256, double L = 10.0, double sigma = 0.0);

    static double callPrice(double S0, double K, double r, double T,
                            const std::function<std::complex<double>(double)>& chf,
                            size_t N = 256, double L = 10.0, double sigma = 0.0);

    static double putPrice(double S0, double K, double r, double T,
                           const std::function<std::complex<double>(double)>& chf,
                           size_t N = 256, double L = 10.0, double sigma = 0.0);

    // vectorized pricing - multiple strikes at same maturity
    // phi(omega_k) depends only on T, not K, so precompute once
    // O(N_strikes * N_terms) CF evals -> O(N_terms)
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

private:
    // chi_k(c,d) = integral[c,d] e^x * cos(k*pi*(x-a)/(b-a)) dx
    // Fang & Oosterlee Eq. (22)
    static double chi(size_t k, double a, double b, double c, double d);

    // psi_k(c,d) = integral[c,d] cos(k*pi*(x-a)/(b-a)) dx
    // Fang & Oosterlee Eq. (23)
    static double psi(size_t k, double a, double b, double c, double d);

    COSPricer() = delete;
};

#endif // CPPFM_COSPRICER_H
