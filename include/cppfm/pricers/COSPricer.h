#ifndef CPPFM_COSPRICER_H
#define CPPFM_COSPRICER_H

#include <complex>
#include <functional>
#include <vector>

// COS Method for European Option Pricing
// Fang & Oosterlee (2008)
class COSPricer
{
public:
    enum class OptionType { Call, Put };

    // --- sigma-hint bounds (existing API, unchanged) ---

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

    // --- cumulant-based bounds ---

    static double price(double S0, double K, double r, double T,
                        const std::function<std::complex<double>(double)>& chf,
                        OptionType type, size_t N, double L,
                        double c1, double c2, double c4);

    static double callPrice(double S0, double K, double r, double T,
                            const std::function<std::complex<double>(double)>& chf,
                            size_t N, double L,
                            double c1, double c2, double c4);

    static double putPrice(double S0, double K, double r, double T,
                           const std::function<std::complex<double>(double)>& chf,
                           size_t N, double L,
                           double c1, double c2, double c4);

    static std::vector<double> prices(double S0, const std::vector<double>& strikes,
                                       double r, double T,
                                       const std::function<std::complex<double>(double)>& chf,
                                       OptionType type, size_t N, double L,
                                       double c1, double c2, double c4);

    static std::vector<double> callPrices(double S0, const std::vector<double>& strikes,
                                           double r, double T,
                                           const std::function<std::complex<double>(double)>& chf,
                                           size_t N, double L,
                                           double c1, double c2, double c4);

    static std::vector<double> putPrices(double S0, const std::vector<double>& strikes,
                                          double r, double T,
                                          const std::function<std::complex<double>(double)>& chf,
                                          size_t N, double L,
                                          double c1, double c2, double c4);

private:
    // inner pricing with explicit bounds
    static double priceWithBounds(double S0, double K, double r, double T,
                                   const std::function<std::complex<double>(double)>& chf,
                                   OptionType type, size_t N, double a, double b);

    static std::vector<double> pricesWithBounds(double S0, const std::vector<double>& strikes,
                                                 double r, double T,
                                                 const std::function<std::complex<double>(double)>& chf,
                                                 OptionType type, size_t N, double a, double b);

    static double chi(size_t k, double a, double b, double c, double d);
    static double psi(size_t k, double a, double b, double c, double d);

    COSPricer() = delete;
};

#endif // CPPFM_COSPRICER_H
