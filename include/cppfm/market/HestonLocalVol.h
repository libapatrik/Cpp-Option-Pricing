#ifndef CPPFM_HESTONLOCALVOL_H
#define CPPFM_HESTONLOCALVOL_H

#include <cstddef>
#include <vector>

// Analytical Local Volatility for Heston Model
// Dupire local variance from Heston CF using COS method + numerical derivatives
// Direct finite differences on Fourier-generated prices (no IV surface intermediary)
//
// sigma^2_LV(K,T) = 2 * (dC/dT + rK*dC/dK) / (K^2 * d2C/dK2)
// Gatheral (2006), "The Volatility Surface", Chapter 3
class HestonLocalVol
{
public:
    HestonLocalVol(double kappa, double vbar, double sigma, double rho,
                   double v0, double r, double S0);

    double localVariance(double K, double T) const;
    double localVol(double K, double T) const;
    double callPrice(double K, double T) const;

    void setCOSParams(size_t N, double L);

    // batch computation (primary implementation, single-point methods delegate here)
    std::vector<double> localVolsAtTime(const std::vector<double>& spots, double T) const;
    std::vector<double> callPrices(const std::vector<double>& strikes, double T) const;

private:
    double _kappa, _vbar, _sigma, _rho, _v0, _r, _S0;
    size_t _N = 256;    // COS terms
    double _L = 12.0;   // Truncation parameter
};

#endif // CPPFM_HESTONLOCALVOL_H
