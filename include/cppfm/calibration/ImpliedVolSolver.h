#ifndef CPPFM_IMPLIEDVOLSOLVER_H
#define CPPFM_IMPLIEDVOLSOLVER_H

#include <cstddef>

// Brent's method + Newton acceleration for IV inversion
class ImpliedVolSolver
{
public:
    static double solve(double price, double S0, double K, double r, double T,
                        bool isCall, double tol = 1e-8, size_t maxIter = 100);

private:
    ImpliedVolSolver() = delete;
};

#endif // CPPFM_IMPLIEDVOLSOLVER_H
