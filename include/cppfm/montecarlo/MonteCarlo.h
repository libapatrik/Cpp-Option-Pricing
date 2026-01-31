/** 
    Main engine for MC

*/


#ifndef CPPFM_MONTECARLO_MONTECARLO_H
#define CPPFM_MONTECARLO_MONTECARLO_H

#include <vector>

struct MonteCarloResult {
    double mean = 0.0;
    double variance = 0.0;
    double stdError = 0.0;
    double ci95Lower = 0.0;
    double ci95Upper = 0.0;
    size_t numPaths = 0;

    static MonteCarloResult compute(const std::vector<double>& payoffs);
    static MonteCarloResult computeAntithetic(const std::vector<double>& payoffs);
};

#endif