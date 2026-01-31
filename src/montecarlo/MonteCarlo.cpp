
#include <cppfm/montecarlo/MonteCarlo.h>
#include <numeric>
#include <cmath>

MonteCarloResult MonteCarloResult::compute(const std::vector<double>& samples) {
    size_t n = samples.size();
    if (n == 0) return {};

    double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
    double mean = sum / static_cast<double>(n);
    
    double sqSum = 0.0;
    for (double s : samples) {
        sqSum += (s - mean) * (s - mean);
    }
    double variance = sqSum / static_cast<double>(n - 1);
    double stdError = std::sqrt(variance / static_cast<double>(n));

    double ci95Lower = mean - 1.96 * stdError;
    double ci95Upper = mean + 1.96 * stdError;

    return {mean, variance, stdError, ci95Lower, ci95Upper, n};
}

MonteCarloResult MonteCarloResult::computeAntithetic(const std::vector<double>& samples) {
    size_t n = samples.size();
    if (n < 2 || n % 2 != 0) return compute(samples);

    // antithetic: average each (original, negated) pair
    size_t numPairs = n / 2;
    std::vector<double> paired(numPairs);
    for (size_t i = 0; i < numPairs; ++i) {
        paired[i] = (samples[2*i] + samples[2*i + 1]) / 2.0;
    }

    // stats on paired values
    double sum = 0.0;
    for (double p : paired) sum += p;
    double mean = sum / static_cast<double>(numPairs);

    double sqSum = 0.0;
    for (double p : paired) sqSum += (p - mean) * (p - mean);
    double variance = sqSum / static_cast<double>(numPairs - 1);
    double stdError = std::sqrt(variance / static_cast<double>(numPairs));

    double ci95Lower = mean - 1.96 * stdError;
    double ci95Upper = mean + 1.96 * stdError;

    return {mean, variance, stdError, ci95Lower, ci95Upper, n};
}