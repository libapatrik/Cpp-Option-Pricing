#include <cppfm/montecarlo/VarianceReduction.h>
#include <cppfm/montecarlo/MonteCarlo.h>
#include <cmath>
#include <stdexcept>

VarianceReduction::VarianceReduction(VarianceReductionMethod method)
    : _method(method)
{
}

double ControlVariateAdjuster::computeOptimalC(
    const std::vector<double>& primary_payoffs,
    const std::vector<double>& control_payoffs)
{
    if (primary_payoffs.size() != control_payoffs.size()) {
        throw std::invalid_argument("payoff vectors must have same size");
    }
    size_t n = primary_payoffs.size();
    if (n < 2) return 0.0;

    // Welford single-pass for mean + covariance
    double mean_p = 0.0, mean_c = 0.0;
    double cov_pc = 0.0, var_c = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double delta_p = primary_payoffs[i] - mean_p;
        double delta_c = control_payoffs[i] - mean_c;
        mean_c += delta_c / (i + 1);
        mean_p += delta_p / (i + 1);
        // online covariance update
        double delta_c2 = control_payoffs[i] - mean_c;
        cov_pc += delta_p * delta_c2;
        var_c += delta_c * delta_c2;
    }

    if (var_c < 1e-15) return 0.0;  // degenerate case
    return cov_pc / var_c;
}

MonteCarloResult ControlVariateAdjuster::computeWithCV(
    const std::vector<double>& primary_payoffs,
    const std::vector<double>& control_payoffs,
    double control_analytical)
{
    if (primary_payoffs.size() != control_payoffs.size()) {
        throw std::invalid_argument("payoff vectors must have same size");
    }
    size_t n = primary_payoffs.size();
    if (n < 2) {
        return MonteCarloResult{};
    }

    double c = computeOptimalC(primary_payoffs, control_payoffs);

    // build adjusted payoffs
    std::vector<double> adjusted(n);
    for (size_t i = 0; i < n; ++i) {
        adjusted[i] = primary_payoffs[i] - c * (control_payoffs[i] - control_analytical);
    }

    return MonteCarloResult::compute(adjusted);
}

MonteCarloResult VarianceReduction::combine(const std::vector<double>& samples) const {
    switch (_method) {
        case VarianceReductionMethod::None:
            return MonteCarloResult::compute(samples);
        case VarianceReductionMethod::Antithetic:
            return MonteCarloResult::computeAntithetic(samples);
        default:
            throw std::runtime_error("Invalid variance reduction method");
    }
}

