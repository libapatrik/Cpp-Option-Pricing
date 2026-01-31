

#ifndef CPPFM_MONTECARLO_VARIANCE_REDUCTION_H

#define CPPFM_MONTECARLO_VARIANCE_REDUCTION_H

#include <cppfm/montecarlo/MonteCarlo.h>
#include <vector>

// Forwards declaration
struct MonteCarloResult;

// enum class based variance reduction switcher
enum class VarianceReductionMethod {
    None = 0,
    Antithetic = 1,
    ControlVariate = 2,
    // ImportanceSampling = 3,
    // StratifiedSampling = 4,
    // QuasiMonteCarlo = 5,
    // MultiLevelMonteCarlo = 6,
};

class VarianceReduction {
public:
    explicit VarianceReduction(VarianceReductionMethod method = VarianceReductionMethod::None);

    virtual ~VarianceReduction() = default;
    virtual VarianceReduction* clone() const = 0;

    // call before each path
    void beforePath(size_t pathIndex, size_t totalPaths) {}
    // transform the random number for the path
    double transformRandom(double z) const { return z; }
    // combine the samples
    MonteCarloResult combine(const std::vector<double>& samples) const;

private:
    VarianceReductionMethod _method;
};

// standalone CV helper - works with any MC that can replay RNG
class ControlVariateAdjuster {
public:
    // c = cov(primary, control) / var(control)
    static double computeOptimalC(
        const std::vector<double>& primary_payoffs,
        const std::vector<double>& control_payoffs);

    // adjusted[i] = primary[i] - c*(control[i] - control_analytical)
    static MonteCarloResult computeWithCV(
        const std::vector<double>& primary_payoffs,
        const std::vector<double>& control_payoffs,
        double control_analytical);
};

#endif