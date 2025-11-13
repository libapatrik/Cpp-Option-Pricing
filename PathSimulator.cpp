#include "PathSimulator.h"

#include <cmath>
#include <vector>
#include <random>

// GBM implementation
//TODO: Make it into one for loop - requires some vectorisation
std::vector<std::vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths) {
    // Seed
    std::default_random_engine generator(1);
    std::normal_distribution<double> distribution(0.0, 1.0); // mean 0, stddev 1

    double dt = T / static_cast<double>(numSteps);     // time step
    double drift = (r - 0.5 * pow(sigma, 2)) * dt;
    double vol_sqrt_dt = sigma * sqrt(dt);
    int numTotal = 2 * numPaths;                       // Antithetic samples
    // std::vector<std::vector<double>> paths(numTotal, std::vector<double>(numSteps + 1));

    // Want to do: np.zeros(); initialise path[0] = S0
    std::vector<std::vector<double>> S(numTotal, std::vector<double>(numSteps + 1, S0));

    // Simulate paths
    for (int i = 0; i < numPaths; ++i) {
        for (int j = 0; j < numSteps; ++j) {
            double Z = distribution(generator);

            // Paths first half
            S[i][j + 1] = S[i][j] * std::exp(drift + vol_sqrt_dt * Z);
            // Paths second half - antithetic
            S[i + numPaths][j + 1] = S[i + numPaths][j] * std::exp(drift + vol_sqrt_dt * (-Z));
        }
    }
    return S;
}

PathSimulator::PathSimulator(const std::vector<double> &timeSteps, const Model &model, size_t randomSeed)
    : _timeSteps(timeSteps), _modelPtr(model.clone()), _randomEngine(randomSeed)
{
    // _modelPtr(&model) IS WRONG!
    // model can be destroyed while a PathSimulator is still alive
    // _modelPtr should point to a COPY of model (use clone() method)
    // But, we cannot call the class of a Model - need define clone!
    // I am looking for a method that return a pointer to a copy of the object "model"

    // Call a method that does SANITY CHECKS on timeSteps
        // TimeSteps init value is zero
        // TimeSteps is strictly increasing sequence

    if (!timeStepsSanityCheck())
        throw std::runtime_error("The Time Steps are not correct!"); // exception
}

PathSimulator::~PathSimulator()
{
    // delete _modelPtr; // the 'new' called pointer in BSM is deleted here - when the PathSimulator object is destroyed
    if (_modelPtr) {  // Check if Null before deleting pointer
        delete _modelPtr;
        _modelPtr = nullptr;  // Prevent double deletion
    }
}


std::vector<double> PathSimulator::path() const
{
    // Initialisation stage : Path[0] = S0
    std::vector<double> path;
    path.push_back(_modelPtr->initValue());               // .append()

    // Iteration stage : Path[i] to Path[i+1]
    for (size_t idx = 0; idx < _timeSteps.size() - 1; ++idx)  // .size() is +1 bigger, so - 1!
        path.push_back(nextStep(idx, path[idx]));

    return path;
}

bool PathSimulator::timeStepsSanityCheck() const
{   // required more logic 
    if (_timeSteps.empty()) return false;
    if (_timeSteps[0] < 0.0) return false;
    for (size_t i = 1; i < _timeSteps.size(); ++i)
    {
        if (!std::isfinite(_timeSteps[i])) return false; // finiteness
        if (_timeSteps[i] <= _timeSteps[i-1]) return false; // increasing
    }
    return true;
}

double EulerPathSimulator::nextStep(size_t timeIndex, double assetPrice) const
{
    double deltaT = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

    // Rnadom Number Generation
    std::normal_distribution<double> distribution(0.0, 1.0);
    double Z = distribution(_randomEngine);

    // Apply Euler scheme
    double nextStep = assetPrice + _modelPtr->drift(_timeSteps[timeIndex], assetPrice) * deltaT
                    + _modelPtr->diffusion(_timeSteps[timeIndex], assetPrice) * std::sqrt(deltaT) * Z;

    return nextStep;
}

double MilsteinPathSimulator::nextStep(size_t timeIndex, double assetPrice) const
{
    double t = _timeSteps[timeIndex];
    double deltaT = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

    // Rnadom Number Generation
    std::normal_distribution<double> distribution(0.0, 1.0);
    double Z = distribution(_randomEngine);

    // 
    double mu = _modelPtr->drift(t, assetPrice);
    double sigma = _modelPtr->diffusion(t, assetPrice);

    // Derivative of diffusion(t, S) wrt S by central difference
    // TODO: Create pure virtual method for diffusion(t,S)/dS
    const double rel_eps = 1e-6; // small perturbation; avoid division by zero
    double sigma_plus = _modelPtr->diffusion(t, assetPrice*(1. + rel_eps));
    double sigma_minus = _modelPtr->diffusion(t, assetPrice * (1. - rel_eps));
    double sigma_prime = (sigma_plus - sigma_minus) / (2.0 * rel_eps * assetPrice);

    double dW = std::sqrt(deltaT) * Z;
    // Milstein scheme: X_{n+1} = X_n + μ(X_n)Δt + σ(X_n)ΔW + 0.5 * σ(X_n) * σ'(X_n) * (ΔW² - Δt)
    double next = assetPrice + mu * deltaT + sigma * dW + 0.5 * sigma * sigma_prime * (std::pow(dW, 2) - deltaT);
    return next;
}
