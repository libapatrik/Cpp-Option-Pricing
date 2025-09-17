#include "PathSimulator.h"

#include <cmath>
#include <vector>
#include <random>

using namespace std;

// GBM implementation
//TODO: Make it into one for loop - requires some vectorisation
vector<vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths) {
    // Seed
    default_random_engine generator(1);
    normal_distribution<double> distribution(0.0, 1.0); // mean 0, stddev 1

    double dt = T / static_cast<double>(numSteps);     // time step
    double drift = r - 0.5 * pow(sigma, 2) * dt;
    double vol_sqrt_dt = sigma * sqrt(dt);
    int numTotal = 2 * numPaths;                       // Antithetic samples
    // vector<vector<double>> paths(numTotal, vector<double>(numSteps + 1));

    // Want to do: np.zeros(); initialise path[0] = S0
    vector<vector<double>> S(numTotal, vector<double>(numSteps + 1, S0));

    // Simulate paths
    for (int i = 0; i < numPaths; ++i) {
        for (int j = 0; j <= numSteps; ++j) {
            double Z = distribution(generator);

            // Paths first half
            S[i][j + 1] = S[i][j] * exp(drift + vol_sqrt_dt * Z);
            // Paths second half - antithetic
            S[i + numPaths][j + 1] = S[i + numPaths][j] * exp(drift + vol_sqrt_dt * (-Z));
        }
    }
    return S;
}


PathSimulator::PathSimulator(const vector<double> &timeSteps, const Model &model, size_t randomSeed)
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
        throw "The Time Steps are not correct!";
}

PathSimulator::~PathSimulator()
{
    delete _modelPtr; // the 'new' called pointer in BSM is deleted here - when the PathSimulator object is destroyed
}


vector<double> PathSimulator::path() const
{
    // Initialisation stage : Path[0] = S0
    vector<double> path;
    path.push_back(_modelPtr->initValue());               // .append()

    // Iteration stage : Path[i] to Path[i+1]
    for (size_t idx = 0; idx < _timeSteps.size() - 1; ++idx)  // .size() is +1 bigger, so - 1!
        path.push_back(nextStep(idx, path[idx]));

    return path;
}

bool PathSimulator::timeStepsSanityCheck() const
{
    return false;
}

double EulerPathSimulator::nextStep(size_t timeIndex, double assetPrice) const
{
    double deltaT = _timeSteps[timeIndex + 1] - _timeSteps[timeIndex];

    // Rnadom Number Generation
    normal_distribution<double> distribution(0.0, 1.0);
    double Z = distribution(_randomEngine);

    // Apply Euler scheme
    double nextStep = assetPrice + _modelPtr->drift(_timeSteps[timeIndex], assetPrice) * deltaT
        + _modelPtr->diffusion(_timeSteps[timeIndex], assetPrice) * sqrt(deltaT) * Z;

    return nextStep;
}
