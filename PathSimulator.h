#ifndef CPPFM_PATHSIMULATOR_H
#define CPPFM_PATHSIMULATOR_H

#include "Model.h"
#include <random>
#include <vector>
using namespace std;

// 2D vector storing the GBM paths
// vector<vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths);
class PathSimulator
{
public:
    PathSimulator(const vector<double>& timeSteps, const Model& model, size_t randomSeed);

    virtual ~PathSimulator();
    // PathSimulator owns the memory pointed to by _modelPtr - it was allocated with new

    vector<double> path() const; // return the path simulated
    virtual double nextStep(size_t timeIndex, double assetPrice) const = 0; // VPM - because I don't know what it does
    // We delegate the implementation of nextStep to the derived classes

protected:
    bool timeStepsSanityCheck() const;

    vector<double> _timeSteps; // information for discretisation; not necessarily equally spaced! {t_0, t_1, ..., t_n}
    const Model* _modelPtr; // Pointer to base class Model pure virtual cannot instantiate - take const pointer
    default_random_engine _randomEngine;
};

class EulerPathSimulator : public PathSimulator
{
public:
    double nextStep(size_t timeIndex, double assetPrice) const override;
};

#endif //CPPFM_PATHSIMULATOR_H;