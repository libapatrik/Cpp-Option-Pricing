#ifndef CPPFM_PATHSIMULATOR_H
#define CPPFM_PATHSIMULATOR_H

#include "Model.h"
#include <random>
#include <vector>

using namespace std;

class PathSimulator
{
public:
    PathSimulator(const vector<double>& timeSteps, const Model& model, size_t randomSeed);

    virtual ~PathSimulator();
    // PathSimulator owns the memory pointed to by _modelPtr - it was allocated with new

    std::vector<double> path() const; // return the path simulated
    virtual double nextStep(size_t timeIndex, double assetPrice) const = 0; // VPM - because I don't know what it does
    // We delegate the implementation of nextStep to the derived classes
    
    // Virtual getter for polymorphic access
    virtual const vector<double>& timeSteps() const 
    { 
        return _timeSteps; 
    }

protected:
    bool timeStepsSanityCheck() const;

    vector<double> _timeSteps; // information for discretisation; not necessarily equally spaced! {t_0, t_1, ..., t_n}
    const Model* _modelPtr; // Pointer to base class Model pure virtual cannot instantiate - take const pointer
    mutable default_random_engine _randomEngine; // had to change to mutable because of the const method path()
};

class EulerPathSimulator : public PathSimulator
{
public:
    EulerPathSimulator(const vector<double>& timeSteps, const Model& model, size_t randomSeed) // constructor with parameters
        : PathSimulator(timeSteps, model, randomSeed) 
        {
        }
    double nextStep(size_t timeIndex, double assetPrice) const override;
};

class MilsteinPathSimulator : public PathSimulator
{
public:
    MilsteinPathSimulator(const vector<double>& timeSteps, const Model& model, size_t randomSeed)
        : PathSimulator(timeSteps, model, randomSeed)
        {
        }
    double nextStep(size_t timeIndex, double assetPrice) const override;
};

#endif //CPPFM_PATHSIMULATOR_H;