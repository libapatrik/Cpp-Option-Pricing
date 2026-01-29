#ifndef CPPFM_PATHSIMULATOR_H
#define CPPFM_PATHSIMULATOR_H

#include <cppfm/models/Model.h>
#include <random>
#include <vector>
#include <cppfm/third_party/pcg/pcg_random.hpp>


class PathSimulator
{
public:
    PathSimulator(const std::vector<double>& timeSteps, const Model& model, size_t randomSeed);

    virtual ~PathSimulator();
    // PathSimulator owns the memory pointed to by _modelPtr - it was allocated with new

    std::vector<double> path() const; // return the path simulated
    virtual double nextStep(size_t timeIndex, double assetPrice) const = 0; // VPM - because I don't know what it does
    // We delegate the implementation of nextStep to the derived classes
    
    // Virtual getter for polymorphic access
    // Other methods will want to use different time steps, so we need to make it virtual ?
    virtual const std::vector<double>& timeSteps() const 
    { 
        return _timeSteps; 
    }

protected:
    bool timeStepsSanityCheck() const;

    std::vector<double> _timeSteps; // information for discretisation; not necessarily equally spaced! {t_0, t_1, ..., t_n}
    const Model* _modelPtr; // Pointer to base class Model pure virtual cannot instantiate - take const pointer
    // mutable std::default_random_engine _randomEngine; // had to change to mutable because of the const method path()
    mutable pcg32 _randomEngine; // PCG32: faster, smaller state (16 bytes vs 2.5KB), better statistical quality
};

class EulerPathSimulator : public PathSimulator
{
public:
    EulerPathSimulator(const std::vector<double>& timeSteps, const Model& model, size_t randomSeed) // constructor with parameters
        : PathSimulator(timeSteps, model, randomSeed) 
        {
        }
    double nextStep(size_t timeIndex, double assetPrice) const override;
};

class MilsteinPathSimulator : public PathSimulator
{
public:
    MilsteinPathSimulator(const std::vector<double>& timeSteps, const Model& model, size_t randomSeed)
        : PathSimulator(timeSteps, model, randomSeed)
        {
        }
    double nextStep(size_t timeIndex, double assetPrice) const override;
};

#endif //CPPFM_PATHSIMULATOR_H;