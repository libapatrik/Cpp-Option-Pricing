//
// Created by Patrik  Liba on 18/09/2025.
//

#ifndef CPPFM_PRICER_H
#define CPPFM_PRICER_H

#include "FinancialInstrument.h" // to access the FinancialInstrument
#include "Model.h" // to access the Model
#include "PathSimulator.h" // to access the discretisation scheme

class Pricer
{
public:
    Pricer() = default;
    Pricer(const Model& model);
    virtual Pricer * clone() const = 0; // Virtual Pure method - clone

    virtual ~Pricer(); // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

    // Common interface for all pricers
    // TODO: What are the other params which I can point to instead?
    virtual double price(const EuropeanOptionPayoff& option,
                        double spot, // Set as default params
                        double volatility,
                        double r,
                        double T,
                        int numPaths) const; // Virtual Pure method

protected:
    const Model* _modelPtr; // Pointer to base class Model pure virtual cannot instantiate - take const pointer
};

class MonteCarloPricer : public Pricer
{
public:
    MonteCarloPricer() = default; // default constructor
    MonteCarloPricer(const Model& model); // constructor with parameter

    // Clone method
    MonteCarloPricer * clone() const override;
    // Destructor
    ~MonteCarloPricer() override;

    double price(const EuropeanOptionPayoff& option,
                        double spot,
                        double volatility,
                        double r,
                        double T,
                        int numPaths
                        ) const override;

protected:


};


#endif //CPPFM_PRICER_H