//
// Created by Patrik  Liba on 18/09/2025.
//

#ifndef CPPFM_PRICER_H
#define CPPFM_PRICER_H

#include "FinancialInstrument.h" // to access the FinancialInstrument
#include "Model.h" // to access the Model
#include "PathSimulator.h" // to access the discretisation scheme
#include "DiscountCurve.h" // discounting

/* Notes
1. will those pointers point to the same model?

*/


class Pricer
{
public:
    Pricer() = default;
    Pricer(const Model& model, const DiscountCurve& discountCurve);
    virtual Pricer * clone() const = 0; // Virtual Pure method - clone

    virtual ~Pricer(); // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

    virtual double price(const EuropeanOptionPayoff& option) const = 0; // Virtual Pure method

protected:
    const Model* _modelPtr; // Pointer to base class Model pure virtual cannot instantiate - take const pointer
    const DiscountCurve* _discountCurvePtr;
    // later multifactor pricing
};

class BlackScholesPricer : public Pricer
{
public:
    BlackScholesPricer(); // default constructor
    BlackScholesPricer(const BlackScholesModel& model, const DiscountCurve& discountCurve); // give the discount curve to  constructor

    // Clone method
    BlackScholesPricer * clone() const override;
    // Destructor
    ~BlackScholesPricer() override;

    double price(const EuropeanOptionPayoff& option) const override;

protected:
    // Analytic formula
    static double d1(double spot, double strike, double r, double sigma, double T);
    static double d2(double d1, double sigma, double T);

    const BlackScholesModel* _bsModelPtr { nullptr };

};

class MonteCarloPricer : public Pricer
{
public:
    MonteCarloPricer() = default; // default constructor
    MonteCarloPricer(const Model &model, const DiscountCurve& discountCurve, const PathSimulator& simulator, size_t numPaths);
    // always keep the base class arguments first, then extras

    // Clone method
    MonteCarloPricer * clone() const override;

    // Destructor
    ~MonteCarloPricer() override;

    double price(const EuropeanOptionPayoff& option) const override;

protected:
    // Sanity check on address of model to be the same as the address in my pricer -> consistent

    // const PathSimulator* _pathSimulatorPtr; // inconsistency
    size_t _numPaths;
    const PathSimulator* _simulatorPtr { nullptr };

};


#endif //CPPFM_PRICER_H