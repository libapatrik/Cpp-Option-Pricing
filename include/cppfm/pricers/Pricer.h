#ifndef CPPFM_PRICER_H
#define CPPFM_PRICER_H

#include <cppfm/market/FinancialInstrument.h> // to access the FinancialInstrument
#include <cppfm/models/Model.h> // to access the Model
#include <cppfm/simulators/PathSimulator.h> // to access the discretisation scheme
#include <cppfm/market/DiscountCurve.h> // discounting
#include <cppfm/pde/Grid.h> // to access the Grid
#include <cppfm/pde/Solver.h> // to access the Solver


class Pricer
{
public:
    Pricer() = default;
    Pricer(const Model1D& model, const DiscountCurve& discountCurve);
    virtual Pricer * clone() const = 0; // Virtual Pure method - clone

    virtual ~Pricer(); // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

    virtual double price(const EuropeanOptionPayoff& option) const = 0; // Virtual Pure method

protected:
    const Model1D* _modelPtr; // Pointer to base class Model1D pure virtual cannot instantiate - take const pointer
    const DiscountCurve* _discountCurvePtr; // Pointer to base class DiscountCurve
    // later multifactor pricing
};

class BlackScholesPricer : public Pricer
{
public:
    BlackScholesPricer(); // default constructor
    BlackScholesPricer(const BlackScholesModel& model, const DiscountCurve& discountCurve);
    BlackScholesPricer(const BlackScholesPricer& other); // deep copy — Pricer owns raw ptrs

    BlackScholesPricer * clone() const override;
    // Destructor
    ~BlackScholesPricer() override;

    double price(const EuropeanOptionPayoff& option) const override;

protected:
    const BlackScholesModel* _bsModelPtr { nullptr };

};

class MonteCarloPricer : public Pricer
{
public:
    MonteCarloPricer() = default;
    MonteCarloPricer(const Model1D &model, const DiscountCurve& discountCurve, const PathSimulator& simulator, size_t numPaths);
    MonteCarloPricer(const MonteCarloPricer& other); // deep copy

    MonteCarloPricer * clone() const override;

    ~MonteCarloPricer() override; // Destructor

    double price(const EuropeanOptionPayoff& option) const override;

protected:
    // Sanity check on address of model to be the same as the address in my pricer -> consistent

    // const PathSimulator* _pathSimulatorPtr; // inconsistency
    size_t _numPaths;
    const PathSimulator* _simulatorPtr { nullptr };

};


class FDPricer : public Pricer
{
public: 
    /**
     * Finite Difference Pricer for European options
     * @param model Pricing model (BlackScholesModel, DupireModel, etc.)
     * @param discountCurve Discount curve for risk-free rate
     * @param S_min Minimum spot price for grid
     * @param S_max Maximum spot price for grid
     * @param N_S Number of spatial grid points
     * @param N_t Number of time steps (per unit time)
     * Note: Time horizon is determined by each option's maturity
     */
    FDPricer(const Model1D& model, const DiscountCurve& discountCurve,
             double S_min, double S_max, size_t N_S, size_t N_t);
    
    FDPricer(const FDPricer& other); // copy constructor, because we make use of unique pointers; which cannot be copied
    
    FDPricer* clone() const override;
    ~FDPricer() override = default;

    double price(const EuropeanOptionPayoff& option) const override;
    
    // Greeks as public
    double delta(const EuropeanOptionPayoff& option) const;
    double gamma(const EuropeanOptionPayoff& option) const;

protected:
    std::unique_ptr<Grid> _grid;

    std::unique_ptr<ThetaMethodSolver> createSolver(const EuropeanOptionPayoff& option) const;

};

#endif //CPPFM_PRICER_H
