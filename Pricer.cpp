
#include <cmath>
#include <algorithm>

#include "Pricer.h"
#include "DiscountCurve.h"
#include "BlackScholesFormulas.h"
#include "FinancialInstrument.h"
#include "PDEs/PDE.h"
#include "PDEs/BoundaryConditions.h"





Pricer::Pricer(const Model1D &model, const DiscountCurve &discountCurve)
    :   _modelPtr(model.clone()), _discountCurvePtr(discountCurve.clone()) // add to memory
{
}

Pricer::~Pricer()
{   // delete the memory allocated with close
    if (_modelPtr) {  // Check if Null before deleting pointer
        delete _modelPtr;
        _modelPtr = nullptr;  // Prevent double deletion; ensuring that the pointer does not point to deallocated member location
    }
    if (_discountCurvePtr) {  // Check if Null before deleting pointer
        delete _discountCurvePtr;
        _discountCurvePtr = nullptr;  // Prevent double deletion
    }
}

BlackScholesPricer::BlackScholesPricer() = default;

BlackScholesPricer::BlackScholesPricer(const BlackScholesModel &model, const DiscountCurve &discountCurve)
    : Pricer(model, discountCurve)
{
    _bsModelPtr = dynamic_cast<const BlackScholesModel*>(_modelPtr); // set it once
    if (!_bsModelPtr)
        throw std::runtime_error("BlackScholesPricer requires BlackScholesModel");
}

BlackScholesPricer::~BlackScholesPricer() = default;
BlackScholesPricer * BlackScholesPricer::clone() const
{
    return new BlackScholesPricer(*this);
}

double BlackScholesPricer::price(const EuropeanOptionPayoff &option) const
{
    if (_bsModelPtr == nullptr)
        throw std::runtime_error("BlackScholesPricer requires BlackScholesModel");

    return BlackScholesFormulas::price(
        _bsModelPtr->initValue(), 
        option.strike(), 
        *_discountCurvePtr,
        _bsModelPtr->volatility(), 
        option.maturity(), 
        option.type()
    );
}


MonteCarloPricer::MonteCarloPricer(const Model1D& model, const DiscountCurve& discountCurve, const PathSimulator& simulator, size_t numPaths) // with params
    : Pricer(model, discountCurve),  _numPaths(numPaths), _simulatorPtr(&simulator) // explicitly passing the model to base class constructor
{  // NOTE: PathSimulator may need clone method for safety
}

MonteCarloPricer* MonteCarloPricer::clone() const
{
    return new MonteCarloPricer(*this); // destroyed by Pricer destructor!
}

MonteCarloPricer::~MonteCarloPricer() = default;

double MonteCarloPricer::price(const EuropeanOptionPayoff &option) const
{   // TODO: Antithetic sampling
    // Simulate paths, sum the payoffs, take average, discount
    double T = option.maturity();
    double r = -std::log(_discountCurvePtr->discount(T)) / T;

    size_t numPaths = _numPaths > 0 ? _numPaths : 10000;
    double sumPayoff = 0.0;
    for (size_t p = 0; p < numPaths; ++p)
    {
        std::vector<double> path = _simulatorPtr->path();
        double ST = path.back(); // last element of the path
        // Compute payoff at maturity
        double payoff = (option.type() == Option::Type::Call) // condition ? value_if_true : value_if_false
            ? std::max(ST - option.strike(), 0.0)
            : std::max(option.strike() - ST, 0.0);
        sumPayoff += payoff; // sum of payoffs
    }

    double avgPayoff = sumPayoff / static_cast<double>(numPaths); // average payoff
    return std::exp(-r * T) * avgPayoff; // discount 
}

FDPricer::FDPricer(const Model1D &model, const DiscountCurve &discountCurve, double S_min, double S_max, size_t N_S, size_t N_t)
    : Pricer(model, discountCurve),
    _grid(std::make_unique<Grid>(S_min, S_max, N_S, 1.0, N_t))  // Placeholder time horizon, overridden per option
{
    // Note: FDPricer works with any Model1D that implements drift() and diffusion()
    // For DupireModel, this will use local volatility surface via localVolatility()
    // For BlackScholesModel, this will use constant volatility
    // Time horizon T is determined by each option's maturity in createSolver()
}

FDPricer::FDPricer(const FDPricer &other)
    : Pricer(*other._modelPtr, *other._discountCurvePtr), 
      _grid(std::unique_ptr<Grid>(other._grid->clone()))
{
}

FDPricer * FDPricer::clone() const
{
    return new FDPricer(*this);
}


/**
 * Create PDE solver for European option pricing
 * 
 * TIME CONVENTION:
 * ================
 * Solver integrates backward in time from maturity to today:
 * - Solver time t ∈ [0, T]: t=0 is maturity, t=T is today
 * - Calendar time from today: τ = T - t
 * - Time remaining to maturity at solver time t: t itself
 * 
 * Model drift(τ,S) and diffusion(τ,S) expect calendar time τ from today,
 * so we transform: τ = T - t before calling model methods.
 * 
 * PDE SOLVED:
 * ===========
 * ∂V/∂t = (1/2)σ²(S,τ)S² ∂²V/∂S² + r(τ)S ∂V/∂S - r(τ)V
 * where σ(S,τ) comes from model.diffusion(τ,S) = σ(S,τ)*S
 */
std::unique_ptr<ThetaMethodSolver> FDPricer::createSolver(const EuropeanOptionPayoff &option) const {
    double strike = option.strike();
    double T = option.maturity();
    Option::Type type = option.type();
    auto a_fn = [this, T](double S, double t) {
        double tau = T - t;  // Transform solver time to time from today
        double diff = _modelPtr->diffusion(tau, S);  // Returns σ(S,t)*S
        return 0.5 * diff * diff;  // 0.5 * σ²(S,t) * S²
    };

    auto b_fn = [this, T](double S, double t){
        double tau = T - t;  // Transform solver time to time from today
        return _modelPtr->drift(tau, S);  // Returns r(t)*S
    };

    auto c_fn = [this, T](double S, double t){
        double tau = T - t;  // Transform solver time to time from today
        double r_t = _discountCurvePtr->instantaneousRate(tau);
        return -r_t;
    };

    auto payoff = [strike, type](double S){
        return type == Option::Type::Call
            ? std::max(S - strike, 0.0)
            : std::max(strike - S, 0.0);
    };
    // Make the PDE
    auto pde = std::make_unique<VariableCoefficientPDE>(a_fn, b_fn, c_fn,
                                                        [](double, double) { return 0.0; },
                                                        payoff);

    double S_max = _grid->spotMax();

    // Boundary conditions: Solver time t goes from 0 (maturity) to T (today)
    // At solver time t, we are at calendar time T-t from today, so remaining time to maturity is t
    auto lower_bc = [strike, type, T, this](double t) {
        // At S=0: Put worth K*exp(-r*t), Call worth 0
        return type == Option::Type::Call
            ? 0.0
            : strike * _discountCurvePtr->discount(t);
    };
    auto upper_bc = [strike, type, T, S_max, this](double t) {
        // At S→∞: Call worth S-K*exp(-r*t), Put worth 0
        return type == Option::Type::Call
            ? S_max - strike * _discountCurvePtr->discount(t)
            : 0.0;
    };

    auto bc = std::make_unique<DirichletBC>(lower_bc, upper_bc);
    
    // Create grid with correct maturity T for the specific option
    Grid grid(_grid->spotMin(), _grid->spotMax(), _grid->numSpotPoints(), T, _grid->numTimeSteps(), _grid->spacingType());
    
    return std::make_unique<ThetaMethodSolver>(*pde, grid, *bc, ThetaMethodSolver::Scheme::CrankNicolson);

}

double FDPricer::price(const EuropeanOptionPayoff &option) const {
    auto solver = createSolver(option);
    solver->solve();
    return solver->valueAt(_modelPtr->initValue());
}

double FDPricer::delta(const EuropeanOptionPayoff &option) const {
    auto solver = createSolver(option);
    solver->solve();
    return solver->derivativeAt(_modelPtr->initValue());
}

double FDPricer::gamma(const EuropeanOptionPayoff &option) const {
    auto solver = createSolver(option);
    solver->solve();
    return solver->secondDerivativeAt(_modelPtr->initValue());
}


// double FDPricer
