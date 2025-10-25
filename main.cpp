#include "Model.h"
#include "PathSimulator.h"
#include "Pricer.h"
#include "FinancialInstrument.h"
#include "DiscountCurve.h"
#include "InterpolationSchemes.h"
#include "VolatilitySurface.h"

#include <iostream>
#include <vector>
#include <iomanip>

// void will be like a procedure -> set of instructions
int main()
{
    // Analytic Black-Scholes pricer
    FlatDiscountCurve flatCurve(0.05);
    BlackScholesModel modelForPricing(100.0, flatCurve, 0.2);
    BlackScholesPricer pricer(modelForPricing, flatCurve);
    EuropeanOptionPayoff callPayoff(Option::Type::Call, 100.0, 1.0);
    double callPrice = pricer.price(callPayoff);
    std::cout << "Call Price under BlackScholesModel: " << callPrice << std::endl;

    // Simulate paths on a daily grid; compare Euler vs Milstein
    size_t steps = 252;
    std::vector<double> timeSteps; timeSteps.reserve(steps + 1);
    for (size_t i = 0; i <= steps; ++i) timeSteps.push_back(1.0 * static_cast<double>(i) / static_cast<double>(steps));
    EulerPathSimulator eulerSim(timeSteps, modelForPricing, 1);
    MilsteinPathSimulator milsteinSim(timeSteps, modelForPricing, 1);

    MonteCarloPricer mcPricer(modelForPricing, flatCurve, eulerSim, 10000);
    double mcPrice = mcPricer.price(callPayoff);
    std::cout << "Call Price under MC (Euler):        " << mcPrice << std::endl;

    MonteCarloPricer mcMil(modelForPricing, flatCurve, milsteinSim, 10000);
    double mcMilPrice = mcMil.price(callPayoff);
    std::cout << "Call Price under MC (Milstein):     " << mcMilPrice << std::endl;


    
    

    return 0;
    // Question is : which drift method is called : the base class one or the derived class one?
    // Virtuality
    // If base class method is declared virtual -> derived class
    // If base class method is not declared virtual -> base class

}