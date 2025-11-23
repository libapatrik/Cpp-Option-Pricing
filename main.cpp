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

    // Finite Difference Pricer


    std::cout << "\nPricing American Options" << std::endl;
    // TODO: Add FD AM Pricer
    

    // ============================================================================
    // Test of InterpolationSchemes
    // ============================================================================
    // Demonstrates enum-based extrapolation design
    std::cout << "\nTesting InterpolationSchemes" << std::endl;
    std::vector<double> xData = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> yData = {2.0, 3.0, 4.0, 5.0, 6.0};
    
    std::cout << "Choose the extrapolation scheme: 1. Flat, 2. Linear, 3. Quadratic: ";
    int choice;
    std::cin >> choice;
    
    ExtrapolationType extraType;
    if (choice == 1) {
        extraType = ExtrapolationType::Flat;
    } else if (choice == 2) {
        extraType = ExtrapolationType::Linear;
    } else {
        // Default extrapolation
        extraType = ExtrapolationType::Quadratic;
        if (choice != 3) {
            std::cout << "Invalid choice. Using Quadratic (default)." << std::endl;
        }
    }
    
    auto interp = std::make_unique<LinearInterpolation>(xData, yData, extraType);
    // auto interp->std::make_unique<CubicSplineInterpolation>()
    // Outputs:
    std::cout << "\nResults:" << std::endl;
    std::cout << "Value at x=0.5 (extrapolation): " << (*interp)(0.5) << std::endl;
    std::cout << "Value at x=3.5 (interpolation): " << (*interp)(3.5) << std::endl;
    std::cout << "Value at x=6.0 (extrapolation): " << (*interp)(6.0) << std::endl;
    
    return 0;
    
    // Question is : which drift method is called : the base class one or the derived class one?
    // Virtuality
    // If base class method is declared virtual -> derived class
    // If base class method is not declared virtual -> base class

}