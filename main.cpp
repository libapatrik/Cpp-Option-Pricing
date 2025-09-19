#include "Model.h"
#include "PathSimulator.h"

#include <iostream>

// void will be like a procedure -> set of instructions
int main()
{
    //BlackScholesModel model;
    // WANT TO TEST: if a copy assignment operator copies correct the data members from one BSM object to another
    BlackScholesModel bs1(100., 0.03, 0.2);  // constructor with parms

    BlackScholesModel bs2(bs1); // The object bs2 does not exist yet

    BlackScholesModel bs3(200., 0.05, 0.1);
    bs3 = bs1;                 // The object bs3 was already constructed -> returns bs3.operator=(bs1)

    BlackScholesModel bs4 = bs1; // Either it calls default constructor -> assigment operator
    // either it calls just the copy constructor


    // Create new object of type BlackScholesModel object on the heap with params
    // and stores its address in the pointer of type Model*. BSM is derived class of Model
    Model* modelPtr = new BlackScholesModel(100., 0.05, 0.2);
    // Call drift method on the object pointed to by modelPtr, pass params, store result in drift variable.
    // Drift is virtual in Model, so the BSM will be called [otherwise, call Base class method]
    double drift = modelPtr->drift(0.5, 120.);
    double volatility = modelPtr->diffusion(0.5, 120.);
    //double drift = (*modelPtr).drift(0.5, 120.);
    std::cout << "Drift: " << drift << std::endl; // 0.05 * 120 = 6
    std::cout << "Volatility: " << volatility << std::endl; // 0.2 * 120 = 2

    delete modelPtr; // delete the memory allocated on the heap when used `new`

    double S0 = 100.0, r = 0.05, sigma = 0.2, T = 1.0;
    int numSteps = 200;
    int numPaths = 500;
    auto paths = GBM_pathSimulator(S0, r, sigma, T, numSteps, numPaths);

    // Print all paths and all time steps - print the whole matrix S of paths from GBM_pathSimulator
    for (int i = 0; i < paths.size(); ++i) {
        std::cout << "Path " << i + 1 << ": ";
        for (int j = 0; j < paths[i].size(); ++j) {
            std::cout << paths[i][j];
            if (j != paths[i].size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }

    return 0;


    // Question is : which drift method is called : the base class one or the derived class one?
    // Virtuality
    // If base class method is declared virtual -> derived class
    // If base class method is not declared virtual -> base class


}