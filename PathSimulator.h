#ifndef CPPFM_PATHSIMULATOR_H
#define CPPFM_PATHSIMULATOR_H

#include <vector>
using namespace std;

// 2D vector storing the GBM paths
vector<vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths);


#endif //CPPFM_PATHSIMULATOR_H;