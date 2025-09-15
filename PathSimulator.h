#ifndef CPPFM_PATHSIMULATOR_H
#define CPPFM_PATHSIMULATOR_H

#include <cmath>
#include <vector>
#include <random>

using namespace std;

// GBM
//TODO: Make it into one for loop - requires some vectorisation
vector<vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths) {
    // Seed
    default_random_engine generator(1);
    normal_distribution<double> distribution(0.0, 1.0); // mean 0, stddev 1

    double dt = T / static_cast<double>(numSteps); // time step
    double drift = r - 0.5 * pow(sigma, 2) * dt;
    double vol_sqrt_dt = sigma * sqrt(dt);


    // Antithetic samples
    int numTotal = 2 * numPaths;
    // vector<vector<double>> paths(numTotal, vector<double>(numSteps + 1));
    vector<vector<double>> S(numTotal, vector<double>(numSteps + 1, S0)); // removes the for loop below

    // // Initialise all paths with S0
    // // Want to do: np.zeros([numSteps, numPaths])
    // for (int i = 0; i < numTotal; i++) {
    //     S[i][0] = S0;
    // } // We can do better right?


    // Simulate paths
    for (int i = 0; i < numTotal; i++) {
        for (int j = 0; j <= numSteps; j++) {
            double Z = distribution(generator);

            // Paths first half
            S[i][j] = S[i][j] * exp(drift + vol_sqrt_dt * Z); // is there +?

            // Paths second half - antithetic
            S[i + numPaths][i + 1] = S[i + numPaths][j] * exp(drift + vol_sqrt_dt * Z);
        }
    }
    return S;
}



#endif //CPPFM_PATHSIMULATOR_H;