#include "PathSimulator.h"

#include <cmath>
#include <vector>
#include <random>

using namespace std;

// GBM implementation
//TODO: Make it into one for loop - requires some vectorisation
vector<vector<double>> GBM_pathSimulator(const double& S0, const double& r, const double& sigma, const double& T, int numSteps, int numPaths) {
    // Seed
    default_random_engine generator(1);
    normal_distribution<double> distribution(0.0, 1.0); // mean 0, stddev 1

    double dt = T / static_cast<double>(numSteps);     // time step
    double drift = r - 0.5 * pow(sigma, 2) * dt;
    double vol_sqrt_dt = sigma * sqrt(dt);
    int numTotal = 2 * numPaths;                       // Antithetic samples
    // vector<vector<double>> paths(numTotal, vector<double>(numSteps + 1));

    // Want to do: np.zeros()
    vector<vector<double>> S(numTotal, vector<double>(numSteps + 1, S0));


    // Simulate paths
    for (int i = 0; i < numTotal; i++) {
        for (int j = 0; j <= numSteps; j++) {
            double Z = distribution(generator);

            // Paths first half
            S[i][j] = S[i][j] * exp(drift + vol_sqrt_dt * Z); // is there +?

            // Paths second half - antithetic
            S[i + numPaths][i + 1] = S[i + numPaths][j] * exp(drift + vol_sqrt_dt * (-Z));
        }
    }
    return S;
}
