/**
    Optimizer
    - Levenberg-Marquardt
    - Trust region

*/

#ifndef CPPFM_OPTIMIZER_H
#define CPPFM_OPTIMIZER_H

#include <cppfm/optimization/LinearAlgebra.h>
#include <functional>
#include <string>
#include <vector>

struct LMResult {
    std::vector<double> params;
    double finalResidual;
    int iterations;
    bool converged;
    std::string message;
};

struct LMOptions {
    double tol = 1e-8;
    double gradTol = 1e-10;
    int maxIter = 100;
    double lambda0 = 1e-3;
    double lambdaUp = 10.0;
    double lambdaDown = 0.1;
    bool verbose = false;
};

using ResidualFunc = std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunc = std::function<Matrix(const std::vector<double>&)>;

class LevenbergMarquardt {
public:
    LMResult solve(ResidualFunc residuals, 
        const std::vector<double>& x0, 
        const std::vector<double>& lb = {}, 
        const std::vector<double>& ub = {}, 
        JacobianFunc jacobian = nullptr, 
        const LMOptions& opts = {});

private:
    // central differences for O(h2) accuracy
    Matrix numericalJacobian(ResidualFunc f, const std::vector<double>& x, double h = 1e-7);
    std::vector<double> clamp(const std::vector<double>& x,
        const std::vector<double>& lb,
        const std::vector<double>& ub);
    void clampInPlace(std::vector<double>& x,
        const std::vector<double>& lb,
        const std::vector<double>& ub);
};


#endif