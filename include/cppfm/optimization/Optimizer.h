/**
    Optimizer
    - Gradient Descent (backtracking Armijo)
    - L-BFGS (Nocedal & Wright Ch. 7)
    - Levenberg-Marquardt
*/

#ifndef CPPFM_OPTIMIZER_H
#define CPPFM_OPTIMIZER_H

#include <cppfm/optimization/LinearAlgebra.h>
#include <functional>
#include <string>
#include <vector>

// ---- function types ----

using ResidualFunc = std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunc = std::function<Matrix(const std::vector<double>&)>;
using ObjectiveFunc = std::function<double(const std::vector<double>&)>;
using GradientFunc = std::function<std::vector<double>(const std::vector<double>&)>;

// ---- result structs ----

struct OptResult {
    std::vector<double> params;
    double finalValue = 0.0;
    int iterations = 0;
    bool converged = false;
    std::string message;
};

struct LMResult {
    std::vector<double> params;
    double finalResidual = 0.0;
    int iterations = 0;
    bool converged = false;
    std::string message;
};

// ---- options ----

struct GDOptions {
    double tol = 1e-8;
    double gradTol = 1e-10;
    int maxIter = 10000;
    bool verbose = false;
};

struct LBFGSOptions {
    double tol = 1e-8;
    double gradTol = 1e-10;
    int maxIter = 1000;
    int memory = 10;
    bool verbose = false;
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

// ---- Gradient Descent ----

class GradientDescent {
public:
    OptResult solve(ObjectiveFunc f, GradientFunc grad,
        const std::vector<double>& x0,
        const GDOptions& opts = {});
};

// ---- L-BFGS ----

class LBFGS {
public:
    OptResult solve(ObjectiveFunc f, GradientFunc grad,
        const std::vector<double>& x0,
        const LBFGSOptions& opts = {});
};

// ---- Levenberg-Marquardt ----

class LevenbergMarquardt {
public:
    LMResult solve(ResidualFunc residuals,
        const std::vector<double>& x0,
        const std::vector<double>& lb = {},
        const std::vector<double>& ub = {},
        JacobianFunc jacobian = nullptr,
        const LMOptions& opts = {});

private:
    Matrix numericalJacobian(ResidualFunc f, const std::vector<double>& x, double h = 1e-7);
    std::vector<double> clamp(const std::vector<double>& x,
        const std::vector<double>& lb,
        const std::vector<double>& ub);
    void clampInPlace(std::vector<double>& x,
        const std::vector<double>& lb,
        const std::vector<double>& ub);
};

#endif