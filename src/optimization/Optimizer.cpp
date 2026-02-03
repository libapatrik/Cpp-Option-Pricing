#include <cppfm/optimization/Optimizer.h>
#include <iostream>

Matrix LevenbergMarquardt::numericalJacobian(ResidualFunc f, const std::vector<double>& x, double h)
{
    size_t n = x.size();
    std::vector<double> xp = x, xm = x;

    // first column to get m
    double h0 = std::max(h, h * std::abs(x[0]));
    xp[0] = x[0] + h0;
    xm[0] = x[0] - h0;
    auto rp = f(xp);
    auto rm = f(xm);
    size_t m = rp.size();

    Matrix J(m, std::vector<double>(n));
    for (size_t i = 0; i < m; ++i)
        J[i][0] = (rp[i] - rm[i]) / (2.0 * h0);
    xp[0] = x[0];
    xm[0] = x[0];

    // remaining columns
    for (size_t j = 1; j < n; ++j) {
        double hj = std::max(h, h * std::abs(x[j]));
        xp[j] = x[j] + hj;
        xm[j] = x[j] - hj;
        rp = f(xp);
        rm = f(xm);
        for (size_t i = 0; i < m; ++i)
            J[i][j] = (rp[i] - rm[i]) / (2.0 * hj);
        xp[j] = x[j];
        xm[j] = x[j];
    }
    return J;
}

std::vector<double> LevenbergMarquardt::clamp(const std::vector<double>& x,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
{
    std::vector<double> xc = x;
    clampInPlace(xc, lb, ub);
    return xc;
}

void LevenbergMarquardt::clampInPlace(std::vector<double>& x,
    const std::vector<double>& lb,
    const std::vector<double>& ub)
{
    for (size_t i = 0; i < x.size(); ++i) {
        if (!lb.empty()) x[i] = std::max(x[i], lb[i]);
        if (!ub.empty()) x[i] = std::min(x[i], ub[i]);
    }
}

LMResult LevenbergMarquardt::solve(ResidualFunc residuals,
    const std::vector<double>& x0,
    const std::vector<double>& lb,
    const std::vector<double>& ub,
    JacobianFunc jacobian,
    const LMOptions& opts)
{
    LMResult result;
    std::vector<double> x = clamp(x0, lb, ub);
    size_t n = x.size();

    auto r = residuals(x);
    double rNorm2 = MatrixOps::dot(r, r);

    double lambda = opts.lambda0;
    int iter = 0;

    // pre-allocate workspace to avoid per-iteration allocations
    std::vector<double> xNew(n);
    std::vector<double> negJtr(n);
    std::vector<double> delta(n);
    double xNorm = MatrixOps::norm2(x);

    while (iter < opts.maxIter) {
        Matrix J = jacobian ? jacobian(x) : numericalJacobian(residuals, x);

        // compute J^T*J and -J^T*r directly without forming J^T
        auto JtJ = MatrixOps::multiplyAtA(J);
        auto Jtr = MatrixOps::multiplyAtb(J, r);

        // (J^T*J + λI) * δ = -J^T*r
        for (size_t i = 0; i < n; ++i)
            JtJ[i][i] += lambda;

        // negate Jtr in place
        for (size_t i = 0; i < n; ++i)
            negJtr[i] = -Jtr[i];

        // try Cholesky first (faster), fall back to SVD if not positive definite
        bool usedSVD = false;
        try {
            auto L = Cholesky::decompose(JtJ);
            delta = Cholesky::solve(L, negJtr);
        } catch (const std::runtime_error&) {
            // Cholesky failed, fall back to SVD
            auto svd = SVD::decompose(JtJ);
            delta = SVD::solve(svd, negJtr);
            usedSVD = true;
        }

        // xNew = x + delta, then clamp in place
        for (size_t i = 0; i < n; ++i)
            xNew[i] = x[i] + delta[i];
        clampInPlace(xNew, lb, ub);

        auto rNew = residuals(xNew);
        double rNewNorm2 = MatrixOps::dot(rNew, rNew);

        if (rNewNorm2 < rNorm2) {
            x = xNew;
            r = rNew;
            rNorm2 = rNewNorm2;
            lambda *= opts.lambdaDown;
            ++iter;

            double stepNorm = MatrixOps::norm2(delta);
            double gradNorm = MatrixOps::norm2(Jtr);
            xNorm = MatrixOps::norm2(x);

            if (opts.verbose) {
                std::cout << "iter " << iter << ": ||r||^2 = " << rNorm2
                          << ", λ = " << lambda
                          << (usedSVD ? " (SVD)" : "") << "\n";
            }

            if (stepNorm < opts.tol * (1.0 + xNorm)) {
                result.converged = true;
                result.message = "converged: step size below tolerance";
                break;
            }
            if (gradNorm < opts.gradTol) {
                result.converged = true;
                result.message = "converged: gradient below tolerance";
                break;
            }
        } else {
            lambda *= opts.lambdaUp;
            if (lambda > 1e16) {
                result.converged = false;
                result.message = "failed: lambda overflow";
                break;
            }
        }

        if (iter >= opts.maxIter && !result.converged) {
            result.message = "stopped: max iterations reached";
        }
    }

    result.params = x;
    result.finalResidual = rNorm2;
    result.iterations = iter;
    return result;
}