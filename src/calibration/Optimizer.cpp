#include <cppfm/calibration/Optimizer.h>
#include <deque>
#include <iostream>

// ============================================================================
// backtracking line search (Armijo condition)
// ============================================================================

static double backtrack(ObjectiveFunc f, const std::vector<double> &x,
						const std::vector<double> &d, double fx,
						const std::vector<double> &g,
						double c1 = 1e-4, double shrink = 0.5)
{
	double alpha = 1.0;
	double slope = MatrixOps::dot(g, d);
	size_t n = x.size();
	std::vector<double> xTry(n);

	for (int i = 0; i < 60; ++i)
	{
		for (size_t j = 0; j < n; ++j)
			xTry[j] = x[j] + alpha * d[j];
		if (f(xTry) <= fx + c1 * alpha * slope)
			return alpha;
		alpha *= shrink;
	}
	return alpha;
}

// ============================================================================
// Gradient Descent
// ============================================================================

OptResult GradientDescent::solve(ObjectiveFunc f, GradientFunc grad,
								 const std::vector<double> &x0,
								 const GDOptions &opts)
{
	OptResult result;
	result.converged = false;
	size_t n = x0.size();
	std::vector<double> x = x0;
	double fx = f(x);

	for (int iter = 0; iter < opts.maxIter; ++iter)
	{
		auto g = grad(x);
		double gNorm = MatrixOps::norm2(g);

		if (gNorm < opts.gradTol)
		{
			result.converged = true;
			result.message = "converged: gradient below tolerance";
			result.iterations = iter;
			break;
		}

		// steepest descent direction
		std::vector<double> d(n);
		for (size_t i = 0; i < n; ++i)
			d[i] = -g[i];

		double alpha = backtrack(f, x, d, fx, g);

		// step
		double stepNorm2 = 0;
		for (size_t i = 0; i < n; ++i)
		{
			double si = alpha * d[i];
			x[i] += si;
			stepNorm2 += si * si;
		}
		fx = f(x);
		result.iterations = iter + 1;

		if (opts.verbose)
			std::cout << "GD iter " << result.iterations
					  << ": f = " << fx << ", ||g|| = " << gNorm << "\n";

		if (std::sqrt(stepNorm2) < opts.tol * (1.0 + MatrixOps::norm2(x)))
		{
			result.converged = true;
			result.message = "converged: step size below tolerance";
			break;
		}
	}

	result.params = x;
	result.finalValue = fx;
	if (!result.converged)
		result.message = "stopped: max iterations reached";
	return result;
}

// ============================================================================
// L-BFGS (two-loop recursion, Nocedal & Wright Algorithm 7.4)
// ============================================================================

OptResult LBFGS::solve(ObjectiveFunc f, GradientFunc grad,
					   const std::vector<double> &x0,
					   const LBFGSOptions &opts)
{
	OptResult result;
	result.converged = false;
	size_t n = x0.size();
	std::vector<double> x = x0;
	double fx = f(x);
	auto g = grad(x);

	// (s,y) pair storage
	std::deque<std::vector<double>> S, Y;
	std::deque<double> rho;

	for (int iter = 0; iter < opts.maxIter; ++iter)
	{
		double gNorm = MatrixOps::norm2(g);
		if (gNorm < opts.gradTol)
		{
			result.converged = true;
			result.message = "converged: gradient below tolerance";
			result.iterations = iter;
			break;
		}

		// --- two-loop recursion ---
		int m = static_cast<int>(S.size());
		std::vector<double> q = g;
		std::vector<double> alphas(m);

		// first loop: newest → oldest
		for (int i = m - 1; i >= 0; --i)
		{
			alphas[i] = rho[i] * MatrixOps::dot(S[i], q);
			MatrixOps::axpy(-alphas[i], Y[i], q);
		}

		// H0 scaling: γ = s^T y / y^T y
		double gamma = 1.0;
		if (m > 0)
			gamma = MatrixOps::dot(S.back(), Y.back()) /
					MatrixOps::dot(Y.back(), Y.back());

		std::vector<double> r = q;
		MatrixOps::scale(gamma, r);

		// second loop: oldest → newest
		for (int i = 0; i < m; ++i)
		{
			double beta = rho[i] * MatrixOps::dot(Y[i], r);
			MatrixOps::axpy(alphas[i] - beta, S[i], r);
		}

		// d = -H*g
		MatrixOps::scale(-1.0, r);

		// line search
		double alpha = backtrack(f, x, r, fx, g);

		// update x
		std::vector<double> s(n), xNew(n);
		for (size_t i = 0; i < n; ++i)
		{
			s[i] = alpha * r[i];
			xNew[i] = x[i] + s[i];
		}
		double fNew = f(xNew);
		auto gNew = grad(xNew);

		// curvature pair
		std::vector<double> y(n);
		for (size_t i = 0; i < n; ++i)
			y[i] = gNew[i] - g[i];

		double ys = MatrixOps::dot(y, s);
		if (ys > 1e-10)
		{
			S.push_back(s);
			Y.push_back(y);
			rho.push_back(1.0 / ys);
			if (static_cast<int>(S.size()) > opts.memory)
			{
				S.pop_front();
				Y.pop_front();
				rho.pop_front();
			}
		}

		double stepNorm = MatrixOps::norm2(s);
		double xNorm = MatrixOps::norm2(xNew);

		x = xNew;
		fx = fNew;
		g = gNew;
		result.iterations = iter + 1;

		if (opts.verbose)
			std::cout << "LBFGS iter " << result.iterations
					  << ": f = " << fx << ", ||g|| = " << MatrixOps::norm2(g) << "\n";

		if (stepNorm < opts.tol * (1.0 + xNorm))
		{
			result.converged = true;
			result.message = "converged: step size below tolerance";
			break;
		}
	}

	result.params = x;
	result.finalValue = fx;
	if (!result.converged)
		result.message = "stopped: max iterations reached";
	return result;
}

// ============================================================================
// Levenberg-Marquardt
// ============================================================================

Matrix LevenbergMarquardt::numericalJacobian(ResidualFunc f, const std::vector<double> &x, double h)
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
	for (size_t j = 1; j < n; ++j)
	{
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

std::vector<double> LevenbergMarquardt::clamp(const std::vector<double> &x,
											  const std::vector<double> &lb,
											  const std::vector<double> &ub)
{
	std::vector<double> xc = x;
	clampInPlace(xc, lb, ub);
	return xc;
}

void LevenbergMarquardt::clampInPlace(std::vector<double> &x,
									  const std::vector<double> &lb,
									  const std::vector<double> &ub)
{
	for (size_t i = 0; i < x.size(); ++i)
	{
		if (!lb.empty() && lb.size() == x.size())
			x[i] = std::max(x[i], lb[i]);
		if (!ub.empty() && ub.size() == x.size())
			x[i] = std::min(x[i], ub[i]);
	}
}

LMResult LevenbergMarquardt::solve(ResidualFunc residuals,
								   const std::vector<double> &x0,
								   const std::vector<double> &lb,
								   const std::vector<double> &ub,
								   JacobianFunc jacobian,
								   const LMOptions &opts)
{
	LMResult result;
	std::vector<double> x = clamp(x0, lb, ub);
	size_t n = x.size();

	auto r = residuals(x);
	double rNorm2 = MatrixOps::dot(r, r);

	double lambda = opts.lambda0;
	int iter = 0;
	int totalIter = 0;

	// pre-allocate workspace
	std::vector<double> xNew(n);
	std::vector<double> negJtr(n);
	std::vector<double> delta(n);
	std::vector<double> actualStep(n);
	double xNorm = MatrixOps::norm2(x);

	// compute initial Jacobian and normal equations
	Matrix J = jacobian ? jacobian(x) : numericalJacobian(residuals, x);
	auto JtJ = MatrixOps::multiplyAtA(J);
	auto Jtr = MatrixOps::multiplyAtb(J, r);

	while (iter < opts.maxIter && totalIter < 10 * opts.maxIter)
	{
		++totalIter;

		// (J^T*J + λI) * δ = -J^T*r
		// work on a copy so we can retry with different λ
		auto JtJ_aug = JtJ;
		for (size_t i = 0; i < n; ++i)
			JtJ_aug[i][i] += lambda;

		for (size_t i = 0; i < n; ++i)
			negJtr[i] = -Jtr[i];

		// Cholesky first, SVD fallback
		bool usedSVD = false;
		try
		{
			auto L = Cholesky::decompose(JtJ_aug);
			delta = Cholesky::solve(L, negJtr);
		}
		catch (const std::runtime_error &)
		{
			auto svd = SVD::decompose(JtJ_aug);
			delta = SVD::solve(svd, negJtr);
			usedSVD = true;
		}

		// trial point
		for (size_t i = 0; i < n; ++i)
			xNew[i] = x[i] + delta[i];
		clampInPlace(xNew, lb, ub);

		// actual step after clamping
		for (size_t i = 0; i < n; ++i)
			actualStep[i] = xNew[i] - x[i];

		auto rNew = residuals(xNew);
		double rNewNorm2 = MatrixOps::dot(rNew, rNew);

		if (rNewNorm2 < rNorm2)
		{
			// accept step
			x = xNew;
			r = rNew;
			rNorm2 = rNewNorm2;
			lambda *= opts.lambdaDown;
			++iter;

			double stepNorm = MatrixOps::norm2(actualStep);
			xNorm = MatrixOps::norm2(x);

			if (opts.verbose)
			{
				std::cout << "iter " << iter << ": ||r||^2 = " << rNorm2
						  << ", λ = " << lambda
						  << (usedSVD ? " (SVD)" : "") << "\n";
			}

			if (stepNorm < opts.tol * (1.0 + xNorm))
			{
				result.converged = true;
				result.message = "converged: step size below tolerance";
				break;
			}

			// recompute Jacobian at new x before gradient check
			J = jacobian ? jacobian(x) : numericalJacobian(residuals, x);
			JtJ = MatrixOps::multiplyAtA(J);
			Jtr = MatrixOps::multiplyAtb(J, r);

			double gradNorm = MatrixOps::norm2(Jtr);
			if (gradNorm < opts.gradTol)
			{
				result.converged = true;
				result.message = "converged: gradient below tolerance";
				break;
			}
		}
		else
		{
			// reject step, increase damping (J unchanged, skip recompute)
			lambda *= opts.lambdaUp;
			if (lambda > 1e16)
			{
				result.converged = false;
				result.message = "failed: lambda overflow";
				break;
			}
		}
	}

	if (!result.converged && result.message.empty())
		result.message = "stopped: max iterations reached";

	result.params = x;
	result.finalResidual = rNorm2;
	result.iterations = iter;
	return result;
}