#include <cmath>
#include <cppfm/calibration/Optimizer.h>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>

// ---------------------------------------------------------------
// basic LM tests
// ---------------------------------------------------------------

TEST(OptimizerTest, LM_SimplePointFinding)
{
	// r = [x-3, y-4], solution at (3,4)
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {x[0] - 3.0, x[1] - 4.0}; };

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {0.0, 0.0});

	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 3.0, 1e-6);
	EXPECT_NEAR(result.params[1], 4.0, 1e-6);
	EXPECT_LT(result.finalResidual, 1e-12);
}

TEST(OptimizerTest, LM_Rosenbrock)
{
	// r = [10(y-x^2), 1-x], minimum at (1,1)
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]}; };

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {-1.0, 1.0});

	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 1.0, 1e-6);
	EXPECT_NEAR(result.params[1], 1.0, 1e-6);
	EXPECT_LT(result.finalResidual, 1e-12);
}

TEST(OptimizerTest, LM_OverdeterminedExact)
{
	// x+y=1, x-y=1, 2x=2 → (1,0)
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {x[0] + x[1] - 1.0, x[0] - x[1] - 1.0, 2.0 * x[0] - 2.0}; };

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {0.0, 0.0});

	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 1.0, 1e-6);
	EXPECT_NEAR(result.params[1], 0.0, 1e-6);
}

TEST(OptimizerTest, LM_AnalyticalJacobian)
{
	// Rosenbrock with explicit J
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]}; };
	auto jacobian = [](const std::vector<double> &x) -> Matrix
	{ return {{-20.0 * x[0], 10.0}, {-1.0, 0.0}}; };

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {-1.0, 1.0}, {}, {}, jacobian);

	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 1.0, 1e-6);
	EXPECT_NEAR(result.params[1], 1.0, 1e-6);
}

TEST(OptimizerTest, LM_UnderdeterminedCircle)
{
	// x^2+y^2=1 — infinite solutions, LM should find *some* point on the circle
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {x[0] * x[0] + x[1] * x[1] - 1.0}; };

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {2.0, 1.0});

	EXPECT_TRUE(result.converged);
	double r2 = result.params[0] * result.params[0] +
				result.params[1] * result.params[1];
	EXPECT_NEAR(r2, 1.0, 1e-6);
	EXPECT_LT(result.finalResidual, 1e-12);
}

TEST(OptimizerTest, LM_BoundedOptimization)
{
	// unconstrained solution at (3,4), bounds [0,2]x[0,2]
	// simple clamping gets to (2,2) but can't formally converge there —
	// every step overshoots past bounds, gets clamped back, no improvement
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {x[0] - 3.0, x[1] - 4.0}; };

	std::vector<double> lb = {0.0, 0.0};
	std::vector<double> ub = {2.0, 2.0};

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {0.5, 0.5}, lb, ub);

	// reaches the constrained optimum even though convergence flag won't
	// trigger
	EXPECT_NEAR(result.params[0], 2.0, 1e-6);
	EXPECT_NEAR(result.params[1], 2.0, 1e-6);
}

TEST(OptimizerTest, LM_CustomOptions)
{
	auto residuals = [](const std::vector<double> &x) -> std::vector<double>
	{ return {x[0] - 1.0, x[1] - 2.0}; };

	LMOptions opts;
	opts.tol = 1e-8;
	opts.maxIter = 200;

	LevenbergMarquardt lm;
	auto result = lm.solve(residuals, {10.0, 10.0}, {}, {}, nullptr, opts);

	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 1.0, 1e-6);
	EXPECT_NEAR(result.params[1], 2.0, 1e-6);
	EXPECT_LT(result.finalResidual, 1e-12);
}

// ============================================================================
// Himmelblau test
// ============================================================================

// ---------------------------------------------------------------
// Himmelblau convergence comparison: GD vs L-BFGS vs LM
//
// f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
//
// four global minima (all f=0):
//   (3.0, 2.0), (-2.805, 3.131), (-3.779, -3.283), (3.584, -1.848)
// ---------------------------------------------------------------

namespace
{

struct Min2D
{
	double x, y;
	const char *label;
};

const Min2D knownMinima[] = {{3.0, 2.0, "( 3.00,  2.00)"},
							 {-2.805118, 3.131312, "(-2.81,  3.13)"},
							 {-3.779310, -3.283186, "(-3.78, -3.28)"},
							 {3.584428, -1.848126, "( 3.58, -1.85)"}};

std::string identifyMin(const std::vector<double> &p)
{
	for (auto &m : knownMinima)
	{
		if (std::hypot(p[0] - m.x, p[1] - m.y) < 0.1)
			return m.label;
	}
	return "  ???";
}

// Himmelblau objective
double himmelblau(const std::vector<double> &x)
{
	double a = x[0] * x[0] + x[1] - 11.0;
	double b = x[0] + x[1] * x[1] - 7.0;
	return a * a + b * b;
}

// gradient
std::vector<double> himmelblauGrad(const std::vector<double> &x)
{
	double a = x[0] * x[0] + x[1] - 11.0;
	double b = x[0] + x[1] * x[1] - 7.0;
	return {4.0 * x[0] * a + 2.0 * b, 2.0 * a + 4.0 * x[1] * b};
}

// residuals for LM: r = [x^2+y-11, x+y^2-7]
std::vector<double> himmelblauResiduals(const std::vector<double> &x)
{
	return {x[0] * x[0] + x[1] - 11.0, x[0] + x[1] * x[1] - 7.0};
}

} // namespace

TEST(HimmelblauTest, ConvergenceComparison)
{
	// starting points that probe different basins of attraction
	std::vector<std::vector<double>> starts = {
		{0.0, 0.0},	 {1.0, 1.0}, {-4.0, 3.0},  {-4.0, -4.0},
		{4.0, -1.0}, {5.0, 5.0}, {-1.0, -1.0}, {2.0, -3.0},
	};

	GradientDescent gd;
	LBFGS lbfgs;
	LevenbergMarquardt lm;

	// header
	std::cout << "\n  Himmelblau convergence: GD vs L-BFGS vs LM\n";
	std::cout << std::string(90, '-') << "\n";
	std::cout << std::setw(14) << "start"
			  << " | " << std::setw(22) << "Gradient Descent"
			  << " | " << std::setw(22) << "L-BFGS"
			  << " | " << std::setw(22) << "LM" << "\n";
	std::cout << std::string(90, '-') << "\n";

	for (auto &x0 : starts)
	{
		auto gdR = gd.solve(himmelblau, himmelblauGrad, x0);
		auto lbR = lbfgs.solve(himmelblau, himmelblauGrad, x0);
		auto lmR = lm.solve(himmelblauResiduals, x0);

		// all should reach f ≈ 0
		EXPECT_TRUE(gdR.converged)
			<< "GD failed from (" << x0[0] << "," << x0[1] << ")";
		EXPECT_TRUE(lbR.converged)
			<< "LBFGS failed from (" << x0[0] << "," << x0[1] << ")";
		EXPECT_TRUE(lmR.converged)
			<< "LM failed from (" << x0[0] << "," << x0[1] << ")";
		EXPECT_LT(gdR.finalValue, 1e-6);
		EXPECT_LT(lbR.finalValue, 1e-6);
		EXPECT_LT(lmR.finalResidual, 1e-6);

		// table row
		std::cout << std::fixed << std::setprecision(1);
		std::cout << "(" << std::setw(5) << x0[0] << "," << std::setw(5)
				  << x0[1] << ")"
				  << " | " << identifyMin(gdR.params) << " " << std::setw(5)
				  << gdR.iterations << "it"
				  << " | " << identifyMin(lbR.params) << " " << std::setw(5)
				  << lbR.iterations << "it"
				  << " | " << identifyMin(lmR.params) << " " << std::setw(5)
				  << lmR.iterations << "it"
				  << "\n";
	}
	std::cout << std::string(90, '-') << "\n\n";
}

// same starting point, all three methods — verify they land on a valid minimum
TEST(HimmelblauTest, AllMethodsReachMinimum)
{
	std::vector<double> x0 = {1.0, 1.0};

	GradientDescent gd;
	auto gdR = gd.solve(himmelblau, himmelblauGrad, x0);
	EXPECT_TRUE(gdR.converged);
	EXPECT_NEAR(himmelblau(gdR.params), 0.0, 1e-8);

	LBFGS lbfgs;
	auto lbR = lbfgs.solve(himmelblau, himmelblauGrad, x0);
	EXPECT_TRUE(lbR.converged);
	EXPECT_NEAR(himmelblau(lbR.params), 0.0, 1e-8);

	LevenbergMarquardt lm;
	auto lmR = lm.solve(himmelblauResiduals, x0);
	EXPECT_TRUE(lmR.converged);
	EXPECT_NEAR(lmR.finalResidual, 0.0, 1e-8);
}

int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
