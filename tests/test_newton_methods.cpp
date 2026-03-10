#include <gtest/gtest.h>
#include <cppfm/cos/COS.h>
#include <cppfm/utils/Utils.h>

// test CDF inversion accuracy against known normal quantiles
TEST(NewtonMethodsTest, InvertCDFAccuracy) {
    double a = -10.0;
    double b = 10.0;
    size_t N = 128;
    double tol = 1e-8;
    size_t max_iter = 100;

    std::vector<double> test_probs = {0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99};

    for (double p : test_probs) {
        auto [quantile, iters] = Transforms::invertCDF(a, b, N, Utils::stdNormChF, p, max_iter, tol);
        double expected = Utils::inverseNormalCDF(p);

        EXPECT_NEAR(quantile, expected, 1e-5)
            << "CDF inversion failed at p=" << p;
        EXPECT_LT(iters, max_iter)
            << "Newton did not converge at p=" << p;
    }
}


