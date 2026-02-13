#include <cmath>
#include <cppfm/calibration/LinearAlgebra.h>
#include <gtest/gtest.h> 
#include <iostream>

TEST(LinearAlgebraTest, MatrixOps_Multiply) {
    Matrix A        = {{1, 2}, {3, 4}};
    Matrix B        = {{5, 6}, {7, 8}};
    Matrix C        = MatrixOps::multiply(A, B);
    Matrix expected = {{19, 22}, {43, 50}};
    EXPECT_EQ(C, expected);
}

TEST(LinearAlgebraTest, MatrixOps_Transpose) {
    Matrix A        = {{1, 2}, {3, 4}};
    Matrix B        = MatrixOps::transpose(A);
    Matrix expected = {{1, 3}, {2, 4}};
    EXPECT_EQ(B, expected);
}

TEST(LinearAlgebraTest, MatrixOps_Norm2) {
    std::vector<double> x        = {1, 2, 3, 4};
    double              expected = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        expected += x[i] * x[i];
    }
    expected = std::sqrt(expected);

    double result = MatrixOps::norm2(x);
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(LinearAlgebraTest, MatrixOps_NormInf) {
    std::vector<double> x = {1, 2, 3, 4};
    // find elem in x with the largest absolute value (the infinity norm)
    // Using a lambda to compare absolute values.
    double expected = *std::max_element(x.begin(), x.end(), [](double a, double b) {
        return std::abs(a) < std::abs(b);
    });
    double result   = MatrixOps::normInf(x);
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(LinearAlgebraTest, MatrixOps_Dot) {
    std::vector<double> x = {1, 3, 5, 7, 9};
    std::vector<double> y = {2, 4, 6, 8, 10};

    double expected = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
        expected += x[i] * y[i];

    double result = MatrixOps::dot(x, y);
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(LinearAlgebraTest, MatrixOps_axpy) {  // update vector y by adding the alpha * x O(n) pass
    std::vector<double> x        = {1, 2, 3, 4, 5};
    std::vector<double> y        = {6, 7, 8, 9, 10};
    double              alpha    = 0.5;
    std::vector<double> expected = y;  // copy y to expected
    for (size_t i = 0; i < x.size(); ++i) {
        expected[i] += alpha * x[i];  // what y should become
    }
    MatrixOps::axpy(alpha, x, y);
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(y[i], expected[i], 1e-10);  // check if it updated and matches
    }
}

TEST(LinearAlgebraTest, MatrixOps_scale) {
    std::vector<double> x        = {1, 2, 3, 4, 5};
    double              alpha    = 2;
    std::vector<double> expected = x;
    for (size_t i = 0; i < x.size(); ++i) {
        expected[i] *= alpha;
    }
    MatrixOps::scale(alpha, x);
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], expected[i], 1e-10);
    }
}

TEST(LinearAlgebraTest, MatrixOps_identity) {
    size_t n = 7;
    // for given dim expand into indentity matrix
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                EXPECT_EQ(MatrixOps::identity(n)[i][j], 1.0);
            } else {
                EXPECT_EQ(MatrixOps::identity(n)[i][j], 0.0);
            }
        }
    }
}

TEST(LinearAlgebraTest, MatrixOps_zeros) {
    size_t m = 3;
    size_t n = 4;
    Matrix A = MatrixOps::zeros(m, n);  // matrix of m rows, n cols
    ASSERT_EQ(A.size(), m);
    for (size_t i = 0; i < m; ++i) {
        ASSERT_EQ(A[i].size(), n);
        for (size_t j = 0; j < n; ++j) {
            EXPECT_EQ(A[i][j], 0.0);
        }
    }
}

TEST(LinearAlgebraTest, MatrixOps_multiplyAtA) {
    Matrix A = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
    Matrix C = MatrixOps::multiplyAtA(A);

    Matrix expected =
        MatrixOps::zeros(A[0].size(), A[0].size());  // for A^T A the output is always square
    for (size_t i = 0; i < A[0].size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            for (size_t k = 0; k < A.size(); ++k) {
                expected[i][j] += A[k][i] * A[k][j];
            }
        }
    }

    ASSERT_EQ(C.size(), expected.size());
    for (size_t i = 0; i < C.size(); ++i) {
        ASSERT_EQ(C[i].size(), expected[i].size());
        for (size_t j = 0; j < C[i].size(); ++j) {
            EXPECT_NEAR(C[i][j], expected[i][j], 1e-10);
        }
    }
}

TEST(LinearAlgebraTest, MatrixOps_multiplyAtb) {
    Matrix              A = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
    std::vector<double> b = {1, 2, 3};
    ASSERT_EQ(A.size(), b.size());
    std::vector<double> y = MatrixOps::multiplyAtb(A, b);  // A^T * b

    // compute expected via explicit formula: y_i = sum_k A_{k,i} * b_k
    std::vector<double> expected(A[0].size(), 0.0);  // size = ncols of A
    for (size_t i = 0; i < A[0].size(); ++i) {
        for (size_t k = 0; k < b.size(); ++k) {
            expected[i] += A[k][i] * b[k];
        }
    }

    ASSERT_EQ(y.size(), expected.size());
    for (size_t i = 0; i < y.size(); ++i) {
        EXPECT_NEAR(y[i], expected[i], 1e-10);
    }
}

/*
    CHOLESKY
*/

TEST(LinearAlgebraTest, Cholesky_solve)  // solve S * x = b via L * L^T
{
    Matrix              S          = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    Matrix              L          = Cholesky::decompose(S);
    std::vector<double> x_expected = {1, 2, 3};  // known solution
    ASSERT_EQ(S.size(), x_expected.size());
    // compute b = A * x
    std::vector<double> b(S.size(), 0.0);
    for (size_t i = 0; i < S.size(); ++i) {
        for (size_t j = 0; j < S[i].size(); ++j) {
            b[i] += S[i][j] * x_expected[j];
        }
    }
    // solve and check
    std::vector<double> x = Cholesky::solve(L, b);

    ASSERT_EQ(x.size(), x_expected.size());
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], x_expected[i], 1e-10);
    }

    // residual check: S * x should equal b
    std::vector<double> Sx = MatrixOps::multiply(S, x);
    for (size_t i = 0; i < b.size(); ++i) {
        EXPECT_NEAR(Sx[i], b[i], 1e-10);
    }
    std::cout << "Cholesky residual = " << std::endl;
    for (size_t i = 0; i < Sx.size(); ++i) {
        std::cout << Sx[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for (size_t i = 0; i < b.size(); ++i) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
}

/*
   QR
*/
TEST(LinearAlgebraTest, QR_solve) {                            // solve the system b = A * x
    Matrix A = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};  // det = 36, non-singular
    // Matrix A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; // singular throws error, correct
    QRResult qr = QR::decompose(A);

    std::vector<double> x_expected = {1, 2, 3};
    std::vector<double> b          = MatrixOps::multiply(A, x_expected);  // b = Ax

    std::vector<double> x = QR::solve(qr, b);

    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], x_expected[i], 1e-10);
    }

    // residual check
    std::vector<double> residual = MatrixOps::multiply(A, x);
    for (size_t i = 0; i < residual.size(); ++i) {
        EXPECT_NEAR(residual[i], b[i], 1e-10);
    }
    std::cout << "QR residual = " << std::endl;
    for (size_t i = 0; i < residual.size(); ++i) {
        std::cout << residual[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for (size_t i = 0; i < b.size(); ++i) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
}

/*
    SVD
*/
TEST(LinearAlgebraTest, SVD_decompose) {
    // start with simple 2x2 to test basic algorithm
    Matrix    A   = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    SVDResult svd = SVD::decompose(A);

    std::vector<double> x_expected = {1, 2, 3};
    std::vector<double> b          = MatrixOps::multiply(A, x_expected);  // b = [3, 14]

    std::vector<double> x = SVD::solve(svd, b);

    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], x_expected[i], 1e-10);
    }

    // residual check
    std::vector<double> residual = MatrixOps::multiply(A, x);
    for (size_t i = 0; i < residual.size(); ++i) {
        EXPECT_NEAR(residual[i], b[i], 1e-10);
    }
    std::cout << "SVD residual = " << std::endl;
    for (size_t i = 0; i < residual.size(); ++i) {
        std::cout << residual[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for (size_t i = 0; i < b.size(); ++i) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
}
