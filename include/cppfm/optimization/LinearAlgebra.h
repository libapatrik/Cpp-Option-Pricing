#ifndef CPPFM_LINEARALGEBRA_H
#define CPPFM_LINEARALGEBRA_H

#include <vector>
#include <cmath>
#include <stdexcept>

// Matrix = row-major 2D vector
using Matrix = std::vector<std::vector<double>>;

// Basic matrix/vector operations
class MatrixOps
{
public:
    // y = A * x
    static std::vector<double> multiply(const Matrix &A, const std::vector<double> &x);
    // C = A * B
    static Matrix multiply(const Matrix &A, const Matrix &B);
    // A^T
    static Matrix transpose(const Matrix &A);
    // ||x||_2
    static double norm2(const std::vector<double> &x);
    // ||x||_inf
    static double normInf(const std::vector<double> &x);
    // dot product
    static double dot(const std::vector<double> &x, const std::vector<double> &y);
    // y = alpha * x + y
    static void axpy(double alpha, const std::vector<double> &x, std::vector<double> &y);
    // x = alpha * x
    static void scale(double alpha, std::vector<double> &x);
    // identity matrix
    static Matrix identity(size_t n);
    // zero matrix
    static Matrix zeros(size_t m, size_t n);
    // A^T * A (symmetric, avoids forming A^T)
    static Matrix multiplyAtA(const Matrix &A);
    // A^T * b (avoids forming A^T)
    static std::vector<double> multiplyAtb(const Matrix &A, const std::vector<double> &b);

private:
    MatrixOps() = delete;
};

// Cholesky: S = L * L^T for symmetric positive definite S
class Cholesky
{
public:
    // returns lower triangular L
    static Matrix decompose(const Matrix &S);
    // solve L * y = b (forward substitution)
    static std::vector<double> solveL(const Matrix &L, const std::vector<double> &b);
    // solve L^T * x = y (back substitution)
    static std::vector<double> solveLT(const Matrix &L, const std::vector<double> &y);
    // solve S * x = b via L * L^T
    static std::vector<double> solve(const Matrix &L, const std::vector<double> &b);

private:
    Cholesky() = delete;
};

// QR decomposition via Householder reflections
// Implicit Q storage: Householder vectors in lower triangle of R, scalars in tau
struct QRResult
{
    Matrix R;                  // upper tri + Householder vectors below diagonal
    std::vector<double> tau;   // Householder scalars
};

class QR
{
public:
    static QRResult decompose(const Matrix &A);
    // solve min ||Ax - b||_2 via QR
    static std::vector<double> solve(const QRResult &qr, const std::vector<double> &b);
    // solve R * x = c (back substitution)
    static std::vector<double> solveR(const Matrix &R, const std::vector<double> &c);
    // apply Q^T to vector without forming Q
    static std::vector<double> applyQT(const QRResult &qr, const std::vector<double> &b);
    // apply Q to vector (for residual computation)
    static std::vector<double> applyQ(const QRResult &qr, const std::vector<double> &x);

private:
    QR() = delete;
};

// SVD: A = U * diag(sigma) * V^T
struct SVDResult
{
    Matrix U;                  // m x min(m,n)
    std::vector<double> sigma; // singular values
    Matrix V;                  // n x min(m,n)
};

class SVD
{
public:
    // Golub-Reinsch algorithm
    static SVDResult decompose(const Matrix &A);
    // solve min ||Ax - b||_2, truncating small singular values
    static std::vector<double> solve(const SVDResult &svd, const std::vector<double> &b, double tol = 1e-10);

private:
    // bidiagonalize A = U1 * B * V1^T
    static void bidiagonalize(const Matrix &A, Matrix &U, std::vector<double> &d,
                              std::vector<double> &e, Matrix &V);
    // diagonalize bidiagonal matrix via implicit QR
    static void diagonalize(Matrix &U, std::vector<double> &d, std::vector<double> &e, Matrix &V);

    // Givens rotation, returns radius r = sqrt(a^2 + b^2)
    static void givens(double a, double b, double &c, double &s, double &r);

    SVD() = delete;
};

// Convenience wrapper for least squares
class LeastSquares
{
public:
    enum class Method
    {
        QR,
        SVD
    };
    // solve min ||Ax - b||_2
    static std::vector<double> solve(const Matrix &A, const std::vector<double> &b,
                                     Method method = Method::SVD, double svdTol = 1e-10);

private:
    LeastSquares() = delete;
};

#endif // CPPFM_LINEARALGEBRA_H
