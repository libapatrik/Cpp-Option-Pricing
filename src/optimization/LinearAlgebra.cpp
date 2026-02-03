#include <cppfm/optimization/LinearAlgebra.h>

// ===========================================================================
// Matrix Operations
// ===========================================================================

std::vector<double> MatrixOps::multiply(const Matrix &A, const std::vector<double> &x)
{
    size_t m = A.size();
    if (m == 0)
        return {};
    size_t n = A[0].size();
    if (n != x.size())
        throw std::invalid_argument("Dimension mismatch");

    std::vector<double> y(m, 0.0);
    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

Matrix MatrixOps::multiply(const Matrix &A, const Matrix &B)
{
    size_t m = A.size();
    if (m == 0)
        return {};
    size_t k = A[0].size();

    if (B.size() != k)
        throw std::invalid_argument("Dimension mismatch");
    size_t n = B[0].size();

    Matrix C = zeros(m, n);
    // i-l-j order for cache efficiency: B[l][j] now contiguous
    for (size_t i = 0; i < m; ++i) {
        for (size_t l = 0; l < k; ++l) {
            double a_il = A[i][l];
            for (size_t j = 0; j < n; ++j)
                C[i][j] += a_il * B[l][j];
        }
    }
    return C;
}

Matrix MatrixOps::transpose(const Matrix &A)
{
    size_t m = A.size();
    if (m == 0)
        return {};
    size_t n = A[0].size();

    Matrix At(n, std::vector<double>(m));

    // blocked transpose for cache efficiency on larger matrices
    constexpr size_t BLOCK = 8;
    for (size_t ii = 0; ii < m; ii += BLOCK) {
        for (size_t jj = 0; jj < n; jj += BLOCK) {
            size_t i_end = std::min(ii + BLOCK, m);
            size_t j_end = std::min(jj + BLOCK, n);
            for (size_t i = ii; i < i_end; ++i)
                for (size_t j = jj; j < j_end; ++j)
                    At[j][i] = A[i][j];
        }
    }
    return At;
}

double MatrixOps::norm2(const std::vector<double> &x)
{
    double sum = 0.0;
    for (double xi : x)
        sum += xi * xi;
    return std::sqrt(sum);
}

double MatrixOps::normInf(const std::vector<double> &x)
{
    double maxVal = 0.0;
    for (double xi : x)
        maxVal = std::max(maxVal, std::abs(xi));
    return maxVal;
}

double MatrixOps::dot(const std::vector<double> &x, const std::vector<double> &y)
{
    if (x.size() != y.size())
        throw std::invalid_argument("Dimension mismatch");
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
        sum += x[i] * y[i];
    return sum;
}

void MatrixOps::axpy(double alpha, const std::vector<double> &x, std::vector<double> &y)
{
    for (size_t i = 0; i < x.size(); ++i)
        y[i] += alpha * x[i];
}

void MatrixOps::scale(double alpha, std::vector<double> &x)
{
    for (double &xi : x)
        xi *= alpha;
}

Matrix MatrixOps::identity(size_t n)
{
    Matrix I = zeros(n, n);
    for (size_t i = 0; i < n; ++i)
        I[i][i] = 1.0;
    return I;
}

Matrix MatrixOps::zeros(size_t m, size_t n)
{
    return Matrix(m, std::vector<double>(n, 0.0));
}

Matrix MatrixOps::multiplyAtA(const Matrix &A)
{
    size_t m = A.size();
    if (m == 0) return {};
    size_t n = A[0].size();

    Matrix C = zeros(n, n);
    // only compute upper triangle, then mirror
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m; ++k)
                sum += A[k][i] * A[k][j];
            C[i][j] = sum;
            C[j][i] = sum;
        }
    }
    return C;
}

std::vector<double> MatrixOps::multiplyAtb(const Matrix &A, const std::vector<double> &b)
{
    size_t m = A.size();
    if (m == 0) return {};
    size_t n = A[0].size();

    std::vector<double> y(n, 0.0);
    for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < m; ++i)
            y[j] += A[i][j] * b[i];
    return y;
}

// ===========================================================================
// Cholesky
// ===========================================================================
Matrix Cholesky::decompose(const Matrix &S)
{
    size_t n = S.size();
    if (n == 0)
        return {};

    // check symmetry
    for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j)
            if (std::abs(S[i][j] - S[j][i]) > 1e-10)
                throw std::invalid_argument("Matrix is not symmetric");

    Matrix L = MatrixOps::zeros(n, n);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double sum = S[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];

            if (i == j)
            {
                if (sum <= 0.0)
                    throw std::runtime_error("Matrix not positive definite");
                L[i][j] = std::sqrt(sum);
            }
            else
            {
                L[i][j] = sum / L[j][j]; // divide diagonal
            }
        }
    }
    return L;
}

std::vector<double> Cholesky::solveL(const Matrix &L, const std::vector<double> &b)
{
    size_t n = L.size();
    std::vector<double> y(n);

    for (size_t i = 0; i < n; ++i)
    {
        double sum = b[i];
        for (size_t j = 0; j < i; ++j)
            sum -= L[i][j] * y[j];
        y[i] = sum / L[i][i];
    }
    return y;
}

std::vector<double> Cholesky::solveLT(const Matrix &L, const std::vector<double> &y)
{
    size_t n = L.size();
    std::vector<double> x(n);

    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = y[i];
        for (size_t j = i + 1; j < n; ++j)
            sum -= L[j][i] * x[j];
        x[i] = sum / L[i][i];
    }
    return x;
}

std::vector<double> Cholesky::solve(const Matrix &L, const std::vector<double> &b)
{
    auto y = solveL(L, b);
    return solveLT(L, y);
}

// ===========================================================================
// QR - Implicit Q via stored Householder vectors
// ===========================================================================
// Householder vectors stored in lower triangle of R, scalars in tau
// Saves O(m²) storage and O(m²) operations vs explicit Q

QRResult QR::decompose(const Matrix &A)
{
    size_t m = A.size();
    if (m == 0)
        return {};
    size_t n = A[0].size();
    size_t k = std::min(m, n);

    Matrix R = A;
    std::vector<double> tau(k);

    for (size_t j = 0; j < k; ++j)
    {
        // Householder vector for col j
        double norm = 0.0;
        for (size_t i = j; i < m; ++i)
            norm += R[i][j] * R[i][j];
        norm = std::sqrt(norm);

        if (norm < 1e-15) {
            tau[j] = 0.0;
            continue;
        }

        double sign = (R[j][j] >= 0) ? 1.0 : -1.0;
        double u1 = R[j][j] + sign * norm;

        // store scaled Householder vector in lower triangle
        for (size_t i = j + 1; i < m; ++i)
            R[i][j] /= u1;

        tau[j] = sign * u1 / norm;

        // apply H to remaining columns of R
        for (size_t c = j + 1; c < n; ++c)
        {
            double dot = R[j][c];
            for (size_t i = j + 1; i < m; ++i)
                dot += R[i][j] * R[i][c];
            dot *= tau[j];
            R[j][c] -= dot;
            for (size_t i = j + 1; i < m; ++i)
                R[i][c] -= R[i][j] * dot;
        }

        R[j][j] = -sign * norm;  // diagonal of R
    }

    return {std::move(R), std::move(tau)};
}

std::vector<double> QR::applyQT(const QRResult &qr, const std::vector<double> &b)
{
    // apply Q^T = H_k * ... * H_1 to b (forward order)
    size_t m = b.size();
    size_t k = qr.tau.size();

    std::vector<double> c = b;

    for (size_t j = 0; j < k; ++j)
    {
        if (std::abs(qr.tau[j]) < 1e-15) continue;

        double dot = c[j];
        for (size_t i = j + 1; i < m; ++i)
            dot += qr.R[i][j] * c[i];
        dot *= qr.tau[j];

        c[j] -= dot;
        for (size_t i = j + 1; i < m; ++i)
            c[i] -= qr.R[i][j] * dot;
    }

    return c;
}

std::vector<double> QR::applyQ(const QRResult &qr, const std::vector<double> &x)
{
    // apply Q = H_1 * ... * H_k to x (reverse order)
    size_t m = x.size();
    size_t k = qr.tau.size();

    std::vector<double> y = x;

    for (int j = static_cast<int>(k) - 1; j >= 0; --j)
    {
        if (std::abs(qr.tau[j]) < 1e-15) continue;

        double dot = y[j];
        for (size_t i = j + 1; i < m; ++i)
            dot += qr.R[i][j] * y[i];
        dot *= qr.tau[j];

        y[j] -= dot;
        for (size_t i = j + 1; i < m; ++i)
            y[i] -= qr.R[i][j] * dot;
    }

    return y;
}

std::vector<double> QR::solve(const QRResult &qr, const std::vector<double> &b)
{
    // O(mn) via implicit Q instead of O(m²)
    auto c = applyQT(qr, b);
    return solveR(qr.R, c);
}

std::vector<double> QR::solveR(const Matrix &R, const std::vector<double> &c)
{
    size_t n = std::min(R.size(), R[0].size());
    std::vector<double> x(n);

    // relative tolerance based on R norm
    double Rnorm = 0.0;
    for (size_t i = 0; i < n; ++i)
        Rnorm = std::max(Rnorm, std::abs(R[i][i]));
    double tol = 1e-14 * Rnorm;

    for (int i = static_cast<int>(n) - 1; i >= 0; --i)
    {
        double sum = c[i];
        for (size_t j = i + 1; j < n; ++j)
            sum -= R[i][j] * x[j];

        if (std::abs(R[i][i]) < tol)
            throw std::runtime_error("R is singular at diagonal " + std::to_string(i));
        x[i] = sum / R[i][i];
    }
    return x;
}

// ===========================================================================
// SVD
// ===========================================================================

SVDResult SVD::decompose(const Matrix &A)
{
    if (A.empty())
        return {};

    size_t m = A.size();
    size_t n = A[0].size();

    Matrix U;
    std::vector<double> d, e;
    Matrix V;

    bidiagonalize(A, U, d, e, V);
    diagonalize(U, d, e, V);

    // sort singular values descending order
    std::vector<size_t> idx(d.size());
    for (size_t i = 0; i < d.size(); ++i)
        idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return std::abs(d[a]) > std::abs(d[b]);
    });

    SVDResult result;
    result.sigma.resize(d.size());
    result.U.resize(m, std::vector<double>(d.size()));
    result.V.resize(n, std::vector<double>(d.size()));

    for (size_t j = 0; j < d.size(); ++j)
    {
        result.sigma[j] = std::abs(d[idx[j]]);
        for (size_t i = 0; i < m; ++i)
            result.U[i][j] = U[i][idx[j]];
        for (size_t i = 0; i < n; ++i)
            result.V[i][j] = V[i][idx[j]];
    }

    return result;
}

std::vector<double> SVD::solve(const SVDResult &svd, const std::vector<double> &b, double tol)
{
    size_t m = svd.U.size();
    size_t n = svd.V.size();
    size_t k = svd.sigma.size();

    double thresh = tol * svd.sigma[0];

    std::vector<double> x(n, 0.0);

    for (size_t j = 0; j < k; ++j)
    {
        if (svd.sigma[j] > thresh)
        {
            double s = 0.0;
            for (size_t i = 0; i < m; ++i)
                s += svd.U[i][j] * b[i];
            s /= svd.sigma[j];
            for (size_t i = 0; i < n; ++i)
                x[i] += s * svd.V[i][j];
        }
    }

    return x;
}

void SVD::bidiagonalize(const Matrix &A, Matrix &U, std::vector<double> &d, std::vector<double> &e, Matrix &V)
{
    size_t m = A.size();
    size_t n = A[0].size();
    size_t k = std::min(m, n);

    U = A;
    V = MatrixOps::identity(n);
    d.resize(k);
    e.resize(k > 0 ? k - 1 : 0);

    for (size_t j = 0; j < k; ++j)
    {
        // left Housholder
        double norm = 0.0;
        for (size_t i = j; i < m; ++i)
            norm += U[i][j] * U[i][j];
        norm = std::sqrt(norm);

        if (norm > 1e-15)
        {
            if (U[j][j] < 0)
                norm = -norm;
            for (size_t i = j; i < m; ++i)
                U[i][j] /= norm;
            U[j][j] += 1.0;

            for (size_t c = j + 1; c < n; ++c)
            {
                double s = 0.0;
                for (size_t i = j; i < m; ++i)
                    s += U[i][j] * U[i][c];
                s /= U[j][j];
                for (size_t i = j; i < m; ++i)
                    U[i][c] -= s * U[i][j];
            }
        }
        d[j] = -norm;

        // right Householder
        if (j < n - 1)
        {
            norm = 0.0;
            for (size_t c = j + 1; c < n; ++c)
                norm += U[j][c] * U[j][c];
            norm = std::sqrt(norm);

            if (norm > 1e-15)
            {
                if (U[j][j + 1] < 0)
                    norm = -norm;
                for (size_t c = j + 1; c < n; ++c)
                    U[j][c] /= norm;
                U[j][j + 1] += 1.0;

                for (size_t r = j + 1; r < m; ++r)
                {
                    double s = 0.0;
                    for (size_t c = j + 1; c < n; ++c)
                        s += U[j][c] * U[r][c];
                    s /= U[j][j + 1];
                    for (size_t c = j + 1; c < n; ++c)
                        U[r][c] -= s * U[j][c];
                }

                for (size_t r = 0; r < n; ++r)
                {
                    double s = 0.0;
                    for (size_t c = j + 1; c < n; ++c)
                        s += U[j][c] * V[r][c];
                    s /= U[j][j + 1];
                    for (size_t c = j + 1; c < n; ++c)
                        V[r][c] -= s * U[j][c];
                }
            }
            e[j] = -norm;
        }
    }

    // accumulate thin U (m × k) instead of full U (m × m)
    // for tall-skinny matrices this is O(mk) instead of O(m²)
    Matrix Uacc = MatrixOps::zeros(m, k);
    for (size_t i = 0; i < k; ++i)
        Uacc[i][i] = 1.0;

    for (int j = static_cast<int>(k) - 1; j >= 0; --j)
    {
        if (std::abs(U[j][j]) > 1e-15)
        {
            for (size_t c = j; c < k; ++c)  // only k columns, not m
            {
                double s = 0.0;
                for (size_t i = j; i < m; ++i)
                    s += U[i][j] * Uacc[i][c];
                s /= U[j][j];
                for (size_t i = j; i < m; ++i)
                    Uacc[i][c] -= s * U[i][j];
            }
        }
    }
    U = Uacc;
}

void SVD::diagonalize(Matrix &U, std::vector<double> &d, std::vector<double> &e, Matrix &V)
{
    size_t n = d.size();
    if (n == 0)
        return;

    const int maxIter = 50 * static_cast<int>(n);
    int iter = 0;

    for (size_t kk = n; kk > 0; )
    {
        size_t k = kk - 1;

        // find block
        size_t l = k;
        while (l > 0 && std::abs(e[l - 1]) > 1e-15 * (std::abs(d[l - 1]) + std::abs(d[l])))
            --l;

        if (l == k)
        {
            // converged
            if (d[k] < 0)
            {
                d[k] = -d[k];
                for (size_t i = 0; i < V.size(); ++i)
                    V[i][k] = -V[i][k];
            }
            --kk;
            continue;
        }

        if (++iter > maxIter)
            throw std::runtime_error("SVD failed to converge after " + std::to_string(maxIter) + " iterations");

        // QR step with Wilkinson shift
        double shift = d[k];
        double g = (d[k - 1] - shift) / (2.0 * e[k - 1]);
        double r = std::sqrt(g * g + 1.0);
        g = d[l] - shift + e[l] / (g + (g >= 0 ? r : -r));

        double c = 1.0;
        double s = 1.0;
        double p = 0.0;

        for (size_t i = l; i < k; ++i)
        {
            double f = s * e[i];
            double b = c * e[i];

            givens(g, f, c, s);
            if (i > l)
                e[i - 1] = std::sqrt(g * g + f * f);

            g = d[i] - p;
            r = (d[i + 1] - g) * s + 2.0 * c * b;
            p = s * r;
            d[i] = g + p;
            g = c * r - b;

            for (size_t j = 0; j < V.size(); ++j)
            {
                double t = V[j][i + 1];
                V[j][i + 1] = s * V[j][i] + c * t;
                V[j][i] = c * V[j][i] - s * t;
            }
            for (size_t j = 0; j < U.size(); ++j)
            {
                double t = U[j][i + 1];
                U[j][i + 1] = s * U[j][i] + c * t;
                U[j][i] = c * U[j][i] - s * t;
            }
        }
        d[k] -= p;
        e[k - 1] = g;
    }
}

void SVD::givens(double a, double b, double &c, double &s)
{
    if (std::abs(b) < 1e-15)
    {
        c = 1.0;
        s = 0.0;
    }
    else if (std::abs(b) > std::abs(a))
    {
        double t = -a / b;
        s = 1.0 / std::sqrt(1.0 + t * t);
        c = s * t;
    }
    else
    {
        double t = -b / a;
        c = 1.0 / std::sqrt(1.0 + t * t);
        s = c * t;
    }
}

// ============================================================================
// LeastSquares
// +===========================================================================

std::vector<double> LeastSquares::solve(const Matrix &A, const std::vector<double> &b, Method method, double svdTol)
{
    if (method == Method::QR) {
        auto qr = QR::decompose(A);
        return QR::solve(qr, b);
    } else {
        auto svd = SVD::decompose(A);
        return SVD::solve(svd, b, svdTol);
    }
}
