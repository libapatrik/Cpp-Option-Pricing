/**
    Linear Algebra
    - Cholesky decomposition
    - QR decomposition
    - SVD decomposition
 

*/

class Cholesky {
public: 
    // S = LL^T; returns L (lower triangular matrix)
    static std::vector<std::vector<double>> decompose(const std::vector<std::vector<double>>& correlationMatrix);
};

class QR {
public: 
    static std::vector<std::vector<double>> decompose(const std::vector<std::vector<double>>& matrix);
};

class SVD {
public: 
    static std::vector<std::vector<double>> decompose(const std::vector<std::vector<double>>& matrix);
};

class LeastSquares {
public: 
    static std::vector<double> solve(const std::vector<std::vector<double>>& matrix, const std::vector<double>& rhs);
};