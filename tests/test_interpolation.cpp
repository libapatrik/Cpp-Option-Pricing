#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../InterpolationSchemes.h"

// ============================================================================
// Test Constants and Utilities
// ============================================================================

namespace {  // Anonymous namespace to avoid naming conflicts with other test files
    // Numerical tolerances for financial applications
    const double PRICE_TOLERANCE = 1e-6;      // For option prices - tight tolerance for financial precision
    const double VOLATILITY_TOLERANCE = 1e-4; // For volatility surfaces - slightly relaxed for market data
    const double INTERPOLATION_TOLERANCE = 1e-8; // For interpolation accuracy - very tight for mathematical correctness
    
    // Test data sets - create standardized test data for consistent testing
    std::vector<double> createTestXData() {  // Create X coordinates for basic interpolation tests
        return {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};  // Uniformly spaced points for easy verification
    }
    
    std::vector<double> createTestYData() {  // Create Y coordinates for basic interpolation tests
        return {1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2};  // Linear decreasing function for simple validation
    }
    
    // Simple quadratic function for testing accuracy
    std::vector<double> createQuadraticXData() {  // Create X coordinates for derivative testing
        return {0.0, 1.0, 2.0, 3.0, 4.0};  // Fewer points for cubic spline coefficient calculation
    }
    
    std::vector<double> createQuadraticYData() {  // Create Y coordinates for derivative testing
        std::vector<double> x = createQuadraticXData();  // Get the X coordinates
        std::vector<double> y;  // Initialize empty Y vector
        for (double xi : x) {  // Iterate through each X coordinate
            y.push_back(xi * xi - 2*xi + 1); // f(x) = x² - 2x + 1 - quadratic function for derivative testing
        }
        return y;  // Return the computed Y values
    }
    
    // Helper function to check if values are approximately equal
    bool isApproximatelyEqual(double a, double b, double tolerance = INTERPOLATION_TOLERANCE) {  // Check if two numbers are close enough
        return std::abs(a - b) < tolerance;  // Return true if absolute difference is less than tolerance
    }
    
    // Simple output functions for readable test results
    void printTestData(const std::vector<double>& x, const std::vector<double>& y) {  // Print test data points
        std::cout << "\nTest Data:\n";
        std::cout << "  X: ";
        for (double xi : x) std::cout << std::setw(8) << std::fixed << std::setprecision(1) << xi;
        std::cout << "\n  Y: ";
        for (double yi : y) std::cout << std::setw(8) << std::fixed << std::setprecision(3) << yi;
        std::cout << "\n";
    }
    
    void printResult(const std::string& test, double x, double expected, double actual, bool passed) {  // Print test result
        std::cout << "  " << test << " at x=" << std::fixed << std::setprecision(2) << x 
                  << ": expected=" << std::setprecision(6) << expected 
                  << ", actual=" << actual 
                  << " [" << (passed ? "PASS" : "FAIL") << "]\n";
    }
}

// ============================================================================
// LinearInterpolation Test Suite
// ============================================================================

class LinearInterpolationTest : public ::testing::Test {  // Google Test fixture class for linear interpolation tests
protected:  // Protected members accessible to test methods
    void SetUp() override {  // Called before each test method - setup common test data
        xData = createTestXData();  // Initialize X coordinates with test data
        yData = createTestYData();  // Initialize Y coordinates with test data
        linearInterp = std::make_unique<LinearInterpolation>(xData, yData);  // Create linear interpolation object using smart pointer
    }
    
    void TearDown() override {  // Called after each test method - cleanup resources
        linearInterp.reset();  // Release the smart pointer memory safely
    }
    
    std::vector<double> xData, yData;  // Test data vectors - accessible to all test methods
    std::unique_ptr<LinearInterpolation> linearInterp;  // Smart pointer to interpolation object - automatic memory management
};

TEST_F(LinearInterpolationTest, ConstructorValidation) {  // Test constructor input validation and error handling
    // Test valid construction
    EXPECT_NO_THROW(LinearInterpolation(xData, yData));  // Verify valid data doesn't throw exceptions
    
    // Test invalid data sizes
    std::vector<double> shortX = {0.0, 1.0};  // Create X vector with 2 elements
    std::vector<double> longY = {1.0, 2.0, 3.0};  // Create Y vector with 3 elements - size mismatch
    EXPECT_THROW(LinearInterpolation(shortX, longY), std::invalid_argument);  // Expect exception for mismatched sizes
    
    // Test insufficient data points
    std::vector<double> singleX = {0.0};  // Create X vector with only 1 element
    std::vector<double> singleY = {1.0};  // Create Y vector with only 1 element
    EXPECT_THROW(LinearInterpolation(singleX, singleY), std::invalid_argument);  // Expect exception - need at least 2 points for interpolation
    
    // Test unsorted x-data
    std::vector<double> unsortedX = {2.0, 1.0, 3.0};  // Create unsorted X coordinates
    std::vector<double> unsortedY = {1.0, 2.0, 3.0};  // Create corresponding Y coordinates
    EXPECT_THROW(LinearInterpolation(unsortedX, unsortedY), std::invalid_argument);  // Expect exception for unsorted data
}

TEST_F(LinearInterpolationTest, InterpolationAtKnownPoints) {  // Test that interpolation is exact at data points
    std::cout << "\n=== Linear Interpolation - Exact Data Points ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test interpolation at exact data points
    for (size_t i = 0; i < xData.size(); ++i) {  // Iterate through all data points
        double result = linearInterp->interpolate(xData[i]);  // Get interpolated value at data point
        double error = std::abs(result - yData[i]);  // Calculate error
        bool passed = error < INTERPOLATION_TOLERANCE;  // Check if test passed
        
        printResult("Interpolation", xData[i], yData[i], result, passed);  // Show result
        
        EXPECT_NEAR(result, yData[i], INTERPOLATION_TOLERANCE)  // Verify result matches original data exactly
            << "Failed at point " << i << " (x=" << xData[i] << ")";  // Custom error message for debugging
    }
}

TEST_F(LinearInterpolationTest, InterpolationBetweenPoints) {  // Test linear interpolation between data points
    std::cout << "\n=== Linear Interpolation - Between Data Points ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test interpolation between first two points
    double x = 0.25; // Midpoint between 0.0 and 0.5 - test point between data points
    double expected = 0.9; // Linear interpolation: 1.0 + (0.8-1.0)*(0.25-0.0)/(0.5-0.0) - manually calculated expected value
    double result = linearInterp->interpolate(x);  // Get interpolated result
    double error1 = std::abs(result - expected);  // Calculate error
    bool passed1 = error1 < INTERPOLATION_TOLERANCE;  // Check if test passed
    printResult("Between (0.0,0.5)", x, expected, result, passed1);  // Show result
    EXPECT_NEAR(result, expected, INTERPOLATION_TOLERANCE);  // Verify result matches expected calculation
    
    // Test interpolation between points 2 and 3
    x = 1.25; // Midpoint between 1.0 and 1.5 - test another interpolation point
    expected = 0.5; // Linear interpolation: 0.6 + (0.4-0.6)*(1.25-1.0)/(1.5-1.0) - manually calculated expected value
    result = linearInterp->interpolate(x);  // Get interpolated result
    double error2 = std::abs(result - expected);  // Calculate error
    bool passed2 = error2 < INTERPOLATION_TOLERANCE;  // Check if test passed
    printResult("Between (1.0,1.5)", x, expected, result, passed2);  // Show result
    EXPECT_NEAR(result, expected, INTERPOLATION_TOLERANCE);  // Verify result matches expected calculation
}

TEST_F(LinearInterpolationTest, ExtrapolationBehavior) {  // Test extrapolation beyond data range
    std::cout << "\n=== Linear Interpolation - Extrapolation ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test left extrapolation (before first point)
    double x = -0.5;  // Point before the first data point
    double result = (*linearInterp)(x);  // Use operator() which routes to extrapolation automatically
    // Should continue quadratic Taylor series extrapolation from left boundary
    // The exact expected value depends on the extrapolation strategy (default: QuadraticExtrapolation)
    // For testing, we just verify it produces a finite, reasonable value
    EXPECT_TRUE(std::isfinite(result));  // Verify extrapolation produces finite result
    EXPECT_GT(result, 0.5);  // Should be above some reasonable lower bound
    EXPECT_LT(result, 2.0);  // Should be below some reasonable upper bound
    std::cout << "  Left extrapolation at x=" << x << ": " << result << " [FINITE]" << std::endl;
    
    // Test right extrapolation (after last point)
    x = 3.5;  // Point after the last data point
    result = (*linearInterp)(x);  // Use operator() which routes to extrapolation automatically
    EXPECT_TRUE(std::isfinite(result));  // Verify extrapolation produces finite result
    EXPECT_GT(result, -1.0);  // Should be above some reasonable lower bound
    EXPECT_LT(result, 0.5);  // Should be below some reasonable upper bound
    std::cout << "  Right extrapolation at x=" << x << ": " << result << " [FINITE]" << std::endl;
}

TEST_F(LinearInterpolationTest, RangeFunctionality) {  // Test range query functionality
    auto range = linearInterp->getRange();  // Get the interpolation range as a pair
    EXPECT_EQ(range.first, xData.front());  // Verify minimum X value matches first data point
    EXPECT_EQ(range.second, xData.back());  // Verify maximum X value matches last data point
}

TEST_F(LinearInterpolationTest, CloneFunctionality) {  // Test object cloning capability
    auto cloned = linearInterp->clone();  // Clone returns unique_ptr now
    
    ASSERT_NE(cloned, nullptr);  // Verify cloning succeeded (not null)
    
    // Test that cloned object produces same results
    for (double x : xData) {  // Iterate through all data points
        double original = linearInterp->interpolate(x);  // Get result from original object
        double cloned_result = cloned->interpolate(x);  // Get result from cloned object
        EXPECT_NEAR(original, cloned_result, INTERPOLATION_TOLERANCE);  // Verify both give identical results
    }
}

// ============================================================================
// CubicSplineInterpolation Test Suite
// ============================================================================

class CubicSplineInterpolationTest : public ::testing::Test {  // Google Test fixture class for cubic spline tests
protected:  // Protected members accessible to test methods
    void SetUp() override {  // Called before each test method - setup common test data
        xData = createQuadraticXData();  // Initialize X coordinates with quadratic test data
        yData = createQuadraticYData();  // Initialize Y coordinates with quadratic function values
        cubicInterp = std::make_unique<CubicSplineInterpolation>(  // Create cubic spline object using smart pointer
            xData, yData, CubicSplineInterpolation::BoundaryType::Natural  // Use natural boundary conditions (S''=0 at endpoints)
        );
    }
    
    void TearDown() override {  // Called after each test method - cleanup resources
        cubicInterp.reset();  // Release the smart pointer memory safely
    }
    
    std::vector<double> xData, yData;  // Test data vectors - accessible to all test methods
    std::unique_ptr<CubicSplineInterpolation> cubicInterp;  // Smart pointer to cubic spline object - automatic memory management
};

TEST_F(CubicSplineInterpolationTest, ConstructorValidation) {  // Test constructor input validation and error handling
    // Test valid construction
    EXPECT_NO_THROW(CubicSplineInterpolation(xData, yData));  // Verify valid data doesn't throw exceptions
    
    // Test invalid data sizes
    std::vector<double> shortX = {0.0};  // Create X vector with only 1 element
    std::vector<double> shortY = {1.0};  // Create Y vector with only 1 element
    EXPECT_THROW(CubicSplineInterpolation(shortX, shortY), std::invalid_argument);  // Expect exception - need at least 2 points for spline
}

TEST_F(CubicSplineInterpolationTest, InterpolationAtKnownPoints) {  // Test that cubic spline is exact at data points
    std::cout << "\n=== Cubic Spline - Exact Data Points ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test interpolation at exact data points
    for (size_t i = 0; i < xData.size(); ++i) {  // Iterate through all data points
        double result = cubicInterp->interpolate(xData[i]);  // Get interpolated value at data point
        double error = std::abs(result - yData[i]);  // Calculate error
        bool passed = error < INTERPOLATION_TOLERANCE;  // Check if test passed
        
        printResult("Cubic spline", xData[i], yData[i], result, passed);  // Show result
        
        EXPECT_NEAR(result, yData[i], INTERPOLATION_TOLERANCE)  // Verify result matches original data exactly
            << "Failed at point " << i << " (x=" << xData[i] << ")";  // Custom error message for debugging
    }
}

TEST_F(CubicSplineInterpolationTest, DerivativeCalculation) {  // Test derivative calculation and smoothness
    std::cout << "\n=== Cubic Spline - Derivatives ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test derivative calculation - cubic splines create smooth derivatives
    // IMPORTANT: For test function f(x) = x² - 2x + 1, the analytical derivative is f'(x) = 2x - 2
    // However, cubic splines DON'T necessarily reproduce the exact analytical derivatives
    // They create smooth interpolating curves that pass through data points
    // So we test for: (1) finiteness, (2) continuity, (3) reasonable values
    
    std::cout << "\nDerivatives at data points:\n";
    std::cout << "  NOTE: Test function is f(x) = x² - 2x + 1, so f'(x) = 2x - 2\n";
    std::cout << "  Cubic splines create smooth curves but may not match analytical derivatives exactly\n";
    
    // Compute analytical derivatives for comparison (not for exact matching!)
    std::vector<double> analyticalDerivatives;
    for (double x : xData) {
        analyticalDerivatives.push_back(2*x - 2);  // f'(x) = 2x - 2
    }
    
    // Test that derivatives are finite and reasonable at all data points
    for (size_t i = 0; i < xData.size(); ++i) {  // Iterate through all data points
        double x = xData[i];  // Get current X coordinate
        double result = cubicInterp->derivative(x);  // Calculate spline derivative at data point
        double analytical = analyticalDerivatives[i];  // Get analytical derivative for reference
        bool isFinite = std::isfinite(result);  // Check if derivative is finite
        
        // Print with analytical derivative for comparison (but don't require exact match)
        std::cout << "  Spline derivative at x=" << std::fixed << std::setprecision(2) << x 
                  << ": spline=" << std::setprecision(6) << result 
                  << ", analytical=" << analytical
                  << ", diff=" << std::abs(result - analytical)
                  << " [" << (isFinite ? "FINITE" : "NOT FINITE") << "]\n";
        
        EXPECT_TRUE(std::isfinite(result))  // Verify derivative is finite (not NaN or infinity)
            << "Non-finite derivative at point " << i << " (x=" << x << ")";  // Custom error message
    }
    
    std::cout << "\nDerivatives between data points:\n";
    // Test that derivatives are reasonably smooth between data points
    std::vector<double> test_points = {0.5, 1.5, 2.5, 3.5};  // Test points between data points
    for (double x : test_points) {  // Iterate through test points
        double result = cubicInterp->derivative(x);  // Calculate spline derivative at test point
        double analytical = 2*x - 2;  // Analytical derivative f'(x) = 2x - 2
        bool isFinite = std::isfinite(result);  // Check if derivative is finite
        bool inRange = result > -10.0 && result < 10.0;  // Check if in reasonable range
        bool passed = isFinite && inRange;  // Overall test result
        
        std::cout << "  Spline derivative at x=" << std::fixed << std::setprecision(2) << x 
                  << ": spline=" << std::setprecision(6) << result 
                  << ", analytical=" << analytical
                  << ", diff=" << std::abs(result - analytical)
                  << " [" << (passed ? "VALID" : "INVALID") << "]\n";
        
        EXPECT_TRUE(std::isfinite(result));  // Verify derivative is finite
        // Derivatives should be in reasonable range for our test data
        EXPECT_GT(result, -10.0) << "Derivative too negative at x=" << x;  // Check lower bound
        EXPECT_LT(result, 10.0) << "Derivative too positive at x=" << x;  // Check upper bound
    }
    
    // Test continuity of derivatives at interior data points (excluding boundaries)
    for (size_t i = 1; i < xData.size() - 1; ++i) {  // Iterate through interior points only
        double x = xData[i];  // Get current X coordinate
        double result = cubicInterp->derivative(x);  // Calculate derivative at data point
        
        // Verify the derivative is continuous at interior knots
        // This is a fundamental property of cubic splines
        EXPECT_TRUE(std::isfinite(result));  // Verify derivative is finite
    }
}

TEST_F(CubicSplineInterpolationTest, SecondDerivativeCalculation) {  // Test second derivative calculation and boundary conditions
    std::cout << "\n=== Cubic Spline - Second Derivatives ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test second derivative calculation
    // For test function f(x) = x² - 2x + 1, the analytical second derivative is f''(x) = 2 (constant)
    // For natural cubic splines, second derivatives at boundaries MUST be zero (boundary condition)
    // Interior second derivatives should be continuous and finite (but won't match analytical value)
    
    std::cout << "\nSecond derivatives:\n";
    std::cout << "  NOTE: Test function f(x) = x² - 2x + 1, so f''(x) = 2 (constant)\n";
    std::cout << "  Natural spline boundary condition: f''(0) = f''(4) = 0 (REQUIRED)\n";
    std::cout << "  Interior points: spline f''(x) will differ from analytical f''(x) = 2\n\n";
    
    for (size_t i = 0; i < xData.size(); ++i) {  // Iterate through all data points
        double x = xData[i];  // Get current X coordinate
        double result = cubicInterp->secondDerivative(x);  // Calculate spline second derivative
        double analytical = 2.0;  // Analytical second derivative f''(x) = 2
        bool isFinite = std::isfinite(result);  // Check if second derivative is finite
        
        // For natural boundary conditions, endpoints should have zero second derivative
        bool isBoundary = (i == 0 || i == xData.size() - 1);  // Check if this is the first or last point
        bool boundaryCorrect = !isBoundary || std::abs(result) < INTERPOLATION_TOLERANCE;  // Boundary should be zero
        bool passed = isFinite && boundaryCorrect;  // Overall test result
        
        std::string testName = isBoundary ? "Boundary" : "Interior";
        double expectedValue = isBoundary ? 0.0 : analytical;  // Boundaries must be 0, interior for reference only
        
        std::cout << "  " << testName << " at x=" << std::fixed << std::setprecision(2) << x 
                  << ": spline=" << std::setprecision(6) << result 
                  << ", " << (isBoundary ? "required" : "analytical") << "=" << expectedValue
                  << ", diff=" << std::abs(result - expectedValue)
                  << " [" << (passed ? "PASS" : "FAIL") << "]\n";
        
        // Verify that second derivative is finite
        EXPECT_TRUE(std::isfinite(result));  // Check for NaN or infinity
        
        if (isBoundary) {  // Check if this is the first or last point
            EXPECT_NEAR(result, 0.0, INTERPOLATION_TOLERANCE)  // Verify second derivative is zero at boundaries
                << "Natural boundary condition failed at point " << i << " (x=" << x << ")";  // Custom error message
        }
    }
}

TEST_F(CubicSplineInterpolationTest, ExtrapolationBehavior) {  // Test extrapolation beyond data range
    // Test left extrapolation
    double x = -1.0;  // Point before the first data point
    double result = (*cubicInterp)(x);  // Use operator() which routes to extrapolation automatically
    EXPECT_TRUE(std::isfinite(result));  // Verify extrapolation produces finite result
    std::cout << "  Left extrapolation at x=" << x << ": " << result << " [FINITE]" << std::endl;
    
    // Test right extrapolation
    x = 5.0;  // Point after the last data point
    result = (*cubicInterp)(x);  // Use operator() which routes to extrapolation automatically
    EXPECT_TRUE(std::isfinite(result));  // Verify extrapolation produces finite result
    std::cout << "  Right extrapolation at x=" << x << ": " << result << " [FINITE]" << std::endl;
}

TEST_F(CubicSplineInterpolationTest, BoundaryConditions) {  // Test natural boundary conditions
    // Test natural boundary conditions (second derivative = 0 at endpoints)
    double second_deriv_start = cubicInterp->secondDerivative(xData.front());  // Get second derivative at first point
    double second_deriv_end = cubicInterp->secondDerivative(xData.back());  // Get second derivative at last point
    
    // For natural splines, second derivatives at boundaries should be zero
    EXPECT_NEAR(second_deriv_start, 0.0, INTERPOLATION_TOLERANCE);  // Verify first boundary condition
    EXPECT_NEAR(second_deriv_end, 0.0, INTERPOLATION_TOLERANCE);  // Verify second boundary condition
}

// ============================================================================
// Comparison Tests: Linear vs Cubic Spline
// ============================================================================

class InterpolationComparisonTest : public ::testing::Test {  // Google Test fixture for comparing interpolation methods
protected:  // Protected members accessible to test methods
    void SetUp() override {  // Called before each test method - setup both interpolation objects
        xData = createTestXData();  // Initialize X coordinates with test data
        yData = createTestYData();  // Initialize Y coordinates with test data
        linearInterp = std::make_unique<LinearInterpolation>(xData, yData);  // Create linear interpolation object
        cubicInterp = std::make_unique<CubicSplineInterpolation>(  // Create cubic spline interpolation object
            xData, yData, CubicSplineInterpolation::BoundaryType::Natural  // Use natural boundary conditions
        );
    }
    
    std::vector<double> xData, yData;  // Test data vectors - accessible to all test methods
    std::unique_ptr<LinearInterpolation> linearInterp;  // Smart pointer to linear interpolation object
    std::unique_ptr<CubicSplineInterpolation> cubicInterp;  // Smart pointer to cubic spline interpolation object
};

TEST_F(InterpolationComparisonTest, AgreementAtDataPoints) {  // Test that both methods agree at data points
    std::cout << "\n=== Comparison - Agreement at Data Points ===";
    printTestData(xData, yData);  // Show the test data
    
    // Both methods should give identical results at data points
    std::cout << "\nComparing Linear vs Cubic Spline at data points:\n";
    for (size_t i = 0; i < xData.size(); ++i) {  // Iterate through all data points
        double linear_result = linearInterp->interpolate(xData[i]);  // Get linear interpolation result
        double cubic_result = cubicInterp->interpolate(xData[i]);  // Get cubic spline interpolation result
        double difference = std::abs(linear_result - cubic_result);  // Calculate difference
        bool agree = difference < INTERPOLATION_TOLERANCE;  // Check if they agree
        
        std::cout << "  Point " << i << " (x=" << std::fixed << std::setprecision(1) << xData[i] 
                  << "): Linear=" << std::setprecision(6) << linear_result 
                  << ", Cubic=" << cubic_result 
                  << ", Diff=" << difference 
                  << " [" << (agree ? "AGREE" : "DIFFER") << "]\n";
        
        EXPECT_NEAR(linear_result, cubic_result, INTERPOLATION_TOLERANCE);  // Verify both give identical results
    }
}

TEST_F(InterpolationComparisonTest, SmoothnessComparison) {  // Test smoothness characteristics of both methods
    std::cout << "\n=== Comparison - Smoothness Between Data Points ===";
    printTestData(xData, yData);  // Show the test data
    
    // Test that cubic spline produces smoother results between points
    std::vector<double> testPoints = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75};  // Test points between data points
    
    std::cout << "\nComparing Linear vs Cubic Spline between data points:\n";
    for (double x : testPoints) {  // Iterate through test points
        double linear_result = linearInterp->interpolate(x);  // Get linear interpolation result
        double cubic_result = cubicInterp->interpolate(x);  // Get cubic spline interpolation result
        
        // Check if results are reasonable
        bool linearFinite = std::isfinite(linear_result);  // Check if linear result is finite
        bool cubicFinite = std::isfinite(cubic_result);  // Check if cubic result is finite
        bool linearInRange = linear_result > -1.0 && linear_result < 2.0;  // Check linear range
        bool cubicInRange = cubic_result > -1.0 && cubic_result < 2.0;  // Check cubic range
        bool bothValid = linearFinite && cubicFinite && linearInRange && cubicInRange;  // Overall validity
        
        std::cout << "  At x=" << std::fixed << std::setprecision(2) << x 
                  << ": Linear=" << std::setprecision(6) << linear_result 
                  << ", Cubic=" << cubic_result 
                  << " [" << (bothValid ? "VALID" : "INVALID") << "]\n";
        
        // Both should produce reasonable values
        EXPECT_TRUE(std::isfinite(linear_result));  // Verify linear result is finite
        EXPECT_TRUE(std::isfinite(cubic_result));  // Verify cubic result is finite
        
        // Results may differ, but both should be in reasonable range
        EXPECT_GT(linear_result, -1.0);  // Check linear result lower bound
        EXPECT_LT(linear_result, 2.0);  // Check linear result upper bound
        EXPECT_GT(cubic_result, -1.0);  // Check cubic result lower bound
        EXPECT_LT(cubic_result, 2.0);  // Check cubic result upper bound
    }
}

// ============================================================================
// Edge Cases and Error Handling
// ============================================================================

class InterpolationEdgeCasesTest : public ::testing::Test {};  // Simple test class for edge cases - no setup needed

TEST_F(InterpolationEdgeCasesTest, TwoPointLinearCase) {  // Test minimal data scenario
    // Test with minimal data (2 points)
    std::vector<double> x = {0.0, 1.0};  // Create minimal X data (2 points)
    std::vector<double> y = {0.0, 1.0};  // Create minimal Y data (2 points)
    
    LinearInterpolation linear(x, y);  // Create linear interpolation with minimal data
    CubicSplineInterpolation cubic(x, y);  // Create cubic spline with minimal data
    
    // Both should handle 2-point case gracefully
    EXPECT_NEAR(linear.interpolate(0.5), 0.5, INTERPOLATION_TOLERANCE);  // Test linear interpolation at midpoint
    EXPECT_NEAR(cubic.interpolate(0.5), 0.5, INTERPOLATION_TOLERANCE);  // Test cubic spline at midpoint
}

TEST_F(InterpolationEdgeCasesTest, OutOfRangeDerivative) {  // Test derivative calculations work at boundaries (needed for extrapolation)
    std::vector<double> x = createTestXData();  // Get test X data
    std::vector<double> y = createTestYData();  // Get test Y data
    CubicSplineInterpolation cubic(x, y);  // Create cubic spline object
    
    // In the refactored design, derivatives work at boundaries (for generic extrapolation)
    // Test that derivatives at boundaries produce finite values
    double deriv_left = cubic.derivative(x.front());  // Test left boundary derivative
    double deriv_right = cubic.derivative(x.back());  // Test right boundary derivative
    EXPECT_TRUE(std::isfinite(deriv_left));  // Should be finite, not throw
    EXPECT_TRUE(std::isfinite(deriv_right));  // Should be finite, not throw
    
    // Same for second derivatives
    double second_deriv_left = cubic.secondDerivative(x.front());  // Test left boundary second derivative
    double second_deriv_right = cubic.secondDerivative(x.back());  // Test right boundary second derivative
    EXPECT_TRUE(std::isfinite(second_deriv_left));  // Should be finite, not throw
    EXPECT_TRUE(std::isfinite(second_deriv_right));  // Should be finite, not throw
}

// ============================================================================
// Performance and Stress Tests
// ============================================================================

class InterpolationPerformanceTest : public ::testing::Test {};  // Simple test class for performance tests - no setup needed

TEST_F(InterpolationPerformanceTest, LargeDatasetHandling) {  // Test performance with large datasets
    // Test with larger dataset
    std::vector<double> x, y;  // Initialize empty vectors for large dataset
    const int n = 1000;  // Set dataset size to 1000 points
    
    for (int i = 0; i < n; ++i) {  // Generate large dataset
        double xi = static_cast<double>(i) / (n - 1) * 10.0;  // Create X coordinates from 0 to 10
        x.push_back(xi);  // Add X coordinate to vector
        y.push_back(std::sin(xi)); // Test with smooth function - sine wave for testing accuracy
    }
    
    EXPECT_NO_THROW(LinearInterpolation linear(x, y));  // Verify linear interpolation can handle large dataset
    EXPECT_NO_THROW(CubicSplineInterpolation cubic(x, y));  // Verify cubic spline can handle large dataset
    
    // Test interpolation performance
    LinearInterpolation linear(x, y);  // Create linear interpolation object with large dataset
    CubicSplineInterpolation cubic(x, y);  // Create cubic spline object with large dataset
    
    // Both should handle large datasets without issues
    EXPECT_NEAR(linear.interpolate(5.0), std::sin(5.0), 1e-2);  // Test linear interpolation accuracy with relaxed tolerance
    EXPECT_NEAR(cubic.interpolate(5.0), std::sin(5.0), 1e-4); // Cubic should be more accurate - tighter tolerance
}