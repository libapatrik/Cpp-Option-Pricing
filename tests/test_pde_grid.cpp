#include <gtest/gtest.h>
#include "../PDEs/Grid.h"
#include <cmath>

// ============================================================================
// Test Tolerances
// ============================================================================

constexpr double GRID_TOL = 1e-10;  // Tight tolerance for grid construction

// ============================================================================
// PDEGrid Construction Tests
// ============================================================================

TEST(PDEGridTest, UniformGridConstruction) {
    // Create uniform grid: [0, 100] with 101 points, T=1 with 100 steps
    Grid grid(0.0, 100.0, 101, 1.0, 100, Grid::SpacingType::Uniform);
    
    // Verify dimensions
    EXPECT_EQ(grid.numSpotPoints(), 101);
    EXPECT_EQ(grid.numTimeSteps(), 100);
    
    // Verify bounds
    EXPECT_NEAR(grid.spotMin(), 0.0, GRID_TOL);
    EXPECT_NEAR(grid.spotMax(), 100.0, GRID_TOL);
    EXPECT_NEAR(grid.timeMax(), 1.0, GRID_TOL);
    
    // Verify time step
    EXPECT_NEAR(grid.dt(), 0.01, GRID_TOL);
    
    // Verify spacing type
    EXPECT_EQ(grid.spacingType(), Grid::SpacingType::Uniform);
}

TEST(PDEGridTest, UniformGridSpacing) {
    Grid grid(0.0, 100.0, 101, 1.0, 100);
    
    // For uniform grid, spacing should be constant
    double expected_dS = 1.0;  // (100 - 0) / (101 - 1) = 1.0
    
    // Check spacing at several interior points
    for (size_t i = 1; i < 100; ++i) {
        double dS = grid.dS(i);
        EXPECT_NEAR(dS, expected_dS, GRID_TOL) 
            << "Spacing not uniform at point " << i;
    }
    
    // Check forward spacing
    for (size_t i = 0; i < 100; ++i) {
        double dS_fwd = grid.dS_forward(i);
        EXPECT_NEAR(dS_fwd, expected_dS, GRID_TOL);
    }
}

TEST(PDEGridTest, UniformGridValues) {
    Grid grid(0.0, 100.0, 11, 1.0, 10);
    
    // Verify spot values at key points
    EXPECT_NEAR(grid.spot(0), 0.0, GRID_TOL);
    EXPECT_NEAR(grid.spot(5), 50.0, GRID_TOL);
    EXPECT_NEAR(grid.spot(10), 100.0, GRID_TOL);
    
    // Verify time values
    EXPECT_NEAR(grid.time(0), 0.0, GRID_TOL);
    EXPECT_NEAR(grid.time(5), 0.5, GRID_TOL);
    EXPECT_NEAR(grid.time(10), 1.0, GRID_TOL);
}

TEST(PDEGridTest, LogSpacedGridConstruction) {
    // Create log-spaced grid
    Grid grid(10.0, 200.0, 101, 1.0, 100, Grid::SpacingType::LogSpaced);
    
    // Verify dimensions (same as uniform)
    EXPECT_EQ(grid.numSpotPoints(), 101);
    EXPECT_EQ(grid.numTimeSteps(), 100);
    
    // Verify bounds (exact)
    EXPECT_NEAR(grid.spot(0), 10.0, GRID_TOL);
    EXPECT_NEAR(grid.spot(100), 200.0, GRID_TOL);
    
    // Verify spacing type
    EXPECT_EQ(grid.spacingType(), Grid::SpacingType::LogSpaced);
}

TEST(PDEGridTest, LogSpacedGridProperties) {
    Grid grid(10.0, 200.0, 101, 1.0, 100, Grid::SpacingType::LogSpaced);
    
    // Log-spaced grid should have:
    // 1. Smaller spacing near S_min (more points near lower boundary)
    // 2. Larger spacing near S_max
    
    double dS_low = grid.spot(1) - grid.spot(0);
    double dS_mid = grid.spot(51) - grid.spot(50);
    double dS_high = grid.spot(100) - grid.spot(99);
    
    EXPECT_LT(dS_low, dS_mid) 
        << "Log-spaced grid should have smaller spacing near lower bound";
    EXPECT_LT(dS_mid, dS_high)
        << "Log-spaced grid should have larger spacing near upper bound";
}

TEST(PDEGridTest, LogSpacedGridMonotonicity) {
    Grid grid(10.0, 200.0, 101, 1.0, 100, Grid::SpacingType::LogSpaced);
    
    // Grid should be strictly increasing
    for (size_t i = 1; i < 101; ++i) {
        EXPECT_GT(grid.spot(i), grid.spot(i-1))
            << "Grid not monotonically increasing at point " << i;
    }
}

// ============================================================================
// Time Grid Tests
// ============================================================================

TEST(PDEGridTest, TimeGridConstruction) {
    Grid grid(0.0, 100.0, 101, 2.0, 200);
    
    // Verify time grid size (should be N_t + 1)
    EXPECT_EQ(grid.times().size(), 201);  // 200 steps + initial point
    
    // Verify time step
    EXPECT_NEAR(grid.dt(), 0.01, GRID_TOL);  // 2.0 / 200 = 0.01
    
    // Verify time values
    EXPECT_NEAR(grid.time(0), 0.0, GRID_TOL);
    EXPECT_NEAR(grid.time(100), 1.0, GRID_TOL);
    EXPECT_NEAR(grid.time(200), 2.0, GRID_TOL);
}

TEST(PDEGridTest, TimeGridUniformity) {
    Grid grid(0.0, 100.0, 101, 1.0, 100);
    
    double dt = grid.dt();
    
    // Check that time steps are uniform
    for (size_t j = 1; j <= 100; ++j) {
        double time_diff = grid.time(j) - grid.time(j-1);
        EXPECT_NEAR(time_diff, dt, GRID_TOL)
            << "Time step not uniform at step " << j;
    }
}

// ============================================================================
// Parameter Validation Tests
// ============================================================================

TEST(PDEGridTest, InvalidSpotRange) {
    // S_min < 0
    EXPECT_THROW(
        Grid grid(-1.0, 100.0, 101, 1.0, 100),
        std::invalid_argument
    );
    
    // S_max <= S_min
    EXPECT_THROW(
        Grid grid(100.0, 50.0, 101, 1.0, 100),
        std::invalid_argument
    );
    
    // S_max == S_min
    EXPECT_THROW(
        Grid grid(50.0, 50.0, 101, 1.0, 100),
        std::invalid_argument
    );
}

TEST(PDEGridTest, InvalidGridDimensions) {
    // Too few spatial points
    EXPECT_THROW(
        Grid grid(0.0, 100.0, 2, 1.0, 100),
        std::invalid_argument
    );
    
    // Zero time steps
    EXPECT_THROW(
        Grid grid(0.0, 100.0, 101, 1.0, 0),
        std::invalid_argument
    );
}

TEST(PDEGridTest, InvalidTimeRange) {
    // Negative T_max
    EXPECT_THROW(
        Grid grid(0.0, 100.0, 101, -1.0, 100),
        std::invalid_argument
    );
    
    // Zero T_max
    EXPECT_THROW(
        Grid grid(0.0, 100.0, 101, 0.0, 100),
        std::invalid_argument
    );
}

// ============================================================================
// Clone and Copy Tests
// ============================================================================

TEST(PDEGridTest, CloneFunctionality) {
    Grid original(0.0, 100.0, 101, 1.0, 100);
    
    // Clone the grid
    Grid* cloned = original.clone();
    
    ASSERT_NE(cloned, nullptr);
    
    // Verify cloned grid has same properties
    EXPECT_EQ(cloned->numSpotPoints(), original.numSpotPoints());
    EXPECT_EQ(cloned->numTimeSteps(), original.numTimeSteps());
    EXPECT_NEAR(cloned->spotMin(), original.spotMin(), GRID_TOL);
    EXPECT_NEAR(cloned->spotMax(), original.spotMax(), GRID_TOL);
    EXPECT_NEAR(cloned->dt(), original.dt(), GRID_TOL);
    
    // Verify grid values match
    for (size_t i = 0; i < original.numSpotPoints(); ++i) {
        EXPECT_NEAR(cloned->spot(i), original.spot(i), GRID_TOL);
    }
    
    // Clean up
    delete cloned;
}

TEST(PDEGridTest, CopyConstructor) {
    Grid original(0.0, 100.0, 101, 1.0, 100);
    Grid copied(original);
    
    // Verify copied grid has same properties
    EXPECT_EQ(copied.numSpotPoints(), original.numSpotPoints());
    EXPECT_EQ(copied.numTimeSteps(), original.numTimeSteps());
    EXPECT_NEAR(copied.dt(), original.dt(), GRID_TOL);
    
    // Verify grid values match
    for (size_t i = 0; i < original.numSpotPoints(); ++i) {
        EXPECT_NEAR(copied.spot(i), original.spot(i), GRID_TOL);
    }
}

// ============================================================================
// Utility Methods Tests
// ============================================================================

TEST(PDEGridTest, GetBounds) {
    Grid grid(10.0, 200.0, 101, 2.0, 100);
    
    auto bounds = grid.getBounds();
    auto spotBounds = bounds.first;
    auto timeBounds = bounds.second;
    
    EXPECT_NEAR(spotBounds.first, 10.0, GRID_TOL);
    EXPECT_NEAR(spotBounds.second, 200.0, GRID_TOL);
    EXPECT_NEAR(timeBounds.first, 0.0, GRID_TOL);
    EXPECT_NEAR(timeBounds.second, 2.0, GRID_TOL);
}

TEST(PDEGridTest, AccessVectors) {
    Grid grid(0.0, 100.0, 11, 1.0, 10);
    
    // Access entire vectors
    const auto& spots = grid.spots();
    const auto& times = grid.times();
    
    EXPECT_EQ(spots.size(), 11);
    EXPECT_EQ(times.size(), 11);  // N_t + 1
    
    // Verify vector contents
    EXPECT_NEAR(spots.front(), 0.0, GRID_TOL);
    EXPECT_NEAR(spots.back(), 100.0, GRID_TOL);
    EXPECT_NEAR(times.front(), 0.0, GRID_TOL);
    EXPECT_NEAR(times.back(), 1.0, GRID_TOL);
}

// ============================================================================
// Edge Cases and Special Scenarios
// ============================================================================

TEST(PDEGridTest, MinimalGrid) {
    // Minimum valid grid: 3 spatial points, 1 time step
    Grid grid(0.0, 100.0, 3, 1.0, 1);
    
    EXPECT_EQ(grid.numSpotPoints(), 3);
    EXPECT_EQ(grid.numTimeSteps(), 1);
    
    EXPECT_NEAR(grid.spot(0), 0.0, GRID_TOL);
    EXPECT_NEAR(grid.spot(1), 50.0, GRID_TOL);
    EXPECT_NEAR(grid.spot(2), 100.0, GRID_TOL);
}

TEST(PDEGridTest, VeryFineGrid) {
    // Very fine grid for accuracy
    Grid grid(0.0, 100.0, 1001, 1.0, 1000);
    
    EXPECT_EQ(grid.numSpotPoints(), 1001);
    EXPECT_EQ(grid.numTimeSteps(), 1000);
    
    // Verify spacing
    double expected_dS = 0.1;  // (100 - 0) / 1000
    EXPECT_NEAR(grid.dS(500), expected_dS, GRID_TOL);
    
    double expected_dt = 0.001;  // 1.0 / 1000
    EXPECT_NEAR(grid.dt(), expected_dt, GRID_TOL);
}

TEST(PDEGridTest, NearZeroLowerBound) {
    // Test with S_min very close to zero (common in practice)
    double tiny = 1e-6;
    Grid grid(tiny, 100.0, 101, 1.0, 100);
    
    EXPECT_NEAR(grid.spotMin(), tiny, GRID_TOL);
    EXPECT_GT(grid.spot(0), 0.0);  // Should be positive
}

// ============================================================================
// Convergence Properties (For PDE Analysis)
// ============================================================================

TEST(PDEGridTest, GridRefinement) {
    // Test grid refinement for convergence studies
    // Create grids with N, 2N, 4N points
    
    Grid grid_N(0.0, 100.0, 51, 1.0, 50);
    Grid grid_2N(0.0, 100.0, 101, 1.0, 100);
    Grid grid_4N(0.0, 100.0, 201, 1.0, 200);
    
    // Verify refinement ratios
    EXPECT_NEAR(grid_N.dt() / grid_2N.dt(), 2.0, 1e-8);
    EXPECT_NEAR(grid_2N.dt() / grid_4N.dt(), 2.0, 1e-8);
}

TEST(PDEGridTest, StabilityCondition) {
    // For explicit Euler: Δt <= ΔS² / (2σ²S_max²)
    // This test verifies we can compute the relevant quantities
    
    double S_max = 100.0;
    double sigma = 0.2;
    
    Grid grid(0.0, S_max, 101, 1.0, 100);
    
    double dS = grid.dS(50);  // Typical spacing
    double dt = grid.dt();
    
    // CFL number for explicit Euler
    double CFL = dt / (dS * dS) * (sigma * sigma * S_max * S_max);
    
    // For stability, we want CFL <= 0.5
    // This is just informational, not a hard requirement for the grid
    EXPECT_TRUE(std::isfinite(CFL));
}

