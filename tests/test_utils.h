//
// Test utilities for VolatilitySurface testing
//

#ifndef CPPFM_TEST_UTILS_H
#define CPPFM_TEST_UTILS_H

#include <gtest/gtest.h>


/*
PURPOSE: Shared utilities for all test files
- Custom assertion macros for financial tolerances
- Market data generators (flat, smile surfaces)
- Test fixtures with common setup
- Numerical verification helpers

TODO:
1. Add includes: gtest, vector, memory, cmath, VolatilitySurface.h, DiscountCurve.h
2. Define tolerances: PRICE_TOLERANCE=1e-6, VOLATILITY_TOLERANCE=1e-4
3. Create macros: EXPECT_PRICE_NEAR, EXPECT_VOLATILITY_NEAR
4. Build data generators: generateFlatSurface(), generateSmileSurface()
5. Add test fixture: VolatilitySurfaceTest with SetUp/TearDown
6. Add helpers: computeNumericalDerivative(), verifyPutCallParity()
*/

// TODO: Add includes and implementation

#endif // CPPFM_TEST_UTILS_H