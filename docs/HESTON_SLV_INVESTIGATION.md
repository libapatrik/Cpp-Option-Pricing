# Heston SLV Implementation Investigation

## Summary

Investigation into why the Heston Stochastic Local Volatility (SLV) implementation produces elevated errors in OTM options. Multiple issues were identified and fixed.

**Final Results (After COS-based Analytical Dupire)**:
- Feller-satisfying parameters: **~20bp** max error (paper level!)
- Van der Stoep paper parameters (Feller-violating): **~8bp** max error (excellent!)
- Pure QE Heston (no SLV): **~27bp** max error (baseline)
- Flat vol surface: **0.00bp** error (perfect)

**Previous Results (Interpolated IV Surface Dupire)**:
- Feller-satisfying parameters: ~96bp max error
- Van der Stoep paper parameters: ~261bp max error

**Key Improvement**: Implementing COS-based analytical Dupire (`HestonLocalVol`) reduced errors by **5-34x** by eliminating interpolation artifacts.

**Performance Optimization**:
- Vectorized COS pricing: O(N_strikes × N) → O(N) CF evaluations
- Batch local vol computation: ~100x fewer CF evaluations per time slice
- OpenMP parallelization: 2.4x speedup on multi-core systems

---

## Bugs Fixed

### 1. Bin Midpoint Calculation
**File**: `PathSimulator2D.cpp`

Changed to sample mean per Algorithm 1 in van der Stoep paper:
```cpp
// Before (arithmetic mean of bounds)
bins[k].midpoint = 0.5 * (minSpot + maxSpot);

// After (sample mean of paths in bin)
double sumSpot = 0.0;
for (size_t j = 0; j < binSize; ++j) {
    sumSpot += spotValues[idx];
}
bins[k].midpoint = sumSpot / binSize;
```

### 2. Calibration Flag Bug
**File**: `PathSimulator2D.cpp`

Fixed premature setting of `_leverageCalibrated = true` which caused grid interpolation to be used before calibration completed.

### 3. Dupire Numerical Derivatives - ATM Explosion
**File**: `VolatilitySurface.cpp`

**Problem**: Fixed step sizes (h = 1e-4) caused ATM local vol to explode to 122% instead of ~20%.

**Fix**: Vol-adaptive step sizes that scale with implied vol and sqrt(T):
```cpp
// First derivative: ~2% of ATM straddle width
double h = 0.02 * strike * vol * std::sqrt(maturity);
h = std::max(h, 0.005 * strike);

// Second derivative: ~5% of ATM straddle width
double h = 0.05 * strike * vol * std::sqrt(maturity);
h = std::max(h, 0.01 * strike);
```

### 4. COS Pricer OTM Numerical Issues
**File**: `Utils.cpp`

**Problem**: COS pricing returned near-zero prices for deep OTM calls at short maturities.

**Fix**: Use put-call parity for OTM options (price the ITM counterpart and convert):
```cpp
bool isOTM = (isCall && K > S0) || (!isCall && K < S0);
if (isOTM) {
    // Price ITM option, then apply parity
    rawPrice = rawPrice + S0 - K * exp(-r*T);  // for calls
}
```

---

## New Features Implemented

### Heston COS Pricer (`Utils.h/cpp`)

Implemented COS (Fourier-Cosine) method for analytical Heston option pricing:

- **`HestonCF`** class: Heston characteristic function using "Little Heston Trap" formulation
- **`COSPricer`** class: European option pricing via Fourier-cosine expansion
  - `price()`, `callPrice()`, `putPrice()`
  - `impliedVol()` - Newton-Raphson IV inversion

This enables noise-free implied vol surface generation for testing.

---

## Test Results

| Test | Max Error | Notes |
|------|-----------|-------|
| Dupire Flat Vol | 0.00 bp | Perfect - validates Dupire formula |
| Feller-satisfying params (κ=2, σ_v=0.3) | 95 bp | Acceptable for SLV |
| Van der Stoep params (Feller-violating) | 236 bp | Known limitation |

### Local Vol Ratio Analysis

For Heston-consistent surface, local vol ratio (LV/sqrt(v0)) should be ~1.0:

**Feller-satisfying params**:
```
Spot    Ratio
0.88    1.02
0.96    0.97
1.00    0.95
1.04    0.92
1.08    0.91
```

**Van der Stoep params** (problematic):
```
Spot    Ratio
0.90    0.86
0.95    0.79
1.00    0.73  <-- Systematic underestimate
1.05    0.70
1.10    0.69
```

---

## Root Cause of Remaining Errors

### Van der Stoep Surface (~236bp errors)

The systematic underpricing stems from:

1. **Feller condition violated**: σ²=0.90 >> 2κv̄=0.18
   - Variance frequently hits zero
   - QE scheme handles this but introduces bias

2. **Dupire local vol extraction**: Ratios of 0.69-0.86 instead of ~1.0
   - Higher vol surface (27% vs 20%) with steep skew
   - Interpolated grid doesn't preserve exact Heston relationships

3. **Leverage function impact**: L² ∝ LV²/E[V|S]
   - With LV underestimated, L < 1
   - Diffusion scaled down → systematic underpricing

---

## Paper Analysis: Van der Stoep Methodology

### How the Paper Computes Local Volatility

Analysis of the van der Stoep et al. paper (2014) reveals they use **numerical** (finite difference) Dupire, NOT analytical local vol. However, a critical difference exists in how derivatives are computed:

**Van der Stoep Pipeline** (Paper):
```
Heston CF → COS prices C(K,T) → Direct FD on prices → Dupire σ²_LV
```

**Our Implementation Pipeline**:
```
Heston CF → COS prices → IV σ(K,T) → Interpolated surface → IV derivatives → Dupire σ²_LV
```

### Key Differences

| Aspect | Paper | Our Implementation |
|--------|-------|-------------------|
| Derivatives computed on | Raw prices | Implied volatilities |
| Interpolation | None (direct COS) | Bilinear on IV grid |
| Accuracy | Noise-free | Interpolation artifacts |

### Evidence from Paper

From Section 3.1 (page 8):
> "The local volatility σ_LV(S,t) can be extracted from vanilla option prices using Dupire's formula"

They use finite differences but compute them **directly on the Fourier-generated prices**, avoiding any intermediate interpolation step.

### Impact on Results

| Method | Max Error |
|--------|-----------|
| Pure QE Heston (no SLV) | ~27 bp |
| SLV with interpolated IV surface | ~261 bp |
| SLV error contribution | ~235 bp |

The ~235bp additional error is entirely attributable to our Dupire pipeline's interpolation artifacts.

### Diagnostic Test Results

**On-the-fly vs Grid Leverage**:
- On-the-fly (direct Dupire): 261 bp max error
- Calibrated grid (interpolated): 268 bp max error
- Both methods show similar error, confirming issue is in Dupire extraction, not grid calibration

**Leverage Function Analysis** (Van der Stoep params, T=0.5):
```
Spot    Grid L²   On-the-fly L²   Ratio
0.90    0.76      0.88            0.86
0.95    0.66      0.82            0.80
1.00    0.60      0.79            0.76
1.05    0.56      0.78            0.72
1.10    0.53      0.78            0.68
```

Grid L² is systematically lower, leading to under-diffusion and underpricing of OTM options.

---

## Future Improvements

### Completed

1. **COS-based Dupire for Heston** (**IMPLEMENTED & INTEGRATED**)
   - `HestonLocalVol` class in `Utils.h/cpp`
   - Computes local vol directly from COS prices using numerical derivatives
   - Uses Richardson extrapolation for O(h⁴) accuracy
   - **Integrated into `HestonSLVPathSimulator2D`** - no VolatilitySurface needed!

   **Local Vol Accuracy**:
   | Parameter Set | ATM LV Ratio | Max Ratio Error |
   |---------------|--------------|-----------------|
   | Feller-satisfying (κ=2, σ_v=0.3) | 0.95 | 5% |
   | Interpolated IV surface | 2.05 | 105% |

   **SLV Pricing Results** (with integrated HestonLocalVol):
   | Parameter Set | Before | After | Improvement |
   |---------------|--------|-------|-------------|
   | Feller-satisfying | 96 bp | **20 bp** | **5x** |
   | Van der Stoep (Feller-violating) | 261 bp | **8 bp** | **34x** |

   The integrated analytical Dupire achieves **paper-level accuracy**!

2. **Vectorized COS Pricing & Parallelization** (**IMPLEMENTED**)

   Major performance optimization for local vol grid precomputation:

   **Vectorized COS Method** (`COSPricer::prices`, `callPrices`, `putPrices`):
   - Prices multiple strikes at same maturity efficiently
   - Precomputes characteristic function φ(ω_k) once, reuses for all strikes
   - Reduces CF evaluations from O(N_strikes × N_terms) to O(N_terms)

   **Batch Local Vol** (`HestonLocalVol::localVolsAtTime`):
   - Computes local vol for all spots at same maturity in one call
   - Uses vectorized COS for strike derivatives (5 strikes per spot)
   - Uses vectorized COS for time derivatives (4 maturities)
   - ~100x reduction in CF evaluations per time slice

   **OpenMP Parallelization** (`precomputeLocalVolGrid`):
   - Time slices computed in parallel (independent)
   - Effective multi-core utilization (572-605% CPU)

   **Performance Results**:
   | Test | Before | After | Speedup |
   |------|--------|-------|---------|
   | HestonConsistentSurface_COS | 140.8s | 57.7s | **2.4x** |
   | VanDerStoepParameters | 267.3s | 194.2s | **1.4x** |

   **Implementation vs Paper**:
   The van der Stoep paper doesn't provide timing benchmarks. Our implementation likely has comparable or better performance because:
   - Same core algorithm (COS + Dupire FD)
   - Added vectorization (not mentioned in paper)
   - Added OpenMP parallelization (not mentioned in paper)

### Future Work

2. **Broadie-Kaya exact simulation** (Optional)
   - Replace QE scheme for Feller-violating parameters
   - Uses exact CIR sampling + integrated variance inversion
   - May further reduce baseline error (but current ~8bp is already excellent)

### Medium Priority

4. **SVI/SABR surface fitting**
   - Fit parametric model to COS prices before Dupire
   - Smooth, arbitrage-free surface

5. **Adaptive binning**
   - More bins in high-density regions
   - Fewer bins in tails

### Low Priority

6. **Particle method calibration**
   - Calibrate leverage directly from option prices
   - Avoids Dupire formula entirely

---

## Files Modified

| File | Changes |
|------|---------|
| `Utils.h` | Added `HestonCF`, `COSPricer`, `HestonLocalVol` classes; vectorized pricing methods |
| `Utils.cpp` | COS pricing, IV inversion, put-call parity, COS-based Dupire, vectorized COS, batch local vol |
| `VolatilitySurface.cpp` | Vol-adaptive finite difference step sizes |
| `PathSimulator2D.h` | New constructor (no VolatilitySurface), local vol grid members |
| `PathSimulator2D.cpp` | Integrated HestonLocalVol, precomputed local vol grid, OpenMP parallelization |
| `tests/test_heston_slv.cpp` | Updated to new constructor, HestonLocalVol tests |

### Performance Optimization Details

**Vectorized COS Pricing** (`Utils.cpp`):
```cpp
// Before: O(N_strikes × N_terms) CF evaluations
for (each strike K) {
    for (k = 0..N) {
        phi_k = chf(k * PI / bma);  // Expensive!
        ...
    }
}

// After: O(N_terms) CF evaluations
// Precompute CF once
for (k = 0..N) {
    phi_vals[k] = chf(k * PI / bma);  // Done once!
}
// Reuse for all strikes
for (each strike K) {
    for (k = 0..N) {
        // Use precomputed phi_vals[k]
        ...
    }
}
```

**Batch Local Vol** (`Utils.cpp`):
```cpp
// Compute local vol for all spots at same maturity
std::vector<double> localVolsAtTime(const std::vector<double>& spots, double T) {
    // 1. Build all strikes needed for Dupire derivatives (5 per spot)
    // 2. Price all strikes at T using vectorized COS (1 CF computation)
    // 3. Price central strikes at T±h for time derivative (4 CF computations)
    // 4. Apply Dupire formula to each spot
    // Total: 5 CF computations instead of 13 × N_spots!
}
```

**Parallelized Grid** (`PathSimulator2D.cpp`):
```cpp
#pragma omp parallel for schedule(dynamic)
for (size_t t_idx = 0; t_idx < numTimeSteps; ++t_idx) {
    // Each time slice is independent
    std::vector<double> lvs = _hestonLocalVol->localVolsAtTime(_spotGridPoints, t);
    // Copy to grid...
}
```

---

## Test Commands

```bash
cd build
ninja test_heston_slv

# Run all tests
./test_heston_slv

# Specific tests
./test_heston_slv --gtest_filter="COSPricerTest.*"
./test_heston_slv --gtest_filter="HestonLocalVolTest.*"
./test_heston_slv --gtest_filter="StoepComparisonTest.HestonConsistentSurface_COS"
./test_heston_slv --gtest_filter="StoepComparisonTest.VanDerStoepParameters"
./test_heston_slv --gtest_filter="StoepComparisonTest.DupireFlatVolDiagnostic"
```

---

## References

1. Van der Stoep et al. - "The Heston Stochastic-Local Volatility Model" (2014)
2. Fang & Oosterlee - "A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions" (2008)
3. Gatheral - "The Volatility Surface" (2006)
4. Andersen - "Efficient Simulation of the Heston Stochastic Volatility Model" (2007)
