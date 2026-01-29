# CppFM Optimization Research Catalog

**Purpose:** Identify all optimization opportunities in the codebase, focusing on linear algebra, vectorization, and numerical methods.

---

## Executive Summary

| Category | Current State | Main Opportunity | Estimated Gain |
|----------|---------------|------------------|----------------|
| **SIMD/Vectorization** | None | Transcendental functions (exp/sin/cos) | 4-8× |
| **Linear Algebra** | Thomas algorithm (sequential) | Cyclic reduction for parallel tridiag | 2-3× (large grids) |
| **Monte Carlo** | OpenMP path-level | Variance reduction methods | 10-100× variance |
| **COS Method** | Sequential sin/cos loops | FFT or vectorized trig | 3-5× |
| **Accelerate/BLAS** | Not used | vvexp, vvlog, vvsqrt | 4-8× |

---

## 1. SIMD / Vectorization Opportunities

### Current State
- **OpenMP**: Used for path-level parallelism in `PathSimulator2D.cpp`
- **TBB**: Used for `parallel_sort` in bin computation
- **SIMD intrinsics**: None
- **Accelerate.h**: Not used (TODO comment exists in Utils.h:20)

### Hot Loops Identified

| Location | Operation | Calls/sec | SIMD Potential |
|----------|-----------|-----------|----------------|
| `PathSimulator2D.cpp:1118-1136` | `std::exp()`, `std::log()` per path | Millions | **High** |
| `Utils.cpp:257-259` | `std::sin()` in COS summation | N×M per pricing | **High** |
| `Utils.cpp:1045-1063` | Complex exp in COSPricer | N per strike | **Medium** |
| `Solver.cpp:247-334` | FD stencil computation | N×Nt per PDE | **Medium** |

### Optimization Options

**Option A: Apple Accelerate (macOS-specific)**
```cpp
#include <Accelerate/Accelerate.h>
// Replace: for (i) y[i] = std::exp(x[i]);
// With:    vvexp(y, x, &n);  // vectorized, uses SIMD
```
- `vvexp`, `vvlog`, `vvsqrt`, `vvsin`, `vvcos`
- 4-8× speedup for transcendental functions
- Zero external dependencies on macOS

**Option B: SLEEF (Cross-platform)**
```cpp
#include <sleef.h>
// Sleef_expd4_u10(x)  // AVX2, 4 doubles
// Sleef_expd8_u10(x)  // AVX-512, 8 doubles
```
- Portable across Linux/Windows/macOS
- Requires external dependency

**Option C: Manual AVX2 intrinsics**
- Most control, most effort
- Use only for critical inner loops

### Recommendation
Start with **Accelerate.h** for macOS builds. Add compile-time switch for SLEEF on other platforms.

---

## 2. Linear Algebra Optimizations

### Current Solvers

| Algorithm | Location | Complexity | Parallelizable? |
|-----------|----------|------------|-----------------|
| Thomas (tridiagonal) | `Utils.cpp:108-191` | O(n) | **No** (sequential dependency) |
| Cubic spline | `InterpolationSchemes.cpp:166-258` | O(n) | No (uses Thomas) |
| θ-method PDE | `Solver.cpp:172-356` | O(n×Nt) | Per-timestep: No |

### Why Thomas Can't Be Parallelized
Forward elimination has recurrence: `c'[i] = f(c'[i-1])`
Back substitution has recurrence: `x[i] = f(x[i+1])`

### Alternative: Cyclic Reduction
For **large grids** (n > 1000), cyclic reduction achieves O(n/p + log p) with p processors.
- Requires restructuring the tridiagonal system
- Overkill for typical PDE grids (n ≤ 500)

### Alternative: Block Tridiagonal (2D PDEs)
For future 2D PDE extension (ADI schemes):
- Use block Thomas with dense block solves
- BLAS `dgetrf`/`dgetrs` for block LU

### Alternative: SPIKE Algorithm
- Intel MKL provides `LAPACKE_dgtsv` (tridiagonal solve)
- Parallel SPIKE for distributed systems
- Again, overkill for current grid sizes

### Recommendation
**Keep Thomas algorithm** for 1D PDEs (n ≤ 500). Focus optimization elsewhere. Consider cyclic reduction only if 2D ADI is added.

---

## 3. Monte Carlo Variance Reduction

### Current State
- Path-level OpenMP parallelism
- No variance reduction techniques

### Techniques to Implement

| Technique | Variance Reduction | Implementation Effort |
|-----------|-------------------|----------------------|
| **Antithetic variates** | 2× | Low |
| **Control variates** | 5-20× | Medium |
| **Importance sampling** | 10-100× (tail events) | High |
| **Stratified sampling** | 2-5× | Medium |
| **Quasi-Monte Carlo** | O(1/N) vs O(1/√N) | Medium |

### Proposed Structure
```
include/cppfm/montecarlo/
├── MonteCarlo.h           # Core engine
├── VarianceReduction.h    # Base + strategies
├── Estimator.h            # Confidence intervals
└── QuasiRandom.h          # Sobol, Halton sequences
```

### Control Variate Example (Heston)
Use BS price with average variance as control:
```cpp
double cv_estimate = mc_price - beta * (mc_control - known_control);
// where beta = Cov(price, control) / Var(control)
```

---

## 4. COS Method Optimizations

### Current Hotspots

**Inner summation loop** (`Utils.cpp:257-259`):
```cpp
for (size_t k = 1; k < N; ++k) {
    cdf += F_k[k] / omega[k] * std::sin(omega[k] * x_shifted);
}
```
- N = 64-512 typically
- Called M times for M evaluation points
- Total: M × N sin() calls

### Optimization Options

**Option A: Vectorized sin/cos (Accelerate)**
```cpp
std::vector<double> args(N), sines(N);
for (k) args[k] = omega[k] * x_shifted;
int n = N;
vvsin(sines.data(), args.data(), &n);  // single vectorized call
for (k) cdf += F_k[k] / omega[k] * sines[k];
```
Expected: 3-5× speedup

**Option B: FFT-based evaluation**
COS method is essentially inverse discrete cosine transform:
```cpp
// Use FFTW or Accelerate's vDSP_DCT
// Evaluate at all grid points simultaneously
```
Expected: 5-10× for large N, but setup overhead

**Option C: Chebyshev approximation**
Pre-compute polynomial approximation of ChF, then evaluate cheaply.

### Recommendation
Start with **vectorized sin/cos** (Option A). Consider FFT only for batch evaluations.

---

## 5. Characteristic Function Caching

### Current Cost
`HestonCF::operator()` (`Utils.cpp:764-803`):
- 4× `std::exp()` + 2× complex `std::sqrt()`
- Called N times per COS pricing (N = 64-512)

### Optimization: Memoization
```cpp
class HestonCF {
    mutable std::unordered_map<double, std::complex<double>> cache_;
public:
    std::complex<double> operator()(double omega) const {
        auto it = cache_.find(omega);
        if (it != cache_.end()) return it->second;
        auto val = compute(omega);
        cache_[omega] = val;
        return val;
    }
};
```

### When Useful
- Same maturity, varying strikes → same ω grid
- Already partially optimized in `COSPricer::prices()` (caches phi values)

---

## 6. Memory Layout Optimizations

### Current Issues
- `std::vector<double>` everywhere (good)
- No explicit cache alignment
- No structure-of-arrays (SoA) for SIMD

### SoA for Path Simulation
```cpp
// Current (AoS-ish):
std::vector<double> spots, variances;

// Better for SIMD:
struct alignas(64) PathData {
    double spots[8];     // AVX-512 width
    double variances[8];
};
std::vector<PathData> paths;
```

### Cache Line Alignment
```cpp
alignas(64) std::vector<double> data;  // 64-byte alignment for cache lines
```

---

## 7. Bessel Function Optimization

### Current Implementation
`modifiedBesselI()` (`Utils.cpp:718-753`):
- Series expansion, 50+ iterations typical
- Called in Broadie-Kaya integrated variance ChF

### Options
1. **Boost.Math**: Already have Boost dependency, `boost::math::cyl_bessel_i`
2. **Lookup table**: Pre-compute for common (order, argument) pairs
3. **Asymptotic expansion**: For large arguments

---

## 8. Priority Ranking

| Priority | Optimization | Impact | Effort | Files |
|----------|-------------|--------|--------|-------|
| **1** | Accelerate.h for transcendentals | 4-8× | Low | Utils.cpp, PathSimulator2D.cpp |
| **2** | Vectorized COS summation | 3-5× | Medium | Utils.cpp |
| **3** | Monte Carlo variance reduction | 10-100× variance | High | New module |
| **4** | ChF memoization | 2-3× | Low | Utils.cpp |
| **5** | SoA memory layout for paths | 2× | Medium | PathSimulator2D.cpp |
| **6** | Cyclic reduction tridiagonal | 2-3× (large n) | High | Utils.cpp |

---

## 9. Implementation Phases

### Phase A: Quick Wins
- [ ] Add Accelerate.h for vvexp/vvlog/vvsqrt in PathSimulator2D
- [ ] Add ChF caching in HestonCF

### Phase B: COS Method
- [ ] Vectorize sin/cos in COS summation
- [ ] Benchmark FFT alternative for batch pricing

### Phase C: Monte Carlo Module
- [ ] Design montecarlo/ directory structure
- [ ] Implement antithetic variates
- [ ] Implement control variates for Heston
- [ ] Add Sobol sequences for QMC

### Phase D: Advanced (future)
- [ ] SIMD path simulation with SoA layout
- [ ] Cross-platform vectorization (SLEEF)
- [ ] 2D PDE with ADI and block solvers

---

## 10. Verification

After each optimization:
1. Run `ctest` - all tests must pass
2. Benchmark before/after with `test_heston_scheme` timing
3. Verify numerical accuracy within tolerance

---

## Key Files Reference

| Component | Header | Implementation |
|-----------|--------|----------------|
| Thomas algorithm | `include/cppfm/utils/Utils.h:62-101` | `src/utils/Utils.cpp:108-191` |
| COS method | `include/cppfm/utils/Utils.h:114-254` | `src/utils/Utils.cpp:218-592` |
| Heston ChF | `include/cppfm/utils/Utils.h` | `src/utils/Utils.cpp:764-803` |
| Path simulation | `include/cppfm/simulators/PathSimulator2D.h` | `src/simulators/PathSimulator2D.cpp` |
| PDE solver | `include/cppfm/pde/Solver.h` | `src/pde/Solver.cpp` |
| Interpolation | `include/cppfm/utils/InterpolationSchemes.h` | `src/utils/InterpolationSchemes.cpp` |
