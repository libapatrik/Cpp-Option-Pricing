

#ifndef CPPFM_RANDOMNUMBERGENERATOR_H
#define CPPFM_RANDOMNUMBERGENERATOR_H


#include <cppfm/third_party/pcg/pcg_random.hpp>
#include <cppfm/utils/Utils.h>
#include <random>
#include <vector>
#include <memory>
#include <cstdint>


// ============================================================================
// BASE CLASS
// ============================================================================

class RandomNumberGenerator {
public:
    virtual ~RandomNumberGenerator() = default;
    virtual RandomNumberGenerator* clone() const = 0;

    // sampling
    virtual double uniform() = 0;
    double normal();

    // seeding
    virtual void seed(uint64_t s) = 0;

    // state management (for antithetic replay)
    virtual void saveState() = 0;
    virtual void restoreState() = 0;

    // parallel streams
    virtual std::vector<std::unique_ptr<RandomNumberGenerator>>
        createParallelStreams(size_t n) const = 0;

protected:
    RandomNumberGenerator() = default;
};

// ============================================================================
// PCG32 GENERATOR
// ============================================================================

class PCG32Generator : public RandomNumberGenerator {
public:
    explicit PCG32Generator(uint64_t seed, uint64_t stream = 0);
    PCG32Generator* clone() const override;

    double uniform() override;
    void seed(uint64_t s) override;
    void saveState() override;
    void restoreState() override;

    std::vector<std::unique_ptr<RandomNumberGenerator>>
        createParallelStreams(size_t n) const override;

private:
    pcg32 _rng;
    pcg32 _savedState;
    uint64_t _seed;
    uint64_t _stream;
};

// ============================================================================
// MERSENNE TWISTER GENERATOR
// ============================================================================

class MersenneTwisterGenerator : public RandomNumberGenerator {
public:
    explicit MersenneTwisterGenerator(uint64_t seed);
    MersenneTwisterGenerator* clone() const override;

    double uniform() override;
    void seed(uint64_t s) override;
    void saveState() override;
    void restoreState() override;

    std::vector<std::unique_ptr<RandomNumberGenerator>>
        createParallelStreams(size_t n) const override;

private:
    std::mt19937_64 _rng;
    std::mt19937_64 _savedState;
    uint64_t _seed;
};

// ============================================================================
// INLINE IMPLEMENTATIONS
// ============================================================================

inline double RandomNumberGenerator::normal() {
    return Utils::inverseNormalCDF(uniform());
}

#endif // CPPFM_RANDOMNUMBERGENERATOR_H
