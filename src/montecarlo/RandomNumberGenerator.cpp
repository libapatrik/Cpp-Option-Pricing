#include <cppfm/montecarlo/RandomNumberGenerator.h>

// ============================================================================
// PCG32 GENERATOR
// ============================================================================

PCG32Generator::PCG32Generator(uint64_t seed, uint64_t stream)
    : _rng(seed, stream)
    , _savedState(seed, stream)
    , _seed(seed)
    , _stream(stream)
{}

PCG32Generator* PCG32Generator::clone() const {
    return new PCG32Generator(*this);
}

double PCG32Generator::uniform() {
    return std::generate_canonical<double, 53>(_rng);
}

void PCG32Generator::seed(uint64_t s) {
    _seed = s;
    _rng.seed(s, _stream);
}

void PCG32Generator::saveState() {
    _savedState = _rng;
}

void PCG32Generator::restoreState() {
    _rng = _savedState;
}

std::vector<std::unique_ptr<RandomNumberGenerator>>
PCG32Generator::createParallelStreams(size_t n) const {
    std::vector<std::unique_ptr<RandomNumberGenerator>> streams;
    streams.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        // each thread gets same seed but different stream
        streams.push_back(std::make_unique<PCG32Generator>(_seed, _stream + i + 1));
    }
    return streams;
}

// ============================================================================
// MERSENNE TWISTER GENERATOR
// ============================================================================

MersenneTwisterGenerator::MersenneTwisterGenerator(uint64_t seed)
    : _rng(seed)
    , _savedState(seed)
    , _seed(seed)
{}

MersenneTwisterGenerator* MersenneTwisterGenerator::clone() const {
    return new MersenneTwisterGenerator(*this);
}

double MersenneTwisterGenerator::uniform() {
    return std::generate_canonical<double, 53>(_rng);
}

void MersenneTwisterGenerator::seed(uint64_t s) {
    _seed = s;
    _rng.seed(s);
}

void MersenneTwisterGenerator::saveState() {
    _savedState = _rng;
}

void MersenneTwisterGenerator::restoreState() {
    _rng = _savedState;
}

std::vector<std::unique_ptr<RandomNumberGenerator>>
MersenneTwisterGenerator::createParallelStreams(size_t n) const {
    std::vector<std::unique_ptr<RandomNumberGenerator>> streams;
    streams.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        // each thread gets different seed offset
        streams.push_back(std::make_unique<MersenneTwisterGenerator>(_seed + i + 1));
    }
    return streams;
}
