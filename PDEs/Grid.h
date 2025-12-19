#ifndef CPPFM_PDEGRID_H
#define CPPFM_PDEGRID_H

#include <vector>
#include <utility>

/**
 * @class Grid
 * @brief Spatial and temporal discretization for PDE finite difference methods
 * 
 * Represents a 1D spatial grid (asset price) and 1D temporal grid (time)
 * for solving partial differential equations in computational finance.
 * 
 * Grid Structure:
 * - Spatial: S in [S_min, S_max] discretised into N_S points
 * - Temporal: t in [0, T] discretised into N_t time steps
 */
class Grid
{
public:
    /**
     * Spacing type for spatial grid
     */
    enum class SpacingType { Uniform, LogSpaced };
    
    /**
     * Constructor with full grid specification
     * @param S_min Minimum spot price (lower boundary, must be > 0)
     * @param S_max Maximum spot price (upper boundary, must be > S_min)
     * @param N_S Number of spatial grid points (must be >= 3) (default: 100)
     * @param T_max Maximum time to maturity (must be > 0) (default: 1.0)
     * @param N_t Number of time steps (must be >= 1) (default: 100)
     * @param spacingType Uniform or LogSpaced spatial grid (default: Uniform)
     */
    Grid(double S_min, double S_max, size_t N_S, double T_max, size_t N_t, SpacingType spacingType = SpacingType::Uniform);
    
    // Default operations
    Grid(const Grid& other) = default;
    Grid& operator=(const Grid& other) = default;
    ~Grid() = default;
    
     /** Clone
     * Create a deep copy of the grid
     * @return Pointer to new PDEGrid instance (caller owns memory)
     */
     Grid* clone() const;

    // Accessors - Grid Parameters
    double spotMin() const { return _S_min; }
    double spotMax() const { return _S_max; }

    double timeMax() const { return _T_max; }

    size_t numSpotPoints() const { return _N_S; }
    size_t numTimeSteps() const { return _N_t; }
    
    SpacingType spacingType() const { return _spacingType; }
    
    // Accessors - Grid Values
    const std::vector<double>& spots() const { return _spotGrid; }
    
    const std::vector<double>& times() const { return _timeGrid; } //  vector of time points [t_0=0, t_1, ..., t_{N_t}=T]
    
    double spot(size_t i) const { return _spotGrid[i]; } // spot at index i = S_i
    
    double time(size_t j) const { return _timeGrid[j]; } // time at index j = t_j
    

    // Grid Spacing Methods
    double dt() const { return _dt; }
    
    double dS(size_t i) const; // by centred difference
    
    /**
     * Get spatial step size between consecutive points
     * @param i Spatial index (0 <= i < N_S-1)
     * @return S_{i+1} - S_i
     */
    double dS_forward(size_t i) const;
    
    /**
     * Get spatial step size for backward difference
     * @param i Spatial index (1 <= i < N_S)
     * @return S_i - S_{i-1}
     */
    double dS_backward(size_t i) const;
    

    // Utility methods
    /**
     * Get grid bounds as pair
     * @return {{S_min, S_max}, {0, T_max}}
     */
    std::pair<std::pair<double, double>, std::pair<double, double>> getBounds() const;

private:
    // Grid parameters
    double _S_min, _S_max;     // Spatial domain [S_min, S_max]
    double _T_max;             // Temporal domain [0, T_max]
    size_t _N_S, _N_t;         // Number of grid points
    double _dt;                // Time step (constant)
    SpacingType _spacingType;  // Uniform or LogSpaced
    
    // Grid storage
    std::vector<double> _spotGrid;  // S_i for i = 0, 1, ..., N_S-1
    std::vector<double> _timeGrid;  // t_j for j = 0, 1, ..., N_t
    
    // Construction helpers
    void constructUniformGrid();
    void constructLogSpacedGrid();
    void constructTimeGrid();
    void validateParameters() const;
};

#endif //CPPFM_PDEGRID_H
