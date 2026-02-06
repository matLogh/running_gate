#pragma once

#include <memory>
#include <vector>

class RunningGateFitter;

/**
 * @struct PeakData
 * @brief Contains fitting function and results for a single peak
 *
 * Stores the fit function, peak count regions, and vectors of fit results
 * including area, errors, and gate center positions.
 */
struct PeakData
{
    std::vector<double> area;                    ///< Fitted peak areas for each gate
    std::vector<double> area_err;                ///< Errors on fitted peak areas
    std::vector<double> running_gate_energy;     ///< Center positions of running gates
    std::vector<double> running_gate_energy_err; ///< Errors on gate center positions

    std::unique_ptr<RunningGateFitter> fitter; ///< Fitting strategy for the peak
    double peak_region_start;
    double peak_region_end;
    double peak_center_energy;
    // PeakRegion region;                         ///< Peak and background region definitions
};