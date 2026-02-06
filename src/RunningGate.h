#include "TF1.h"
#include "TFile.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include <iostream>
#include <vector>

#include "Fitters.h"
#include "data_structs.h"

// using namespace std;

#ifndef RUNNING_GATE
#define RUNNING_GATE

class RunningGate;

/**
 * @class RunningGate
 * @brief Performs running gate analysis on 2D histograms
 *
 * This class analyzes a 2D histogram by creating a series of 1D projections
 * (gates) along one axis and fitting peaks in each projection. It supports
 * multiple fitting methods and tracks peak evolution across gates.
 */
class RunningGate
{

  public:
    /**
     * @brief Constructor that initializes the running gate analysis
     * @param matrix Pointer to the 2D histogram to analyze
     */
    RunningGate(TH2 *matrix);

    /**
     * @brief Destructor
     */
    ~RunningGate();

    /**
     * @brief Executes the running gate analysis
     *
     * Loops through all gate positions, creates projections, fits peaks,
     * and stores results for each configured peak.
     */
    void RunGate(const double energy_gate_start, const double energy_gate_end, const unsigned int bin_gate_width);

    void SetCoincidencePeak_countEvents(const double peak_center, const double peak_en_start, const double peak_en_end,
                                        const double bcgleft_en_start, const double bcgleft_en_end,
                                        const double bcgright_en_start, const double bcgright_en_end,
                                        const unsigned short background_polynomial_order = 2);

    void SetCoincidencePeak_gausFit(const double peak_center, const double peak_en_start, const double peak_en_end,
                                    const double bcgleft_en_start, const double bcgleft_en_end,
                                    const double bcgright_en_start, const double bcgright_en_end,
                                    const unsigned short background_polynomial_order = 2);

    std::shared_ptr<TMultiGraph> GetMultigraph()
    {
        return fMultiGraph;
    };
    void SetSaveFileName(const std::string &filename)
    {
        fSaveFileName = filename;
    };

  private:
    TF1 *fcn_peak1; ///< First peak fitting function (legacy)
    TF1 *fcn_peak2; ///< Second peak fitting function (legacy)

    std::vector<PeakData> fCoincPeaks; ///< Vector of peak data objects for all configured peaks

    TH2 *fMatrix; ///< clone of the input 2D histogram
    // TH1 *fGate_slice; ///< Temporary 1D histogram for current gate projection

    /**
     * @brief Creates TGraphErrors objects for all configured peaks
     *
     * Initializes the graph vector with TGraphErrors objects that will
     * store the evolution of peak areas vs gate position.
     */
    std::vector<TGraph *> fGraphs;
    std::shared_ptr<TMultiGraph> fMultiGraph;
    std::string fSaveFileName{};

    std::shared_ptr<TMultiGraph> CreateGraphs(std::shared_ptr<TFile> write_file = nullptr);

  public:
    static TH1F *GraphErrorsToHistogram(TGraphErrors *graph, const char *histName = "hist");
    static double GetXAtMaxY(const TGraph *g, int *outIndex = nullptr);
};

#endif