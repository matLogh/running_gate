#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <array>
#include <utility>

#include "TheuerkaufPeak.h"
#include "data_structs.h"

class RunningGateFitter
{
  public:
    virtual std::pair<double, double> Eval(TH1 *histo, const PeakData &peak) = 0;
    virtual int Initialize(TH1 *total_gate_projection, const PeakData &peak,
                           const std::vector<std::pair<TH1 *, double>> &gate_projections) = 0;
    virtual void SaveFitToFile(TFile *file, TH1 *histo, const PeakData &peak) = 0;

  protected:
    std::string GetFuncUniqueName(const char *prefix, void *ptr);
};

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////

class EventCounter : public RunningGateFitter
{
  public:
    EventCounter(std::array<double, 2> bcg_left, std::array<double, 2> bcg_right, int bcg_order = 1)
        : fBcgLeft(bcg_left), fBcgRight(bcg_right), fBcgPolyOrder(bcg_order){};

    std::pair<double, double> Eval(TH1 *histo, const PeakData &peak) override;
    int Initialize(TH1 *total_gate_projection, const PeakData &peak,
                   const std::vector<std::pair<TH1 *, double>> &gate_projections) override;
    void SaveFitToFile(TFile *file, TH1 *histo, const PeakData &peak) override;

  private:
    double EvalBcg(double *x, double *p);

  private:
    unsigned int fBcgPolyOrder;
    std::array<double, 2> fBcgLeft;
    std::array<double, 2> fBcgRight;

    TF1 fBcgFcn;
};

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////

class GausFitter : public RunningGateFitter
{
  public:
    GausFitter(std::array<double, 2> bcg_left, std::array<double, 2> bcg_right, int bcg_order = 1);

    std::pair<double, double> Eval(TH1 *histo, const PeakData &peak) override;
    int Initialize(TH1 *total_gate_projection, const PeakData &peak,
                   const std::vector<std::pair<TH1 *, double>> &gate_projections) override;
    void SaveFitToFile(TFile *file, TH1 *histo, const PeakData &peak) override;

  private:
    double EvalFcn(double *x, double *p);

  private:
    unsigned int fBcgPolyOrder;
    std::array<double, 2> fBcgLeft;
    std::array<double, 2> fBcgRight;
    std::array<double, 2> fPeakRegion;

    bool fIsInitialized;
    TF1 fPeakFcn;
    // TF1 fBcgFcn;
};

////////////////////////////////////////
////////////////////////////////////////
////////////////////////////////////////

// class TailedGaussFitter : public RunningGateFitter
// {
//   public:
//     TailedGaussFitter(int bcg_order = 2, bool has_step = false);

//     std::pair<double, double> Eval(TH1 *histo, const PeakData &peak) override;
//     int Initialize(TH1 *total_gate_projection, const PeakData &peak) override;
//     void SaveFitToFile(TFile *file) override;

//   private:
//     TheuerkaufFitter fTheuFitter;
// };
