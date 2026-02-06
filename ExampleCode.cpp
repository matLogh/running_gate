#include "RunningGate.h"

#include "TApplication.h"
#include "TBrowser.h"

void SampleCode()
{
    TApplication app("app", 0, 0);

    // pass symetric matrix to the running gate
    // TFile file("example_data/RawDataBin02.root", "READ");
    TFile file("example_data/example_matrix.root", "READ");
    // TH2D *mat = (TH2D *)file.Get("GammaGamma/IndividualDetectors/mGammaGammaRaw_RestvsDet0");
    TH2D *mat = (TH2D *)file.Get("mGammaGammaRaw_RestvsDet0");
    // create file for storing results

    RunningGate gate(mat);
    gate.SetSaveFileName("running_gate_output.root");

    // run the gate
    gate.SetCoincidencePeak_gausFit(543.7, 542, 546, 530, 540, 550, 560);
    gate.SetCoincidencePeak_countEvents(613.6, 612, 616, 600, 610, 620, 630);
    gate.RunGate(150.5, 156.5, 1);

    new TBrowser();
    app.Run();
}

int main(int argc, char **argv)
{
    // Run the sample code
    SampleCode();

    return 0;
}