#include "RunningGate.h"
#include <cmath>
#include <limits>

void RunningGate::RunGate(const double energy_gate_start, const double energy_gate_end,
                          const unsigned int bin_gate_width)
{

    // supress ROOT fit info
    gPrintViaErrorHandler = kTRUE;
    gErrorIgnoreLevel = 6001;

    if (fCoincPeaks.size() == 0)
    {
        std::cerr << "No function/fit region was set!" << std::endl;
        std::cerr << "Nothing to be done..." << std::endl;
        return;
    }

    const int running_gate_start_bin = fMatrix->GetXaxis()->FindBin(energy_gate_start);
    const int running_gate_end_bin = fMatrix->GetXaxis()->FindBin(energy_gate_end);
    const double bin_width = fMatrix->GetXaxis()->GetBinWidth(1);

    std::shared_ptr<TFile> fout{nullptr};
    if (!fSaveFileName.empty())
    {
        fout = std::shared_ptr<TFile>(new TFile(fSaveFileName.c_str(), "RECREATE"));
    }

    // make full gate projection
    auto h_proj_total = fMatrix->ProjectionY("proj_gate_total", running_gate_start_bin, running_gate_end_bin, "e");
    if (fout)
    {
        fout->cd();
        h_proj_total->Write();
    }

    // create projections
    std::vector<std::pair<TH1 *, double>> gate_projections;
    for (int bin_i = running_gate_start_bin; bin_i <= running_gate_end_bin; bin_i++)
    {

        std::cout << "STATUS: "
                  << "Slicing gg-matrix "
                  << "                       \r" << std::flush;

        auto h_slice = fMatrix->ProjectionY("proj_gate_", bin_i, bin_i + bin_gate_width);

        const double center_energy =
            (h_slice->GetXaxis()->GetBinLowEdge(bin_i) + h_slice->GetXaxis()->GetBinUpEdge(bin_i + bin_gate_width)) /
            2.0;

        h_slice->SetName(Form("proj_gate_%d", center_energy));
        h_slice->SetTitle(
            Form("Projection of bins %i-%i with center energy %.2f", bin_i, bin_i + bin_gate_width, center_energy));

        gate_projections.push_back(std::make_pair(h_slice, center_energy));
    }

    for (auto &peak : fCoincPeaks)
    {
        std::cout << "STATUS: "
                  << "Initializing fitters "
                  << "                       \r" << std::flush;
        peak.fitter->Initialize(h_proj_total, peak, gate_projections);
    }

    int counter = 0;
    for (auto &p : gate_projections)
    {

        std::cout << "STATUS: "
                  << "Calculating projection volume " << (double)(counter++) / (double)gate_projections.size() * 100
                  << "% done    \r" << std::flush;

        TH1 *h_slice = p.first;
        const double center_energy = p.second;

        for (auto &peak : fCoincPeaks)
        {
            auto res = peak.fitter->Eval(h_slice, peak);
            peak.area.push_back(res.first);
            peak.area_err.push_back(res.second);
            peak.running_gate_energy.push_back(center_energy);
            peak.running_gate_energy_err.push_back(bin_width / 2.);

            if (fout)
            {
                peak.fitter->SaveFitToFile(fout.get(), h_slice, peak);
            }
        }

        if (fout)
        {
            std::string folder_name{"projections"};
            auto *dir = fout->GetDirectory(folder_name.c_str());
            if (!dir)
            {
                dir = fout->mkdir(folder_name.c_str());
            }
            dir->cd();
            h_slice->Write();
            fout->cd();
        }
    }

    for (auto &gp : gate_projections)
    {
        delete gp.first;
    }

    std::cout << "STATUS: "
              << "Drawing and saving to file (if selected)"
              << "                         \r" << std::flush;

    CreateGraphs(fout);
    if (fout)
    {
        fout->Close();
    }
}

void RunningGate::SetCoincidencePeak_countEvents(const double peak_center, const double peak_en_start,
                                                 const double peak_en_end, const double bcgleft_en_start,
                                                 const double bcgleft_en_end, const double bcgright_en_start,
                                                 const double bcgright_en_end,
                                                 const unsigned short background_polynomial_order)
{
    PeakData pd;
    pd.peak_center_energy = peak_center;
    pd.fitter = std::make_unique<EventCounter>(std::array<double, 2>{bcgleft_en_start, bcgleft_en_end},
                                               std::array<double, 2>{bcgright_en_start, bcgright_en_end},
                                               background_polynomial_order);

    pd.peak_region_start = peak_en_start;
    pd.peak_region_end = peak_en_end;

    fCoincPeaks.push_back(std::move(pd));
}

void RunningGate::SetCoincidencePeak_gausFit(const double peak_center, const double peak_en_start,
                                             const double peak_en_end, const double bcgleft_en_start,
                                             const double bcgleft_en_end, const double bcgright_en_start,
                                             const double bcgright_en_end,
                                             const unsigned short background_polynomial_order)
{
    PeakData pd;
    pd.peak_center_energy = peak_center;
    pd.fitter = std::make_unique<GausFitter>(std::array<double, 2>{bcgleft_en_start, bcgleft_en_end},
                                             std::array<double, 2>{bcgright_en_start, bcgright_en_end},
                                             background_polynomial_order);

    pd.peak_region_start = peak_en_start;
    pd.peak_region_end = peak_en_end;

    fCoincPeaks.push_back(std::move(pd));
}

RunningGate::RunningGate(TH2 *m)
{
    fMatrix = dynamic_cast<TH2 *>(m->Clone(Form("%s_clone", m->GetName())));
    fMatrix->SetDirectory(0);
}

RunningGate::~RunningGate()
{
    delete fMatrix;
}

std::shared_ptr<TMultiGraph> RunningGate::CreateGraphs(std::shared_ptr<TFile> write_file)
{
    fMultiGraph = std::make_shared<TMultiGraph>();
    fMultiGraph->SetName("running_gate_multigraph");
    fMultiGraph->SetTitle("Running Gate Analysis;Gate center energy [keV];Peak area [counts]");

    static unsigned int ___color_index = 0;
    static std::vector<std::string> ___color_array{
        "#2196f3", "#f44336", "#4caf50", "#ff9800", "#000000", "#673ab7", "#795548", "#e91e63", "#8bc34a", "#ff5722",
        "#607d8b", "#9c27b0", "#00bcd4", "#ffeb3b", "#3f51b5", "#03a9f4", "#009688", "#cddc39", "#ffc107", "#9e9e9e"};

    for (auto &peak : fCoincPeaks)
    {
        // data graph
        TGraphErrors *g_data =
            new TGraphErrors(peak.running_gate_energy.size(), &(peak.running_gate_energy[0]), &(peak.area[0]),
                             &(peak.running_gate_energy_err[0]), &(peak.area_err[0]));
        g_data->SetName(Form("data_graph_peak_%.2f", peak.peak_center_energy));
        {
            Color_t c = TColor::GetColor(___color_array[___color_index++].c_str());
            g_data->SetLineColor(c);
            g_data->SetMarkerColor(c);
            g_data->SetTitle(Form("Data peak %.2f", peak.peak_center_energy));
            if (___color_index >= ___color_array.size())
            {
                ___color_index = 0;
            }
            g_data->SetMarkerStyle(20); // Small filled circle instead of open circle
            g_data->SetMarkerSize(0.8);
            fGraphs.emplace_back(g_data);
        }
        fMultiGraph->Add(g_data, "EP");

        // fit the data and make fit graph with confidence intervals
        {
            TheuerkaufFitter tfit(g_data->GetXaxis()->GetXmin(), g_data->GetXaxis()->GetXmax());
            tfit.AddPeak(this->GetXAtMaxY(g_data), true, true, false);

            tfit.SetBackgroundPoly(2);

            // create temporary histogram from graph to fit the data, this is needed because TheuerkaufFitter works with
            // histograms
            auto ghist = RunningGate::GraphErrorsToHistogram(g_data, "rg_fit_hist");
            tfit.Fit(ghist, "");

            if (write_file)
            {
                gROOT->SetBatch(kTRUE);
                auto c = tfit.Analyze(ghist);
                c->cd(1);
                gPad->SetLogy(false);
                gPad->Modified();
                c->Update();
                c->Write(Form("fit_coinc_%.2f", peak.peak_center_energy));
                gROOT->SetBatch(kFALSE);
            }
            // delete ghist;

            // now we need to make a tgrapherror from the fit curve to have a nice plot
            const int npoints = g_data->GetN() * 10;
            const double step = (g_data->GetPointX(g_data->GetN() - 1) - g_data->GetPointX(0)) / (double)(npoints);
            double x[npoints], y[npoints], ex[npoints], ey[npoints];
            for (int i = 0; i < npoints; i++)
            {
                x[i] = g_data->GetPointX(0) + i * step;
                ex[i] = 0;
                y[i] = tfit.GetTotalFunction()->Eval(x[i]);
            }
            tfit.GetConfidenceIntervals(npoints, x, ey, 0.95, false);

            TGraphErrors *fit_graph = new TGraphErrors(npoints, x, y, ex, ey);
            fit_graph->SetName(Form("fit_graph_peak_%.2f", peak.peak_center_energy));
            fit_graph->SetTitle(Form("Fit & 95%% CI for peak %.2f", peak.peak_center_energy));
            fit_graph->SetLineColor(g_data->GetLineColor());
            fit_graph->SetFillStyle(3001);
            fit_graph->SetFillColor(g_data->GetLineColor());
            fGraphs.emplace_back(fit_graph);

            fMultiGraph->Add(fit_graph, "4C");

            {
                double x[2] = {tfit.GetPeak(0)->GetPos(), tfit.GetPeak(0)->GetPos()};
                double y[2] = {0, tfit.GetTotalFunction()->Eval(
                                      tfit.GetPeak(0)->GetPos())}; // arbitrary y values to make a vertical line
                TGraph *g_centroid_marker = new TGraph(2, x, y);
                g_centroid_marker->SetName(Form("centroid_marker_peak_%.2f", peak.peak_center_energy));
                g_centroid_marker->SetTitle(Form("Fit centroid for peak %.2f", peak.peak_center_energy));
                g_centroid_marker->SetLineColor(g_data->GetLineColor());
                g_centroid_marker->SetLineWidth(2);
                g_centroid_marker->SetLineStyle(7); // dashed line
                fGraphs.emplace_back(g_centroid_marker);

                fMultiGraph->Add(g_centroid_marker, "L");
            }
        }
    }

    auto *c_mg = new TCanvas("canvas_multigraph", "Running Gate results", 1400, 800);
    fMultiGraph->Draw("AEP");
    c_mg->BuildLegend();

    if (write_file)
    {
        write_file->cd();
        // write_file->mkdir("graphs");
        // write_file->cd("graphs");
        for (auto &g : fGraphs)
        {
            g->Write();
        }
        // write_file->cd("..");
        c_mg->Write();
        fMultiGraph->Write();
    }

    fGraphs.clear(); // we don't want to keep the graphs in memory, they are owned by the multigraph and will be deleted
                     // with it

    return fMultiGraph;
}

double RunningGate::GetXAtMaxY(const TGraph *g, int *outIndex)
{
    if (!g || g->GetN() == 0)
        return std::numeric_limits<double>::quiet_NaN();

    double x = 0.0, y = 0.0;
    double maxY = -std::numeric_limits<double>::infinity();
    double xAtMax = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < g->GetN(); ++i)
    {
        g->GetPoint(i, x, y);
        if (std::isnan(y))
            continue;
        if (y > maxY)
        {
            maxY = y;
            xAtMax = x;
            if (outIndex)
                *outIndex = i;
        }
    }
    return xAtMax;
}

TH1F *RunningGate::GraphErrorsToHistogram(TGraphErrors *graph, const char *histName)
{
    int nPoints = graph->GetN();
    if (nPoints < 2)
    {
        std::cerr << "Graph must have at least 2 points" << std::endl;
        return nullptr;
    }

    double *xValues = graph->GetX();
    double *yValues = graph->GetY();
    double *yErrors = graph->GetEY();

    // Calculate bin width from first two points
    double binWidth = xValues[1] - xValues[0];

    // Create bin edges (centered on the point values)
    std::vector<double> binEdges;
    binEdges.push_back(xValues[0] - binWidth / 2.0);
    for (int i = 1; i < nPoints; i++)
    {
        binEdges.push_back(xValues[i] - binWidth / 2.0);
    }
    binEdges.push_back(xValues[nPoints - 1] + binWidth / 2.0);

    // Create histogram with variable bin edges
    TH1F *hist =
        new TH1F(histName, histName, nPoints, xValues[0] - binWidth / 2.0, xValues[nPoints - 1] + binWidth / 2.0);
    hist->SetDirectory(0);
    // Fill histogram with Y values and errors
    for (int i = 0; i < nPoints; i++)
    {
        hist->SetBinContent(i + 1, yValues[i]);
        hist->SetBinError(i + 1, yErrors[i]);
    }

    return hist;
}