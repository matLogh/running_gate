#include "Fitters.h"

#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TROOT.h>

std::string RunningGateFitter::GetFuncUniqueName(const char *prefix, void *ptr)
{
    static int ___num = 0;
    static std::mutex __getfunctuniquename__mtx;
    std::lock_guard<std::mutex> lock(__getfunctuniquename__mtx);

    std::ostringstream name;
    name << prefix << "_" << ptr << "_" << ++___num;

    return name.str();
}

double EventCounter::EvalBcg(double *x, double *p)
{
    if (!(x[0] >= fBcgLeft[0] && x[0] <= fBcgLeft[1] || (x[0] >= fBcgRight[0] && x[0] <= fBcgRight[1])))
    {
        TF1::RejectPoint();
        return 0;
    }

    double res = 0;
    for (int i = 0; i <= fBcgPolyOrder; i++)
    {
        res += p[i] * TMath::Power(x[0], i);
    }

    return res;
}

int EventCounter::EventCounter::Initialize(TH1 *total_gate_projection, const PeakData &peak,
                                           const std::vector<std::pair<TH1 *, double>> &gate_projections)
{
    // round up the background boundaries to nearest bin edges
    fBcgLeft[0] =
        total_gate_projection->GetXaxis()->GetBinLowEdge(total_gate_projection->GetXaxis()->FindBin(fBcgLeft[0]));
    fBcgLeft[1] =
        total_gate_projection->GetXaxis()->GetBinUpEdge(total_gate_projection->GetXaxis()->FindBin(fBcgLeft[1]));
    fBcgRight[0] =
        total_gate_projection->GetXaxis()->GetBinLowEdge(total_gate_projection->GetXaxis()->FindBin(fBcgRight[0]));
    fBcgRight[1] =
        total_gate_projection->GetXaxis()->GetBinUpEdge(total_gate_projection->GetXaxis()->FindBin(fBcgRight[1]));

    assert((fBcgLeft[0] < fBcgLeft[1]) && (fBcgRight[0] < fBcgRight[1]) && (fBcgLeft[1] < fBcgRight[0]) &&
           "Background regions are not properly defined or are overlapping!");

    fBcgFcn = TF1("bcg_fit", this, &EventCounter::EvalBcg, fBcgLeft[0], fBcgRight[1], fBcgPolyOrder + 1);
    auto fname = this->GetFuncUniqueName(fBcgFcn.GetName(), &fBcgFcn);
    fBcgFcn.SetName(fname.c_str());
    return 0;
}

std::pair<double, double> EventCounter::Eval(TH1 *histo, const PeakData &peak)
{
    // Calculate background using polynomial fit
    TFitResultPtr fit_res = histo->Fit(&fBcgFcn, "RSQNFC");
    // std::cout << "bcg fit end" << std::endl;

    const int peak_start_bin = histo->GetXaxis()->FindBin(peak.peak_region_start);
    const int peak_end_bin = histo->GetXaxis()->FindBin(peak.peak_region_end);

    const double peak_start = histo->GetXaxis()->GetBinLowEdge(peak_start_bin);
    const double peak_end = histo->GetXaxis()->GetBinUpEdge(peak_end_bin);
    const double bin_width = histo->GetXaxis()->GetBinWidth(1);

    const double total_bcg = fBcgFcn.Integral(peak_start, peak_end) / bin_width;
    const double total_bcg_err = fBcgFcn.IntegralError(peak_start, peak_end, fBcgFcn.GetParameters(),
                                                       fit_res->GetCovarianceMatrix().GetMatrixArray()) /
                                 bin_width;

    double area = 0;
    double area_err = 0;

    for (int i = peak_start_bin; i <= peak_end_bin; i++)
    {
        area += histo->GetBinContent(i);
    }

    area_err = TMath::Sqrt(area);

    // subtract background
    area = area - total_bcg;
    area_err = TMath::Sqrt(area_err * area_err + total_bcg_err * total_bcg_err);

    return std::pair<double, double>(area, area_err);
}

void EventCounter::SaveFitToFile(TFile *file, TH1 *histo, const PeakData &peak)
{
    // prepare
    std::string folder_name = Form("peak_%0.2f", peak.peak_center_energy);
    auto *dir = file->GetDirectory(folder_name.c_str());
    if (!dir)
    {
        dir = file->mkdir(folder_name.c_str());
    }
    dir->cd();
    gROOT->SetBatch(true);

    // divide bcg function into 2 for nice drawing
    TF1 fBcgFcnTemp(Form("%s_left", fBcgFcn.GetName()), Form("pol%i", fBcgPolyOrder + 1),
                    std::min(fBcgLeft[0], fBcgRight[0]), std::max(fBcgLeft[1], fBcgRight[1]));
    fBcgFcnTemp.SetParameters(fBcgFcn.GetParameters());

    auto *htemp1 = dynamic_cast<TH1 *>(histo->Clone(Form("%s_clone1", histo->GetName())));
    auto *htemp2 = dynamic_cast<TH1 *>(histo->Clone(Form("%s_clone2", histo->GetName())));
    htemp1->Reset();
    htemp2->Reset();
    for (int bin = htemp1->GetXaxis()->FindBin(peak.peak_region_start);
         bin <= htemp1->GetXaxis()->FindBin(peak.peak_region_end); bin++)
    {
        htemp1->SetBinContent(bin, histo->GetBinContent(bin));
    }

    htemp1->SetName(Form("%s_peakRegion", histo->GetName()));
    htemp1->SetFillColor(kRed);
    htemp1->SetFillStyle(3003);
    htemp1->SetLineWidth(0);

    for (int bin = htemp2->GetXaxis()->FindBin(fBcgLeft[0]); bin <= htemp2->GetXaxis()->FindBin(fBcgLeft[1]); bin++)
    {
        htemp2->SetBinContent(bin, histo->GetBinContent(bin));
    }
    for (int bin = htemp2->GetXaxis()->FindBin(fBcgRight[0]); bin <= htemp2->GetXaxis()->FindBin(fBcgRight[1]); bin++)
    {
        htemp2->SetBinContent(bin, histo->GetBinContent(bin));
    }
    htemp2->SetName(Form("%s_peakRegion", histo->GetName()));
    htemp2->SetFillColor(kBlack);
    htemp2->SetFillStyle(3003);
    htemp2->SetLineWidth(0);

    histo->GetXaxis()->SetRangeUser(std::min(peak.peak_region_start, fBcgLeft[0]),
                                    std::max(peak.peak_region_end, fBcgRight[1]));
    TCanvas c(Form("c_%s", histo->GetName()), histo->GetTitle(), 800, 600);
    histo->Draw();
    fBcgFcnTemp.Draw("same");
    htemp1->Draw("same");
    htemp2->Draw("same");
    c.Write();
    // cleanup
    histo->GetXaxis()->SetRangeUser(histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
    gROOT->SetBatch(false);
    file->cd("..");

    delete htemp1;
    delete htemp2;
}

///////////////////////////////////////// GausFitter /////////////////////////////////////////

double GausFitter::EvalFcn(double *x, double *p)
{

    if (!(x[0] >= fBcgLeft[0] && x[0] <= fBcgLeft[1] || (x[0] >= fBcgRight[0] && x[0] <= fBcgRight[1]) ||
          (x[0] >= fPeakRegion[0] && x[0] <= fPeakRegion[1])))

    {
        TF1::RejectPoint();
        return 0;
    }
    static const double sqrt_2pi = TMath::Sqrt(2. * TMath::Pi());
    double volume_normalization = p[0] / (sqrt_2pi * p[2]);
    double peak =
        volume_normalization * exp(-0.5 * TMath::Power((x[0] - p[1]) / p[2], 2.) / (sqrt(2 * TMath::Pi()) * p[2]));

    double bcg = 0;
    for (int i = 0; i <= fBcgPolyOrder; i++)
    {
        bcg += p[3 + i] * TMath::Power(x[0], i);
    }
    return peak + bcg;
}

GausFitter::GausFitter(std::array<double, 2> bcg_left, std::array<double, 2> bcg_right, int bcg_order)
    : fBcgLeft(bcg_left), fBcgRight(bcg_right), fBcgPolyOrder(bcg_order), fIsInitialized(false)
{
    fPeakFcn = TF1("peak_fcn", this, &GausFitter::EvalFcn, 0, 1000, 3 + fBcgPolyOrder + 1);
    fPeakFcn.SetParName(0, "Volume");
    fPeakFcn.SetParName(1, "Mean");
    fPeakFcn.SetParName(2, "Sigma");
    for (int i = 0; i <= fBcgPolyOrder; i++)
    {
        fPeakFcn.SetParName(3 + i, Form("pol_%i", i));
    }
}

std::pair<double, double> GausFitter::Eval(TH1 *histo, const PeakData &peak)
{
    assert(("GausFitter not initialized!" && fIsInitialized));

    TFitResultPtr fit_res = histo->Fit(&fPeakFcn, "RSQMEN");
    if (!fit_res->IsValid())
    {
        std::cerr << "Fit failed for GausFitter" << std::endl;
        return std::pair<double, double>(0, 0);
    }

    // bin width modulates the volume!
    const double bin_width_factor = 1. / histo->GetXaxis()->GetBinWidth(1);

    const double area = fPeakFcn.GetParameter(0) * bin_width_factor;
    const double area_err = fPeakFcn.GetParError(0) * bin_width_factor;

    return std::pair<double, double>(area, area_err);
}

int GausFitter::Initialize(TH1 *total_gate_projection, const PeakData &peak,
                           const std::vector<std::pair<TH1 *, double>> &gate_projections)
{
    fPeakRegion = {peak.peak_region_start, peak.peak_region_end};
    assert(gate_projections.size() > 0 && "GausFitter Initialize requires at least one gate projection!");

    const int bin_start = total_gate_projection->FindBin(peak.peak_region_start);
    const int bin_end = total_gate_projection->FindBin(peak.peak_region_end);
    const double bin_width = total_gate_projection->GetXaxis()->GetBinWidth(1);

    // we need to estimate initial fit parameters for the peak in projection - so we need to find projection with the
    // most counts in the peak region

    auto it = std::max_element(gate_projections.begin(), gate_projections.end(),
                               [&](const std::pair<TH1 *, double> &a, const std::pair<TH1 *, double> &b) {
                                   const double a_counts = a.first->Integral(bin_start, bin_end);
                                   const double b_counts = b.first->Integral(bin_start, bin_end);
                                   return a_counts < b_counts;
                               });
    auto h_proj = it->first;

    const double fit_range_start = std::min({peak.peak_region_start, fBcgLeft[0], fBcgRight[0]});
    const double fit_range_end = std::max({peak.peak_region_end, fBcgLeft[1], fBcgRight[1]});

    fPeakFcn.SetRange(fit_range_start, fit_range_end);
    fPeakFcn.SetNpx((fit_range_end - fit_range_start) / bin_width * 10.);

    // estimate parameters of the peak
    // estimate bcg
    const int est_bcg = (h_proj->GetBinContent(bin_start) + h_proj->GetBinContent(bin_end)) / 2.;

    fPeakFcn.SetParameter(3, est_bcg); // set background level
    for (int i = 1; i <= fBcgPolyOrder; i++)
    {
        fPeakFcn.SetParameter(3 + i, 0); // set higher order background terms to 0
    }

    // estimate area
    const double est_area = h_proj->Integral(bin_start, bin_end) - est_bcg * (bin_end - bin_start + 1);
    fPeakFcn.SetParameter(0, est_area);
    fPeakFcn.SetParLimits(0, 0, est_area * 100);

    // find max bin for mean estimation
    int max_bin = bin_start;

    for (int i = bin_start; i <= bin_end; i++)
    {
        if (h_proj->GetBinContent(i) > h_proj->GetBinContent(max_bin))
        {
            max_bin = i;
        }
    }
    // fPeakFcn.SetParameter(1, total_gate_projection->GetXaxis()->GetBinCenter(max_bin));
    fPeakFcn.SetParameter(1, peak.peak_center_energy);
    fPeakFcn.SetParLimits(1, peak.peak_region_start, peak.peak_region_end);

    // estimate sigma
    const double peak_max = h_proj->GetBinContent(max_bin) - est_bcg;
    static const double sqrt_2pi = TMath::Sqrt(2. * TMath::Pi());
    const double est_sigma = est_area / (peak_max * sqrt_2pi) * bin_width;
    fPeakFcn.SetParameter(2, est_sigma);
    fPeakFcn.SetParLimits(2, 0, est_sigma * 2.);

    TFitResultPtr fit_res = h_proj->Fit(&fPeakFcn, "RSQMEN");
    if (!fit_res->IsValid())
    {
        std::cerr << "Fit failed for GausFitter for peak " << peak.peak_region_start << " - " << peak.peak_region_end
                  << std::endl;
        std::runtime_error("GausFitter initialization fit failed! Check your parameters.");
        return -1;
    }

    // now fix peak parameters - mean and sigma
    fPeakFcn.SetParLimits(0, 0, fPeakFcn.GetParameter(0) * 2);
    fPeakFcn.FixParameter(1, fPeakFcn.GetParameter(1));
    fPeakFcn.FixParameter(2, fPeakFcn.GetParameter(2));

    fIsInitialized = true;
    return 0;
}

void GausFitter::SaveFitToFile(TFile *file, TH1 *histo, const PeakData &peak)
{
    // prepare
    std::string folder_name = Form("peak_%0.2f", peak.peak_center_energy);
    auto *dir = file->GetDirectory(folder_name.c_str());
    if (!dir)
    {
        dir = file->mkdir(folder_name.c_str());
    }
    dir->cd();
    gROOT->SetBatch(true);

    // for nice drawing
    TF1 fPeakFcnTemp(Form("%s_peaktmp", fPeakFcn.GetName()), Form("gausn(0) + pol%i(3)", fBcgPolyOrder),
                     std::min({fPeakRegion[0], fBcgLeft[0], fBcgRight[0]}),
                     std::max({fPeakRegion[1], fBcgLeft[1], fBcgRight[1]}));
    fPeakFcnTemp.SetParameters(fPeakFcn.GetParameters());
    fPeakFcnTemp.SetLineColor(kRed);

    auto *htemp1 = dynamic_cast<TH1 *>(histo->Clone(Form("%s_clone1", histo->GetName())));
    auto *htemp2 = dynamic_cast<TH1 *>(histo->Clone(Form("%s_clone2", histo->GetName())));
    htemp1->Reset();
    htemp2->Reset();
    for (int bin = htemp1->GetXaxis()->FindBin(peak.peak_region_start);
         bin <= htemp1->GetXaxis()->FindBin(peak.peak_region_end); bin++)
    {
        htemp1->SetBinContent(bin, histo->GetBinContent(bin));
    }

    htemp1->SetName(Form("%s_peakRegion", histo->GetName()));
    htemp1->SetFillColor(kRed);
    htemp1->SetFillStyle(3003);
    htemp1->SetLineWidth(0);

    for (int bin = htemp2->GetXaxis()->FindBin(fBcgLeft[0]); bin <= htemp2->GetXaxis()->FindBin(fBcgLeft[1]); bin++)
    {
        htemp2->SetBinContent(bin, histo->GetBinContent(bin));
    }
    for (int bin = htemp2->GetXaxis()->FindBin(fBcgRight[0]); bin <= htemp2->GetXaxis()->FindBin(fBcgRight[1]); bin++)
    {
        htemp2->SetBinContent(bin, histo->GetBinContent(bin));
    }
    htemp2->SetName(Form("%s_peakRegion", histo->GetName()));
    htemp2->SetFillColor(kBlack);
    htemp2->SetFillStyle(3003);
    htemp2->SetLineWidth(0);

    histo->GetXaxis()->SetRangeUser(std::min(peak.peak_region_start, fBcgLeft[0]),
                                    std::max(peak.peak_region_end, fBcgRight[1]));
    TCanvas c(Form("c_%s", histo->GetName()), histo->GetTitle(), 800, 600);
    histo->Draw();
    fPeakFcnTemp.Draw("same");
    htemp1->Draw("same");
    htemp2->Draw("same");
    c.Write();
    // cleanup
    histo->GetXaxis()->SetRangeUser(histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
    gROOT->SetBatch(false);
    file->cd("..");

    delete htemp1;
    delete htemp2;
}

///////////////////////////////////////// TheuerkaufFitter /////////////////////////////////////////
