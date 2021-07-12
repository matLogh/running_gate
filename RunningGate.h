#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFormula.h"
#include <vector>
#include <iostream>
#include "TMath.h"

using namespace std;

#ifndef RUNNING_GATE
#define RUNNING_GATE

struct PeakCounts
{
	int peak_bin_start;
	int peak_bin_end;
	int bcgl_bin_start;
	int bcgl_bin_end;
	int bcgr_bin_start;
	int bcgr_bin_end;

	double peak_en_start;
	double peak_en_end;
	double bcgl_en_start;
	double bcgl_en_end;
	double bcgr_en_start;
	double bcgr_en_end;
};

struct PeakData
{
	TF1* fcn;
	PeakCounts* pc;
	int fcn_area_parameter;
	vector< double > fit_area;
	vector< double > fit_area_err;
	vector< double > running_gate_center;
	vector< double > running_gate_center_err;
};



class RunningGate
{
	
public:
	RunningGate(TH2D* m);
	void SetRunningGate(double g_start, double g_end, double g_width);
	~RunningGate();
	void RunGate();
	TGraphErrors* GetFcn(int i);
	vector<TGraphErrors*> GetFcn();
	void SetFcn(double en_centroid, double e_start, double e_end);
	void SetFcn(TF1* f, int area_parameter);
	void SetFcn(double pe_start, double pe_end,
				double bl_start, double bl_end,
				double br_start, double br_end);
private:
	TF1* fcn_peak1;
	TF1* fcn_peak2;

	vector< PeakData* > peak_object;

	double running_gate_start_energy;
	double running_gate_end_energy;
	double running_gate_width_energy;
	int running_gate_start_bin = -999;
	int running_gate_end_bin = -999;
	int running_gate_width_bin = -999;
	TH2D* matrix;
	TH1D* gate_slice;
	double bin_width;

	void CreateGraphs();
   // double GausLin(double* x, double* p);

	vector <TGraphErrors*> graph;
};

#endif