#include "RunningGate.h"

void RunningGate::RunGate()
{
	if(peak_object.size() == 0)
	{
		cout << "No function/fit region was set!" << endl;
		cout << "Nothing to be done..." << endl;
		return;
	}
	if(running_gate_start_bin == -999)
	{
		running_gate_start_bin = matrix->GetXaxis()->FindBin(running_gate_start_energy);
		running_gate_end_bin = matrix->GetXaxis()->FindBin(running_gate_end_energy);
		
		running_gate_width_bin = running_gate_width_energy/(double)
									matrix->GetXaxis()->GetBinWidth(running_gate_start_bin);
	}
	cout << peak_object.size() << endl;

	//get the data
	TFile fx("slice.root","RECREATE");

	for(int bin_start = running_gate_start_bin; bin_start < running_gate_end_bin + running_gate_width_bin; bin_start++)
	{
		gate_slice = matrix->ProjectionY("p",bin_start, bin_start+running_gate_width_bin);
		double center = (bin_start+running_gate_width_bin/(double)2.)*bin_width;
		for(unsigned int i=0; i<peak_object.size(); i++)
		{
			PeakData *peak = peak_object[i];
			//if we are going to fit distribution
			if(peak->fcn != NULL)
			{
				gate_slice->Fit(peak->fcn, "RLQ");

				char ccc[50];
				sprintf(ccc,"slice_%i_peak_%i",bin_start, i);
				gate_slice->SetName(ccc);
				fx.cd();
				gate_slice->Write();

				double area = peak->fcn->GetParameter(peak->fcn_area_parameter);
				double area_err = peak->fcn->GetParError(peak->fcn_area_parameter);
				cout << area << " " << area_err << endl;

				peak->fit_area.push_back(area);
				peak->fit_area_err.push_back(area_err);
				peak->running_gate_center.push_back(center);
				peak->running_gate_center_err.push_back(bin_width/(double)2.);
			}
			//or if we are going to count events
			else
			{
				int area=0;
				int bcg =0;
				int area_err=0;
				int bcg_err=0;
				for(int i=peak->pc->peak_bin_start; i<=peak->pc->peak_bin_end; i++)
				{
					area += gate_slice->GetBinContent(i);
					area_err += gate_slice->GetBinError(i);
				}
				for(int i=peak->pc->bcgl_bin_start; i<=peak->pc->bcgl_bin_end; i++)
				{
					bcg += gate_slice->GetBinContent(i);
					bcg_err += gate_slice->GetBinError(i);
				}
				for(int i=peak->pc->bcgr_bin_start; i<=peak->pc->bcgr_bin_end; i++)
				{
					bcg += gate_slice->GetBinContent(i);
					bcg_err += gate_slice->GetBinError(i);
				}

				double norm = (peak->pc->peak_bin_end-peak->pc->peak_bin_start)/(double)
							(peak->pc->bcgl_bin_end-peak->pc->bcgl_bin_start +
							peak->pc->bcgr_bin_end-peak->pc->bcgr_bin_start);

				bcg_err = bcg_err*(norm*norm);
				area = area - bcg*norm;
				area_err = TMath::Sqrt(area_err+bcg);


				peak->fit_area.push_back(area);
				peak->fit_area_err.push_back(area_err);
				peak->running_gate_center.push_back(center);
				peak->running_gate_center_err.push_back(bin_width/(double)2.);
			}
			
		}
	}
	CreateGraphs();
}


void RunningGate::SetFcn(double pe_start, double pe_end,
					 	 double bl_start, double bl_end,
					 	 double br_start, double br_end)
{
	PeakCounts* ps = new PeakCounts;
	ps->peak_en_start  = pe_start;
	ps->peak_en_end  = pe_end;
	ps->bcgl_en_start  = bl_start;
	ps->bcgl_en_end  = bl_end;
	ps->bcgr_en_start  = br_start;
	ps->bcgr_en_end  = br_end;

	ps->peak_bin_start = pe_start/(double)bin_width;
	ps->peak_bin_end = pe_end/(double)bin_width;
	ps->bcgl_bin_start = bl_start/(double)bin_width;
	ps->bcgl_bin_end = bl_end/(double)bin_width;
	ps->bcgr_bin_start = br_start/(double)bin_width;
	ps->bcgr_bin_end = br_end/(double)bin_width;

	PeakData* pd = new PeakData;
	pd->pc = ps;
	pd->fcn = NULL;
	pd->fcn_area_parameter = -1;
	peak_object.push_back(pd);
}

void RunningGate::SetRunningGate(double g_start, double g_end, double g_width)
{
	running_gate_start_energy = g_start;
	running_gate_end_energy   = g_end;
	running_gate_width_energy = g_width;
}



void RunningGate::SetFcn(TF1* f, int area_parameter)
{
	PeakData *pd = new PeakData;
	pd->fcn = f;
	pd->fcn_area_parameter = area_parameter;
	peak_object.push_back(pd);
}

void RunningGate::SetFcn(double en_centroid, double e_start, double e_end)
{
	char str[50];
	sprintf(str, "fcn_%zu",peak_object.size());
	TF1* f = new TF1(str,"GausLin",e_start, e_end,5);
	f->SetParameter(0,50); //TODO ???
	f->SetParLimits(0,0,1000000); //TODO ???
	f->SetParameter(1,en_centroid);
	f->SetParLimits(1,en_centroid-bin_width, en_centroid+bin_width);
	f->SetParameter(2,(e_end-e_start)/4);
	f->SetParLimits(2,0,100000);

	//estimate linear parameters
	int b_start = matrix->GetXaxis()->FindBin(e_start);
	int b_end = matrix->GetXaxis()->FindBin(e_end);
	double p = (matrix->GetBinContent(b_start) - matrix->GetBinContent(b_end)) 
				/(double) (b_start - b_end);
	double q = matrix->GetBinContent(b_start) - p*b_start;
	f->SetParameter(3, q);
	f->SetParameter(4, p);

	//fcn_peak.push_back(f);
	//fcn_area_parameter.push_back(0);

	PeakData *pd = new PeakData;
	pd->fcn = f;
	pd->fcn_area_parameter = 0;
	peak_object.push_back(pd);
}


RunningGate::RunningGate(TH2D* m)
{
	matrix = m;
	bin_width = matrix->GetXaxis()->GetBinWidth(1);
}

RunningGate::~RunningGate()
{
	//no idea what can i delete due to ROOT inheritance rules
	//so do nothing
}


void RunningGate::CreateGraphs()
{

	for(unsigned int i =0; i<peak_object.size(); i++)
	{
		char str[50];
		PeakData *peak = peak_object[i];
		int en=0;
		if(peak->fcn == NULL)
		{
			en = (peak->pc->peak_en_end-peak->pc->peak_en_end)/2. + peak->pc->peak_en_end;
		}
		for(int j=0; j<peak->running_gate_center.size();j++)
		{
			cout << peak->running_gate_center.at(j) << " " << peak->fit_area.at(j) << endl;
		}
		sprintf(str, "g_%zu_en_%i",peak_object.size(), en);
		TGraphErrors *g = new TGraphErrors(peak->running_gate_center.size(),
							&(peak->running_gate_center[0]),
							&(peak->fit_area[0]),
							&(peak->running_gate_center_err[0]),
							&(peak->fit_area_err[0])
							);
		g->SetName(str);
		g->SetTitle(str);
		g->SetMarkerColor(4);
   		g->SetMarkerStyle(21);

   		graph.push_back(g);
	}
}

TGraphErrors* RunningGate::GetFcn(int i)
{
	return graph.at(i);
}

vector<TGraphErrors*> RunningGate::GetFcn()
{
	return graph;
}