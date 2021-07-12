#include "RunningGate.h"

double GausLin(double* x, double* p)
{
	return p[0]*TMath::Exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]))/(TMath::Sqrt(2*TMath::Pi())*p[2]) + p[3] + p[4]*x[0];
}

double GausGausLin(double* x, double* p)
{
	return p[0]*TMath::Exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]))/(TMath::Sqrt(2*TMath::Pi())*p[2]) +
	       p[3]*TMath::Exp(-0.5*((x[0]-p[4])/p[5])*((x[0]-p[4])/p[5]))/(TMath::Sqrt(2*TMath::Pi())*p[5]) +
		 p[6] + p[7]*x[0];
}

double GausGausGausLin(double* x, double* p)
{
	return p[0]*TMath::Exp(-0.5*((x[0]-p[1])/p[2])*((x[0]-p[1])/p[2]))/(TMath::Sqrt(2*TMath::Pi())*p[2]) +
	       p[3]*TMath::Exp(-0.5*((x[0]-p[4])/p[5])*((x[0]-p[4])/p[5]))/(TMath::Sqrt(2*TMath::Pi())*p[5]) +
	       p[6]*TMath::Exp(-0.5*((x[0]-p[7])/p[8])*((x[0]-p[7])/p[8]))/(TMath::Sqrt(2*TMath::Pi())*p[8]) +
		 p[9] + p[10]*x[0];
}

void RunningGate_test()
{
	TFile file("matrix.root","READ");
	TH2D* mat = (TH2D*)file.Get("m");

	TFile ffile("testFile.root","RECREATE");

	RunningGate gate(mat);
	gate.SetRunningGate(120,160,1);
	gate.SetFcn(240,245,233,238,251,255);
	gate.SetFcn(218,223,197,204,233,238);

	gate.RunGate();

	TGraphErrors* f1 = gate.GetFcn(0);
	TGraphErrors* f2 = gate.GetFcn(1);

	TF1* fgg = new TF1("fgg",GausGausLin,130, 160,8);
	fgg->SetParNames("area1", "mu1","sig1","area2", "mu2","sig2", "offset", "lin");
	fgg->SetParameter(0,1000); //TODO ???
	fgg->SetParLimits(0,0,10000); //TODO ???
	fgg->SetParameter(1,146);
	fgg->SetParLimits(1,144, 147);
	fgg->SetParameter(2,1.4);
	fgg->SetParLimits(2,0.5,2);

	fgg->SetParameter(3,1000); //TODO ???
	fgg->SetParLimits(3,0,10000); //TODO ???
	fgg->SetParameter(4,152.5);
	fgg->SetParLimits(4,151, 154);
	fgg->SetParameter(5,1.4);
	fgg->SetParLimits(5,0.5,2);

	fgg->SetParameter(6,50);
	fgg->SetParameter(7,0);

	TF1* fggg = new TF1("fggg",GausGausGausLin,130, 160,11);
	fggg->SetParNames("area1", "mu1","sig1","area2", "mu2","sig2","area3", "mu3","sig3", "offset", "lin");
	fggg->SetParameter(0,1000); //TODO ???
	fggg->SetParLimits(0,0,10000); //TODO ???
	fggg->SetParameter(1,146);
	fggg->SetParLimits(1,144, 147);
	fggg->SetParameter(2,1.4);
	fggg->SetParLimits(2,0.5,2);

	fggg->SetParameter(3,1000); //TODO ???
	fggg->SetParLimits(3,0,10000); //TODO ???
	fggg->SetParameter(4,152.5);
	fggg->SetParLimits(4,151, 154);
	fggg->SetParameter(5,1.4);
	fggg->SetParLimits(5,0.5,2);

	fggg->SetParameter(6,1000); //TODO ???
	fggg->SetParLimits(6,0,10000); //TODO ???
	fggg->SetParameter(7,140);
	fggg->SetParLimits(7,138, 144);
	fggg->SetParameter(8,1.4);
	fggg->SetParLimits(8,0.5,5);

	fggg->SetParameter(9,50);
	fggg->SetParameter(10,0);



	TF1* ff = new TF1("ff",GausLin,130, 160,5);
	ff->SetParameter(0,50); //TODO ???
	ff->SetParLimits(0,0,10000); //TODO ???
	ff->SetParameter(1,146);
	ff->SetParLimits(1,145, 147);
	ff->SetParameter(2,1.4);
	ff->SetParLimits(2,0.5,2);
	ff->SetParameter(3,10);
	ff->SetParameter(4,0);

	cout << "245:" << endl;
	f1->Fit(fggg,"REM");
	cout << "221:" << endl;
	f2->Fit(fgg,"REM");

	ffile.cd();
	f1->Write();
	f2->Write();
	file.Close();
	ffile.Close();
}