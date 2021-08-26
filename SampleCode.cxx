#include "RunningGate.h"

//if you run this in ROOT:
//first compile and load the RunningGate library by
// .L RunningGate.cxx+
//then you can start this macro by
//.x SampleCode.cxx

//the RunningGate_test.cxx has in addition also some fitting
//of the resulting TGraphs

void SampleCode()
{
  //pass symetric matrix to the running gate
	TFile file("matrix.root","READ");
	TH2D* mat = (TH2D*)file.Get("m");
  //create file for storing results
	TFile ffile("runningGate.root","RECREATE");

	RunningGate gate(mat);
	gate.SetRunningGate(120, //start running the gate at this energy
                      160, //stop running the gate at this energy
                      1);  //width of the running gate in energy
                           //you should make this larger in case of poor statistics
	//setting the gate1
  gate.SetFcn(240,245, //start and stop energy of the gate 
              233,238, //start and stop energy of the background1
              251,255);//start and stop energy of the background2
	//setting gate2
  gate.SetFcn(218,223,197,204,233,238);
  //run the gate
	gate.RunGate();
  //get the resulting TGraphs
	TGraphErrors* f1 = gate.GetFcn(0);
	TGraphErrors* f2 = gate.GetFcn(1);
  //save results to the file
	ffile.cd();
	f1->Write();
	f2->Write();
	file.Close();
	ffile.Close();
}
