//----------------------------------------------------------------------------------------------------------------------------------//
//     //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` LikCut_MinimumEvolution.C -o LikCut_MinimumEvolution    //
//               ./LikCut_MinimumEvolution                                                                                          //
//----------------------------------------------------------------------------------------------------------------------------------//

#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[]){

  //Given input parameters are the following:
  // 1) Number of root files
  // 2) - ...) Path of root file
  int nrRootFiles = 0;
  std::vector<std::string> rootFilePath;
  rootFilePath.clear();
  double CutValue[100];
  if( argc >= 1){
    nrRootFiles = atoi(string(argv[1]).c_str());
    for(int ii = 0; ii < nrRootFiles; ii++){
      rootFilePath.push_back(string(argv[2+ii]));
      CutValue[ii] = atof(string(argv[2+nrRootFiles+ii]).c_str());
    }
  }
  else{
    std::cout << " ERROR: script can only work when root files are given as input !! " << std::endl;
    std::cout << "   --> This script will create the minimum position for different likelihood cuts! " << std::endl;
    exit(0);
  }

  //Now loop over all the different ROOT files, open the histSum histogram, get the minimum and store this in a sort of calibration-curve TGraph.
  // --> This will show how the position of the minimum moves over the considered range when changing the likelihood cut!
  //TH2F* minimumCalibration = new TH2F("minimumCalib","Determined minimum for different likelihood cuts",

  double Minimum[100] = {0.};
  double MinError[100] = {0.};
  for(int iFile = 0; iFile < nrRootFiles; iFile++){

    TFile* file = new TFile(rootFilePath[iFile].c_str(),"READ");
    TH1F* h_SummedWeights = (TH1F*) file->Get("FitSums/SummedHist_Acc");
    TF1* pol_Sum = h_SummedWeights->GetFunction("polFit_Acc_SummedHist");

    Minimum[iFile] = pol_Sum->GetMinimumX();
    MinError[iFile] = pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, pol_Sum->GetMinimumX(), 0.4) - pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, -0.4, pol_Sum->GetMinimumX());
    std::cout << " Cut value is : " << CutValue[iFile] << std::endl;
  }

  TGraphErrors *gr_HistSum = new TGraphErrors(nrRootFiles, CutValue, Minimum, 0, MinError);
  gr_HistSum->SetMarkerStyle(1); gr_HistSum->SetMarkerColor(1); gr_HistSum->SetLineColor(1);
  TCanvas *canv = new TCanvas("canv","canv");
  canv->cd();
  gr_HistSum->Draw("AP");
  canv->SaveAs("MinimumPosition_Reco.pdf");
  gr_HistSum->SaveAs("MinimumPosition_Reco.root");
}
