//---------------------------------------------------------------------------------------------------------------------------//
//  Script which will compare the measured minimum position with the variable used for creating the MadGraph sample          //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` CalibCurve_Parabola3.C -o CalibCurve_Parabola3   //
//               ./CalibCurve_Parabola3                                                                                      //
//---------------------------------------------------------------------------------------------------------------------------//

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

using namespace std;

double getVarNr(std::string consVar){  //Update on 02/01/2016: Wrong minus-sign in Lagrangian, so invert all gR values!
  double ConsVar = 999;

  if( consVar == "Pos04")      ConsVar = -0.4; 
  else if( consVar == "Pos03") ConsVar = -0.3; 
  else if( consVar == "Pos02") ConsVar = -0.2; 
  else if( consVar == "Pos015") ConsVar = -0.15; 
  else if( consVar == "Pos01") ConsVar = -0.1; 
  else if( consVar == "Pos005") ConsVar = -0.05; 
  else if( consVar == "SM") ConsVar = 0.0; 
  else if( consVar == "Neg005") ConsVar = 0.05; 
  else if( consVar == "Neg01") ConsVar = 0.1; 
  else if( consVar == "Neg015") ConsVar = 0.15; 
  else if( consVar == "Neg02") ConsVar = 0.2; 
  else if( consVar == "Neg03") ConsVar = 0.3; 
  else if( consVar == "Neg04") ConsVar = 0.4; 
  
  return ConsVar;
}

std::string getVarStr(double consVar){  //Update on 02/01/2016: Wrong minus-sign in Lagrangian, so invert all gR values!
  std::string ConsVar = "";

  if( consVar == -0.4)        ConsVar = "Pos04";
  else if( consVar == -0.3)   ConsVar = "Pos03";
  else if( consVar == -0.2)   ConsVar = "Pos02";
  else if( consVar == -0.15)  ConsVar = "Pos015";
  else if( consVar == -0.1)   ConsVar = "Pos01";
  else if( consVar == -0.05)  ConsVar = "Pos005";
  else if( consVar == 0.0)   ConsVar = "SM";
  else if( consVar == 0.05) ConsVar = "Neg005";
  else if( consVar == 0.1)  ConsVar = "Neg01";
  else if( consVar == 0.15) ConsVar = "Neg015";
  else if( consVar == 0.2)  ConsVar = "Neg02";
  else if( consVar == 0.3)  ConsVar = "Neg03";
  else if( consVar == 0.4)  ConsVar = "Neg04";
  
  return ConsVar;
}

int main(int argc, char *argv[]){

  // ------------------------------------------------------- //
  // --- Use the input variables to specify some choices --- //
  // ---     * Are cuts applied (y/n)                    --- //
  // ---     * Should outer bins be excluded for fit?    --- //
  // ---     * Should a cut-value be applied?            --- //
  // ------------------------------------------------------- //
  const int NrSamples = 11;
  std::string sample[NrSamples] = {"Neg03","Neg02","Neg015","Neg01","Neg005","SM","Pos005","Pos01","Pos015","Pos02","Pos03"};
  double VarGraph[NrSamples] = {-0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3};
  double VarGraphSign[NrSamples] = {0.3, 0.2, 0.15, 0.1, 0.05, 0.0, -0.05, -0.1, -0.15, -0.2, -0.3};  //Update on 02/01/2016: Wrong minus-sign in Lagrangian, so invert all gR values!
  const int NrConfigs = 13;
  double Var[NrConfigs] = {-0.4, -0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
  double VarSign[NrConfigs] = {0.4, 0.3, 0.2, 0.15, 0.1, 0.05, -0.0, -0.05, -0.1, -0.15, -0.2, -0.3, -0.4};  //Update on 02/01/2016: Wrong minus-sign in Lagrangian, so invert all gR values!
  bool savePlots = true;

  std::string whichCut = "";
  if( argc >= 2){
    if( string(argv[1]) == "y") whichCut = "CutsAlsoOnMET";
    else                        whichCut = "NoCuts";
  }

  std::string outerBins = "";
  double FitMin = -0.4, FitMax = 0.4;
  std::cout << " Number of arguments = " << argc << std::endl;
  if( argc >= 3){
    std::cout << " Value of 2nd argument = " << string(argv[2]) << std::endl;
    if( string(argv[2]) == "y"){ outerBins = "_OuterBinsExclForFit"; FitMin = -0.3; FitMax = 0.3;}
  }

  std::string whichLikCut = "";
  if(argc >= 4){
    whichLikCut = ("_LikelihoodCut"+string(argv[3])).c_str();
  }

  //Store the minimum values and their uncertainties
  double Min_PolFit[NrSamples] = {0.}, Err_PolFit[NrSamples] = {0.};
  double Min_Parab3[NrSamples] = {0.}, Err_Parab3[NrSamples] = {0.};
  double Min_Parab5[NrSamples] = {0.}, Err_Parab5[NrSamples] = {0.};

  //Loop over the different gR-values and open the correct ROOT file
  for(int iSample = 0; iSample < NrSamples; iSample++){

    TFile* file = new TFile( ("Events_CalibCurve/CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample"+sample[iSample]+"_20000Evts_"+whichCut+"/FitDistributions_CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample"+sample[iSample]+"_20000Evts_"+whichCut+""+whichLikCut+""+outerBins+"_20000Evts.root").c_str() ,"READ");

    //Now get the relevant histogram from this file:
    TH1F* h_SummedWeights = (TH1F*) file->Get("FitSums/SummedHist_Acc");
  
    // 1) Access the TF1 over the TH1F using GetFunction:
    TF1* polFit_Sum = h_SummedWeights->GetFunction("polFit_Acc_SummedHist");
    polFit_Sum->SetParameter(1,polFit_Sum->GetParameter(1)*-1);               //Mirroring a quadratic function corresponds to changing its second parameter (~bx) to -b!!
    Min_PolFit[iSample] = polFit_Sum->GetMinimumX();
    std::cout << " Looking at variable " << getVarStr(VarGraph[iSample]) << " which has minimum = " << polFit_Sum->GetMinimumX() << std::endl;
    Err_PolFit[iSample] = (polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, polFit_Sum->GetMinimumX(), FitMax) - polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, FitMin, polFit_Sum->GetMinimumX()))/2.0;
    
    // 2) Get the output from a parabola through the inner 3 points!
    double Entries[NrSamples];
    double PosConfig = VarSign[0], NegConfig = VarSign[NrConfigs-1], PosPosConfig = VarSign[0], NegNegConfig = VarSign[NrConfigs-1];
    for(int iConfig = 1; iConfig < NrConfigs-1; iConfig++){    //Only want to run from -0.3 to 0.3!
      if(Var[iConfig+2] == getVarNr(sample[iSample]) ){ NegNegConfig = Var[iConfig];}
      if(Var[iConfig+1] == getVarNr(sample[iSample]) ){ NegConfig = Var[iConfig];}
      if(Var[iConfig-1] == getVarNr(sample[iSample]) ){ PosConfig = Var[iConfig];}
      if(Var[iConfig-2] == getVarNr(sample[iSample]) ){ PosPosConfig = Var[iConfig];}

      Entries[iConfig-1] = h_SummedWeights->GetBinContent(h_SummedWeights->FindBin(VarSign[iConfig]));
    }
    std::cout << " For sample " << getVarStr(VarGraph[iSample]) << " the +1 and -1 configs are : " << PosConfig << " and " << NegConfig << std::endl;
    TGraph* gr_Summed = new TGraph(NrSamples, VarGraph, Entries);
    gr_Summed->GetXaxis()->SetTitle("g_{R} coefficient");
    gr_Summed->GetYaxis()->SetTitle("-ln(L_{MEM}) output"); gr_Summed->GetYaxis()->SetTitleOffset(1.5);
    gr_Summed->SetMarkerStyle(22);

    TF1* polFitParabola_Sum = new TF1(("polFitParabola3_AroundMin_Acc_Var"+sample[iSample]).c_str(),"pol2",FitMin,FitMax);
    polFitParabola_Sum->SetLineColor(6);
    gr_Summed->Fit(polFitParabola_Sum,"Q","",NegConfig, PosConfig);
    std::cout << " Trying to find minimum between " << NegConfig << " and " << PosConfig << std::endl;
    Min_Parab3[iSample] = polFitParabola_Sum->GetMinimumX();
    Err_Parab3[iSample] = (polFitParabola_Sum->GetX(polFitParabola_Sum->GetMinimum()+0.5, polFitParabola_Sum->GetMinimumX(), 0.4) - polFitParabola_Sum->GetX(polFitParabola_Sum->GetMinimum()+0.5, -0.4, polFitParabola_Sum->GetMinimumX()))/2.0;

    TF1* polFitParabola5_Sum = new TF1(("polFitParabola5_AroundMin_Acc_Var"+sample[iSample]).c_str(),"pol2",NegNegConfig,PosPosConfig);
    polFitParabola5_Sum->SetLineColor(7);
    gr_Summed->Fit(polFitParabola5_Sum,"Q","",NegNegConfig, PosPosConfig);
    Min_Parab5[iSample] = polFitParabola5_Sum->GetMinimumX();
    Err_Parab5[iSample] = (polFitParabola5_Sum->GetX(polFitParabola5_Sum->GetMinimum()+0.5, polFitParabola5_Sum->GetMinimumX(), 0.4) - polFitParabola5_Sum->GetX(polFitParabola5_Sum->GetMinimum()+0.5, -0.4, polFitParabola5_Sum->GetMinimumX()))/2.0;

    TCanvas* canv_CompMinFit = new TCanvas(("Comparison_FitMethod_Var"+getVarStr(VarGraph[iSample])).c_str(), ("Comparison_FitMethod_Var"+getVarStr(VarGraph[iSample])).c_str());
    canv_CompMinFit->cd();
    gr_Summed->Draw("AP");
    polFit_Sum->Draw("same");
    polFitParabola_Sum->Draw("same");
    polFitParabola5_Sum->Draw("same");
    //polParabolaSM_Sum->Draw("same");
    canv_CompMinFit->SaveAs(("Events_CalibCurve/Comparison_FitMethod_Var"+getVarStr(VarGraph[iSample])+".pdf").c_str());
  }

  //Once all minima are stored, plot them into a TGraphErrors
  TGraphErrors* gr_FitResult = new TGraphErrors(NrSamples, VarGraphSign, Min_PolFit, 0, Err_PolFit);
  TGraphErrors* gr_Parabola3Result = new TGraphErrors(NrSamples, VarGraphSign, Min_Parab3, 0, Err_Parab3);
  TGraphErrors* gr_Parabola5Result = new TGraphErrors(NrSamples, VarGraphSign, Min_Parab5, 0, Err_Parab5);

  TGraph *gr_StraightLine = new TGraph(NrSamples, VarGraphSign, VarGraphSign);
  gr_StraightLine->SetLineStyle(2); gr_StraightLine->SetLineColor(13);  
  gr_StraightLine->GetXaxis()->SetTitle("Input value of the g_{R} coefficient");
  gr_StraightLine->GetYaxis()->SetTitle("Measured g_{R} coefficient");

  //Fit the calibration curves with a straight line (slope should be as close as possible to 1)
  TF1* polFit_Line = new TF1("polFitLine","pol1",-0.16,0.16); 

  TCanvas* canv_FitResult = new TCanvas("LinearityTest_FitResult","LinearityTest_FitResult");
  canv_FitResult->cd();
  gr_FitResult->Fit(polFit_Line,"Q","",-0.16, 0.16);
  gr_StraightLine->SetTitle("");
  gr_StraightLine->Draw("AC");
  gr_FitResult->Draw("P");
  gr_FitResult->SetMarkerStyle(20);
  gr_FitResult->SetMarkerSize(0.35);
  canv_FitResult->SaveAs("Events_CalibCurve/LinearityTest_FitResult.pdf");
  std::cout << " Slope obtained from FitResult is : " << polFit_Line->GetParameter(1) << std::endl;
  std::cout << " Bias obtained from FitResult is  : " << polFit_Line->GetParameter(0) << std::endl;

  TCanvas* canv_Parabola3Result = new TCanvas("LinearityTest_Parabola3Result","LinearityTest_Parabola3Result");
  canv_Parabola3Result->cd();
  gr_Parabola3Result->Fit(polFit_Line,"Q","",-0.16, 0.16);
  gr_StraightLine->SetTitle("Linearity test using the parabola through the 3 points of the expected minimum");
  gr_StraightLine->Draw("AC");
  gr_Parabola3Result->Draw("P");
  canv_Parabola3Result->SaveAs("Events_CalibCurve/LinearityTest_Parabola3Result.pdf");
  std::cout << " Slope obtained from Parabola3Result is : " << polFit_Line->GetParameter(1) << std::endl;
  std::cout << " Bias obtained from Parabola3Result is  : " << polFit_Line->GetParameter(0) << std::endl;
  
  TCanvas* canv_Parabola5Result = new TCanvas("LinearityTest_Parabola5Result","LinearityTest_Parabola5Result");
  canv_Parabola5Result->cd();
  gr_Parabola5Result->Fit(polFit_Line,"Q","",-0.16, 0.16);
  gr_StraightLine->SetTitle("Linearity test using the parabola through the 5 points of the expected minimum");
  gr_StraightLine->Draw("AC");
  gr_Parabola5Result->Draw("P");
  canv_Parabola5Result->SaveAs("Events_CalibCurve/LinearityTest_Parabola5Result.pdf");
  std::cout << " Slope obtained from Parabola5Result is : " << polFit_Line->GetParameter(1) << std::endl;
  std::cout << " Bias obtained from Parabola5Result is  : " << polFit_Line->GetParameter(0) << std::endl;
  
  return 0;
}
