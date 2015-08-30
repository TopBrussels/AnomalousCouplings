//-----------------------------------------------------------------------------------------------------------------------------//
//                                                                                                                             //
//  Create a script which will compare the measured minimum position with the variable used for creating the MadGraph sample   //
//                                                                                                                             //
//-----------------------------------------------------------------------------------------------------------------------------//

#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"

void getMinimum(TH1F* h_Fst, TH1F* h_Scd, float rangeVal, int range, double MinFst[], double MinScd[], double eMinFst[], double eMinScd[]){
  
  TF1* polFit_Fst = new TF1("polFit_Fst","pol2",-0.5,0.5); 
  h_Fst->Fit(polFit_Fst,"Q");

  MinFst[range] = polFit_Fst->GetMinimumX();
  double Neg1Sigma_Fst = polFit_Fst->GetX(h_Fst->GetBinContent(h_Fst->FindBin(polFit_Fst->GetMinimumX()))+0.5,-0.5,polFit_Fst->GetMinimumX());
  double Pos1Sigma_Fst = polFit_Fst->GetX(h_Fst->GetBinContent(h_Fst->FindBin(polFit_Fst->GetMinimumX()))+0.5,polFit_Fst->GetMinimumX(), 0.5);
  eMinFst[range] = Pos1Sigma_Fst - Neg1Sigma_Fst;
  
  TF1* polFit_Scd = new TF1("polFit_Scd","pol2",-0.5,0.5); 
  h_Scd->Fit(polFit_Scd,"Q");

  MinScd[range] = polFit_Scd->GetMinimumX();
  double Neg1Sigma_Scd = polFit_Scd->GetX(h_Scd->GetBinContent(h_Scd->FindBin(polFit_Scd->GetMinimumX()))+0.5,-0.3,                      polFit_Scd->GetMinimumX());
  double Pos1Sigma_Scd = polFit_Scd->GetX(h_Scd->GetBinContent(h_Scd->FindBin(polFit_Scd->GetMinimumX()))+0.5,polFit_Scd->GetMinimumX(), 0.3);
  eMinScd[range] = Pos1Sigma_Scd - Neg1Sigma_Scd;
}

void CalibrationCurve(){

  //------------------------------------------------------//
  //  Set the input correct for the considered directory  //	
  //------------------------------------------------------//
  bool saveAsPDF = false;

  int choice = 1;
  std::string whichDir_[] = {"Events_RgRScan_CalibrationCurve","Events_RgRScan_CalibrationCurve","Events_RgRScan_AllDeltaTF_EDep"};  std::string whichDir = whichDir_[choice];
  std::string fileTitle_[] = {"RgR_CalibrationCurve_AllDeltaTF","RgR_CalibrationCurve_AllDeltaTF","AllDeltaTF"};                     std::string FileTitle = fileTitle_[choice];

  double rangeZero[] = {-0.4, -0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3};
  double rangeOne[]  = {-0.4, -0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
  double rangeTwo[]  = {-0.3, -0.2, -0.1, -0.05,  0.0,  0.05, 0.1, 0.2,  0.3};
  std::string RNameZero[] = {"Neg04","Neg03","Neg02","Neg015","Neg01","Neg005","SM","Pos005","Pos01","Pos015","Pos02","Pos03"};
  std::string RNameOne[] = {"Neg04","Neg03","Neg02","Neg015","Neg01","Neg005","SM","Pos005","Pos01","Pos015","Pos02","Pos03","Pos04"};
  std::string RNameTwo[] = {"Neg03","Neg02","Neg01","Neg005","SM","Pos005","Pos01","Pos02","Pos03"};
  int rangeSize_[] = {sizeof(rangeZero)/sizeof(rangeZero[0]), sizeof(rangeOne)/sizeof(rangeOne[1]), sizeof(rangeTwo)/sizeof(rangeTwo[1])}; const int rangeSize = rangeSize_[choice];

  double Range[rangeSize];
  std::string RangeName[rangeSize];
  switch (choice){
    case 0: for(int i = 0; i < rangeSize; i++){ Range[i] = rangeZero[i]; RangeName[i] = RNameZero[i];} break;
    case 1: for(int i = 0; i < rangeSize; i++){ Range[i] = rangeOne[i];  RangeName[i] = RNameOne[i];}  break;
    case 2: for(int i = 0; i < rangeSize; i++){ Range[i] = rangeTwo[i];  RangeName[i] = RNameTwo[i];}  break;
  }

  //Set the cut information!
  int cutChoice = 1;
  std::string cuts_[] = {"NoCuts", "CutsOnMET", "AllCuts"}; std::string cuts = cuts_[cutChoice];

  //CosThetaCorrection or not??
  int cosThetaChoice = 0;
  std::string cosTh_[] ={"","_CosThetaReweightingApplied"}; std::string cosTh = cosTh_[cosThetaChoice];

  //Initialize the arrays for the TGraphErrors:
  double min_Fst[rangeSize], min_Scd[rangeSize], cal_Fst[rangeSize], cal_Scd[rangeSize], calRel_Fst[rangeSize], calRel_Scd[rangeSize];
  double Emin_Fst[rangeSize], Emin_Scd[rangeSize], Ecal_Fst[rangeSize], Ecal_Scd[rangeSize], EcalRel_Fst[rangeSize], EcalRel_Scd[rangeSize];

  //Now loop over all the configurations that need to be considered and get the correct distributions
  for( int iRange = 0; iRange < rangeSize; iRange++){
    TFile* file = new TFile((whichDir+"/"+FileTitle+"_RgRMGSample"+RangeName[iRange]+"_25000Evts_"+cuts+"/FitOptimizations_"+FileTitle+"_RgRMGSample"+RangeName[iRange]+"_25000Evts_"+cuts+""+cosTh+"_25000Evts.root").c_str(),"READ");
  
    //Now get the relevant histograms from this file:
    TH1F* h_FirstPol = (TH1F*) file->Get("OriginalDistributions/FirstPolAcc_Summed");
    TH1F* h_SecondPol = (TH1F*) file->Get("OriginalDistributions/SecondPolAcc_Summed");

    //Now send this to the class which gets the minimum
    getMinimum(h_FirstPol, h_SecondPol, Range[iRange], iRange, min_Fst, min_Scd, Emin_Fst, Emin_Scd);
  
    //Now fill the other arrays:
    cal_Fst[iRange] = min_Fst[iRange] - Range[iRange];  Ecal_Fst[iRange] = Emin_Fst[iRange];
    cal_Scd[iRange] = min_Scd[iRange] - Range[iRange];  Ecal_Scd[iRange] = Emin_Scd[iRange];
    calRel_Fst[iRange] = cal_Fst[iRange]/min_Fst[iRange]; EcalRel_Fst[iRange] = (Emin_Fst[iRange]*abs(Range[iRange]))/(Emin_Fst[iRange]*Emin_Fst[iRange]);
    calRel_Scd[iRange] = cal_Scd[iRange]/min_Fst[iRange]; EcalRel_Scd[iRange] = (Emin_Scd[iRange]*abs(Range[iRange]))/(Emin_Scd[iRange]*Emin_Scd[iRange]);
  }
  
  //ROOT file where everything will be stored!  
  TFile* outputFile = new TFile((whichDir+"/"+FileTitle+"_25000Evts_"+cuts+""+cosTh+"_Calibration.root").c_str(),"RECREATE");
  outputFile->cd();

  //Make a straight line which indicates where you expect the points!
  //double x[rangeSize] = Range;
  //double yMin[rangeSize] = Range;
  //double yCal[] = {0.0};
  TGraph* gr_LineMin = new TGraph(rangeSize, Range, Range); 
  //TGraph* gr_LineCal = new TGraph(rangeSize, Range, yCal);

  //------------------------//
  //   Fill the TGraphs !   //
  //------------------------//

  //Make the TGraphError!
  TGraphErrors *gr_MinFst = new TGraphErrors(rangeSize, Range, min_Fst, 0, Emin_Fst);
  TGraphErrors *gr_CalFst = new TGraphErrors(rangeSize, Range, cal_Fst, 0, Ecal_Fst);
  TGraphErrors *gr_RelCalFst = new TGraphErrors(rangeSize, Range, calRel_Fst, 0, EcalRel_Fst);
  TGraphErrors *gr_MinScd = new TGraphErrors(rangeSize, Range, min_Scd, 0, Emin_Scd);
  TGraphErrors *gr_CalScd = new TGraphErrors(rangeSize, Range, cal_Scd, 0, Ecal_Scd);
  TGraphErrors *gr_RelCalScd = new TGraphErrors(rangeSize, Range, calRel_Scd, 0, EcalRel_Scd);

  //------  FirstPol Minimum  ------//
  TCanvas* canv_MinFst = new TCanvas(("MinimumCurve_FirstPol_"+cuts+""+cosTh).c_str(),"MinimumCurve_FirstFit"); canv_MinFst->cd();
  gr_MinFst->Draw("AP");
  gr_MinFst->SetMarkerStyle(1);
  gr_MinFst->SetTitle("Measured minimum versus used value of MG sample (using pol2 on all points)");
  gr_MinFst->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_MinFst->GetYaxis()->SetTitle("Minimum obtained from fit");
  gr_LineMin->SetLineColor(3); gr_LineMin->Draw("C"); 

  //Fit the minimum curve with a straight line!
  TF1* polFit_Line = new TF1("polFitLine","pol1",-0.4,0.4); 
  gr_MinFst->Fit(polFit_Line,"","",-0.25,0.25);

  //Add the fit information to this canvas
  stringstream ssSlope;        ssSlope << polFit_Line->GetParameter(1);       string sSlope = ssSlope.str();
  stringstream ssSlopeErr;     ssSlopeErr << polFit_Line->GetParError(1);     string sSlopeErr = ssSlopeErr.str();
  stringstream ssIntercept;    ssIntercept << polFit_Line->GetParameter(0);   string sIntercept = ssIntercept.str();
  stringstream ssInterceptErr; ssInterceptErr << polFit_Line->GetParError(0); string sInterceptErr = ssInterceptErr.str();
  stringstream ssChiSq;        ssChiSq << polFit_Line->GetChisquare();        string sChiSq = ssChiSq.str();
  stringstream ssNDF;          ssNDF << polFit_Line->GetNDF();                string sNDF = ssNDF.str();
  TLatex t;
  t.SetTextSize(0.025);
  t.DrawLatex(0.,-0.3,("Fit output (for "+cuts+"):").c_str());
  t.DrawLatex(0.,-0.35,(" - slope = "+sSlope+" #pm "+sSlopeErr).c_str());
  t.DrawLatex(0.,-0.4,(" - intercept = "+sIntercept+" #pm "+sInterceptErr).c_str());
  t.DrawLatex(0., -0.45,(" - #chi^{2}/ndf = "+sChiSq+"/"+sNDF).c_str());
  //t.Draw();
  if(saveAsPDF) canv_MinFst->SaveAs((whichDir+"/MinimumCurveFirstPol_"+cuts+""+cosTh+".pdf").c_str());
  canv_MinFst->Write();

  //------  FirstPol Calibration  ------//
  TCanvas* canv_CalFst = new TCanvas("CalibrationCurve_FirstFit","CalibrationCurve_FirstFit"); canv_CalFst->cd();
  gr_CalFst->Draw("A*");
  gr_CalFst->SetTitle("Deviation from expected minimum versus used value of MG sample (using pol2 on all points)");
  gr_CalFst->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_CalFst->GetYaxis()->SetTitle("Fit Minimum - expected minimum");
  //gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  if(saveAsPDF) canv_CalFst->SaveAs((whichDir+"/CalibrationCurveFirstPol_"+cuts+""+cosTh+".pdf").c_str());	
  canv_CalFst->Write();

  //------  FirstPol Relative calibration  ------//
  TCanvas* canv_RelCalFst = new TCanvas("RelativeCalibrationCurve_FirstFit","RelativeCalibrationCurve_FirstFit"); canv_RelCalFst->cd();
  gr_RelCalFst->Draw("A*");
  gr_RelCalFst->SetTitle("Relative deviation from expected minimum versus used value of MG sample (using pol2 on all points)");
  gr_RelCalFst->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_RelCalFst->GetYaxis()->SetTitle("Fit Min - expected min relative to Fit min");
  //gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  if(saveAsPDF) canv_RelCalFst->SaveAs((whichDir+"/RelativeCalibrationCurveFirstPol_"+cuts+""+cosTh+".pdf").c_str());	
  canv_RelCalFst->Write();

  //------  SecondPol Minimum  ------//
  TCanvas* canv_MinScd = new TCanvas("MinimumCurve_SecondFit","MinimumCurve_SecondFit"); canv_MinScd->cd();
  gr_MinScd->Draw("A*");
  gr_MinScd->SetTitle("Measured minimum versus used value of MG sample (using pol2 on reduced points)");
  gr_MinScd->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_MinScd->GetYaxis()->SetTitle("Minimum obtained from fit");
  gr_LineMin->SetLineColor(3); gr_LineMin->Draw("C"); 
  if(saveAsPDF) canv_MinScd->SaveAs((whichDir+"/MinimumCurveSecondPol_"+cuts+""+cosTh+".pdf").c_str());
  canv_MinScd->Write();

  //------  SecondPol Calibration  ------//
  TCanvas* canv_CalScd = new TCanvas("CalibrationCurve_SecondFit","CalibrationCurve_SecondFit"); canv_CalScd->cd();
  gr_CalScd->Draw("A*");
  gr_CalScd->SetTitle("Deviation from expected minimum versus used value of MG sample (using pol2 on reduced points)");
  gr_CalScd->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_CalScd->GetYaxis()->SetTitle("Fit Minimum - expected minimum");
  //gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  if(saveAsPDF) canv_CalScd->SaveAs((whichDir+"/CalibrationCurveSecondPol_"+cuts+""+cosTh+".pdf").c_str());
  canv_CalScd->Write();

  //------  SecondPol Relative Calibration  ------//
  TCanvas* canv_RelCalScd = new TCanvas("RelativeCalibrationCurve_SecondFit","RelativeCalibrationCurve_SecondFit"); canv_RelCalScd->cd();
  gr_RelCalScd->Draw("A*");
  gr_RelCalScd->SetTitle("Relative deviation from expected minimum versus used value of MG sample (using pol2 on reduced points)");
  gr_RelCalScd->GetXaxis()->SetTitle("gR value used for MG sample generation");
  gr_RelCalScd->GetYaxis()->SetTitle("Fit Min - expected min relative to Fit min");
  //gr_LineCal->SetLineColor(3); gr_LineCal->Draw("C"); 
  if(saveAsPDF) canv_RelCalScd->SaveAs((whichDir+"/RelativeCalibrationCurveSecondPol_"+cuts+""+cosTh+".pdf").c_str());
  canv_RelCalScd->Write();
}
