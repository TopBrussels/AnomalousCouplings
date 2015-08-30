#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>

/////////////////////////////////////////////////////////////
// Specify whether the stacked canvasses have to be stored //
bool storeSplittedCanvas = true; 
std::string SplittedDir = "Events/MTop_MGSampleCreatedWith172_SingleGausTF_10000Evts_Narrow/SplittedCanvasses"; 
/////////////////////////////////////////////////////////////

std::string VarValues[] = {"m_{top} = 171","m_{top} = 172","m_{top} = 173","m_{top} = 174","m_{top} = 175"}; 
double Var[] = {171.0,172.0,173.0,174.0,175.0}; 
double MGXS[] = {10.70485,10.8257,10.96469,11.08428,11.22448}; 
double MGXSCut[] = {10.70485,10.8257,10.96469,11.08428,11.22448}; 
int xBin = 5; 
float xLow = 170.5; 
float xHigh = 175.5; 
int xMinValue[] = {5,4,2}; 
std::string KinVar = "m_{top}"; 
int VarWindow = 3; 
int xPos[] = {3,4}; 
int xNeg[] = {1,0}; 
std::string title = "MGSample_MTop"; 
TFile* Tfile = new TFile("Events/MTop_MGSampleCreatedWith172_SingleGausTF_10000Evts_Narrow/FitDeviation_MGSample_MTop_1000Evts.root","RECREATE"); 

//ROOT file to store the Fit functions --> Will fasten the study of the cut-influences ...
TFile* file_FitDist = new TFile("Events/MTop_MGSampleCreatedWith172_SingleGausTF_10000Evts_Narrow/FitDistributions_MGSample_MTop_1000Evts.root","RECREATE"); 
TDirectory *dir_OriginalLL = file_FitDist->mkdir("OriginalLL"),        *dir_OriginalLLXS = file_FitDist->mkdir("OriginalLL_XS"),        *dir_OriginalLLAcc = file_FitDist->mkdir("OriginalLL_Acc");
TDirectory *dir_FirstFit = file_FitDist->mkdir("FirstPolynomialFit"),  *dir_FirstFitXS = file_FitDist->mkdir("FirstPolynomialFit_XS"),  *dir_FirstFitAcc = file_FitDist->mkdir("FirstPolynomialFit_Acc");
TDirectory *dir_SecondFit = file_FitDist->mkdir("SecondPolynomialFit"),*dir_SecondFitXS = file_FitDist->mkdir("SecondPolynomialFit_XS"),*dir_SecondFitAcc = file_FitDist->mkdir("SecondPolynomialFit_Acc");

const int NrConfigs = 5; 
const int nEvts = 1000; 
const int NrToDel = 1; 
int NrRemaining = NrConfigs-NrToDel;
std::string sNrCanvas ="0";
std::string sNrRemaining = ""; std::stringstream ssNrRemaining; 

//Information for the stackedCanvas division!
int NrCanvas = 0, xDivide = 5, yDivide = 4;

TF1 *polFit_AllPoints, *polFit_ReducedPoints;
TF1 *polFit_AllPointsPos, *polFit_AllPointsNeg;
TH1F *h_FitDeviation[NrConfigs], *h_FitDeviationRel[NrConfigs];
TH1F *h_PointsDelByFitDev[3], *h_PointsDelByFitDevRel[3];
TDirectory *dir_FitDevDelete = Tfile->mkdir("PointsDeletedByFitDev"), *dir_RelFitDevDelete = Tfile->mkdir("PointsDeletedByRelFitDev");
TH1F *h_ChiSquaredFirstFit[3], *h_ChiSquaredSecondFit[3];

//Method to sort a pair based on the second element!
struct sort_pred {
  bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
  }
};

void PaintOverflow(TH1F *h, TFile *FileToWrite, std::string dirName){ //TDirectory *dir){// const char* path = 0){     // This function draws the histogram h with an extra bin for overflows
  Int_t nx    = h->GetNbinsX()+1;
  Double_t x1 = h->GetBinLowEdge(1), bw = h->GetBinWidth(nx), x2 = h->GetBinLowEdge(nx)+bw;

  //Define a temporary histogram having an extra bin for overflows
  char newTitle[100], newName[100];
  strcpy(newTitle,h->GetTitle()); strcat(newTitle," (under- and overflow added)" );
  strcpy(newName,h->GetName());  strcat(newName, "_Flow");
  TH1F *h_tmp = new TH1F(newName, newTitle, nx, x1, x2);

  // Fill the new histogram including the extra bin for overflows
  for (Int_t i=1; i<=nx; i++)
    h_tmp->Fill(h_tmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows
  h_tmp->Fill(x1-1, h->GetBinContent(0));

  // Restore the number of entries
  h_tmp->SetEntries(h->GetEntries());

  //Set the correct path to save the file
  const char* path = 0;
  if(dirName != ""){
    TDirectory *dir = FileToWrite->GetDirectory(dirName.c_str());
    if (!dir)
      dir = FileToWrite->mkdir(dirName.c_str());
    path = dir->GetName();
  }

  FileToWrite->cd(path);
  h_tmp->Write();  
  FileToWrite->cd();  //Reset to general directory! 
}                  

void calculateFit(TH1F *h_LogLik, string EvtNumber, std::string Type, int evtCounter, TCanvas *canv_SplittedLL){
  file_FitDist->cd();
  if(evtCounter == 1 && sNrRemaining == "") ssNrRemaining << NrRemaining; sNrRemaining = ssNrRemaining.str();

  double LogLikelihood[NrConfigs];
  int TypeNr;
  std::string TypeName[3] = {"no norm","XS norm","Acc norm"};
  for(int ii = 0; ii < NrConfigs; ii++)
    LogLikelihood[ii] = h_LogLik->GetBinContent(h_LogLik->FindBin(Var[ii]));

  if(Type == ""){         dir_OriginalLL->cd();    TypeNr = 0;}
  else if(Type == "XS"){  dir_OriginalLLXS->cd();  TypeNr = 1;}
  else if(Type == "Acc"){ dir_OriginalLLAcc->cd(); TypeNr = 2;}
  h_LogLik->Write();
  std::string YAxisTitle = "-ln(L) value ("+TypeName[TypeNr]+" -- evt "+EvtNumber+")";

  //Set name of chisquared distributions!
  if(evtCounter == 1){
    h_ChiSquaredFirstFit[TypeNr]  = new TH1F(("ChiSquared_"+Type+"FirstFit").c_str(), ("Distribution of the chi-squared after the fit on all the points (norm = "+Type+")").c_str(),200,0,0.005);
    h_ChiSquaredSecondFit[TypeNr] = new TH1F(("ChiSquared_"+Type+"SecondFit").c_str(),("Distribution of the chi-squared after the fit on the reduced points (norm = "+Type+")").c_str(),200,0,0.005);
    h_PointsDelByFitDev[TypeNr]    = new TH1F(("PointsDelBy"+Type+"FitDev").c_str(),   ("Overview of deleted points due to largest FitDeviation (norm = "+Type+")").c_str(),xBin,xLow,xHigh);
    h_PointsDelByFitDevRel[TypeNr] = new TH1F(("PointsDelBy"+Type+"FitDevRel").c_str(),("Overview of deleted points due to largest relative FitDeviation (norm = "+Type+")").c_str(),xBin,xLow,xHigh);
  }
 
  polFit_AllPoints = new TF1(("polFit"+Type+"_AllPoints_Evt"+EvtNumber).c_str(),"pol2",Var[0],Var[NrConfigs-1]); 
  TGraph* gr_LnLik = new TGraph(NrConfigs,Var, LogLikelihood);
  gr_LnLik->Fit(polFit_AllPoints,"Q");
  h_ChiSquaredFirstFit[TypeNr]->Fill(polFit_AllPoints->GetChisquare());
  if(Type == "")         dir_FirstFit->cd();
  else if(Type == "XS")  dir_FirstFitXS->cd();
  else if(Type == "Acc") dir_FirstFitAcc->cd();
  polFit_AllPoints->Write();

  std::vector<std::pair<int, double> > FitDeviation, FitDeviationRel;
  double LogLikFit[NrConfigs] = {-9999}; 
  for(int iConfig = 0; iConfig < NrConfigs; iConfig++){
    if(polFit_AllPoints->GetNpar() == 3)
      LogLikFit[iConfig] = polFit_AllPoints->GetParameter(0)+polFit_AllPoints->GetParameter(1)*Var[iConfig]+polFit_AllPoints->GetParameter(2)*Var[iConfig]*Var[iConfig];
    else if(polFit_AllPoints->GetNpar() == 5)
      LogLikFit[iConfig] = polFit_AllPoints->GetParameter(0)+polFit_AllPoints->GetParameter(1)*Var[iConfig]+polFit_AllPoints->GetParameter(2)*Var[iConfig]*Var[iConfig]+polFit_AllPoints->GetParameter(3)*pow(Var[iConfig],3)+polFit_AllPoints->GetParameter(4)*pow(Var[iConfig],4)+polFit_AllPoints->GetParameter(5)*pow(Var[iConfig],5);
    FitDeviation.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig]) ) );
    FitDeviationRel.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig])/LogLikelihood[iConfig] ) );
  }
  //Sort the fitdeviation values depending on the second value!
  std::sort(FitDeviation.begin(), FitDeviation.end(), sort_pred() );
  std::sort(FitDeviationRel.begin(), FitDeviationRel.end(), sort_pred() );
  
  //Now loop again over all configurations, plot the FitDeviation in sorted order and save the configNr's which should be excluded 
  std::vector<int> FitDevPointsToDel, FitDevRelPointsToDel;
  for(int itSortedConfig = NrConfigs-1; itSortedConfig >= 0 ; itSortedConfig--){  //Looping from high to low values of the deviation!
    h_FitDeviation[itSortedConfig]->Fill(FitDeviation[itSortedConfig].second);
    h_FitDeviationRel[itSortedConfig]->Fill(FitDeviationRel[itSortedConfig].second);
    
    //Store the 'NrToDel' points which need to be excluded from the TGraph!
    if(FitDevPointsToDel.size() < NrToDel){ FitDevPointsToDel.push_back(FitDeviation[itSortedConfig].first); h_PointsDelByFitDev[TypeNr]->Fill(Var[FitDeviation[itSortedConfig].first]);}
    if(FitDevRelPointsToDel.size() < NrToDel){FitDevRelPointsToDel.push_back(FitDeviationRel[itSortedConfig].first);h_PointsDelByFitDevRel[TypeNr]->Fill(Var[FitDeviationRel[itSortedConfig].first]);}
  }

  //Create new arrays with the reduced information!
  double ReducedVar[NrConfigs-NrToDel], ReducedLogLik[NrConfigs-NrToDel];
  int iCounter = 0;
  for(int iConf = 0; iConf < NrConfigs; iConf++){
    if(!(std::find(FitDevPointsToDel.begin(), FitDevPointsToDel.end(), iConf) != FitDevPointsToDel.end()) ){
      ReducedVar[iCounter] = Var[iConf];
      ReducedLogLik[iCounter] = LogLikelihood[iConf];
      iCounter++;
    }
  }

  //Define new TGraph and fit again
  TGraph* gr_ReducedLnLik = new TGraph(NrConfigs-NrToDel, ReducedVar, ReducedLogLik);
  polFit_ReducedPoints = new TF1(("polFit"+Type+"_"+sNrRemaining+"ReducedPoints_Evt"+EvtNumber).c_str(),"pol2",Var[0],Var[NrConfigs-1]); 
  gr_ReducedLnLik->Fit(polFit_ReducedPoints,"Q"); 
  h_ChiSquaredSecondFit[TypeNr]->Fill(polFit_ReducedPoints->GetChisquare());   //As expected NDF is always equal to NrConfigs-NrToDel-3 (= nr params needed to define a parabola)
  if(Type == "")         dir_SecondFit->cd();
  else if(Type == "XS")  dir_SecondFitXS->cd();
  else if(Type == "Acc") dir_SecondFitAcc->cd();
  stringstream ssChiSq; ssChiSq << polFit_ReducedPoints->GetChisquare(); string sChiSq = ssChiSq.str();
  h_LogLik->SetTitle(("Polynomial 2nd fit ("+TypeName[TypeNr]+" -- "+sNrRemaining+" used points -- #chi^{2} = "+sChiSq+")").c_str());
  polFit_ReducedPoints->Write();
  Tfile->cd();

  if( storeSplittedCanvas == true){
    h_LogLik->GetYaxis()->SetTitle(YAxisTitle.c_str());
    h_LogLik->GetYaxis()->SetTitleOffset(1.4);
    h_LogLik->GetXaxis()->SetTitle(KinVar.c_str());
    canv_SplittedLL->cd(evtCounter - (xDivide*yDivide*NrCanvas) ); h_LogLik->Draw("p"); polFit_ReducedPoints->Draw("same"); canv_SplittedLL->Update();
  }

  delete gr_ReducedLnLik;
  delete gr_LnLik;
}

void fitDeviationMacro(){
  
  int xMin = xMinValue[VarWindow-1];
  float xStep[] = {Var[int(xNeg[0])]-Var[int(xNeg[1])], Var[int(xMin)]-Var[int(xNeg[0])], Var[int(xPos[0])]-Var[int(xMin)], Var[int(xPos[1])]-Var[int(xPos[0])] };
  std::cout << " Value of xStep is : " << xStep[0] << std::endl;

  //--- Initialize histograms ---//
  TH1F* LnLikAll    = new TH1F("LnLikAll",   ("Distribution of -ln(L) using all events (no norm -- "+title+" evts)").c_str(), xBin,xLow,xHigh);
  TH1F* LnLikXSAll  = new TH1F("LnLikXSAll", ("Distribution of -ln(L) using all events (XS norm -- "+title+" evts)").c_str(), xBin,xLow,xHigh);
  TH1F* LnLikAccAll = new TH1F("LnLikAccAll",("Distribution of -ln(L) using all events (Acc norm -- "+title+" evts)").c_str(),xBin,xLow,xHigh);

  TH1F* YPlusGausTest   =new TH1F("YPlusGausTest",   ("Comparison of fit deviation in "+VarValues[xPos[1]]+" and "+VarValues[xPos[0]]+" (no norm -- "+title+" evts)").c_str(), 250,-0.15,0.15);
  TH1F* YPlusGausTestXS =new TH1F("YPlusGausTestXS", ("Comparison of fit deviation in "+VarValues[xPos[1]]+" and "+VarValues[xPos[0]]+" (XS norm -- "+title+" evts)").c_str(), 250,-0.15,0.15);
  TH1F* YPlusGausTestAcc=new TH1F("YPlusGausTestAcc",("Comparison of fit deviation in "+VarValues[xPos[1]]+" and "+VarValues[xPos[0]]+" (Acc norm -- "+title+" evts)").c_str(),250,-0.15,0.15);
  TH1F* YPlusGausTestPosScdDer    = new TH1F("YPlusGausTestPosScdDer",   ("Gaussian fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestXSPosScdDer  = new TH1F("YPlusGausTestXSPosScdDer", ("Gaussian fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestAccPosScdDer = new TH1F("YPlusGausTestAccPosScdDer",("Gaussian fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),250,-0.1,0.1);
  TH1F* YPlusGausTestNegScdDer    = new TH1F("YPlusGausTestNegScdDer",   ("Gaussian fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestXSNegScdDer  = new TH1F("YPlusGausTestXSNegScdDer", ("Gaussing fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(), 250,-0.1,0.1);
  TH1F* YPlusGausTestAccNegScdDer = new TH1F("YPlusGausTestAccNegScdDer",("Gaussian fit deviation using "+VarValues[xPos[1]]+" & "+VarValues[xPos[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),250,-0.1,0.1);

  TH1F* YPlus    = new TH1F("YPlus",   ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YPlusXS  = new TH1F("YPlusXS", ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YPlusAcc = new TH1F("YPlusAcc",("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.25,0.25);
  TH1F* YPlusPosScdDer    = new TH1F("YPlusPosScdDer",   ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusXSPosScdDer  = new TH1F("YPlusXSPosScdDer", ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusAccPosScdDer = new TH1F("YPlusAccPosScdDer",("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  TH1F* YPlusNegScdDer    = new TH1F("YPlusNegScdDer",   ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusXSNegScdDer  = new TH1F("YPlusXSNegScdDer", ("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YPlusAccNegScdDer = new TH1F("YPlusAccNegScdDer",("Deviation from parabolic fit for "+VarValues[xPos[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);

  TH1F* YPlusPlus    = new TH1F("YPlusPlus",   ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YPlusPlusXS  = new TH1F("YPlusPlusXS", ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YPlusPlusAcc = new TH1F("YPlusPlusAcc",("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.5,0.5);
  TH1F* YPlusPlusPosScdDer    = new TH1F("YPlusPlusPosScdDer",   ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YPlusPlusXSPosScdDer  = new TH1F("YPlusPlusXSPosScdDer", ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YPlusPlusAccPosScdDer = new TH1F("YPlusPlusAccPosScdDer",("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  TH1F* YPlusPlusNegScdDer    = new TH1F("YPlusPlusNegScdDer",   ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YPlusPlusXSNegScdDer  = new TH1F("YPlusPlusXSNegScdDer", ("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YPlusPlusAccNegScdDer = new TH1F("YPlusPlusAccNegScdDer",("Deviation from parabolic fit for "+VarValues[xPos[1]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);

  TH1F* YMin    = new TH1F("YMin",   ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YMinXS  = new TH1F("YMinXS", ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.25,0.25);
  TH1F* YMinAcc = new TH1F("YMinAcc",("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.25,0.25);
  TH1F* YMinPosScdDer    = new TH1F("YMinPosScdDer",   ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinXSPosScdDer  = new TH1F("YMinXSPosScdDer", ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinAccPosScdDer = new TH1F("YMinAccPosScdDer",("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  TH1F* YMinNegScdDer    = new TH1F("YMinNegScdDer",   ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinXSNegScdDer  = new TH1F("YMinXSNegScdDer", ("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.2,0.2);
  TH1F* YMinAccNegScdDer = new TH1F("YMinAccNegScdDer",("Deviation from parabolic fit for "+VarValues[xNeg[0]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.2,0.2);
  
  TH1F* YMinMin    = new TH1F("YMinMin",   ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YMinMinXS  = new TH1F("YMinMinXS", ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.5,0.5);
  TH1F* YMinMinAcc = new TH1F("YMinMinAcc",("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.5,0.5);
  TH1F* YMinMinPosScdDer    = new TH1F("YMinMinPosScdDer",   ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (no norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinXSPosScdDer  = new TH1F("YMinMinXSPosScdDer", ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (XS norm -- outer 2nd der > 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinAccPosScdDer = new TH1F("YMinMinAccPosScdDer",("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (Acc norm -- outer 2nd der > 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  TH1F* YMinMinNegScdDer    = new TH1F("YMinMinNegScdDer",   ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (no norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinXSNegScdDer  = new TH1F("YMinMinXSNegScdDer", ("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (XS norm -- outer 2nd der < 0 -- "+title+" evts)").c_str() ,150,-0.3,0.3);
  TH1F* YMinMinAccNegScdDer = new TH1F("YMinMinAccNegScdDer",("Deviation from parabolic fit for "+VarValues[xNeg[1]]+" (Acc norm -- outer 2nd der < 0 -- "+title+" evts)").c_str(),150,-0.3,0.3);
  
  TH1F* YRelPlus    = new TH1F("YRelPlus",   ("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,250,-0.1,0.1);
  TH1F* YRelPlusXS  = new TH1F("YRelPlusXS", ("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,250,-0.1,0.1);
  TH1F* YRelPlusAcc = new TH1F("YRelPlusAcc",("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),250,-0.1,0.1);
  
  TH1F* YRelMin    = new TH1F("YRelMin",   ("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (no norm -- "+title+" evts)").c_str() ,150,-0.1,0.1);
  TH1F* YRelMinXS  = new TH1F("YRelMinXS", ("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (XS norm -- "+title+" evts)").c_str() ,150,-0.1,0.1);
  TH1F* YRelMinAcc = new TH1F("YRelMinAcc",("Relative deviation from parabolic fit for "+VarValues[xNeg[0]]+" (Acc norm -- "+title+" evts)").c_str(),150,-0.1,0.1);
  
  TH1F* FstDer    = new TH1F("FirstDer",   ("First derivative of -ln(likelihood) distribution -- "+title+" evts").c_str(), 4,0,4);
  std::ostringstream ssStepSize; ssStepSize << xStep[0]; std::string sStepSize(ssStepSize.str());
  FstDer->GetXaxis()->SetBinLabel(1,("(LL("+VarValues[xNeg[1]]+")-LL("+VarValues[xNeg[0]]+"))/"+sStepSize).c_str());
  FstDer->GetXaxis()->SetBinLabel(2,("(LL("+VarValues[xNeg[0]]+")-LL("+VarValues[xMin]+"))/"+sStepSize).c_str());
  FstDer->GetXaxis()->SetBinLabel(3,("(LL("+VarValues[xMin]+")-LL("+VarValues[xPos[0]]+"))/"+sStepSize).c_str());
  FstDer->GetXaxis()->SetBinLabel(4,("(LL("+VarValues[xPos[1]]+")-LL("+VarValues[xPos[1]]+"))/"+sStepSize).c_str());
  //TH1F* FstDerXS  = new TH1F("FirstDerivativeXS", ("First derivative of -ln(likelihood) distribution (XS norm -- "+title+" evts)").c_str(), 5,-0.25,0.25);
  //TH1F* FstDerAcc = new TH1F("FirstDerivativeAcc",("First derivative of -ln(likelihood) distribution (Acc norm -- "+title+" evts)").c_str(),5,-0.25,0.25);
  TH1F* ScdDerInner    = new TH1F("SecondDerivativeInner",   ("Second derivative of -ln(likelihood) distribution (no norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSInner  = new TH1F("SecondDerivativeXSInner", ("Second derivative of -ln(likelihood) distribution (XS norm -- using inner points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccInner = new TH1F("SecondDerivativeAccInner",("Second derivative of -ln(likelihood) distribution (Acc norm -- using inner points -- "+title+" evts)").c_str(),250,-5,5);
  TH1F* ScdDerOuter    = new TH1F("SecondDerivativeOuter",   ("Second derivative of -ln(likelihood) distribution (no norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerXSOuter  = new TH1F("SecondDerivativeXSOuter", ("Second derivative of -ln(likelihood) distribution (XS norm -- using outer points -- "+title+" evts)").c_str(), 250,-5,5);
  TH1F* ScdDerAccOuter = new TH1F("SecondDerivativeAccOuter",("Second derivative of -ln(likelihood) distribution (Acc norm -- using outer points -- "+title+" evts)").c_str(),250,-5,5);
  TH2F* ScdDerScatter    = new TH2F("ScdDerScatterPlot",   ("Second derivative of -ln(L) using inner points versus using outer points (no norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerXSScatter  = new TH2F("ScdDerXSScatterPlot", ("Second derivative of -ln(L) using inner points versus using outer points (XS norm -- "+title+" evts)").c_str(), 250,-5,5,250,-5,5);
  TH2F* ScdDerAccScatter = new TH2F("ScdDerAccScatterPlot",("Second derivative of -ln(L) using inner points versus using outer points (Acc norm -- "+title+" evts)").c_str(),250,-5,5,250,-5,5);
  
  //TH1F* ProcContrDist = new TH1F("ProcentualContribution","Procentual contribution of each of the events to the total likelihood",250,0,1);
  
  //TH1F* TotalFitDevDist    = new TH1F("TotalFitDeviation",   ("Sum of difference between likelihood and fit value for each point in range (no norm -- "+title+" events)").c_str(), 250,0,5);
  //TH1F* TotalFitDevXSDist  = new TH1F("TotalFitXSDeviation", ("Sum of difference between likelihood and fit value for each point in range (XS norm -- "+title+" events)").c_str(), 250,0,5);
  //TH1F* TotalFitDevAccDist = new TH1F("TotalFitAccDeviation",("Sum of difference between likelihood and fit value for each point in range (Acc norm -- "+title+" events)").c_str(),250,0,5);
  //TH1F* TotalFctDevDist    = new TH1F("TotalFctDeviation",   ("Sum of difference between likelihood and function value for each point in range (no norm -- "+title+" events)").c_str(), 250,0,5);
  //TH1F* TotalFctDevXSDist  = new TH1F("TotalFctXSDeviation", ("Sum of difference between likelihood and function value for each point in range (XS norm -- "+title+" events)").c_str(), 250,0,5);
  //TH1F* TotalFctDevAccDist = new TH1F("TotalFctAccDeviation",("Sum of difference between likelihood and function value for each point in range (Acc norm -- "+title+" events)").c_str(),250,0,5);

  TH1F *LnLikDist = 0, *LnLikXSDist = 0, *LnLikAccDist = 0; 
  TCanvas *LnLikCanv = 0, *LnLikXSCanv = 0, *LnLikAccCanv = 0;
  TGraph *LnLikFctOuter = 0, *LnLikXSFctOuter = 0, *LnLikAccFctOuter = 0;
  
  //TDirectory* FitComp = Tfile->mkdir("FitComparison");
  //TDirectory* LnLikAccDirVarVsUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanAvgUnc");
  //TDirectory* LnLikAccDirVarVsDUnc = Tfile->mkdir("LnLikAccDist_VarLargerThanTwiceAvgUnc");
  TDirectory* LnLikDir = Tfile->mkdir("LnLikDist");       TDirectory* LnLikStackDir = LnLikDir->mkdir("LnLikStackDist");
  TDirectory* LnLikXSDir = Tfile->mkdir("LnLikXSDist");   TDirectory* LnLikXSStackDir = LnLikXSDir->mkdir("LnLikXSStackDist");
  TDirectory* LnLikAccDir = Tfile->mkdir("LnLikAccDist"); TDirectory* LnLikAccStackDir = LnLikAccDir->mkdir("LnLikAccStackDist");
  TDirectory* FstDerDir = Tfile->mkdir("FirstDerivativeDist");

  //-----------------------------------//
  //----  Specific cut histograms  ----//
  //-----------------------------------//
  TH1F* LLSmallFctDevDist = new TH1F("LLSmallFctDeviation","-ln(L) when sum of deviations between -ln(L)- and function-value < 0.5 (acc norm -- outer points used)",xBin,xLow,xHigh);

  TH1F* LLPosScdDerDistBoth     = new TH1F("LLPosScdDerBoth",    "-ln(L) when both 2nd derivatives > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistBoth   = new TH1F("LLXSPosScdDerBoth",  "-ln(L) when both 2nd derivatives > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistBoth  = new TH1F("LLAccPosScdDerBoth", "-ln(L) when both 2nd derivatives > 0 (Acc norm -- ",xBin,xLow,xHigh);
  TH1F* LLPosScdDerDistInner    = new TH1F("LLPosScdDerInner",   "-ln(L) when inner 2nd derivative > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistInner  = new TH1F("LLXSPosScdDerInner", "-ln(L) when inner 2nd derivative > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistInner = new TH1F("LLAccPosScdDerInner","-ln(L) when inner 2nd derivative > 0 (Acc norm -- ",xBin,xLow,xHigh);
  TH1F* LLPosScdDerDistOuter    = new TH1F("LLPosScdDerOuter",   "-ln(L) when outer 2nd derivative > 0 (no norm -- ", xBin,xLow,xHigh);
  TH1F* LLXSPosScdDerDistOuter  = new TH1F("LLXSPosScdDerOuter", "-ln(L) when outer 2nd derivative > 0 (XS norm -- ", xBin,xLow,xHigh);
  TH1F* LLAccPosScdDerDistOuter = new TH1F("LLAccPosScdDerOuter","-ln(L) when outer 2nd derivative > 0 (Acc norm -- ",xBin,xLow,xHigh);

  //Initialize the event counters:
  vector<int> EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc;
  int EvtsWithPosScdDerInner = 0, EvtsWithPosScdDerXSInner = 0, EvtsWithPosScdDerAccInner = 0;
  int EvtsWithPosScdDerOuter = 0, EvtsWithPosScdDerXSOuter = 0, EvtsWithPosScdDerAccOuter = 0;
  int EvtsWithPosScdDerBoth  = 0, EvtsWithPosScdDerXSBoth  = 0, EvtsWithPosScdDerAccBoth  = 0;
  //vector<int> EvtsPosScdDerInner, EvtsPosScdDerXSInner, EvtsPosScdDerAccInner;
  //vector<int> EvtsPosScdDerOuter, EvtsPosScdDerXSOuter, EvtsPosScdDerAccOuter;
  vector<int> EvtsWithSmallFctDev;
  //double weightsValue[nEvts][NrConfigs];
  int consEvts = 0;
  TCanvas *StackCanvasLL = 0, *StackCanvasLLXS = 0, *StackCanvasLLAcc = 0; 

  double aHat[2] = {0.0}, aHatXS[2] = {0.0}, aHatAcc[2] = {0.0}, bHat[2] = {0.0}, bHatXS[2] = {0.0}, bHatAcc[2] = {0.0}, cHat[2] = {0.0}, cHatXS[2] = {0.0}, cHatAcc[2] = {0.0};
  double LnLik[NrConfigs] = {0.0}, LnLikXS[NrConfigs] = {0.0}, LnLikAcc[NrConfigs] = {0.0};        

  //--- Read all likelihood values ! ---//
  std::ifstream ifs ("Events/MTop_MGSampleCreatedWith172_SingleGausTF_10000Evts_Narrow/weights_NoZero.out", std::ifstream::in); 
  std::cout << " Value of ifs : " << ifs.eof() << std::endl;
  std::string line;
  int evt,config,tf;
  double weight, weightUnc;
  while( std::getline(ifs,line) && consEvts < nEvts){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc){
      if(config == 1 && ((consEvts+1) % 50 == 0) ) std::cout << " Looking at event : " << consEvts+1 << std::endl;
      stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();
      stringstream ssConfig; ssConfig << config; string sConfig = ssConfig.str();

      //Initialize the fitDeviation histograms (array of TH1F, one for each configuration ...)
      if(consEvts == 0){
        h_FitDeviation[config-1] = new TH1F(("FitDeviation_LowestConfig"+sConfig).c_str(),("FitDeviation histogram for "+sConfig+"st lowest value").c_str(), 200, 0, 0.1);
        h_FitDeviationRel[config-1] = new TH1F(("FitDeviationRelative_LowestConfig"+sConfig).c_str(), ("FitDeviationRelative histogram for "+sConfig+"st lowest value").c_str(),200,0,0.005);
      }

      //--- Initialize the event-per-event variables! ---//
      if( config == 1){

        LnLikDist = new TH1F(("LnLik_Evt"+sEvt).c_str(),("LnLik distribution for event "+sEvt+" -- "+title+" evts").c_str(),xBin,xLow,xHigh);
        LnLikDist->SetMarkerStyle(20); LnLikDist->SetLineColor(1); LnLikDist->SetMarkerColor(1); LnLikDist->SetMarkerSize(1.2);
        LnLikXSDist = new TH1F(("LnLikXS_Evt"+sEvt).c_str(),("LnLikXS distribution for event "+sEvt+" -- "+title+" evts").c_str(),xBin,xLow,xHigh);
        LnLikXSDist->SetMarkerStyle(21); LnLikXSDist->SetLineColor(3); LnLikXSDist->SetMarkerColor(3); LnLikXSDist->SetMarkerSize(1.2);
        LnLikAccDist = new TH1F(("LnLikAcc_Evt"+sEvt).c_str(),("LnLikAcc distribution for event "+sEvt+" -- "+title+" evts").c_str(),xBin,xLow,xHigh);
        LnLikAccDist->SetMarkerStyle(22); LnLikAccDist->SetLineColor(4); LnLikAccDist->SetMarkerColor(4); LnLikAccDist->SetMarkerSize(1.2);
      }
      //Lik[config-1] = weight;         LikXS[config-1] = weight/MGXS[config-1];              LikAcc[config-1] = weight/MGXSCut[config-1];
      LnLik[config-1] = -log(weight); LnLikXS[config-1] = -log(weight)+log(MGXS[config-1]); LnLikAcc[config-1] = -log(weight)+log(MGXSCut[config-1]);
      //weightsValue[evt-1][config-1] = weight;

      //---  Fill the LnLik histograms for each event and for all events together  ---//
      LnLikDist->SetBinContent(   LnLikDist->FindBin(Var[config-1]),    LnLik[config-1]);
      LnLikXSDist->SetBinContent( LnLikXSDist->FindBin(Var[config-1]),  LnLikXS[config-1]);
      LnLikAccDist->SetBinContent(LnLikAccDist->FindBin(Var[config-1]), LnLikAcc[config-1]);
      LnLikAll->SetBinContent(   LnLikAll->FindBin(Var[config-1]),    LnLikAll->GetBinContent(LnLikAll->FindBin(Var[config-1]))       + LnLik[config-1]   );
      LnLikXSAll->SetBinContent( LnLikXSAll->FindBin(Var[config-1]),  LnLikXSAll->GetBinContent(LnLikXSAll->FindBin(Var[config-1]))   + LnLikXS[config-1] );
      LnLikAccAll->SetBinContent(LnLikAccAll->FindBin(Var[config-1]), LnLikAccAll->GetBinContent(LnLikAccAll->FindBin(Var[config-1])) + LnLikAcc[config-1]);

      //---  Only perform the fit after all configurations are considered!  ---//
      if( config == NrConfigs){
        consEvts++;   //Count the number of full events!
        double LnLikFunction[2][NrConfigs], LnLikXSFunction[2][NrConfigs], LnLikAccFunction[2][NrConfigs];
        vector<double> FctDevOuter, RelFctDevOuter,FctDevXSOuter, RelFctDevXSOuter , FctDevAccOuter, RelFctDevAccOuter;
        double TotalFctDevOuter = 0, TotalRelFctDevOuter = 0, TotalFctDevXSOuter = 0, TotalRelFctDevXSOuter = 0, TotalFctDevAccOuter = 0, TotalRelFctDevAccOuter = 0;

        //Save xDivide*yDivide of these histograms in one TCanvas!
        if( storeSplittedCanvas == true){
          if(consEvts == 1){
            StackCanvasLL    = new TCanvas("StackCanvasLL_Nr0",   "StackCanvasLL");    StackCanvasLL->Divide(xDivide,yDivide);    
            StackCanvasLLXS  = new TCanvas("StackCanvasLLXS_Nr0", "StackCanvasLLXS");  StackCanvasLLXS->Divide(xDivide,yDivide);  
            StackCanvasLLAcc = new TCanvas("StackCanvasLLAcc_Nr0","StackCanvasLLAcc"); StackCanvasLLAcc->Divide(xDivide,yDivide); 
          }
        }

        //-- Send the array containing the log(weights) to the predefined function to define the TGraph, fit this, detect the deviation points and fit again! --//
	calculateFit(LnLikDist,   sEvt,"",    consEvts, StackCanvasLL);
	calculateFit(LnLikXSDist, sEvt,"XS",  consEvts, StackCanvasLLXS);
        calculateFit(LnLikAccDist,sEvt,"Acc", consEvts, StackCanvasLLAcc);

        if( storeSplittedCanvas == true){
          if( consEvts == (xDivide*yDivide*(NrCanvas+1)) || consEvts == nEvts){
            //std::string NameTest = (StackCanvasLLAcc->GetTitle()+"_Nr"+sNrCanvas).c_str();   //Why can't the Title or name be used ... ?
            StackCanvasLL->Print((SplittedDir+"/StackCanvasLL_Nr"+sNrCanvas+".pdf").c_str());       LnLikStackDir->cd();    StackCanvasLL->Write();
            StackCanvasLLXS->Print((SplittedDir+"/StackCanvasLLXS_Nr"+sNrCanvas+".pdf").c_str());   LnLikXSStackDir->cd();  StackCanvasLLXS->Write();
            StackCanvasLLAcc->Print((SplittedDir+"/StackCanvasLLAcc_Nr"+sNrCanvas+".pdf").c_str()); LnLikAccStackDir->cd(); StackCanvasLLAcc->Write();
            if( consEvts != nEvts){
              NrCanvas++; stringstream ssNrCanvas; ssNrCanvas << NrCanvas; sNrCanvas = ssNrCanvas.str();
              StackCanvasLL    = new TCanvas(("StackCanvasLL_Nr"+sNrCanvas).c_str(),   "SplittedCanvasLL");    StackCanvasLL->Divide(xDivide,yDivide);   
              StackCanvasLLXS  = new TCanvas(("StackCanvasLLXS_Nr"+sNrCanvas).c_str(), "SplittedCanvasLLXS");  StackCanvasLLXS->Divide(xDivide,yDivide); 
              StackCanvasLLAcc = new TCanvas(("StackCanvasLLAcc_Nr"+sNrCanvas).c_str(),"SplittedCanvasLLAcc"); StackCanvasLLAcc->Divide(xDivide,yDivide);
            }
          }
        }//Draw the stacked canvasses!        

        for(int ii=0; ii < 2; ii++){
          cHat[ii] = LnLik[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLik[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatXS[ii] = LnLikXS[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikXS[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          cHatAcc[ii] = LnLikAcc[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikAcc[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          bHat[ii] = LnLik[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLik[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatXS[ii] = LnLikXS[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikXS[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));
          bHatAcc[ii] = LnLikAcc[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikAcc[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]));

          aHat[ii] = LnLik[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatXS[ii] = LnLikXS[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));
          aHatAcc[ii] = LnLikAcc[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]));

          for( int iConfig = 0; iConfig < NrConfigs; iConfig++){
            LnLikFunction[ii][iConfig] = aHat[ii]*Var[iConfig]*Var[iConfig]+bHat[ii]*Var[iConfig]+cHat[ii];
            LnLikXSFunction[ii][iConfig] = aHatXS[ii]*Var[iConfig]*Var[iConfig]+bHatXS[ii]*Var[iConfig]+cHatXS[ii];
            LnLikAccFunction[ii][iConfig] = aHatAcc[ii]*Var[iConfig]*Var[iConfig]+bHatAcc[ii]*Var[iConfig]+cHatAcc[ii];

            if(ii == 1){
              FctDevOuter.push_back(abs(LnLikFunction[1][iConfig]-LnLik[iConfig]));       RelFctDevOuter.push_back(FctDevOuter[iConfig]/LnLik[iConfig]);
              FctDevXSOuter.push_back(abs(LnLikXSFunction[1][iConfig]-LnLikXS[iConfig]));    RelFctDevXSOuter.push_back(FctDevXSOuter[iConfig]/LnLikXS[iConfig]);
              FctDevAccOuter.push_back(abs(LnLikAccFunction[1][iConfig]-LnLikAcc[iConfig])); RelFctDevAccOuter.push_back(FctDevAccOuter[iConfig]/LnLikAcc[iConfig]);
              TotalFctDevOuter    = TotalFctDevOuter+FctDevOuter[iConfig];       TotalRelFctDevOuter    = TotalRelFctDevOuter+RelFctDevOuter[iConfig];
              TotalFctDevXSOuter  = TotalFctDevXSOuter+FctDevXSOuter[iConfig];   TotalRelFctDevXSOuter  = TotalRelFctDevXSOuter+RelFctDevXSOuter[iConfig];
              TotalFctDevAccOuter = TotalFctDevAccOuter+FctDevAccOuter[iConfig]; TotalRelFctDevAccOuter = TotalRelFctDevAccOuter+RelFctDevAccOuter[iConfig];
            }
          }
        }//Looping over inner and outer parabolic function!

        //---  Save the LnLik distributions (event-per-event) together with the function "fit" ---//
        //-- Transformation needed in order to work with TGraphs --//
        double yLnLik[NrConfigs], yLnLikXS[NrConfigs], yLnLikAcc[NrConfigs];
        for( int i = 0; i < NrConfigs; i++){ yLnLik[i] = LnLikFunction[1][i]; yLnLikXS[i] = LnLikXSFunction[1][i]; yLnLikAcc[i] = LnLikAccFunction[1][i]; }

        stringstream ssTotalFctDevOuter;       ssTotalFctDevOuter << TotalFctDevOuter;             string sTotalFctDevOuter       = ssTotalFctDevOuter.str();
        stringstream ssTotalRelFctDevOuter;    ssTotalRelFctDevOuter << TotalRelFctDevOuter;       string sTotalRelFctDevOuter    = ssTotalRelFctDevOuter.str();
        stringstream ssTotalFctDevXSOuter;     ssTotalFctDevXSOuter << TotalFctDevXSOuter;         string sTotalFctDevXSOuter     = ssTotalFctDevXSOuter.str();
        stringstream ssTotalRelFctDevXSOuter;  ssTotalRelFctDevXSOuter << TotalRelFctDevXSOuter;   string sTotalRelFctDevXSOuter  = ssTotalRelFctDevXSOuter.str();
        stringstream ssTotalFctDevAccOuter;    ssTotalFctDevAccOuter << TotalFctDevAccOuter;       string sTotalFctDevAccOuter    = ssTotalFctDevAccOuter.str();
        stringstream ssTotalRelFctDevAccOuter; ssTotalRelFctDevAccOuter << TotalRelFctDevAccOuter; string sTotalRelFctDevAccOuter = ssTotalRelFctDevAccOuter.str();

        LnLikFctOuter    = new TGraph(NrConfigs,Var,yLnLik);    LnLikFctOuter->SetMarkerColor(2);    LnLikFctOuter->SetLineColor(2);
        LnLikXSFctOuter  = new TGraph(NrConfigs,Var,yLnLikXS);  LnLikXSFctOuter->SetMarkerColor(2);  LnLikXSFctOuter->SetLineColor(2);
        LnLikAccFctOuter = new TGraph(NrConfigs,Var,yLnLikAcc); LnLikAccFctOuter->SetMarkerColor(2); LnLikAccFctOuter->SetLineColor(2);

        LnLikFctOuter->SetTitle(("Outer fct dev = "+sTotalFctDevOuter+" & "+sTotalRelFctDevOuter).c_str());
        LnLikFctOuter->GetYaxis()->SetTitle(("-ln(L) value (no norm -- evt "+sEvt+")").c_str()); LnLikFctOuter->GetYaxis()->SetTitleOffset(1.5);
        LnLikXSFctOuter->SetTitle(("Outer fct dev = "+sTotalFctDevOuter+" & "+sTotalRelFctDevOuter).c_str());
        LnLikXSFctOuter->GetYaxis()->SetTitle(("-ln(L) value (XS norm -- evt "+sEvt+")").c_str()); LnLikXSFctOuter->GetYaxis()->SetTitleOffset(1.5);
        LnLikAccFctOuter->SetTitle(("Outer fct dev = "+sTotalFctDevOuter+" & "+sTotalRelFctDevOuter).c_str());
        LnLikAccFctOuter->GetYaxis()->SetTitle(("-ln(L) value (Acc norm -- evt "+sEvt+")").c_str()); LnLikAccFctOuter->GetYaxis()->SetTitleOffset(1.5);
        LnLikCanv =    new TCanvas(("LnLikCanv_"+sEvt).c_str(),"LnLik");      LnLikCanv->cd();   LnLikFctOuter->Draw("AC*");   LnLikDist->Draw("samep");   LnLikDir->cd();   LnLikCanv->Write();
        LnLikXSCanv =  new TCanvas(("LnLikXSCanv_"+sEvt).c_str(),"LnLikXS");  LnLikXSCanv->cd(); LnLikXSFctOuter->Draw("AC*"); LnLikXSDist->Draw("samep"); LnLikXSDir->cd(); LnLikXSCanv->Write();
        LnLikAccCanv = new TCanvas(("LnLikAccCanv_"+sEvt).c_str(),"LnLikAcc");LnLikAccCanv->cd();LnLikAccFctOuter->Draw("AC*");LnLikAccDist->Draw("samep");LnLikAccDir->cd();LnLikAccCanv->Write();

  
        //---  Calculate the first derivative distribution and save them in event-by-event plot  ---//
        FstDer->SetName(("FirstDerivative_Evt"+sEvt).c_str());
        FstDer->SetTitle(("First derivative of -ln(likelihood) distribution for event"+sEvt+" (no normalisation) -- "+title+" evts").c_str());
        for(int ii = 0; ii < FstDer->GetNbinsX(); ii++)
          FstDer->SetBinContent(ii+1, (LnLik[ii+2] - LnLik[ii+3])/xStep[ii]);
        FstDerDir->cd(); FstDer->Write();

        //---  Calculate the positive and negative deviation  ---//
        double yPlusPlus[3] = {LnLik[xPos[1]] - LnLikFunction[0][xPos[1]], LnLikXS[xPos[1]] - LnLikXSFunction[0][xPos[1]], LnLikAcc[xPos[1]] - LnLikAccFunction[0][xPos[1]]};
        double yMinMin[3]   = {LnLik[xNeg[1]] - LnLikFunction[0][xNeg[1]], LnLikXS[xNeg[1]] - LnLikXSFunction[0][xNeg[1]], LnLikAcc[xNeg[1]] - LnLikAccFunction[0][xNeg[1]]};
        double yPlus[3] = {LnLik[xPos[0]] - LnLikFunction[1][xPos[0]], LnLikXS[xPos[0]] - LnLikXSFunction[1][xPos[0]], LnLikAcc[xPos[0]] - LnLikAccFunction[1][xPos[0]]};
        double yMin[3]  = {LnLik[xNeg[0]] - LnLikFunction[1][xNeg[0]], LnLikXS[xNeg[0]] - LnLikXSFunction[1][xNeg[0]], LnLikAcc[xNeg[0]] - LnLikAccFunction[1][xNeg[0]]};
        double scdDerInner[3] = {LnLik[xNeg[0]]-2*LnLik[xMin]+LnLik[xPos[0]], LnLikXS[xNeg[0]]-2*LnLikXS[xMin]+LnLikXS[xPos[0]], LnLikAcc[xNeg[0]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[0]]};
        double scdDerOuter[3]   = {LnLik[xNeg[1]]-2*LnLik[xMin]+LnLik[xPos[1]], LnLikXS[xNeg[1]]-2*LnLikXS[xMin]+LnLikXS[xPos[1]], LnLikAcc[xNeg[1]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[1]]};

        //--- Fill the histograms  ---//
        YPlus->Fill(yPlus[0]);                                     YPlusXS->Fill(yPlus[1]);                                     YPlusAcc->Fill(yPlus[2]);
        YPlusPlus->Fill(yPlusPlus[0]);                             YPlusPlusXS->Fill(yPlusPlus[1]);                             YPlusPlusAcc->Fill(yPlusPlus[2]);
        YPlusGausTest->Fill(yPlus[0] + yPlusPlus[0]/4);            YPlusGausTestXS->Fill(yPlus[1] + yPlusPlus[1]/4);            YPlusGausTestAcc->Fill(yPlus[2] + yPlusPlus[2]/4);
        //YRelPlus->Fill(yPlus[0]/LnLik[6]);                         YRelPlusXS->Fill(yPlus[1]/LnLikXS[6]);                       YRelPlusAcc->Fill(yPlus[2]/LnLikAcc[6]);
        YMin->Fill(yMin[0]);                                       YMinXS->Fill(yMin[1]);                                       YMinAcc->Fill(yMin[2]);
        YMinMin->Fill(yMinMin[0]);                                 YMinMinXS->Fill(yMinMin[1]);                                 YMinMinAcc->Fill(yMinMin[2]);
        //YRelMin->Fill(yMin[0]/LnLik[2]);                           YRelMinXS->Fill(yMin[1]/LnLikXS[2]);                         YRelMinAcc->Fill(yMin[2]/LnLikAcc[2]);
        ScdDerInner->Fill(scdDerInner[0]);                         ScdDerXSInner->Fill(scdDerInner[1]);                         ScdDerAccInner->Fill(scdDerInner[2]);
        ScdDerOuter->Fill(scdDerOuter[0]);                         ScdDerXSOuter->Fill(scdDerOuter[1]);                         ScdDerAccOuter->Fill(scdDerOuter[2]);
        ScdDerScatter->Fill(scdDerOuter[0], scdDerInner[0]);       ScdDerXSScatter->Fill(scdDerOuter[1], scdDerInner[1]);       ScdDerAccScatter->Fill(scdDerOuter[2], scdDerInner[2]);
        //if( TotalFctDevAccOuter > 5) std::cout << "Overflow found for TotalFctDevDist : " << TotalFctDevAccOuter << std::endl;

        //-- Apply cut on YPlusGausTest --//
        if( (yPlus[0] + yPlusPlus[0]/4) <= 0.025 && (yPlus[0] + yPlusPlus[0]/4) >= -0.025) EvtsWithYPlusGausSmall.push_back(evt);
        if( (yPlus[1] + yPlusPlus[1]/4) <= 0.025 && (yPlus[1] + yPlusPlus[1]/4) >= -0.025) EvtsWithYPlusGausSmallXS.push_back(evt);
        if( (yPlus[2] + yPlusPlus[2]/4) <= 0.025 && (yPlus[2] + yPlusPlus[2]/4) >= -0.025) EvtsWithYPlusGausSmallAcc.push_back(evt);

        //------------------------------------------------------------//
        //---  Now apply cuts and save the interesting histograms  ---//
        //------------------------------------------------------------//
        for(int iConf = 0; iConf < NrConfigs; iConf++){

          //-- Apply cut on TotalFctDeviation --//
          if( TotalFctDevAccOuter <= 0.1){
            if(iConf == 0) EvtsWithSmallFctDev.push_back(evt);  //Still an interesting variable??
            LLSmallFctDevDist->SetBinContent(LLSmallFctDevDist->FindBin(Var[iConf]), LLSmallFctDevDist->GetBinContent(LLSmallFctDevDist->FindBin(Var[iConf])) + LnLikAcc[iConf]);
          }
          //-- Apply cut on ScdDer (using inner Var points) --//
          if( scdDerInner[0] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerInner++; 
            LLPosScdDerDistInner->SetBinContent(LLPosScdDerDistInner->FindBin(Var[iConf]), LLPosScdDerDistInner->GetBinContent(LLPosScdDerDistInner->FindBin(Var[iConf])) + LnLik[iConf] );          
          }
          if( scdDerInner[1] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerXSInner++;
            LLXSPosScdDerDistInner->SetBinContent(LLXSPosScdDerDistInner->FindBin(Var[iConf]), LLXSPosScdDerDistInner->GetBinContent(LLXSPosScdDerDistInner->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          if( scdDerInner[2] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerAccInner++;
            LLAccPosScdDerDistInner->SetBinContent(LLAccPosScdDerDistInner->FindBin(Var[iConf]), LLAccPosScdDerDistInner->GetBinContent(LLAccPosScdDerDistInner->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }

          //-- Apply cut on both ScdDer --//        
          if( scdDerInner[0] > 0.0 && scdDerOuter[0] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerBoth++;
            LLPosScdDerDistBoth->SetBinContent(LLPosScdDerDistBoth->FindBin(Var[iConf]), LLPosScdDerDistBoth->GetBinContent(LLPosScdDerDistBoth->FindBin(Var[iConf])) + LnLik[iConf] );
          }
          if( scdDerInner[1] > 0.0 && scdDerOuter[1] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerXSBoth++;
            LLXSPosScdDerDistBoth->SetBinContent(LLXSPosScdDerDistBoth->FindBin(Var[iConf]), LLXSPosScdDerDistBoth->GetBinContent(LLXSPosScdDerDistBoth->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          if( scdDerInner[2] > 0.0 && scdDerOuter[2] > 0.0){
            if(iConf == 0) EvtsWithPosScdDerAccBoth++;
            LLAccPosScdDerDistBoth->SetBinContent(LLAccPosScdDerDistBoth->FindBin(Var[iConf]), LLAccPosScdDerDistBoth->GetBinContent(LLAccPosScdDerDistBoth->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }

          //-- Apply cut on ScdDer (using outer Var points) --//
          if( scdDerOuter[0] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerOuter++;
              YPlusGausTestPosScdDer->Fill(yPlus[0] + yPlusPlus[0]/4);
              YPlusPosScdDer->Fill(yPlus[0]); YPlusPlusPosScdDer->Fill(yPlusPlus[0]); YMinPosScdDer->Fill(yMin[0]); YMinMinPosScdDer->Fill(yMinMin[0]);
            }
            LLPosScdDerDistOuter->SetBinContent(LLPosScdDerDistOuter->FindBin(Var[iConf]), LLPosScdDerDistOuter->GetBinContent(LLPosScdDerDistOuter->FindBin(Var[iConf])) + LnLik[iConf] );
          }
          else{  //Are these distributions for negative scdDer still interesting?
            if(iConf == 0){
              YPlusGausTestNegScdDer->Fill(yPlus[0] + yPlusPlus[0]/4);
              YPlusNegScdDer->Fill(yPlus[0]); YPlusPlusNegScdDer->Fill(yPlusPlus[0]); YMinNegScdDer->Fill(yMin[0]); YMinMinNegScdDer->Fill(yMinMin[0]);
            }
          }

          if( scdDerOuter[1] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerXSOuter++;
              YPlusGausTestXSPosScdDer->Fill(yPlus[1] + yPlusPlus[1]/4);
              YPlusXSPosScdDer->Fill(yPlus[1]); YPlusPlusXSPosScdDer->Fill(yPlusPlus[1]); YMinXSPosScdDer->Fill(yMin[1]); YMinMinXSPosScdDer->Fill(yMinMin[1]);
            }
            LLXSPosScdDerDistOuter->SetBinContent(LLXSPosScdDerDistOuter->FindBin(Var[iConf]), LLXSPosScdDerDistOuter->GetBinContent(LLXSPosScdDerDistOuter->FindBin(Var[iConf])) + LnLikXS[iConf] );
          }
          else{
            if(iConf == 0){
              YPlusGausTestXSNegScdDer->Fill(yPlus[1] + yPlusPlus[1]/4);
              YPlusXSNegScdDer->Fill(yPlus[1]); YPlusPlusXSNegScdDer->Fill(yPlusPlus[1]); YMinXSNegScdDer->Fill(yMin[1]); YMinMinXSNegScdDer->Fill(yMinMin[1]);
            }
          }

          if( scdDerOuter[2] > 0.0){
            if(iConf == 0){
              EvtsWithPosScdDerAccOuter++;
              YPlusGausTestAccPosScdDer->Fill(yPlus[2] + yPlusPlus[2]/4);
              YPlusAccPosScdDer->Fill(yPlus[2]); YPlusPlusAccPosScdDer->Fill(yPlusPlus[2]); YMinAccPosScdDer->Fill(yMin[2]); YMinMinAccPosScdDer->Fill(yMinMin[2]);
            }
            LLAccPosScdDerDistOuter->SetBinContent(LLAccPosScdDerDistOuter->FindBin(Var[iConf]), LLAccPosScdDerDistOuter->GetBinContent(LLAccPosScdDerDistOuter->FindBin(Var[iConf])) + LnLikAcc[iConf] );
          }
          else{
            if(iConf == 0){  
              YPlusGausTestAccNegScdDer->Fill(yPlus[2] + yPlusPlus[2]/4);
              YPlusAccNegScdDer->Fill(yPlus[2]); YPlusPlusAccNegScdDer->Fill(yPlusPlus[2]); YMinAccNegScdDer->Fill(yMin[2]); YMinMinAccNegScdDer->Fill(yMinMin[2]);
            }
          }
        }

      }
    }
  }
  ifs.close();

  //--- Save all the histograms containing information about all the events! ---//
  Tfile->cd();
  LnLikAll->Write();                LnLikXSAll->Write();                LnLikAccAll->Write();
  if(storeSplittedCanvas == true){
    TCanvas* LnLikAllCanv    = new TCanvas("LnLikAllCanv",   "LnLikAllCanv");    LnLikAllCanv->cd();    LnLikAll->Draw();    LnLikAllCanv->Print((SplittedDir+"/TotalLnLik.pdf").c_str());
    TCanvas* LnLikAllXSCanv  = new TCanvas("LnLikAllXSCanv", "LnLikAllXSCanv");  LnLikAllXSCanv->cd();  LnLikXSAll->Draw();  LnLikAllXSCanv->Print((SplittedDir+"/TotalLnLikXS.pdf").c_str());
    TCanvas* LnLikAllAccCanv = new TCanvas("LnLikAllAccCanv","LnLikAllAccCanv"); LnLikAllAccCanv->cd(); LnLikAccAll->Draw(); LnLikAllAccCanv->Print((SplittedDir+"/TotalLnLikAcc.pdf").c_str());
  }
  YPlusGausTest->Write();           YPlusGausTestXS->Write();           YPlusGausTestAcc->Write();
  YPlus->Write();                   YPlusXS->Write();                   YPlusAcc->Write();
  YPlusPlus->Write();               YPlusPlusXS->Write();               YPlusPlusAcc->Write();
  YRelPlus->Write();                YRelPlusXS->Write();                YRelPlusAcc->Write();
  YMin->Write();                    YMinXS->Write();                    YMinAcc->Write();
  YMinMin->Write();                 YMinMinXS->Write();                 YMinMinAcc->Write();
  YRelMin->Write();                 YRelMinXS->Write();                 YRelMinAcc->Write();
  ScdDerInner->Write();             ScdDerXSInner->Write();             ScdDerAccInner->Write();
  ScdDerOuter->Write();             ScdDerXSOuter->Write();             ScdDerAccOuter->Write();
  ScdDerScatter->Write();           ScdDerXSScatter->Write();           ScdDerAccScatter->Write();
  //-- Save the histograms for which oveflow information is needed! --//
  for(int iConf = 0; iConf < NrConfigs; iConf++){
    PaintOverflow(h_FitDeviation[iConf], Tfile, "FitDeviation");
    PaintOverflow(h_FitDeviationRel[iConf], Tfile, "RelativeFitDeviation");
  }
  for(int ii = 0; ii < 3; ii++){
    dir_FitDevDelete->cd(); h_PointsDelByFitDev[ii]->Write();
    dir_RelFitDevDelete->cd(); h_PointsDelByFitDevRel[ii]->Write();
    PaintOverflow(h_ChiSquaredFirstFit[ii], Tfile, "ChiSquaredDist");
    PaintOverflow(h_ChiSquaredSecondFit[ii],Tfile, "ChiSquaredDist");
  }

  //---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---//
  std::cout << "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = " << Var[xNeg[0]] << "/" << Var[xMin] << "/" << Var[xPos[0]] << ") : " << EvtsWithPosScdDerInner << ", " << EvtsWithPosScdDerXSInner << " & " << EvtsWithPosScdDerAccInner << std::endl;
  std::cout << "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = " << Var[xNeg[1]] << "/" << Var[xMin] << "/" << Var[xPos[1]] << ") : " << EvtsWithPosScdDerOuter << ", " << EvtsWithPosScdDerXSOuter << " & " << EvtsWithPosScdDerAccOuter << std::endl;
  std::cout << "Nr of events with both 2nd derivatives > 0 (LnLik, LnLikXS & LnLikAcc) : " << EvtsWithPosScdDerOuter << ", " << EvtsWithPosScdDerXSOuter << " & " << EvtsWithPosScdDerAccOuter << std::endl;
  std::cout << "Nr of events with small Gaussian-test deviation for plus-case (LnLik, LnLikXS & LnLikAcc) " << EvtsWithYPlusGausSmall.size() << ", " << EvtsWithYPlusGausSmallXS.size() << " & " << EvtsWithYPlusGausSmallAcc.size() << std::endl;
  std::cout << "Nr of events with total function deviation < 0.5 : " << EvtsWithSmallFctDev.size() << std::endl;

  TDirectory* AppliedCutsDir = Tfile->mkdir("LikelihoodAfterCuts");
  TDirectory* SmallFctDevDir = AppliedCutsDir->mkdir("SmallFunctionDeviation");
  SmallFctDevDir->cd(); LLSmallFctDevDist->Write();

  stringstream ssEvtsWithPosScdDerOuter;    ssEvtsWithPosScdDerOuter << EvtsWithPosScdDerOuter;    string sEvtsWithPosScdDerOuter    = ssEvtsWithPosScdDerOuter.str();
  stringstream ssEvtsWithPosScdDerXSOuter;  ssEvtsWithPosScdDerXSOuter << EvtsWithPosScdDerOuter;  string sEvtsWithPosScdDerXSOuter  = ssEvtsWithPosScdDerXSOuter.str();
  stringstream ssEvtsWithPosScdDerAccOuter; ssEvtsWithPosScdDerAccOuter << EvtsWithPosScdDerOuter; string sEvtsWithPosScdDerAccOuter = ssEvtsWithPosScdDerAccOuter.str();
  stringstream ssEvtsWithPosScdDerInner;    ssEvtsWithPosScdDerInner << EvtsWithPosScdDerInner;    string sEvtsWithPosScdDerInner    = ssEvtsWithPosScdDerInner.str();
  stringstream ssEvtsWithPosScdDerXSInner;  ssEvtsWithPosScdDerXSInner << EvtsWithPosScdDerInner;  string sEvtsWithPosScdDerXSInner  = ssEvtsWithPosScdDerXSInner.str();
  stringstream ssEvtsWithPosScdDerAccInner; ssEvtsWithPosScdDerAccInner << EvtsWithPosScdDerInner; string sEvtsWithPosScdDerAccInner = ssEvtsWithPosScdDerAccInner.str();
  stringstream ssEvtsWithPosScdDerBoth;     ssEvtsWithPosScdDerBoth << EvtsWithPosScdDerBoth;      string sEvtsWithPosScdDerBoth     = ssEvtsWithPosScdDerBoth.str();
  stringstream ssEvtsWithPosScdDerXSBoth;   ssEvtsWithPosScdDerXSBoth << EvtsWithPosScdDerBoth;    string sEvtsWithPosScdDerXSBoth   = ssEvtsWithPosScdDerXSBoth.str();
  stringstream ssEvtsWithPosScdDerAccBoth;  ssEvtsWithPosScdDerAccBoth << EvtsWithPosScdDerBoth;   string sEvtsWithPosScdDerAccBoth  = ssEvtsWithPosScdDerAccBoth.str();
  stringstream ssNEvts; ssNEvts << nEvts; string sNEvts = ssNEvts.str();

  LLPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (no norm -- "+sEvtsWithPosScdDerBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLXSPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (XS norm -- "+sEvtsWithPosScdDerXSBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLAccPosScdDerDistBoth->SetTitle( ("-ln(L) when both 2nd derivatives > 0 (Acc norm -- "+sEvtsWithPosScdDerAccBoth+"/"+sNEvts+" evts -- "+title+")").c_str());
  LLPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (no norm -- "+sEvtsWithPosScdDerInner+"/"+sNEvts+" evts -- "+title+")").c_str()    );
  LLXSPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (XS norm -- "+sEvtsWithPosScdDerXSInner+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLAccPosScdDerDistInner->SetTitle( ("-ln(L) when inner 2nd derivative > 0 (Acc norm -- "+sEvtsWithPosScdDerAccInner+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (no norm -- "+sEvtsWithPosScdDerOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLXSPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (XS norm -- "+sEvtsWithPosScdDerXSOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );
  LLAccPosScdDerDistOuter->SetTitle( ("-ln(L) when outer 2nd derivative > 0 (Acc norm -- "+sEvtsWithPosScdDerAccOuter+"/"+sNEvts+" evts -- "+title+")").c_str() );

  TDirectory* SignScdDerDir = AppliedCutsDir->mkdir("SignSecondDerivative"); SignScdDerDir->cd();
  LLPosScdDerDistInner->Write(); LLXSPosScdDerDistInner->Write(); LLAccPosScdDerDistInner->Write();
  LLPosScdDerDistOuter->Write(); LLXSPosScdDerDistOuter->Write(); LLAccPosScdDerDistOuter->Write();
  LLPosScdDerDistBoth->Write();  LLXSPosScdDerDistBoth->Write();  LLAccPosScdDerDistBoth->Write();

  Tfile->Close();
  file_FitDist->Close();
  /*for(int iEvt = 0; iEvt < nEvts; iEvt++){
    if(std::find(EvtsWithSmallFctDev.begin(), EvtsWithSmallFctDev.end(), iEvt) != EvtsWithSmallFctDev.end() ){
      std::cout << iEvt << ") Event number is found : " << EvtsWithSmallFctDev[iEvt]  << std::endl;
    }
  }*/

}


