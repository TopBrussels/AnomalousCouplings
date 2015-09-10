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
bool storeSplittedCanvas = false; 
std::string SplittedDir = "Events_RecoTest/Reco_UnmatchedEvts_DblGausTF_8Sept/SplittedCanvasses"; 
/////////////////////////////////////////////////////////////

//std::string VarValues[] = {"Re(g_{R}) = -0.3","Re(g_{R}) = -0.2","Re(g_{R}) = -0.1","Re(g_{R}) = -0.05","Re(g_{R}) = 0.0","Re(g_{R}) = 0.05","Re(g_{R}) = 0.1","Re(g_{R}) = 0.2","Re(g_{R}) = 0.3"};
double Var[] = {-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2};
double MGXS[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double MGXSCut[] = {0.948253,1.13713,1.36514,1.63988,1.96892,2.35987,2.8203,3.3578,3.97995};
int xBin = 9; 
float xLow = -0.225; 
float xHigh = 0.225; 
std::string KinVar = "Re(g_{R})"; 
int xMin = 4; 
std::string title = "Reco_UnmatchedEvts_DblGausTF_8Sept"; 
std::string NormTypeName[2] = {"","_Acc"}; 
std::string NormType[2] = {"no","acceptance"}; 
const int nrNorms = sizeof(NormType)/sizeof(NormType[0]);

//ROOT file to store the Fit functions --> Will fasten the study of the cut-influences ...
TFile* file_FitDist = new TFile("Events_RecoTest/Reco_UnmatchedEvts_DblGausTF_8Sept/FitDistributions_Reco_UnmatchedEvts_DblGausTF_8Sept_64268Evts.root","RECREATE"); 
TDirectory *dir_OriginalLL[nrNorms] = {0}, *dir_FirstFit[nrNorms] = {0}, *dir_SecondFit[nrNorms] = {0};

const int NrConfigs = 9; 
const int nEvts = 64268; 
const unsigned int NrToDel = 2; 
int NrRemaining = NrConfigs-NrToDel;
std::string sNrCanvas ="0";
std::string sNrRemaining = ""; std::stringstream ssNrRemaining; 
std::string sNEvts = ""; std::stringstream ssNEvts;
std::string sNrConfigs = ""; std::stringstream ssNrConfigs;

//Information for the stackedCanvas division!
int NrCanvas = 0, xDivide = 4, yDivide = 4;

TF1 *polFit_AllPoints, *polFit_ReducedPoints;
TH1F *h_FitDeviation[NrConfigs], *h_FitDeviationRel[NrConfigs];
TH1F *h_PointsDelByFitDev[3], *h_PointsDelByFitDevRel[3];
TH1F *h_ChiSquaredFirstFit[3], *h_ChiSquaredSecondFit[3];
TH2F* h_TotalFitDevVSChiSq = new TH2F("TotalFitDevVSChiSq","Total fit deviation versus chi-squared",200,0,0.000005, 200, 0, 0.0005);
TH1F* h_TotalRelFitDeviationReduced = new TH1F("TotalRelFitDeviationReduced","TotalRelFitDeviationReduced",200,0,0.001);
TH1F* h_TotalRelFitDeviation = new TH1F("TotalRelFitDeviation","TotalRelFitDeviation",200,0,0.01);

//Store all the fit parameters into a vector of doubles
vector<double> FitParams_FirstFit[nrNorms];
vector<double> FitParams_SecondFit[nrNorms];

//Method to sort a pair based on the second element!
struct sort_pred {
  bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
  }
};

void PaintOverflow(TH1F *h, TFile *FileToWrite, std::string dirName){  // This function draws the histogram h with an extra bin for overflows
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
  FileToWrite->cd();
  if(dirName != ""){
    if(dirName.find("_") < 999){

      std::string firstDirName = dirName.substr(0,dirName.find("_"));
      TDirectory *firstDir = FileToWrite->GetDirectory(firstDirName.c_str());
      if (!firstDir) firstDir = FileToWrite->mkdir(firstDirName.c_str());
      firstDir->cd();

      std::string secondDirName = dirName.substr(dirName.find("_")+1);
      TDirectory *dir = firstDir->GetDirectory(secondDirName.c_str());
      if (!dir) dir = firstDir->mkdir(secondDirName.c_str());
      dir->cd();
    }
    else{
      TDirectory *dir = FileToWrite->GetDirectory(dirName.c_str());
      if (!dir) dir = FileToWrite->mkdir(dirName.c_str());
      dir->cd();
      //delete dir;       -->Cannot delete this directory since it should exist for the following histograms ...
    }
  }

  h_tmp->Write();  
  FileToWrite->cd();  //Reset to general directory! 

  delete h_tmp;
}                  

void calculateFit(TH1F *h_LogLik, string EvtNumber, int normType, int evtCounter, TCanvas *canv_SplittedLL){
  file_FitDist->cd();

  double LogLikelihood[NrConfigs];
  for(int ii = 0; ii < NrConfigs; ii++)
    LogLikelihood[ii] = h_LogLik->GetBinContent(h_LogLik->FindBin(Var[ii]));

  dir_OriginalLL[normType]->cd();
  h_LogLik->Write();
  std::string YAxisTitle = "-ln(L) value ("+NormType[normType]+" norm -- evt "+EvtNumber+")";

  //Set name of chisquared distributions!
  if(evtCounter == 1){
    h_ChiSquaredFirstFit[normType]  = new TH1F(("ChiSquared"+NormTypeName[normType]+"_FirstFit").c_str(), ("Distribution of the #chi^{2} after the fit on all the points ("+NormType[normType]+" norm)").c_str(),200,0,0.05);
    h_ChiSquaredSecondFit[normType] = new TH1F(("ChiSquared"+NormTypeName[normType]+"_SecondFit").c_str(),("Distribution of the #chi^{2} after the fit on the reduced points ("+NormType[normType]+" norm)").c_str(),200,0,0.005);
    h_PointsDelByFitDev[normType]    = new TH1F(("DeletedPointsFitDev"+NormTypeName[normType]).c_str(),   ("Overview of deleted points due to largest FitDeviation ("+NormType[normType]+" norm)").c_str(),xBin,xLow,xHigh);
    h_PointsDelByFitDevRel[normType] = new TH1F(("DeletedPointsFitDevRel"+NormTypeName[normType]).c_str(),("Overview of deleted points due to largest relative FitDeviation ("+NormType[normType]+" norm)").c_str(),xBin,xLow,xHigh);
  }
 
  polFit_AllPoints = new TF1(("polFit"+NormTypeName[normType]+"_AllPoints_Evt"+EvtNumber).c_str(),"pol2",Var[1],Var[NrConfigs-2]); 
  TGraph* gr_LnLik = new TGraph(NrConfigs,Var, LogLikelihood);
  gr_LnLik->Fit(polFit_AllPoints,"Q","",polFit_AllPoints->GetXmin(), polFit_AllPoints->GetXmax());
  h_ChiSquaredFirstFit[normType]->Fill(polFit_AllPoints->GetChisquare());

  if(evtCounter == 1){
    FitParams_FirstFit[normType].clear();
    for(int ipar = 0; ipar < polFit_AllPoints->GetNpar(); ipar++)
      FitParams_FirstFit[normType].push_back(polFit_AllPoints->GetParameter(ipar));
  }
  else{
    for(int ipar = 0; ipar < polFit_AllPoints->GetNpar(); ipar++)
      FitParams_FirstFit[normType][ipar] += polFit_AllPoints->GetParameter(ipar);
  }

  dir_FirstFit[normType]->cd();
  polFit_AllPoints->Write();

  //Calculate the lowest value for the LnLik!
  float LnLikMin = LogLikelihood[0];
  for(int iConf = 0; iConf < NrConfigs; iConf++){
    if(LnLikMin > LogLikelihood[iConf]){
      LnLikMin = LogLikelihood[iConf];	
    }
  }

  std::vector<std::pair<int, double> > FitDeviation, FitDeviationRel;
  double LogLikFit[NrConfigs] = {-9999};
  double TotalRelFitDeviation = 0;
  for(int iConfig = 0; iConfig < NrConfigs; iConfig++){
    if(polFit_AllPoints->GetNpar() == 3)
      LogLikFit[iConfig] = polFit_AllPoints->GetParameter(0)+polFit_AllPoints->GetParameter(1)*Var[iConfig]+polFit_AllPoints->GetParameter(2)*Var[iConfig]*Var[iConfig];
    else if(polFit_AllPoints->GetNpar() == 5)
      LogLikFit[iConfig] = polFit_AllPoints->GetParameter(0)+polFit_AllPoints->GetParameter(1)*Var[iConfig]+polFit_AllPoints->GetParameter(2)*Var[iConfig]*Var[iConfig]+polFit_AllPoints->GetParameter(3)*pow(Var[iConfig],3)+polFit_AllPoints->GetParameter(4)*pow(Var[iConfig],4)+polFit_AllPoints->GetParameter(5)*pow(Var[iConfig],5);
    FitDeviation.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig]) ) );
    FitDeviationRel.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig])/LogLikelihood[iConfig] ) );
    TotalRelFitDeviation += abs(LogLikelihood[iConfig]-LogLikFit[iConfig])/LogLikelihood[iConfig];
  }
  //Sort the fitdeviation values depending on the second value!
  std::sort(FitDeviation.begin(), FitDeviation.end(), sort_pred() );
  std::sort(FitDeviationRel.begin(), FitDeviationRel.end(), sort_pred() );
  if(NormTypeName[normType] == "_Acc") h_TotalRelFitDeviation->Fill(TotalRelFitDeviation);

  //Now loop again over all configurations, plot the FitDeviation in sorted order and save the configNr's which should be excluded 
  std::vector<int> FitDevPointsToDel, FitDevRelPointsToDel;
  for(int itSortedConfig = NrConfigs-1; itSortedConfig >= 0 ; itSortedConfig--){  //Looping from high to low values of the deviation!
    if(NormTypeName[normType] == "_Acc") h_FitDeviation[itSortedConfig]->Fill(FitDeviation[itSortedConfig].second);
    if(NormTypeName[normType] == "_Acc") h_FitDeviationRel[itSortedConfig]->Fill(FitDeviationRel[itSortedConfig].second);
    
    //Store the 'NrToDel' points which need to be excluded from the TGraph!
    if(FitDevPointsToDel.size() < NrToDel){ FitDevPointsToDel.push_back(FitDeviation[itSortedConfig].first); h_PointsDelByFitDev[normType]->Fill(Var[FitDeviation[itSortedConfig].first]);}
    if(FitDevRelPointsToDel.size() < NrToDel){FitDevRelPointsToDel.push_back(FitDeviationRel[itSortedConfig].first);h_PointsDelByFitDevRel[normType]->Fill(Var[FitDeviationRel[itSortedConfig].first]);}
  }

  //Create new arrays with the reduced information!
  double ReducedVar[NrConfigs-NrToDel], ReducedLogLik[NrConfigs-NrToDel];
  double TotalRelFitDeviationReduced = 0;
  int iCounter = 0;
  for(int iConf = 0; iConf < NrConfigs; iConf++){
    if(!(std::find(FitDevPointsToDel.begin(), FitDevPointsToDel.end(), iConf) != FitDevPointsToDel.end()) ){
      ReducedVar[iCounter] = Var[iConf];
      ReducedLogLik[iCounter] = LogLikelihood[iConf];
      TotalRelFitDeviationReduced += (abs(LogLikelihood[iConf]-LogLikFit[iConf])/LogLikelihood[iConf]);
      iCounter++;
    }
  }
  h_TotalRelFitDeviationReduced->Fill(TotalRelFitDeviationReduced);

  //Define new TGraph and fit again
  TGraph* gr_ReducedLnLik = new TGraph(NrConfigs-NrToDel, ReducedVar, ReducedLogLik);
  polFit_ReducedPoints = new TF1(("polFit"+NormTypeName[normType]+"_"+sNrRemaining+"ReducedPoints_Evt"+EvtNumber).c_str(),"pol2",Var[1],Var[NrConfigs-2]); 
  gr_ReducedLnLik->Fit(polFit_ReducedPoints,"Q","",polFit_AllPoints->GetXmin(), polFit_AllPoints->GetXmax()); 
  h_ChiSquaredSecondFit[normType]->Fill(polFit_ReducedPoints->GetChisquare());   //As expected NDF is always equal to NrConfigs-NrToDel-3 (= nr params needed to define a parabola)

  //Add the different parameters of the individual fits in order to form a final fit!
  if(evtCounter == 1){
    FitParams_SecondFit[normType].clear();
    for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
      FitParams_SecondFit[normType].push_back(polFit_ReducedPoints->GetParameter(ipar));
  }
  else{
    for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
      FitParams_SecondFit[normType][ipar] += polFit_ReducedPoints->GetParameter(ipar);
  }

  dir_SecondFit[normType]->cd();
  //Create a 2D-plot which contains the deviation of the expected minimum wrt the chi-squared of the fit
  if(NormTypeName[normType] == "_Acc"){
    h_TotalFitDevVSChiSq->Fill(polFit_ReducedPoints->GetChisquare(), TotalRelFitDeviationReduced);
  }

  polFit_ReducedPoints->Write();
  stringstream ssChiSq; ssChiSq << polFit_AllPoints->GetChisquare(); string sChiSq = ssChiSq.str();
  std::ostringstream osParam0; osParam0 << polFit_AllPoints->GetParameter(0); std::string sParam0 = osParam0.str();
  std::ostringstream osParam1; osParam1 << polFit_AllPoints->GetParameter(1); std::string sParam1 = osParam1.str();
  std::ostringstream osParam2; osParam2 << polFit_AllPoints->GetParameter(2); std::string sParam2 = osParam2.str();
  std::ostringstream osParam3; osParam3 << polFit_AllPoints->GetParameter(3); std::string sParam3 = osParam3.str();
  std::ostringstream osParam4; osParam4 << polFit_AllPoints->GetParameter(4); std::string sParam4 = osParam4.str();
  h_LogLik->SetTitle(("FitInfo (slope req): "+sParam0+", "+sParam1+", "+sParam2+", "+sParam3+" & "+sParam4).c_str());

  if( storeSplittedCanvas == true){
    h_LogLik->GetYaxis()->SetTitle(YAxisTitle.c_str());
    h_LogLik->GetYaxis()->SetTitleOffset(1.4);
    h_LogLik->GetXaxis()->SetTitle(KinVar.c_str());
    h_LogLik->SetMaximum(polFit_ReducedPoints->GetMaximum(Var[0], Var[NrConfigs-1]));
    h_LogLik->SetMinimum(polFit_ReducedPoints->GetMinimum(Var[0], Var[NrConfigs-1]));
    polFit_ReducedPoints->SetLineColor(7);
    canv_SplittedLL->cd(evtCounter - (xDivide*yDivide*NrCanvas) ); h_LogLik->Draw("p"); polFit_AllPoints->Draw("same"); polFit_ReducedPoints->Draw("same"); canv_SplittedLL->Update();
    //canv_SplittedLL->cd(evtCounter - (xDivide*yDivide*NrCanvas) ); h_LogLik->Draw("p"); polFit_ReducedPoints->Draw("same"); canv_SplittedLL->Update();
  }
  else{
    delete polFit_ReducedPoints;
    delete polFit_AllPoints;
  }

  delete gr_ReducedLnLik;
  delete gr_LnLik;
}

void doublePolFitMacro(){

  clock_t start = clock();
  
  ssNrRemaining << NrRemaining; sNrRemaining = ssNrRemaining.str();
  ssNEvts << nEvts; sNEvts = ssNEvts.str();
  ssNrConfigs << NrConfigs; sNrConfigs = ssNrConfigs.str();

  TH1F *histSum[nrNorms];
  for(int iNorm = 0; iNorm < nrNorms; iNorm++){
    histSum[iNorm] = new TH1F(("SummedHist"+NormTypeName[iNorm]).c_str(),("Sum of all individual histograms ("+NormType[iNorm]+" norm)").c_str(),xBin,xLow,xHigh); 
  }

  TH1F *h_LnLik[nrNorms] = {0};
  TDirectory *dir_SplitCanv = 0, *dir_LLSplit[nrNorms] = {0};
  if(storeSplittedCanvas)
    dir_SplitCanv = file_FitDist->mkdir("SplitCanvasses");

  //Name the directories according to the number of normalisations that should be considered! 
  for(int iNorm = 0; iNorm < nrNorms; iNorm++){
    dir_OriginalLL[iNorm] = file_FitDist->mkdir(("OriginalLL"+NormTypeName[iNorm]).c_str());
    dir_FirstFit[iNorm] = file_FitDist->mkdir(("FirstPolynomialFit"+NormTypeName[iNorm]).c_str());
    dir_SecondFit[iNorm] = file_FitDist->mkdir(("SecondPolynomialFit"+NormTypeName[iNorm]).c_str());
  
    if(storeSplittedCanvas)
      dir_LLSplit[iNorm] = dir_SplitCanv->mkdir(("LnLik"+NormTypeName[iNorm]).c_str());
  }
 
  int consEvts = 0;
  TCanvas *canv_SplitLL[nrNorms] = {0};
  double LnLik[nrNorms][NrConfigs] = {{0.0}}; //, LnLikXS[NrConfigs] = {0.0}, LnLikAcc[NrConfigs] = {0.0};        

  //--- Read all likelihood values ! ---//
  std::ifstream ifs ("Events_RecoTest/Reco_UnmatchedEvts_DblGausTF_8Sept/weights.out", std::ifstream::in); 
  std::cout << " Value of ifs : " << ifs.eof() << std::endl;
  std::string line;
  int evt,config,tf;
  double weight, weightUnc;
  double CosThetaCorr = 1;
  while( std::getline(ifs,line) && consEvts < nEvts){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc){ 
      if(config == 1 && ((consEvts+1) % 2000 == 0) ) std::cout << " Looking at event : " << consEvts+1 << " (" << (double)(consEvts+1)*100/(double)nEvts << "%)" << flush<<"\r";
      stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();
      stringstream ssConfig; ssConfig << config; string sConfig = ssConfig.str();

      //Initialize the fitDeviation histograms (array of TH1F, one for each configuration ...)
      if(consEvts == 0){
        h_FitDeviation[config-1] = new TH1F(("FitDeviation_LowestConfig"+sConfig).c_str(),("FitDeviation histogram for "+sConfig+"st lowest value").c_str(), 200, 0, 0.1);
        h_FitDeviationRel[config-1] = new TH1F(("FitDeviationRelative_LowestConfig"+sConfig).c_str(), ("FitDeviationRelative histogram for "+sConfig+"st lowest value").c_str(),200,0,0.005);
      }

      //--- Initialize the event-per-event variables! ---//
      // --    Loop over the different allowed normalisations
      for(int iNorm = 0; iNorm < nrNorms; iNorm++){
  
        if( config == 1){
          //std::cout << " \n Looking at event : " << sEvt << " with cos theta* weight = " << CosThetaCorr  << std::endl;
          h_LnLik[iNorm] = new TH1F(("LnLik"+NormTypeName[iNorm]+"_Evt"+sEvt).c_str(),("LnLik"+NormTypeName[iNorm]+" distribution for event "+sEvt+" -- "+title+" evts").c_str(),xBin,xLow,xHigh);
          h_LnLik[iNorm]->SetMarkerStyle(20+iNorm); h_LnLik[iNorm]->SetLineColor(3+iNorm); h_LnLik[iNorm]->SetMarkerColor(1); h_LnLik[iNorm]->SetMarkerSize(1.2);
        }
        //std::cout << " Changing LnLik[2] from " << (-log(weight)+log(MGXSCut[config-1])) << " to " << (-log(weight)+log(MGXSCut[config-1]))*CosThetaCorr << std::endl;	
        if(iNorm == 0)                                                        LnLik[iNorm][config-1] = -log(weight)*CosThetaCorr;
        else if(iNorm == 1 && nrNorms == 3)                                   LnLik[iNorm][config-1] = (-log(weight)+log(MGXS[config-1]))*CosThetaCorr;
        else if((nrNorms == 3 && iNorm == 2) || (nrNorms == 2 && iNorm == 1)) LnLik[iNorm][config-1] = (-log(weight)+log(MGXSCut[config-1]))*CosThetaCorr;
        //std::cout << " Checking whether cos theta* correction changes likelihood : " << (-log(weight)+log(MGXSCut[config-1])) << " --> " << LnLik[2][config-1] << " ==> A difference of " << LnLik[2][config-1] - (-log(weight)+log(MGXSCut[config-1])) << std::endl;

        //---  Fill the LnLik histograms for each event and for all events together  ---//
        h_LnLik[iNorm]->SetBinContent(h_LnLik[iNorm]->FindBin(Var[config-1]), LnLik[iNorm][config-1]);

        //---  Only perform the fit after all configurations are considered!  ---//
        if( config == NrConfigs){
          histSum[iNorm]->Add( h_LnLik[iNorm] );
          if(iNorm == 0) consEvts++;   //Count the number of full events!

          //Save xDivide*yDivide of these histograms in one TCanvas!
          if( storeSplittedCanvas == true){
            if(consEvts == 1){
              canv_SplitLL[iNorm] = new TCanvas(("SplitCanvasLL"+NormTypeName[iNorm]+"_Nr0").c_str(),("SplitCanvasLL"+NormTypeName[iNorm]).c_str()); canv_SplitLL[iNorm]->Divide(xDivide,yDivide);    
            }
          }

          //-- Send the array containing the ln(weights) to the predefined function to define the TGraph, fit this, detect the deviation points and fit again! --//
	  calculateFit(h_LnLik[iNorm], sEvt, iNorm, consEvts, canv_SplitLL[iNorm]);

          if( storeSplittedCanvas == true){
            if( consEvts == (xDivide*yDivide*(NrCanvas+1)) || consEvts == nEvts){
              canv_SplitLL[iNorm]->Print((SplittedDir+"/SplitCanvasLL"+NormTypeName[iNorm]+"_Nr"+sNrCanvas+".pdf").c_str()); dir_LLSplit[iNorm]->cd(); canv_SplitLL[iNorm]->Write(); delete canv_SplitLL[iNorm];
              delete h_LnLik[iNorm];
              if( consEvts != nEvts){
                if(iNorm == nrNorms-1) NrCanvas++; 
                stringstream ssNrCanvas; ssNrCanvas << NrCanvas; sNrCanvas = ssNrCanvas.str();
                canv_SplitLL[iNorm] = new TCanvas(("SplitCanvasLL"+NormTypeName[iNorm]+"_Nr"+sNrCanvas).c_str(), ("SplitCanvasLL"+NormTypeName[iNorm]).c_str()); canv_SplitLL[iNorm]->Divide(xDivide,yDivide);   
              }
            }
          }//Draw the stacked canvasses!
          else{
            delete h_LnLik[iNorm];
          }
        }
      }//End of loop over nrNorms!
    }
  }
  ifs.close();


  //-- Save the histograms for which oveflow information is needed! --//
  for(int iConf = 0; iConf < NrConfigs; iConf++){
    PaintOverflow(h_FitDeviation[iConf],    file_FitDist, "FitResults_FitDeviation");         delete h_FitDeviation[iConf];
    PaintOverflow(h_FitDeviationRel[iConf], file_FitDist, "FitResults_RelativeFitDeviation"); delete h_FitDeviationRel[iConf];
  }
  TDirectory *dir_FitResults = file_FitDist->GetDirectory("FitResults");
  if (!dir_FitResults) dir_FitResults = file_FitDist->mkdir("FitResults");
  dir_FitResults->cd();

  h_TotalFitDevVSChiSq->Write();
  h_TotalRelFitDeviationReduced->Write();
  h_TotalRelFitDeviation->Write();

  TDirectory *dir_FitDevDelete = dir_FitResults->GetDirectory("PointsDeletedByFitDev");
  if(!dir_FitDevDelete) dir_FitDevDelete = dir_FitResults->mkdir("PointsDeletedByFitDev");

  TDirectory *dir_RelFitDevDelete = dir_FitResults->GetDirectory("PointsDeletedByRelFitDev");
  if(!dir_RelFitDevDelete) dir_RelFitDevDelete = dir_FitResults->mkdir("PointsDeletedByRelFitDev");

  for(int ii = 0; ii < nrNorms; ii++){    
    dir_FitDevDelete->cd(); h_PointsDelByFitDev[ii]->Write();               delete h_PointsDelByFitDev[ii];
    dir_RelFitDevDelete->cd(); h_PointsDelByFitDevRel[ii]->Write();         delete h_PointsDelByFitDevRel[ii];
    PaintOverflow(h_ChiSquaredFirstFit[ii], file_FitDist, "FitResults_ChiSquaredFit"); delete h_ChiSquaredFirstFit[ii];
    PaintOverflow(h_ChiSquaredSecondFit[ii],file_FitDist, "FitResults_ChiSquaredFit"); delete h_ChiSquaredSecondFit[ii];
  }
  
  TDirectory *dir_FitSums = file_FitDist->GetDirectory("FitSums");
  if (!dir_FitSums) dir_FitSums = file_FitDist->mkdir("FitSums");
  dir_FitSums->cd();
  
  //Now create the final fit from the individually summed fit parameters!
  std::cout << " " << endl;
  for(int iNorm = 0; iNorm < nrNorms; iNorm++){
    histSum[iNorm]->Write();

    TF1* FitSum_FirstFit = new TF1(("SummedFit_FirstFit"+NormTypeName[iNorm]).c_str(),"pol2",xLow,xHigh);
    TF1* FitSum_SecondFit = new TF1(("SummedFit_SecondFit"+NormTypeName[iNorm]).c_str(),"pol2",xLow,xHigh);
    FitSum_FirstFit->SetTitle(("Distribution of first fit after summing over "+sNrConfigs+" points ("+sNEvts+" events -- "+NormType[iNorm]+" norm)").c_str());
    FitSum_SecondFit->SetTitle(("Distribution of second fit after summing over "+sNrConfigs+" points ("+sNEvts+" events -- "+NormType[iNorm]+" norm)").c_str());

    for(int ipar = 0; ipar < FitSum_FirstFit->GetNpar(); ipar++) FitSum_FirstFit->SetParameter(ipar, FitParams_FirstFit[iNorm][ipar]);
    FitSum_FirstFit->Write();
    std::cout << " MinimumX value for " << FitSum_FirstFit->GetName() << " is : " << FitSum_FirstFit->GetMinimumX() << " +- " << FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, FitSum_FirstFit->GetMinimumX(), 0.2) - FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, -0.2, FitSum_FirstFit->GetMinimumX()) << endl;
    delete FitSum_FirstFit;

    for(int ipar = 0; ipar < FitSum_SecondFit->GetNpar(); ipar++) FitSum_SecondFit->SetParameter(ipar, FitParams_SecondFit[iNorm][ipar]);
    FitSum_SecondFit->Write();
    std::cout << " MinimumX value for " << FitSum_SecondFit->GetName() << " is : " << FitSum_SecondFit->GetMinimumX() << " +- " << FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, FitSum_SecondFit->GetMinimumX(), 0.2) - FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, -0.2, FitSum_SecondFit->GetMinimumX()) << endl;
    delete FitSum_SecondFit;
  }

  file_FitDist->Close();

  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
}
