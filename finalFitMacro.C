//------------------------------------------------------------------------------------------------------------//
//  Scripts which compares the different root files and directories and directly saves the histogram!         //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` finalFitMacr.C -o finalFitMacro   //
//               ./finalFitMacro                                                                              //
//------------------------------------------------------------------------------------------------------------//

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

using namespace std;

int main(int argc, char *argv[]){

  //Input variables will be:
  // 1) Likelihood cut value!
  // 2-..) Weight files which will be considered

  double LikCut = 100;
  if( argc >= 2)
    LikCut = atoi(argv[1]);
  std::cout << " LikCut value is : " << LikCut << std::endl;

  vector<string> inputFiles;
  if( argc >= 3){
    for(int iFile = 2; iFile < argc; iFile++){
      inputFiles.push_back(string(argv[iFile]).c_str());
      std::cout << " Stored file name is : " << inputFiles[inputFiles.size()-1] << std::endl;
    }
  }

  //Considered values and corresponding XS-values
  const int NrConfigs = 9; 
  double Var[NrConfigs]     = {-0.2,     -0.15,   -0.1,    -0.05,   0.0,     0.05,    0.1,     0.15,    0.2    };
  double MGXSCut[NrConfigs] = {0.947244, 1.13624, 1.36448, 1.63952, 1.96892, 2.36027, 2.82111, 3.35903, 3.98157};

  return 0;
}

/*

  //Information for the histograms
  int NrBins = 8; 
  float xLow, xHigh;
  double FitMin = -0.15, FitMax = 0.15; 
  const unsigned int NrToDel = 2; 
  int NrRemaining = NrConfigs-NrToDel;


//ROOT file to store the Fit functions --> Will fasten the study of the cut-influences ...
TFile* file_FitDist = new TFile("Events_RecoTest/Reco_CorrectEvts_DblGausTF_LeptDelta_NonBinned_AllEvts_ISR1_28Oct/FitDistributions_Reco_CorrectEvts_DblGausTF_LeptDelta_NonBinned_AllEvts_ISR1_28Oct_SFAdded_LikelihoodCut66_OuterBinsExclForFit_117658Evts.root","RECREATE"); 
TDirectory *dir_OriginalLL = 0, *dir_FirstFit = 0, *dir_SecondFit = 0;

std::string sNrRemaining = ""; std::stringstream ssNrRemaining; 
std::string sNEvts = ""; std::stringstream ssNEvts;
std::string sNrConfigs = ""; std::stringstream ssNrConfigs;

TF1 *polFit_AllPoints, *polFit_ReducedPoints;
TH1F *h_FitDeviation[NrConfigs], *h_FitDeviationRel[NrConfigs];
TH1F *h_PointsDelByFitDev[3], *h_PointsDelByFitDevRel[3];
TH1F *h_ChiSquaredFirstFit[3], *h_ChiSquaredSecondFit[3];
TH2F* h_TotalFitDevVSChiSq = new TH2F("TotalFitDevVSChiSq","Total fit deviation versus chi-squared",200,0,0.000005, 200, 0, 0.0005);
TH1F* h_TotalRelFitDeviationReduced = new TH1F("TotalRelFitDeviationReduced","TotalRelFitDeviationReduced",200,0,0.001);
TH1F* h_TotalRelFitDeviation = new TH1F("TotalRelFitDeviation","TotalRelFitDeviation",200,0,0.01);
TH1F* h_SMLikelihoodValue = new TH1F("SMLikelihoodValue","Distribution of likelihood value at gR = 0.0",500,30,90);
TH1F* h_SMLikValue_AfterCut = new TH1F("SMLikValue_AfterCut","Distribution of likelihood value at gR = 0.0 (after cut)",500,30,90);
TH2F* h_SMLikelihoodValue_vs_DeltaLikelihood = new TH2F("SMLikelihoodValue_vs_DeltaLikelihood","Likelihood value at gR = 0 versus difference in likelihood",500,30,90,100,0,5);

//Store all the fit parameters into a vector of doubles
vector<double> FitParams_FirstFit;
vector<double> FitParams_SecondFit;

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

void calculateFit(TH1F *h_LogLik, string EvtNumber, int evtCounter){
  file_FitDist->cd();

  double LogLikelihood[NrConfigs];
  for(int ii = 0; ii < NrConfigs; ii++)
    LogLikelihood[ii] = h_LogLik->GetBinContent(h_LogLik->FindBin(Var[ii]));

  dir_OriginalLL->cd();
  h_LogLik->Write();
  std::string YAxisTitle = "-ln(L) value (evt "+EvtNumber+")";

  //Set name of chisquared distributions!
  if(evtCounter == 1){
    h_ChiSquaredFirstFit  = new TH1F("ChiSquared_FirstFit", "Distribution of the #chi^{2} after the fit on all the points",200,0,0.05);
    h_ChiSquaredSecondFit = new TH1F("ChiSquared_SecondFit","Distribution of the #chi^{2} after the fit on the reduced points",200,0,0.005);
    h_PointsDelByFitDev    = new TH1F("DeletedPointsFitDev",   "Overview of deleted points due to largest FitDeviation",NrBins+1,xLow,xHigh);
    h_PointsDelByFitDevRel = new TH1F("DeletedPointsFitDevRel","Overview of deleted points due to largest relative FitDeviation",NrBins+1,xLow,xHigh);
  }
 
  polFit_AllPoints = new TF1(("polFit_AllPoints_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
  TGraph* gr_LnLik = new TGraph(NrConfigs,Var, LogLikelihood);
  gr_LnLik->Fit(polFit_AllPoints,"Q","",polFit_AllPoints->GetXmin(), polFit_AllPoints->GetXmax());
  h_ChiSquaredFirstFit->Fill(polFit_AllPoints->GetChisquare());

  if(evtCounter == 1){
    FitParams_FirstFit.clear();
    for(int ipar = 0; ipar < polFit_AllPoints->GetNpar(); ipar++)
      FitParams_FirstFit.push_back(polFit_AllPoints->GetParameter(ipar));
  }
  else{
    for(int ipar = 0; ipar < polFit_AllPoints->GetNpar(); ipar++)
      FitParams_FirstFit[ipar] += polFit_AllPoints->GetParameter(ipar);
  }

  dir_FirstFit->cd();
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
  if(NormTypeName == "_Acc") h_TotalRelFitDeviation->Fill(TotalRelFitDeviation);

  //Now loop again over all configurations, plot the FitDeviation in sorted order and save the configNr's which should be excluded 
  std::vector<int> FitDevPointsToDel, FitDevRelPointsToDel;
  for(int itSortedConfig = NrConfigs-1; itSortedConfig >= 0 ; itSortedConfig--){  //Looping from high to low values of the deviation!
    h_FitDeviation[itSortedConfig]->Fill(FitDeviation[itSortedConfig].second);
    h_FitDeviationRel[itSortedConfig]->Fill(FitDeviationRel[itSortedConfig].second);
    
    //Store the 'NrToDel' points which need to be excluded from the TGraph!
    if(FitDevPointsToDel.size() < NrToDel){ FitDevPointsToDel.push_back(FitDeviation[itSortedConfig].first); h_PointsDelByFitDev->Fill(Var[FitDeviation[itSortedConfig].first]);}
    if(FitDevRelPointsToDel.size() < NrToDel){FitDevRelPointsToDel.push_back(FitDeviationRel[itSortedConfig].first);h_PointsDelByFitDevRel->Fill(Var[FitDeviationRel[itSortedConfig].first]);}
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
  polFit_ReducedPoints = new TF1(("polFit_"+sNrRemaining+"ReducedPoints_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
  gr_ReducedLnLik->Fit(polFit_ReducedPoints,"Q","",polFit_AllPoints->GetXmin(), polFit_AllPoints->GetXmax()); 
  h_ChiSquaredSecondFit->Fill(polFit_ReducedPoints->GetChisquare());   //As expected NDF is always equal to NrConfigs-NrToDel-3 (= nr params needed to define a parabola)

  //Add the different parameters of the individual fits in order to form a final fit!
  if(evtCounter == 1){
    FitParams_SecondFit.clear();
    for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
      FitParams_SecondFit.push_back(polFit_ReducedPoints->GetParameter(ipar));
  }
  else{
    for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
      FitParams_SecondFit[ipar] += polFit_ReducedPoints->GetParameter(ipar);
  }

  dir_SecondFit->cd();
  //Create a 2D-plot which contains the deviation of the expected minimum wrt the chi-squared of the fit
  h_TotalFitDevVSChiSq->Fill(polFit_ReducedPoints->GetChisquare(), TotalRelFitDeviationReduced);

  polFit_ReducedPoints->Write();
  stringstream ssChiSqAll; ssChiSqAll << polFit_AllPoints->GetChisquare(); string sChiSqAll = ssChiSqAll.str();
  stringstream ssChiSqRed; ssChiSqRed << polFit_ReducedPoints->GetChisquare(); string sChiSqRed = ssChiSqRed.str();
  std::ostringstream osParam0; osParam0 << polFit_AllPoints->GetParameter(0); std::string sParam0 = osParam0.str();
  std::ostringstream osParam1; osParam1 << polFit_AllPoints->GetParameter(1); std::string sParam1 = osParam1.str();
  std::ostringstream osParam2; osParam2 << polFit_AllPoints->GetParameter(2); std::string sParam2 = osParam2.str();
  h_LogLik->SetTitle(("FitInfo: #chi^{2} = "+sChiSqAll+" & "+sChiSqRed).c_str());

  delete polFit_ReducedPoints;
  delete polFit_AllPoints;

  delete gr_ReducedLnLik;
  delete gr_LnLik;
}

void doublePolFitMacro(){

  clock_t start = clock();
  
  double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
  xLow = Var[0] - halfBinWidth; xHigh = Var[NrConfigs-1] + halfBinWidth; 
  
  ssNrRemaining << NrRemaining; sNrRemaining = ssNrRemaining.str();
  ssNEvts << nEvts; sNEvts = ssNEvts.str();
  ssNrConfigs << NrConfigs; sNrConfigs = ssNrConfigs.str();

  TH1F *histSum = new TH1F("SummedHist","Sum of all individual histograms",NrBins+1,xLow,xHigh); 
  TH1F *h_LnLik = 0;

  int consEvts = 0;
  TCanvas *canv_SplitLL = 0;
  double LnLik[NrConfigs] = {0.0}; //, LnLikXS[NrConfigs] = {0.0}, LnLikAcc[NrConfigs] = {0.0};        

  //--- Read all likelihood values ! ---//
  std::ifstream ifs ("Events_RecoTest/Reco_CorrectEvts_DblGausTF_LeptDelta_NonBinned_AllEvts_ISR1_28Oct/weights_SFAdded.out", std::ifstream::in); 
  std::cout << " Value of ifs : " << ifs.eof() << std::endl;
  std::string line;
  int evt,config,tf;
  double weight, weightUnc;
  double CosThetaCorr = 1;
  double MCScaleFactor = 1, Luminosity = 1, NormFactor = 1;
  int SMConfig = 99; 
  while( std::getline(ifs,line) && consEvts < nEvts){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc >> MCScaleFactor >> Luminosity >> NormFactor ){ 
      if(config == 1 && ((consEvts+1) % 2000 == 0) ) std::cout << " Looking at event : " << consEvts+1 << " (" << (double)(consEvts+1)*100/(double)nEvts << "%)" << flush<<"\r";
      stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

      //--- Initialize the event-per-event variables! ---//          
      double MaxLik = 0, MinLik = 9999;
      if( config == 1){
        h_LnLik = new TH1F(("LnLik_Evt"+sEvt).c_str(),("LnLik distribution for event "+sEvt+" -- "+title).c_str(),NrBins+1,xLow,xHigh);
        h_LnLik->SetMarkerStyle(20); h_LnLik->SetLineColor(3); h_LnLik->SetMarkerColor(1); h_LnLik->SetMarkerSize(1.2);
      }
      LnLik[config-1] = (-log(weight)+log(MGXSCut[config-1]));

      //Get the maximum and minimum likelihood value
      if(LnLik[config-1] > MaxLik) MaxLik = LnLik[config-1];
      if(LnLik[config-1] < MinLik) MinLik = LnLik[config-1];

      //Set the SMConfig:
      if( Var[config-1] == 0.0) SMConfig = config-1;

      //---  Fill the LnLik histograms for each event and for all events together  ---//
      h_LnLik->SetBinContent(h_LnLik->FindBin(Var[config-1]), LnLik[config-1]*MCScaleFactor*Luminosity*NormFactor);

      //---  Only perform the fit after all configurations are considered!  ---//
      if( config == NrConfigs){

        //Plot the SM likelihood value:
        h_SMLikelihoodValue->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);

        //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
        //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
        if( LnLik[SMConfig] > LikCut){ if(LikCut == 100) std::cout << " ******** \n ERROR: Should not reject any event when cut-value = 100 !!! \n ******** \n " << std::endl; delete h_LnLik; continue; }
        consEvts++;   //Count the number of full events!

        //Initialize the fitDeviation histograms (array of TH1F, one for each configuration ...)
        if(consEvts == 1){
          for(int ii = 0; ii < NrConfigs; ii++){
            stringstream ssConfig; ssConfig << ii+1; string sConfig = ssConfig.str();
            h_FitDeviation[ii] = new TH1F(("FitDeviation_LowestConfig"+sConfig).c_str(),("FitDeviation histogram for "+sConfig+"st lowest value").c_str(), 200, 0, 0.1);
            h_FitDeviationRel[ii] = new TH1F(("FitDeviationRelative_LowestConfig"+sConfig).c_str(), ("FitDeviationRelative histogram for "+sConfig+"st lowest value").c_str(),200,0,0.005);
          }
        }

        h_SMLikValue_AfterCut->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
        histSum->Add( h_LnLik);
        h_SMLikelihoodValue_vs_DeltaLikelihood->Fill(LnLik[SMConfig], MaxLik-MinLik);

        //-- Send the array containing the ln(weights) to the predefined function to define the TGraph, fit this, detect the deviation points and fit again! --//
        calculateFit(h_LnLik, sEvt, consEvts);

        delete h_LnLik;
      }
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
  h_SMLikelihoodValue->GetXaxis()->SetTitle("Likelihood value at gR = 0.0"); h_SMLikelihoodValue->GetYaxis()->SetTitle("# events");
  h_SMLikelihoodValue->Write();
  h_SMLikValue_AfterCut->Write();
  h_SMLikelihoodValue_vs_DeltaLikelihood->Write();

  TDirectory *dir_FitDevDelete = dir_FitResults->GetDirectory("PointsDeletedByFitDev");
  if(!dir_FitDevDelete) dir_FitDevDelete = dir_FitResults->mkdir("PointsDeletedByFitDev");

  TDirectory *dir_RelFitDevDelete = dir_FitResults->GetDirectory("PointsDeletedByRelFitDev");
  if(!dir_RelFitDevDelete) dir_RelFitDevDelete = dir_FitResults->mkdir("PointsDeletedByRelFitDev");

  dir_FitDevDelete->cd(); h_PointsDelByFitDev->Write();               delete h_PointsDelByFitDev;
  dir_RelFitDevDelete->cd(); h_PointsDelByFitDevRel->Write();         delete h_PointsDelByFitDevRel;
  PaintOverflow(h_ChiSquaredFirstFit, file_FitDist, "FitResults_ChiSquaredFit"); delete h_ChiSquaredFirstFit;
  PaintOverflow(h_ChiSquaredSecondFit,file_FitDist, "FitResults_ChiSquaredFit"); delete h_ChiSquaredSecondFit;
  
  TDirectory *dir_FitSums = file_FitDist->GetDirectory("FitSums");
  if (!dir_FitSums) dir_FitSums = file_FitDist->mkdir("FitSums");
  dir_FitSums->cd();
  
  //Now create the final fit from the individually summed fit parameters!
  std::cout << " " << endl;

  TF1* polFit_histSum = new TF1("polFit_SummedHist","pol2",FitMin, FitMax); 
  histSum->Fit(polFit_histSum,"Q","",polFit_histSum->GetXmin(), polFit_histSum->GetXmax());
  std::cout << " Minimum for " << polFit_histSum->GetName() << " is : " << polFit_histSum->GetMinimumX() << " +- " << polFit_histSum->GetX(polFit_histSum->GetMinimum()+0.5, polFit_histSum->GetMinimumX(),0.2) - polFit_histSum->GetX(polFit_histSum->GetMinimum()+0.5, -0.2, polFit_histSum->GetMinimumX()) << endl;
  histSum->Write();

  TF1* FitSum_FirstFit = new TF1("SummedFit_FirstFit","pol2",FitMin,FitMax);
  TF1* FitSum_SecondFit = new TF1("SummedFit_SecondFit","pol2",FitMin,FitMax);
  FitSum_FirstFit->SetTitle(("Distribution of first fit after summing over "+sNrConfigs+" points ("+sNEvts+" events)").c_str());
  FitSum_SecondFit->SetTitle(("Distribution of second fit after summing over "+sNrConfigs+" points ("+sNEvts+" events)").c_str());

  for(int ipar = 0; ipar < FitSum_FirstFit->GetNpar(); ipar++) FitSum_FirstFit->SetParameter(ipar, FitParams_FirstFit[ipar]);
  FitSum_FirstFit->Write();
  std::cout << " MinimumX value for " << FitSum_FirstFit->GetName() << " is : " << FitSum_FirstFit->GetMinimumX() << " +- " << FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, FitSum_FirstFit->GetMinimumX(), 0.2) - FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, -0.2, FitSum_FirstFit->GetMinimumX()) << endl;
  delete FitSum_FirstFit;

  for(int ipar = 0; ipar < FitSum_SecondFit->GetNpar(); ipar++) FitSum_SecondFit->SetParameter(ipar, FitParams_SecondFit[ipar]);
  FitSum_SecondFit->Write();
  std::cout << " MinimumX value for " << FitSum_SecondFit->GetName() << " is : " << FitSum_SecondFit->GetMinimumX() << " +- " << FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, FitSum_SecondFit->GetMinimumX(), 0.2) - FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, -0.2, FitSum_SecondFit->GetMinimumX()) << endl;
  delete FitSum_SecondFit;

  file_FitDist->Close();

  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
}
*/
