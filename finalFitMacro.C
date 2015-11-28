//------------------------------------------------------------------------------------------------------------//
//  Scripts which compares the different root files and directories and directly saves the histogram!         //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` finalFitMacro.C -o finalFitMacro  //
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

//Convert integer to string!
template <typename T>
string NumberToString ( T Number ){
    stringstream ss; ss << Number;
    return ss.str();
}

//Add a timestamp!
std::string timestamp(){
    time_t ltime;
    struct tm *Tm;
 
    ltime=time(NULL);
    Tm=localtime(&ltime);

    std::string TIME = (NumberToString(Tm->tm_mday)+NumberToString(Tm->tm_mon+1)+NumberToString(Tm->tm_year+1900)+"_"+NumberToString(Tm->tm_hour)+NumberToString(Tm->tm_min)+NumberToString(Tm->tm_sec)).c_str();
    return TIME;
}

//Method to sort a pair based on the second element!
struct sort_pred {
  bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
    return left.second < right.second;
  }
};

void calculateFit(double LogLikelihood[], string EvtNumber, int evtCounter, bool consAllPoints, const int NConfigs, vector<double> FitParamsFirstFit, vector<double> FitParamsSecondFit, double Var[], double FitMin, double FitMax){
  const unsigned int NrToDel = 2; 

  std::string whichFit = "AllPoints";
  if(!consAllPoints) whichFit = "ReducedPoints";

  //Set name of chisquared distributions!
  //if(evtCounter == 1){
  //  if(consAllPoints) h_ChiSquaredFirstFit  = new TH1F("ChiSquared_FirstFit", "Distribution of the #chi^{2} after the fit on all the points",200,0,0.05);
  //  if(!consAllPoints){
  //    h_ChiSquaredSecondFit = new TH1F("ChiSquared_SecondFit","Distribution of the #chi^{2} after the fit on the reduced points",200,0,0.005);
  //    h_PointsDelByFitDev    = new TH1F("DeletedPointsFitDev",   "Overview of deleted points due to largest FitDeviation",NrBins+1,xLow,xHigh);
  //    h_PointsDelByFitDevRel = new TH1F("DeletedPointsFitDevRel","Overview of deleted points due to largest relative FitDeviation",NrBins+1,xLow,xHigh);
  //  }
  //}
 
  TF1* polFit = new TF1(("polFit_"+whichFit+"_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
  TGraph* gr_LnLik = new TGraph(NConfigs,Var, LogLikelihood);
  gr_LnLik->Fit(polFit,"Q","",polFit->GetXmin(), polFit->GetXmax());
  //if(consAllPoints) h_ChiSquaredFirstFit->Fill(polFit->GetChisquare());

  if(consAllPoints && evtCounter == 1){
    for(int ipar = 0; ipar < polFit->GetNpar(); ipar++)
      FitParamsFirstFit.push_back(polFit->GetParameter(ipar));
  }
  else if(consAllPoints && evtCounter != 1){
    for(int ipar = 0; ipar < polFit->GetNpar(); ipar++)
      FitParamsFirstFit[ipar] += polFit->GetParameter(ipar);
  }

  if(!consAllPoints){
    std::vector<std::pair<int, double> > FitDeviation, FitDeviationRel;
    double LogLikFit[] = {-9999};
    //double TotalRelFitDeviation = 0;
    for(int iConfig = 0; iConfig < NConfigs; iConfig++){
      if(polFit->GetNpar() == 3)
        LogLikFit[iConfig] = polFit->GetParameter(0)+polFit->GetParameter(1)*Var[iConfig]+polFit->GetParameter(2)*Var[iConfig]*Var[iConfig];
      else if(polFit->GetNpar() == 5)
        LogLikFit[iConfig] = polFit->GetParameter(0)+polFit->GetParameter(1)*Var[iConfig]+polFit->GetParameter(2)*Var[iConfig]*Var[iConfig]+polFit->GetParameter(3)*pow(Var[iConfig],3)+polFit->GetParameter(4)*pow(Var[iConfig],4)+polFit->GetParameter(5)*pow(Var[iConfig],5);
      FitDeviation.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig]) ) );
      FitDeviationRel.push_back( std::make_pair(iConfig, abs(LogLikelihood[iConfig]-LogLikFit[iConfig])/LogLikelihood[iConfig] ) );
      //TotalRelFitDeviation += abs(LogLikelihood[iConfig]-LogLikFit[iConfig])/LogLikelihood[iConfig];
    }
    //Sort the fitdeviation values depending on the second value!
    std::sort(FitDeviation.begin(), FitDeviation.end(), sort_pred() );
    std::sort(FitDeviationRel.begin(), FitDeviationRel.end(), sort_pred() );
    //h_TotalRelFitDeviation->Fill(TotalRelFitDeviation);

    //Now loop again over all configurations, plot the FitDeviation in sorted order and save the configNr's which should be excluded 
    std::vector<int> FitDevPointsToDel, FitDevRelPointsToDel;
    for(int itSortedConfig = NConfigs-1; itSortedConfig >= 0 ; itSortedConfig--){  //Looping from high to low values of the deviation!
      //h_FitDeviation[itSortedConfig]->Fill(FitDeviation[itSortedConfig].second);
      //h_FitDeviationRel[itSortedConfig]->Fill(FitDeviationRel[itSortedConfig].second);
    
      //Store the 'NrToDel' points which need to be excluded from the TGraph!
      if(FitDevPointsToDel.size() < NrToDel){ FitDevPointsToDel.push_back(FitDeviation[itSortedConfig].first); }//h_PointsDelByFitDev->Fill(Var[FitDeviation[itSortedConfig].first]);
      if(FitDevRelPointsToDel.size() < NrToDel){FitDevRelPointsToDel.push_back(FitDeviationRel[itSortedConfig].first); }//h_PointsDelByFitDevRel->Fill(Var[FitDeviationRel[itSortedConfig].first]);
    }

    //Create new arrays with the reduced information!
    double ReducedVar[NConfigs-NrToDel], ReducedLogLik[NConfigs-NrToDel];
    //double TotalRelFitDeviationReduced = 0;
    int iCounter = 0;
    for(int iConf = 0; iConf < NConfigs; iConf++){
      if(!(std::find(FitDevPointsToDel.begin(), FitDevPointsToDel.end(), iConf) != FitDevPointsToDel.end()) ){
        ReducedVar[iCounter] = Var[iConf];
        ReducedLogLik[iCounter] = LogLikelihood[iConf];
        //TotalRelFitDeviationReduced += (abs(LogLikelihood[iConf]-LogLikFit[iConf])/LogLikelihood[iConf]);
        iCounter++;
      }
    }
    //h_TotalRelFitDeviationReduced->Fill(TotalRelFitDeviationReduced);

    //Define new TGraph and fit again
    TGraph* gr_ReducedLnLik = new TGraph(NConfigs-NrToDel, ReducedVar, ReducedLogLik);
    TF1* polFit_ReducedPoints = new TF1(("polFit_ReducedPoints_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
    gr_ReducedLnLik->Fit(polFit_ReducedPoints,"Q","",polFit->GetXmin(), polFit->GetXmax()); 
    //h_ChiSquaredSecondFit->Fill(polFit_ReducedPoints->GetChisquare());   //As expected NDF is always equal to NConfigs-NrToDel-3 (= nr params needed to define a parabola)

    //Add the different parameters of the individual fits in order to form a final fit!
    if(evtCounter == 1){
      for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
        FitParamsSecondFit.push_back(polFit_ReducedPoints->GetParameter(ipar));
    }
    else{
      for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
        FitParamsSecondFit[ipar] += polFit_ReducedPoints->GetParameter(ipar);
    }

    //Create a 2D-plot which contains the deviation of the expected minimum wrt the chi-squared of the fit
    //h_TotalFitDevVSChiSq->Fill(polFit_ReducedPoints->GetChisquare(), TotalRelFitDeviationReduced);

    //polFit_ReducedPoints->Write();
    delete polFit_ReducedPoints;
    delete gr_ReducedLnLik;
  }
  
  delete polFit;
  delete gr_LnLik;
}

//------------------------------------------------------------------//
//---------         Start of the main analysis !!!         ---------//
//------------------------------------------------------------------//
int main(int argc, char *argv[]){
  clock_t start = clock();

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

  TFile* outputFile = new TFile(("Events_MC/OutFile_"+timestamp()).c_str(),"RECREATE");   //So what if also Data is added?
  outputFile->cd();

  //Considered values and corresponding XS-values
  const int NrConfigs = 9; 
  std::stringstream ssNrConfigs; ssNrConfigs << NrConfigs; std::string sNrConfigs = ssNrConfigs.str();
  double Var[NrConfigs]     = {-0.2,     -0.15,   -0.1,    -0.05,   0.0,     0.05,    0.1,     0.15,    0.2    };
  double MGXSCut[NrConfigs] = {0.947244, 1.13624, 1.36448, 1.63952, 1.96892, 2.36027, 2.82111, 3.35903, 3.98157};

  double FitMin = -0.15, FitMax = 0.15;
  double LnLik[NrConfigs] = {0.0};

  //** Loop over the different input files and read out the necessary info **//
  double Luminosity = 19646.8;
  for(int iWeightFile = 0; iWeightFile < inputFiles.size(); iWeightFile++){

    //--- Read the likelihood values ! ---//
    std::ifstream ifs (inputFiles[iWeightFile].c_str(), std::ifstream::in); 
    std::string line;
    int evt,config,tf;
    double weight, weightUnc;
    double CosThetaCorr = 1;
    double MCScaleFactor = 1, Luminosity = 1, NormFactor = 1;
    double SampleName;
    std::string sSampleName;
    bool sampleNameSet = false;
    int SMConfig = 99, consEvts = 0;
    double MaxLik = 0, MinLik = 9999;

    vector<double> FitParametersFirstFit, FitParametersSecondFit;
    FitParametersFirstFit.clear(); FitParametersSecondFit.clear();
    bool doFits = true, allPointsInFit = true; 
  
    TH1F *h_LnLik = 0;
    int NrBins = 8; 
    double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
    float xLow = Var[0] - halfBinWidth, xHigh = Var[NrConfigs-1] + halfBinWidth; 
    TH1F *histSum = new TH1F("SummedHist","Sum of all individual histograms",NrBins+1,xLow,xHigh); 
    TH1F* h_SMLikelihoodValue = new TH1F("SMLikelihoodValue","Distribution of likelihood value at gR = 0.0",500,30,90);
    TH1F* h_SMLikValue_AfterCut = new TH1F("SMLikelihoodValue_AfterCut","Distribution of likelihood value at gR = 0.0 (after cut of ..)",500,30,90);

    //while( std::getline(ifs,line) && consEvts < nEvts){
    while( std::getline(ifs,line)){
      std::istringstream iss(line);

      if( iss >> evt >> config >> tf >> weight >> weightUnc >> MCScaleFactor ){
        if(!sampleNameSet) std::cout << " Already went into if .... " << std::endl;

        //if(config == 1 && ((consEvts+1) % 2000 == 0) ) std::cout << " Looking at event : " << consEvts+1 << " (" << (double)(consEvts+1)*100/(double)nEvts << "%)" << flush<<"\r";
        if(config == 1 && ((consEvts+1) % 2000 == 0) ) std::cout << " Looking at event : " << consEvts+1 << flush<<"\r";
        stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

        //--- Initialize the event-per-event variables! ---//          
        if( config == 1){
          MaxLik = 0; MinLik = 9999; //Reset this variable at the beginning of each event!
          SMConfig = 99;

          h_LnLik = new TH1F(("LnLik_Evt"+sEvt+"_"+sSampleName).c_str(),("LnLik distribution for event "+sEvt+" -- "+sSampleName).c_str(),NrBins+1,xLow,xHigh);
          h_LnLik->SetMarkerStyle(20); h_LnLik->SetLineColor(3); h_LnLik->SetMarkerColor(1); h_LnLik->SetMarkerSize(1.2);
        }

        //Set the SMConfig:
        if( Var[config-1] == 0.0) SMConfig = config-1;

        //---  Fill the LnLik histograms for each event and for all events together  ---//
        LnLik[config-1] = (-log(weight)+log(MGXSCut[config-1]));
        h_LnLik->SetBinContent(h_LnLik->FindBin(Var[config-1]), LnLik[config-1]*MCScaleFactor*Luminosity*NormFactor);

        //Get the maximum and minimum likelihood value
        if(LnLik[config-1] > MaxLik) MaxLik = LnLik[config-1];
        if(LnLik[config-1] < MinLik) MinLik = LnLik[config-1];

        //---  Only perform the fit after all configurations are considered!  ---//
        if( config == NrConfigs){

          //Plot the SM likelihood value:
          h_SMLikelihoodValue->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);

          //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
          //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
          if( LnLik[SMConfig] > LikCut){ if(LikCut == 100) std::cout << " ******** \n ERROR: Should not reject any event when cut-value = 100 !!! \n ******** \n " << std::endl; delete h_LnLik; continue; }
          consEvts++;   //Count the number of full events!

          //Initialize the fitDeviation histograms (array of TH1F, one for each configuration ...)
          //if(consEvts == 1){
          //  for(int ii = 0; ii < NrConfigs; ii++){
          //    stringstream ssConfig; ssConfig << ii+1; string sConfig = ssConfig.str();
          //    h_FitDeviation[ii] = new TH1F(("FitDeviation_LowestConfig"+sConfig).c_str(),("FitDeviation histogram for "+sConfig+"st lowest value").c_str(), 200, 0, 0.1);
          //    h_FitDeviationRel[ii] = new TH1F(("FitDeviationRelative_LowestConfig"+sConfig).c_str(), ("FitDeviationRelative histogram for "+sConfig+"st lowest value").c_str(),200,0,0.005);
          //  }
          //}

          h_SMLikValue_AfterCut->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
          histSum->Add(h_LnLik);
          //h_SMLikelihoodValue_vs_DeltaLikelihood->Fill(LnLik[SMConfig], MaxLik-MinLik);

          //-- Send the array containing the ln(weights) to the predefined function to define the TGraph, fit this, detect the deviation points and fit again! --//
          if(doFits)
            calculateFit(LnLik, sEvt, consEvts, allPointsInFit, NrConfigs, FitParametersFirstFit, FitParametersSecondFit, Var, FitMin, FitMax);

          delete h_LnLik;
        }
      }
      /*if( iss >> evt >> config >> tf >> weight ){ 
        std::cout << " In else if line : " << line << std::endl;
        //stringstream ssSampleName; ssSampleName << SampleName; sSampleName = ssSampleName.str();
        std::cout << " Looking at line : " << line << std::endl;
      }*/
      else{
        sampleNameSet = true;
        iss >> sSampleName >> NormFactor;
        //stringstream ssSampleName; ssSampleName << SampleName; sSampleName = ssSampleName.str();
        std::cout << " What is sample name : " << sSampleName.c_str() << std::endl; 
        std::cout << " What is norm-factor : " << NormFactor << std::endl;
        std::cout << " ELSE -> Line = " << line << std::endl;
      }
    }
    ifs.close();

    TF1* polFit_histSum = new TF1("polFit_SummedHist","pol2",FitMin, FitMax); 
    histSum->Fit(polFit_histSum,"Q","",polFit_histSum->GetXmin(), polFit_histSum->GetXmax());
    std::cout << " Minimum for " << polFit_histSum->GetName() << " is : " << polFit_histSum->GetMinimumX() << " +- " << polFit_histSum->GetX(polFit_histSum->GetMinimum()+0.5, polFit_histSum->GetMinimumX(),0.2) - polFit_histSum->GetX(polFit_histSum->GetMinimum()+0.5, -0.2, polFit_histSum->GetMinimumX()) << endl;
    histSum->Write();

    if(doFits){
      TF1* FitSum_FirstFit = new TF1("SummedFit_FirstFit","pol2",FitMin,FitMax);
      TF1* FitSum_SecondFit = new TF1("SummedFit_SecondFit","pol2",FitMin,FitMax);
      FitSum_FirstFit->SetTitle(("Distribution of first fit after summing over "+sNrConfigs+" points").c_str());
      FitSum_SecondFit->SetTitle(("Distribution of second fit after summing over "+sNrConfigs+" points").c_str());

      for(int ipar = 0; ipar < FitSum_FirstFit->GetNpar(); ipar++) FitSum_FirstFit->SetParameter(ipar, FitParametersFirstFit[ipar]);
      FitSum_FirstFit->Write();
      std::cout << " MinimumX value for " << FitSum_FirstFit->GetName() << " is : " << FitSum_FirstFit->GetMinimumX() << " +- " << FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, FitSum_FirstFit->GetMinimumX(), 0.2) - FitSum_FirstFit->GetX(FitSum_FirstFit->GetMinimum()+0.5, -0.2, FitSum_FirstFit->GetMinimumX()) << endl;
      delete FitSum_FirstFit;

      for(int ipar = 0; ipar < FitSum_SecondFit->GetNpar(); ipar++) FitSum_SecondFit->SetParameter(ipar, FitParametersSecondFit[ipar]);
      FitSum_SecondFit->Write();
      std::cout << " MinimumX value for " << FitSum_SecondFit->GetName() << " is : " << FitSum_SecondFit->GetMinimumX() << " +- " << FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, FitSum_SecondFit->GetMinimumX(), 0.2) - FitSum_SecondFit->GetX(FitSum_SecondFit->GetMinimum()+0.5, -0.2, FitSum_SecondFit->GetMinimumX()) << endl;
      delete FitSum_SecondFit;
    }

  }//End of looping over the different weight files

  outputFile->Close();
  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
  return 0;
}

/*

  //Information for the histograms
  int NrBins = 8; 
  float xLow, xHigh;
  const unsigned int NrToDel = 2; 
  int NrRemaining = NrConfigs-NrToDel;

std::string sNrRemaining = ""; std::stringstream ssNrRemaining; 
std::string sNEvts = ""; std::stringstream ssNEvts;
std::string sNrConfigs = ""; std::stringstream ssNrConfigs;

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
    }
  }

  h_tmp->Write();  
  FileToWrite->cd();  //Reset to general directory! 

  delete h_tmp;
}                  


void doublePolFitMacro(){

  
  double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
  xLow = Var[0] - halfBinWidth; xHigh = Var[NrConfigs-1] + halfBinWidth; 
  
  ssNrRemaining << NrRemaining; sNrRemaining = ssNrRemaining.str();
  ssNEvts << nEvts; sNEvts = ssNEvts.str();
  ssNrConfigs << NrConfigs; sNrConfigs = ssNrConfigs.str();

  TH1F *histSum = new TH1F("SummedHist","Sum of all individual histograms",NrBins+1,xLow,xHigh); 
  TH1F *h_LnLik = 0;

  TCanvas *canv_SplitLL = 0;
  double LnLik[NrConfigs] = {0.0}; //, LnLikXS[NrConfigs] = {0.0}, LnLikAcc[NrConfigs] = {0.0};        


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

}
*/
