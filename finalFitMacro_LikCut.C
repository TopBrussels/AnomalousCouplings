//--------------------------------------------------------------------------------------------------------------------------//
//  Scripts which compares the different root files and directories and directly saves the histogram!                       //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` finalFitMacro_LikCut.C -o finalFitMacro_LikCut  //
//               ./finalFitMacro_LikCut                                                                                     //
//--------------------------------------------------------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <map>
#include "TGraphErrors.h"
#include "TLegend.h"

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

void getMinimum(vector< vector< vector<double> > > LnLikArray, vector< vector<double> > SF, vector< vector<double> > likVal, vector<std::string> name, vector<double> norm, double Lumi, TFile* outFile, double Var[], double FitMin, double FitMax, std::string GenOrReco){

  TH1D* h_MinComp = new TH1D("MinimumComparison","Comparison of obtained minimum for the different considered samples",LnLikArray.size()+1, -0.5, LnLikArray.size()+0.5);
  double LikCutOpt[10]= {0};
  double LikCutOptReco[10] = {60, 62, 63, 64, 65, 66, 67, 68, 69, 120};
  double LikCutOptGen[10] =  {55, 57, 58, 59, 60, 61, 62, 63, 64, 120};
  double likFitMin = 59, likFitMax = 68.5;
  if(GenOrReco == "RECO"){for(int i = 0; i < 10; i++) LikCutOpt[i] = LikCutOptReco[i];}
  else                   {for(int i = 0; i < 10; i++) LikCutOpt[i] = LikCutOptGen[i]; likFitMin = 54.5; likFitMax = 64;}

  const int NrCuts = 10;
  //TLegend* leg = new TLegend(0.7, 0.1, 0.9, 0.3);
  TCanvas* canv = new TCanvas("canv","canv");
  //TLegend* legUnw = new TLegend(0.7, 0.1, 0.9, 0.3);
  TCanvas* canvUnw = new TCanvas("canvUnw","canvUnw");

  std::cout << "\n Considered a total of : " << LnLikArray.size() << " samples" << std::endl;
  vector<double> summedEntries; summedEntries.clear();
  vector<double> summedEntriesUnw; summedEntriesUnw.clear();
  vector<double> minimum, error; minimum.clear(); error.clear();
  vector<double> minimumUnw, errorUnw; minimumUnw.clear(); errorUnw.clear();
  for(int iCut = 0; iCut < 10; iCut++){

    for(int iFile = 0; iFile < LnLikArray.size(); iFile++){
      vector<double> summedSampleEntries; summedSampleEntries.clear();
      vector<double> summedSampleEntriesUnw; summedSampleEntriesUnw.clear();
      for(int iEvt = 0; iEvt < LnLikArray[iFile].size(); iEvt++){
        for(int iConf = 0; iConf < LnLikArray[iFile][iEvt].size(); iConf++){
          if(iEvt == 0){ 
            summedSampleEntries.push_back(0); 
            summedSampleEntriesUnw.push_back(0);
            if(iFile == 0){
              summedEntries.push_back(0);
              summedEntriesUnw.push_back(0);
            }
          }
          if(likVal[iFile][iEvt] < LikCutOpt[iCut]){
            summedSampleEntries[iConf] += LnLikArray[iFile][iEvt][iConf]*SF[iFile][iEvt]*Lumi*norm[iFile];
            summedSampleEntriesUnw[iConf] += LnLikArray[iFile][iEvt][iConf]*SF[iFile][iEvt];
          }  
          if(iEvt == LnLikArray[iFile].size() -1){
            if(iFile == 0){ summedEntries[iConf] = summedSampleEntries[iConf];  summedEntriesUnw[iConf] = summedSampleEntriesUnw[iConf]; }
            else{           summedEntries[iConf] += summedSampleEntries[iConf]; summedEntriesUnw[iConf] += summedSampleEntriesUnw[iConf];}
          }
        }
      }
    }

    //Do the calculations once all events have been processed:
    double Entries[50], EntriesUnw[50];
    for(int i = 0; i < LnLikArray[0][0].size(); i++){   //Just loop over the (9) configurations!
      Entries[i] = summedEntries[i];
      EntriesUnw[i] = summedEntriesUnw[i];
    }

    TGraph* gr_Sum = new TGraph(LnLikArray[0][0].size(), Var, Entries);
    TF1* polFit_Sum = new TF1("polFit_SummedGraph","pol2",FitMin, FitMax);
    gr_Sum->Fit(polFit_Sum,"Q","",polFit_Sum->GetXmin(), polFit_Sum->GetXmax());
    double Minimum = polFit_Sum->GetMinimumX();
    double Error = (polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, polFit_Sum->GetMinimumX(),0.2) - polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, -0.2, polFit_Sum->GetMinimumX()))/2.0;
    delete gr_Sum;
    delete polFit_Sum;

    minimum.push_back(Minimum);
    error.push_back(Error);

    TGraph* gr_SumUnw = new TGraph(LnLikArray[0][0].size(), Var, EntriesUnw);
    TF1* polFit_SumUnw = new TF1("polFit_SummedGraphUnw","pol2",FitMin, FitMax);
    gr_SumUnw->Fit(polFit_SumUnw,"Q","",polFit_SumUnw->GetXmin(), polFit_SumUnw->GetXmax());
    double MinimumUnw = polFit_SumUnw->GetMinimumX();
    double ErrorUnw = (polFit_SumUnw->GetX(polFit_SumUnw->GetMinimum()+0.5, polFit_SumUnw->GetMinimumX(),0.2) - polFit_SumUnw->GetX(polFit_SumUnw->GetMinimum()+0.5,-0.2,polFit_SumUnw->GetMinimumX()))/2.0;
    delete gr_SumUnw;
    delete polFit_SumUnw;

    minimumUnw.push_back(MinimumUnw);
    errorUnw.push_back(ErrorUnw);
  }

  //Now create the TGraphErrors for all the samples and cuts!
  double Min[50], Err[50], MinUnw[50], ErrUnw[50];
  for(int i = 0; i < minimum.size(); i++){
    Min[i] = minimum[i]; MinUnw[i] = minimumUnw[i];
    Err[i] = error[i];   ErrUnw[i] = errorUnw[i];
  }
  
  outFile->cd();
  TGraphErrors* gr_MinComp = new TGraphErrors(NrCuts, LikCutOpt, Min, 0, Err);
  gr_MinComp->SetName("MinComp"); gr_MinComp->SetTitle("Minimum comparison for all samples");
  gr_MinComp->SetLineColor(1);    gr_MinComp->SetMarkerStyle(22);    gr_MinComp->SetMarkerColor(1);
  gr_MinComp->GetXaxis()->SetTitle("MadWeight -ln(L) cut-value");    gr_MinComp->GetYaxis()->SetTitle("Obtained minimum for gR");
  gr_MinComp->Write();
  //leg->AddEntry(gr_MinComp,(name[iFile]+" events").c_str(),"p");
   
  TGraphErrors* gr_MinCompUnw = new TGraphErrors(NrCuts, LikCutOpt, MinUnw, 0, ErrUnw);
  gr_MinCompUnw->SetName("MinCompUnw"); gr_MinCompUnw->SetTitle("Minimum comparison for all samples");
  gr_MinCompUnw->SetLineColor(1); gr_MinCompUnw->SetMarkerStyle(22); gr_MinCompUnw->SetMarkerColor(1);
  gr_MinCompUnw->GetXaxis()->SetTitle("MadWeight -ln(L) cut-value"); gr_MinCompUnw->GetYaxis()->SetTitle("Obtained minimum for gR");
  gr_MinCompUnw->Write();
  //legUnw->AddEntry(gr_MinCompUnw,(name[iFile]+" events (unw)").c_str(),"p");

  canv->cd();    gr_MinComp->Draw("AP");
  canvUnw->cd(); gr_MinCompUnw->Draw("AP");
  canv->Write();
  canvUnw->Write();
  
  //Add a fit to decide on the optimal cut-value!
  std::cout << " Will be doing a fit to determine the optimal cut value between " << likFitMin << " and " << likFitMax << std::endl;
  TF1* minLikCut = new TF1("minLikCut_pol","pol3");  minLikCut->SetLineColor(LnLikArray.size()+1);
  gr_MinComp->Fit(minLikCut,"Q","",likFitMin,likFitMax);
  std::cout << " Optimal cut position is : " << minLikCut->GetX(0, likFitMin, likFitMax) << std::endl;
  gr_MinComp->SetName("LikCutFit");
  gr_MinComp->Write();

}

//------------------------------------------------------------------//
//---------         Start of the main analysis !!!         ---------//
//------------------------------------------------------------------//
int main(int argc, char *argv[]){
  clock_t start = clock();

  //Input variables will be:
  // 1) MG or Reco
  // 2-..) Weight files which will be considered

  std::string MGorRECO = "RECO";
  if( argc >= 2)
    MGorRECO = string(argv[1]);

  vector<string> inputFiles;
  if( argc >= 3 ){
    for(int iFile = 2; iFile < argc; iFile++){
      inputFiles.push_back(string(argv[iFile]).c_str());
      std::cout << " - Stored file name is : " << inputFiles[inputFiles.size()-1] << std::endl;
    }
  }
  std::string consSamples;

  //Considered values and corresponding XS-values
  const int NrConfigsReco = 9;
  double VarReco[NrConfigsReco]     = {-0.2,     -0.15,   -0.1,    -0.05,   0.0,     0.05,    0.1,     0.15,    0.2    };
  double MGXSCutReco[NrConfigsReco] = {0.450287, 0.540131, 0.648626, 0.77937, 0.935959, 1.12199, 1.34106, 1.59677, 1.8927};  //Updated with fine-tuning cuts!
  //double MGXSCutReco[NrConfigsReco] = {0.947244, 1.13624,  1.36448,  1.63952, 1.96892,  2.36027, 2.82111, 3.35903, 3.98157};
  const int NrConfigsGen = 13;
  double VarGen[NrConfigsGen]     = {-0.4,    -0.3,    -0.2,     -0.15,   -0.1,    -0.05,   0.0,     0.05,    0.1,     0.15,    0.2,    0.3,      0.4     };
  double MGXSCutGen[NrConfigsGen] = {0.93159, 1.27966, 1.825208, 2.194079, 2.6393, 3.17698, 3.80921, 4.5645, 5.45665, 6.47791, 7.66805, 10.63243, 14.46786};
  
  //Now set the correct ones (Gen or Reco)
  int NrConfigs = 0; 
  vector<double> MGXSCut; MGXSCut.clear();
  double Var[13] = {0};
  if(MGorRECO == "RECO"){    NrConfigs = NrConfigsReco; for(int i = 0; i < NrConfigs; i++){ MGXSCut.push_back(MGXSCutReco[i]); Var[i] = VarReco[i];}}
  else if(MGorRECO == "MG"){
    NrConfigs = NrConfigsGen; 
    for(int i = 0; i < NrConfigs; i++){ 
      MGXSCut.push_back(MGXSCutGen[i]);
      Var[i] = VarGen[i];
    }
  }
  else{ std::cout << " ERROR: should specify whether it is MG or RECO ! " << std::endl; exit(-1);}
  std::stringstream ssNrConfigs; ssNrConfigs << NrConfigs; std::string sNrConfigs = ssNrConfigs.str();

  double FitMin = -0.15, FitMax = 0.15;
  double LnLik[13] = {0.0};    //Should find a way to initialize this ... maybe use vector ...

  //vector< vector<TH1D> > indivLnLik;
  vector< vector< vector<double> > > indivLnLik;
  vector< vector<double> > LnLikCutValue;
  vector< std::string > sampleName; sampleName.clear();
  vector< double > normFactor;     normFactor.clear();
  vector< vector<double> > scaleFactor;

  //Store the TH1Ds in a map such that they can easily be looped over!
  std::map<string,TH1D*> histo1D;  

  //** Loop over the different input files and read out the necessary info **//
  double Luminosity = 19646.8;
  for(int iWeightFile = 0; iWeightFile < inputFiles.size(); iWeightFile++){

    //Make sure that for each dataset the vector storing the individual histograms is empty
    vector< vector<double> > indivSampleLnLik; indivSampleLnLik.clear();
    vector<double> sampleSF; sampleSF.clear();
    vector<double> LnLikCut;

    //--- Read the likelihood values ! ---//
    std::ifstream ifs (inputFiles[iWeightFile].c_str(), std::ifstream::in); 
    std::string line;
    int evt,config,tf;
    double weight, weightUnc;
    double MCScaleFactor = 1, NormFactor = 1;
    double SampleName;
    std::string sSampleName, sNormFactor;
    bool sampleNameSet = false;
    int SMConfig = 99, consEvts = 0;
    double MaxLik = 0, MinLik = 9999;

    vector<double> FitParametersFirstFit, FitParametersSecondFit;
    FitParametersFirstFit.clear(); FitParametersSecondFit.clear();
    bool doFits = false, allPointsInFit = false;       //Before changing this to true the issue of passing on the fit parameters should be resolved!!
  
    int NrBins = 8; 
    double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
    double xLow = Var[0] - halfBinWidth, xHigh = Var[NrConfigs-1] + halfBinWidth; 

    while( std::getline(ifs,line)){
      std::istringstream iss(line);

      if( iss >> evt >> config >> tf >> weight >> weightUnc){
        if(MGorRECO == "RECO"){ iss >> MCScaleFactor;}

        if(config == 1 && (consEvts+1) % 5000 == 0) std::cout << " Looking at event : " << consEvts+1 << flush<<"\r";
        stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

        //--- Initialize the event-per-event variables! ---//          
        if(config == 1){
          MaxLik = 0; MinLik = 9999; //Reset this variable at the beginning of each event!
          SMConfig = 99;

          if(consEvts == 0){
            //Loop over all histo1Ds and add the sampleName!
            if(histo1D.size() > 0){
              for(std::map<std::string,TH1D*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
                TH1D *temp = it->second;
                temp->SetName((string(temp->GetName())+"_"+sampleName[iWeightFile]).c_str());
              }
            }
          }
        }

        //Set the SMConfig:
        if( Var[config-1] == 0.0) SMConfig = config-1;

        //---  Fill the LnLik histograms for each event and for all events together  ---//
        LnLik[config-1] = (-log(weight)+log(MGXSCut[config-1]));

        //Get the maximum and minimum likelihood value
        if(LnLik[config-1] > MaxLik) MaxLik = LnLik[config-1];
        if(LnLik[config-1] < MinLik) MinLik = LnLik[config-1];

        //---  Only perform the fit after all configurations are considered!  ---//
        if( config == NrConfigs){

          //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
          //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
          LnLikCut.push_back(LnLik[SMConfig]);
          consEvts++;   //Count the number of full events!

          vector<double> indivEvtLnLik;
          indivEvtLnLik.clear();

          for(int i=0; i< NrConfigs; i++)
            indivEvtLnLik.push_back(LnLik[i]);
          indivSampleLnLik.push_back(indivEvtLnLik);
          sampleSF.push_back(MCScaleFactor);
        }
      }
      else{
        sampleNameSet = true;
        std::cout << " In ELSE -> Line = " << line << std::endl;
        sSampleName = line.substr(0, line.find(" "));
        sNormFactor = line.substr(sSampleName.length(), line.length());
        NormFactor = atof(sNormFactor.c_str());
        std::cout << "  * What is sample name : " << sSampleName << std::endl;
        consSamples += sSampleName;
        if(iWeightFile != inputFiles.size()-1) consSamples += "AND"; 
        sampleName.push_back(sSampleName);
        std::cout << "  * What is norm-factor : " << NormFactor << std::endl;    //To notice: 1/normFactor = lumi of sample!
        normFactor.push_back(NormFactor);
      }
    }
    ifs.close();

    //Push back all the vectors containing event-per-event information!
    indivLnLik.push_back(indivSampleLnLik);
    scaleFactor.push_back(sampleSF);
    LnLikCutValue.push_back(LnLikCut);

  }//End of looping over the different weight files
  
  //Store all information into a ROOT file:
  TFile* outputFile = new TFile(("Events_Nom/OutFile_LikelihoodCutComparison_"+consSamples+".root").c_str(),"RECREATE");   //So what if also Data is added?

  getMinimum(indivLnLik, scaleFactor, LnLikCutValue, sampleName, normFactor, Luminosity, outputFile, Var, FitMin, FitMax, MGorRECO);

  outputFile->Close();
  cout << " Info stored in : " << outputFile->GetName() << std::endl;

  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
  return 0;
}
