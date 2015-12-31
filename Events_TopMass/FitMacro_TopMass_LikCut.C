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
#include "THStack.h"
#include "TLine.h"
#include "TStyle.h"

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

TGraphErrors* getMinimum(vector< vector< vector<double> > > LnLikArray, vector< vector<double> > SF, vector< vector<double> > likVal, vector<std::string> name, vector<double> norm, double Lumi, TFile* outFile, double Var[]){

  TH1D* h_MinComp = new TH1D("MinimumComparison","Comparison of obtained minimum for the different considered samples",LnLikArray.size()+1, -0.5, LnLikArray.size()+0.5);
  //double LikCutOpt[10]= {0};
  double LikCutOpt[10] = {53, 54, 55, 56, 57, 58, 60, 62, 65, 68};
  //double LikCutOptGen[10] =  {55, 56, 57, 58, 59, 60, 61, 62, 63, 64};
  double likFitMin = 52.5, likFitMax = 68.5;
  //if(GenOrReco == "RECO"){for(int i = 0; i < 10; i++) LikCutOpt[i] = LikCutOptReco[i];}
  //else                   {for(int i = 0; i < 10; i++) LikCutOpt[i] = LikCutOptGen[i]; likFitMin = 54.5; likFitMax = 64;}

  const int NrCuts = 10;

  vector<double> summedEntries; summedEntries.clear();
  vector<double> minimum, error; minimum.clear(); error.clear();
  int nrSamples = 0;
  for(int iCut = 0; iCut < 10; iCut++){

    for(int iFile = 0; iFile < LnLikArray.size(); iFile++){
      if(iFile == 0) nrSamples = 0;
      //if(MCorData == "MC" && name[iFile] == "Data") continue;
      //if(MCorData == "Data" && name[iFile] != "Data") continue;

      nrSamples++;
      vector<double> summedSampleEntries; summedSampleEntries.clear();
      for(int iEvt = 0; iEvt < LnLikArray[iFile].size(); iEvt++){
        for(int iConf = 0; iConf < LnLikArray[iFile][iEvt].size(); iConf++){
          if(iEvt == 0){ 
            summedSampleEntries.push_back(0); 
            if(nrSamples == 1){
              summedEntries.push_back(0);
            }
          }
          if(likVal[iFile][iEvt] < LikCutOpt[iCut]){
            summedSampleEntries[iConf]    += LnLikArray[iFile][iEvt][iConf]*SF[iFile][iEvt]*Lumi*norm[iFile];
          }  
          if(iEvt == LnLikArray[iFile].size()-1){
            if(nrSamples == 1){ summedEntries[iConf] = summedSampleEntries[iConf]; }
            else{               summedEntries[iConf] += summedSampleEntries[iConf];}
          }
        }
      }
    }

    //Do the calculations once all events have been processed:
    double Entries[50];
    for(int i = 0; i < LnLikArray[0][0].size(); i++){   //Just loop over the (9) configurations!
      Entries[i] = summedEntries[i];
    }

    //Minimum obtained using the standard fitting approach
    TGraph* gr_Sum = new TGraph(LnLikArray[0][0].size(), Var, Entries);
    TF1* polFit_Sum = new TF1("polFit_SummedGraph","pol2",170, 175);
    gr_Sum->Fit(polFit_Sum,"Q","",polFit_Sum->GetXmin(), polFit_Sum->GetXmax());
    double Minimum = polFit_Sum->GetMinimumX();
    double Error = (polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, polFit_Sum->GetMinimumX(),175) - polFit_Sum->GetX(polFit_Sum->GetMinimum()+0.5, 170, polFit_Sum->GetMinimumX()))/2.0;
    delete gr_Sum;
    delete polFit_Sum;

    minimum.push_back(Minimum);
    error.push_back(Error);
  }

  std::cout << "\n Considered a total of : " << nrSamples << " samples" << std::endl;

  //Now create the TGraphErrors for all the samples and cuts!
  double Min[20], Err[20];
  for(int i = 0; i < minimum.size(); i++){
    Min[i] = minimum[i]; Err[i] = error[i];  
  }
  
  outFile->cd();
  TGraphErrors* gr_MinComp = new TGraphErrors(NrCuts, LikCutOpt, Min, 0, Err);
  gr_MinComp->SetName("MinComp"); gr_MinComp->SetTitle("Minimum comparison for all samples");
  //if(MCorData == "Data"){    gr_MinComp->SetLineColor(1); gr_MinComp->SetMarkerStyle(20); gr_MinComp->SetMarkerColor(1); }
  //else if(MCorData == "MC"){ 
  gr_MinComp->SetLineColor(4); gr_MinComp->SetMarkerStyle(22); gr_MinComp->SetMarkerColor(4); //}
  gr_MinComp->GetXaxis()->SetTitle("MadWeight -ln(L) cut-value");    gr_MinComp->GetYaxis()->SetTitle("Obtained minimum for gR");
  gr_MinComp->Write();
  //leg->AddEntry(gr_MinComp,(name[iFile]+" events").c_str(),"p");
   
  TCanvas* canv = new TCanvas("canv","canv"); 
  canv->cd();    gr_MinComp->Draw("AP");    canv->Write();
  
  //Add a fit to decide on the optimal cut-value!
  std::cout << " Will be doing a fit to determine the optimal cut value between " << likFitMin << " and " << likFitMax << std::endl;
  TF1* minLikCut = new TF1("minLikCut_pol","pol3");  
  //if(MCorData == "Data")    minLikCut->SetLineColor(LnLikArray.size()+1);
  //else if(MCorData == "MC") 
  minLikCut->SetLineColor(4);
  gr_MinComp->Fit(minLikCut,"Q","",likFitMin,likFitMax);
  std::cout << " Optimal cut position is : " << minLikCut->GetX(172.5, likFitMin, likFitMax) << std::endl;
  gr_MinComp->SetName("LikCutFit");
  gr_MinComp->Write();
  
  return gr_MinComp;
}

//------------------------------------------------------------------//
//---------         Start of the main analysis !!!         ---------//
//------------------------------------------------------------------//
int main(int argc, char *argv[]){
  clock_t start = clock();

  gStyle->SetOptStat(0);

  vector<string> inputFiles;
  inputFiles.push_back("weights_CheckedEvts_SFAdded.out");
  //inputFiles.push_back("weights_CheckedEvts_SFAdded.out");

  std::string consSamples;
  //Considered values and corresponding XS-values
  const int NrConfigs = 6; 
  std::stringstream ssNrConfigs; ssNrConfigs << NrConfigs; std::string sNrConfigs = ssNrConfigs.str();
  double Var[NrConfigs]     = {170,      171,      172,     173,      174,      175     };
  double MGXSCut[NrConfigs] = {0.877965, 0.901163, 0.92436, 0.947558, 0.970755, 0.993953};  //Updated with fine-tuning cuts!

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
  double Luminosity = 37.14979;  
  TH1D* h_SMLikValue_Sum = new TH1D("SMLikelihoodValue_Sum","Distribution of likelihood value at gR = 0.0 (all samples)",85,40,85);
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

    vector<double> FitParametersFirstFit, FitParametersSecondFit;
    FitParametersFirstFit.clear(); FitParametersSecondFit.clear();
    bool doFits = false, allPointsInFit = false;       //Before changing this to true the issue of passing on the fit parameters should be resolved!!
  
    int NrBins = 8; 
    double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
    double xLow = Var[0] - halfBinWidth, xHigh = Var[NrConfigs-1] + halfBinWidth; 

    while( std::getline(ifs,line)){
      std::istringstream iss(line);

      if( iss >> evt >> config >> tf >> weight >> weightUnc){

        if(config == 1 && (consEvts+1) % 5000 == 0) std::cout << " Looking at event : " << consEvts+1 << flush<<"\r";
        stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

        //--- Initialize the event-per-event variables! ---//          
        if(config == 1){
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
        if( Var[config-1] == 173) SMConfig = config-1;

        //---  Fill the LnLik histograms for each event and for all events together  ---//
        LnLik[config-1] = (-log(weight)+log(MGXSCut[config-1]));

        //---  Only perform the fit after all configurations are considered!  ---//
        if( config == NrConfigs){

          //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
          //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
          h_SMLikValue_Sum->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
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
  TFile* outputFile = new TFile("OutFile_LikelihoodCutComparison.root","RECREATE"); 
  
  TCanvas* canv = new TCanvas("SMLikCanvas","SMLikCanvas");
  canv->cd();
  h_SMLikValue_Sum->SetTitle("");
  h_SMLikValue_Sum->GetXaxis()->SetTitle("-ln(L_{MEM}) value evaluated at g_{R} = 0"); h_SMLikValue_Sum->GetXaxis()->SetTitleOffset(1.2);
  h_SMLikValue_Sum->GetYaxis()->SetTitle("# normalised events");                       h_SMLikValue_Sum->GetYaxis()->SetTitleOffset(1.15);
  h_SMLikValue_Sum->DrawNormalized();
  h_SMLikValue_Sum->Write();
  canv->SaveAs("SMLikelihoodValue.pdf");

  TGraphErrors* MC = getMinimum(indivLnLik, scaleFactor, LnLikCutValue, sampleName, normFactor, Luminosity, outputFile, Var);
  TCanvas* canv_MC = new TCanvas("MinComp_MCOnly","MinComp_MCOnly");
  canv_MC->cd();
  MC->SetTitle("");
  MC->Draw("AP");
  TF1 *line = new TF1("line","172.5",58,70);
  line->SetLineColor(2);
  line->Draw("same");
  MC->Draw("P");
  canv_MC->Write();
  canv_MC->SaveAs("MinComp_MCOnly_StraightLineAtZero.pdf"); 

  outputFile->Close();
  cout << "\n Info stored in : " << outputFile->GetName() << std::endl;

  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
  return 0;
}
