//------------------------------------------------------------------------------------------------------------------//
//  Scripts which compares the different root files and directories and directly saves the histogram!               //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` FitMacro_TopMass.C -o FitMacro_TopMass  //
//               ./finalFitMacro                                                                                    //
//------------------------------------------------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
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
#include "THStack.h"
#include "TLegend.h"

using namespace std;

//Convert integer to string!
template <typename T>
string NumberToString ( T Number ){
    stringstream ss; ss << Number;
    return ss.str();
}

void getMinimum(vector< vector< vector<double> > > LnLikArray, vector< vector<double> > SF, vector<std::string> name, vector<double> norm, double Lumi, TFile* outFile, double Var[], std::string dir){

  TH1D* h_MinComp = 0;
  h_MinComp = new TH1D("MinimumComparison","Comparison of obtained minimum for the different considered samples",LnLikArray.size()+1, -0.5, LnLikArray.size()+0.5);

  std::cout << "\n\n ** Calculating minimum for " << LnLikArray.size() << " samples ** " << std::endl;
  vector<double> summedEntries; summedEntries.clear();
  for(int iFile = 0; iFile < LnLikArray.size(); iFile++){
    std::cout << "  - Considering : " << name[iFile] << " with " << LnLikArray[iFile].size() << " stored events " << std::endl;

    vector<double> summedSampleEntries; summedSampleEntries.clear();
    for(int iEvt = 0; iEvt < LnLikArray[iFile].size(); iEvt++){

      for(int iConf = 0; iConf < LnLikArray[iFile][iEvt].size(); iConf++){
        if(iEvt == 0){
          summedSampleEntries.push_back(0);
          if(iFile == 0){
            summedEntries.push_back(0);
          }
        }
        summedSampleEntries[iConf] += LnLikArray[iFile][iEvt][iConf]*SF[iFile][iEvt]*Lumi*norm[iFile];
        if(iEvt == LnLikArray[iFile].size() -1){
          if(iFile == 0) summedEntries[iConf] = summedSampleEntries[iConf];
          else           summedEntries[iConf] += summedSampleEntries[iConf];
        }
      }
    }

    double SampleEntries[15];
    for(int i = 0; i < LnLikArray[iFile][0].size(); i++){
      SampleEntries[i] = summedSampleEntries[i];
    }
  }
  
  double Entries[15];
  double MinEntries = summedEntries[0]; for(int i = 0; i < LnLikArray[0][0].size(); i++){ if(MinEntries > summedEntries[i]) MinEntries = summedEntries[i]; }
  double ySM = 0, yNeg = 0, yPos = 0;
  for(int iEnt = 0; iEnt < LnLikArray[0][0].size(); iEnt++){
    Entries[iEnt] = 2*(summedEntries[iEnt] - MinEntries);
  }
  TGraph* gr_totalSum_NoFit = new TGraph(LnLikArray[0][0].size(), Var, Entries);
  gr_totalSum_NoFit->SetMarkerStyle(22);
  gr_totalSum_NoFit->SetMarkerSize(1.3);
  gr_totalSum_NoFit->SetMinimum(-10);
  gr_totalSum_NoFit->GetXaxis()->SetTitle("top-quark mass");
  gr_totalSum_NoFit->GetYaxis()->SetTitle("#Delta #chi^{2}_{MEM} value"); gr_totalSum_NoFit->GetYaxis()->SetTitleOffset(1.5);
  gr_totalSum_NoFit->SetTitle("");

  TCanvas* canv_totalSum_NoFit = new TCanvas("canv_totalSum_NoFit","canv_totalSum_NoFit");
  canv_totalSum_NoFit->cd();
  gr_totalSum_NoFit->Draw("AP");
  canv_totalSum_NoFit->SaveAs((dir+"/OverallLikelihoodCurve_NoFit_"+dir+".pdf").c_str());

  TGraph* gr_totalSum = new TGraph(LnLikArray[0][0].size(), Var, Entries);
  gr_totalSum->SetMarkerStyle(22);
  gr_totalSum->SetMarkerSize(1.3);
  gr_totalSum->SetMinimum(-10);
  gr_totalSum->GetXaxis()->SetTitle("top-quark mass (GeV)");
  gr_totalSum->GetYaxis()->SetTitle("#Delta #chi^{2}_{MEM} value");
  gr_totalSum->GetYaxis()->SetTitleOffset(1.5);

  TDirectory* grdir = outFile->GetDirectory("graphs");   //Check whether directory already exists ..
  if(!grdir) grdir = outFile->mkdir("graphs");          // .. and otherwise create it!
  grdir->cd();

  //Fit through the reduced range!
  TF1* polFit_totalSum = new TF1("polFit_SummedTotalGraph","pol2",171, 175);
  gr_totalSum->Fit(polFit_totalSum,"Q","",polFit_totalSum->GetXmin(), polFit_totalSum->GetXmax());
  double MinimumComb = polFit_totalSum->GetMinimumX();
  //double ErrorComb = (polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+0.5, polFit_totalSum->GetMinimumX(),0.2) - polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+0.5, -0.2, polFit_totalSum->GetMinimumX()))/2.0;
  double ErrorComb = (polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+1.0, polFit_totalSum->GetMinimumX(),175) - polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+1.0, 170, polFit_totalSum->GetMinimumX()))/2.0;

  std::cout << " * Minimum for " << polFit_totalSum->GetName() << " is : \n            " << MinimumComb << " \\pm " << ErrorComb << endl;

  TCanvas *canv_totalSum = new TCanvas("canv_FittedSummedGraph","canv_FittedSummedGraph");
  canv_totalSum->cd();
  gr_totalSum->SetTitle("");
  gr_totalSum->Draw("AP");
  canv_totalSum->Write();
  canv_totalSum->SaveAs((dir+"/OverallLikelihoodCurve_"+dir+".pdf").c_str());

  gr_totalSum->SetName("FittedGraph_AllSamples");
  gr_totalSum->SetTitle("Graph containing the summed events for all samples");
  gr_totalSum->Write();
}

//------------------------------------------------------------------//
//---------         Start of the main analysis !!!         ---------//
//------------------------------------------------------------------//
int main(int argc, char *argv[]){
  clock_t start = clock();

  //Input variables will be:
  // 1) Likelihood cut value!
  // 2) Which combination of files to considerd (DataMC, MC, TTbarJets, file)
  // 3-..) Weight files which will be considered (only if 3rd argument is equal to file!!)

  double LikCut = 120;
  if( argc >= 2){
    if( string(argv[1]) == "opt") LikCut = 63.869;
    else                          LikCut = atof(argv[1]);
  }
  std::cout << " - Applied likelihood cut value is : " << LikCut << std::endl;
  std::string sLikCut = string(NumberToString(LikCut));

  //Also give the luminosity from the command line (easy to change between ttbar (semi-lept) and data lumi!
  double Luminosity = 37.14979;
  if(argc >= 3){
    if(string(argv[2]) == "MCLumi")        Luminosity = 226504.2973;
    else if(string(argv[2]) == "EvtLumi")  Luminosity = 37.14979; //corresponds to 4000 TT events!
    else if(string(argv[2]) == "DataLumi") Luminosity = 19646.8;
    else if(string(argv[2]) == "NoLumi")   Luminosity = 1;
    else                                   Luminosity = atof(argv[3]);
  }
  std::cout << " - Applied luminosity is : " << Luminosity << std::endl;

  std::string whichDir = "";
  if( argc >= 3)
    whichDir = string(argv[3]);
  vector<string> inputFiles;
  if(whichDir.find("TTSemiLept") < whichDir.size())
    inputFiles.push_back((whichDir+"/weights_CheckedEvts_SFAdded.out").c_str());
  else
    inputFiles.push_back((whichDir+"/weights_CheckedEvts.out").c_str());
  //inputFiles.push_back("weights_CheckedEvts.out");

  //Store all information into a ROOT file:
  TFile* outputFile = new TFile((whichDir+"/OutFile_TopMassExample.root").c_str(),"RECREATE");
  outputFile->cd();

  //Considered values and corresponding XS-values
  const int NrConfigs = 6; 
  std::stringstream ssNrConfigs; ssNrConfigs << NrConfigs; std::string sNrConfigs = ssNrConfigs.str();
  double Var[NrConfigs]     = {170,      171,      172,      173,      174,      175     };
  double MGXSCut[NrConfigs] = {0.877635, 0.900965, 0.924294, 0.947624, 0.970953, 0.994283};

  double LnLik[NrConfigs] = {0.0};

  //vector< vector<TH1D> > indivLnLik;
  vector< vector< vector<double> > > indivLnLik;
  vector< std::string > sampleName; sampleName.clear();
  vector< double > normFactor;     normFactor.clear();
  vector< vector<double> > scaleFactor;

  //Store the TH1Ds in a map such that they can easily be looped over!
  std::map<string,TH1D*> histo1D;  
  std::map<string,TH2D*> histo2D;  
  
  for(int iWeightFile = 0; iWeightFile < inputFiles.size(); iWeightFile++){

    //Initialise histos
    bool fileInitialised = false;

    //Make sure that for each dataset the vector storing the individual histograms is empty
    vector< vector<double> > indivSampleLnLik; indivSampleLnLik.clear();
    vector<double> sampleSF; sampleSF.clear();

    //--- Read the likelihood values ! ---//
    std::ifstream ifs (inputFiles[iWeightFile].c_str(), std::ifstream::in); 
    std::string line;
    int evt,config,tf;
    double weight, weightUnc;
    double CosThetaCorr = 1;
    double MCScaleFactor = 1, NormFactor = 1;
    double SampleName;
    std::string sSampleName, sNormFactor;
    bool sampleNameSet = false;
    int SMConfig = 99, xOuterL = 99, xOuterR = 99;
    int consEvts = 0;
    double MaxLik = 0, MinLik = 9999;
    double MingR = 99;
    double ScdDer = 0;
    int xOuterLFine = 99, xOuterRFine = 99;

    vector<double> FitParametersFirstFit, FitParametersSecondFit;
    FitParametersFirstFit.clear(); FitParametersSecondFit.clear();
    bool doFits = false, allPointsInFit = false;       //Before changing this to true the issue of passing on the fit parameters should be resolved!!
  
    histo1D["SMLikValue"] = new TH1D("SMLikelihoodValue","Distribution of likelihood value at gR = 0.0",85,40,85);
    histo1D["OuterLikValue"] = new TH1D("OuterLikelihoodValue","Distribution of likelihood value at gR = 0.1",85,40,85);
    histo1D["SMLikValue_AfterCut"] = new TH1D("SMLikelihoodValue_AfterCut","Distribution of likelihood value at gR = 0.0 (after cut of ..)",85,40,85);
    histo1D["MCScaleFactor"] = new TH1D("MCScaleFactor","Scale factor for MC sample", 100,0,2);
    histo1D["Luminosity"] = new TH1D("Luminosity","Luminosity used", 100,15000,23000);
    histo1D["NormFactor"] = new TH1D("NormFactor","Norm factor for MC sample", 100,0,0.00001);

    while( std::getline(ifs,line)){
      std::istringstream iss(line);

      if( iss >> evt >> config >> tf >> weight >> weightUnc){ // >> MCScaleFactor ){
        if(whichDir.find("TTSemiLept") < whichDir.size())
          iss >> MCScaleFactor;

        //if(config == 1 && (consEvts+1) % 2000 == 0) std::cout << " Looking at event : " << consEvts+1 << flush<<"\r";
        stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

        //--- Initialize the event-per-event variables! ---//          
        if(config == 1){
          MaxLik = 0; MinLik = 9999; //Reset this variable at the beginning of each event!
          MingR = 99;
          SMConfig = 99; xOuterL = 99; xOuterR = 99;
	  ScdDer = 0; xOuterLFine = 99; xOuterRFine = 99;

          if(consEvts == 0 && fileInitialised == false){
            //Keep track of the number of considerd samples of the same type!
            fileInitialised = true;

            //Loop over all histo1Ds and add the sampleName!
            if(histo1D.size() > 0){
              for(std::map<std::string,TH1D*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
                TH1D *temp = it->second;
                temp->SetName((string(temp->GetName())+"_"+sampleName[iWeightFile]).c_str());
              }
            }
            //Loop over all histo2Ds and add the sampleName!
            if(histo2D.size() > 0){
              for(std::map<std::string,TH2D*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
                TH2D *temp = it->second;
                temp->SetName((string(temp->GetName())+"_"+sampleName[iWeightFile]).c_str());
              }
            }
          }
        }

        //Set the SMConfig:
        if(Var[config-1] == 173) SMConfig = config-1;

        //---  Fill the LnLik histograms for each event and for all events together  ---//
        LnLik[config-1] = (-log(weight)+log(MGXSCut[config-1]));

        //Get the maximum and minimum likelihood value
        if(LnLik[config-1] > MaxLik) MaxLik = LnLik[config-1];
        if(LnLik[config-1] < MinLik){ MinLik = LnLik[config-1]; MingR = Var[config-1];}

        //---  Only perform the fit after all configurations are considered!  ---//
        if(config == NrConfigs){

          //Plot the SM likelihood value:
          histo1D["SMLikValue"]->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
          histo1D["MCScaleFactor"]->Fill(MCScaleFactor);
          histo1D["Luminosity"]->Fill(Luminosity);
          histo1D["NormFactor"]->Fill(NormFactor);

          //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
          //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
          //if(LnLik[SMConfig] > LikCut || MaxLik - MinLik < 1.0){if(LikCut == 120) std::cout << " *****\n ERROR: Should not reject any event when cut-value = 120 !!! \n *****\n --> Value is : " << LnLik[SMConfig] << "\n " << std::endl; continue;}
          if(LnLik[SMConfig] > LikCut){if(consEvts == 0) std::cout << " Rejecting the first element ! " << std::endl; if(LikCut == 120) std::cout << " *****\n ERROR: Should not reject any event when cut-value = 120 !!! \n *****\n --> Value is : " << LnLik[SMConfig] << "\n " << std::endl; continue;}
          consEvts++;   //Count the number of full events!

          vector<double> indivEvtLnLik;
          indivEvtLnLik.clear();

          for(int i=0; i< NrConfigs; i++)
            indivEvtLnLik.push_back(LnLik[i]);
          indivSampleLnLik.push_back(indivEvtLnLik);
          sampleSF.push_back(MCScaleFactor);

          histo1D["SMLikValue_AfterCut"]->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
          
          //Save the first six LnLik curves
          if(consEvts % 25 == 0 && consEvts <= 4000){
            double lnLikArr[NrConfigs] = {0};
            for(int i = 0; i < NrConfigs; i++)
              lnLikArr[i] = LnLik[i];

            TGraph* gr_lnLik = new TGraph(NrConfigs, Var, lnLikArr);
            gr_lnLik->SetMarkerStyle(22);
            gr_lnLik->SetTitle("");
            gr_lnLik->GetXaxis()->SetTitle("top-quark mass (GeV)");
            gr_lnLik->GetYaxis()->SetTitle("-ln(L) value");

            //TF1* fit = new TF1("IndivFit","pol2",170,175);
            //gr_lnLik->Fit(fit,"Q","",170,175);
            
            TCanvas* canv = new TCanvas(("canv_"+NumberToString(consEvts)).c_str(),("canv_"+NumberToString(consEvts)).c_str());
            canv->cd();
            gr_lnLik->Draw("AP");
            canv->SaveAs((whichDir+"/IndivLikelihoods/EventLikelihood_"+whichDir+"_"+NumberToString(consEvts)+".pdf").c_str());
          }

        }
      }
      else{
        sampleNameSet = true;
        std::cout << "\n In ELSE -> Line = " << line << std::endl;
        sSampleName = line.substr(0, line.find(" "));
        sNormFactor = line.substr(sSampleName.length(), line.length());
        NormFactor = atof(sNormFactor.c_str());
        sampleName.push_back(sSampleName);
        normFactor.push_back(NormFactor);
        std::cout << "  * Recovered sampleName = " << sSampleName << " and norm-factor = " << NormFactor << std::endl;
      }
    }
    ifs.close();

    //Push back all the vectors containing event-per-event information!
    indivLnLik.push_back(indivSampleLnLik);
    scaleFactor.push_back(sampleSF);

    if(histo1D.size() > 0){
      TDirectory* th1dir = outputFile->GetDirectory(("1D_histo_"+sampleName[iWeightFile]).c_str());   //Check whether directory already exists ..
      if(!th1dir) th1dir = outputFile->mkdir(("1D_histo_"+sampleName[iWeightFile]).c_str());          // .. and otherwise create it!
      th1dir->cd();
      for(std::map<std::string,TH1D*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
        TH1D *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
      }
    }

    if(histo2D.size() > 0){
      TDirectory* th2dir = outputFile->GetDirectory(("2D_histo_"+sampleName[iWeightFile]).c_str());
      if(!th2dir) th2dir = outputFile->mkdir(("2D_histo_"+sampleName[iWeightFile]).c_str());
      th2dir->cd();
      for(std::map<std::string,TH2D*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
        TH2D *temp = it->second;
        temp->Write();
      }
    }
  }//End of looping over the different weight files
 
  outputFile->cd();
  //Get the minimum!!
  getMinimum(indivLnLik, scaleFactor, sampleName, normFactor, Luminosity, outputFile, Var, whichDir);

  outputFile->Close();
  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
  return 0;
}
