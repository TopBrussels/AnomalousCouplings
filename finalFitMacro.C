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
#include <cstdlib>
#include <map>

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

  TF1* polFit = new TF1(("polFit_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
  TGraph* gr_LnLik = new TGraph(NConfigs,Var, LogLikelihood);
  gr_LnLik->Fit(polFit,"Q","",polFit->GetXmin(), polFit->GetXmax());

  if(evtCounter == 1){
    for(int ipar = 0; ipar < polFit->GetNpar(); ipar++)
      FitParamsFirstFit.push_back(polFit->GetParameter(ipar));
  }
  else{
    for(int ipar = 0; ipar < polFit->GetNpar(); ipar++)
      FitParamsFirstFit[ipar] += polFit->GetParameter(ipar);
  }

  if(consAllPoints){
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
    }
    //Sort the fitdeviation values depending on the second value!
    std::sort(FitDeviation.begin(), FitDeviation.end(), sort_pred() );
    std::sort(FitDeviationRel.begin(), FitDeviationRel.end(), sort_pred() );

    //Now loop again over all configurations, plot the FitDeviation in sorted order and save the configNr's which should be excluded 
    std::vector<int> FitDevPointsToDel, FitDevRelPointsToDel;
    for(int itSortedConfig = NConfigs-1; itSortedConfig >= 0 ; itSortedConfig--){  //Looping from high to low values of the deviation!
      //Store the 'NrToDel' points which need to be excluded from the TGraph!
      if(FitDevPointsToDel.size() < NrToDel){ FitDevPointsToDel.push_back(FitDeviation[itSortedConfig].first); }//h_PointsDelByFitDev->Fill(Var[FitDeviation[itSortedConfig].first]);
      if(FitDevRelPointsToDel.size() < NrToDel){FitDevRelPointsToDel.push_back(FitDeviationRel[itSortedConfig].first); }//h_PointsDelByFitDevRel->Fill(Var[FitDeviationRel[itSortedConfig].first]);
    }

    //Create new arrays with the reduced information!
    double ReducedVar[NConfigs-NrToDel], ReducedLogLik[NConfigs-NrToDel];
    int iCounter = 0;
    for(int iConf = 0; iConf < NConfigs; iConf++){
      if(!(std::find(FitDevPointsToDel.begin(), FitDevPointsToDel.end(), iConf) != FitDevPointsToDel.end()) ){
        ReducedVar[iCounter] = Var[iConf];
        ReducedLogLik[iCounter] = LogLikelihood[iConf];
        iCounter++;
      }
    }

    //Define new TGraph and fit again
    TGraph* gr_ReducedLnLik = new TGraph(NConfigs-NrToDel, ReducedVar, ReducedLogLik);
    TF1* polFit_ReducedPoints = new TF1(("polFit_ReducedPoints_Evt"+EvtNumber).c_str(),"pol2",FitMin, FitMax); 
    gr_ReducedLnLik->Fit(polFit_ReducedPoints,"Q","",polFit->GetXmin(), polFit->GetXmax()); 

    //Add the different parameters of the individual fits in order to form a final fit!
    if(evtCounter == 1){
      for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
        FitParamsSecondFit.push_back(polFit_ReducedPoints->GetParameter(ipar));
    }
    else{
      for(int ipar = 0; ipar < polFit_ReducedPoints->GetNpar(); ipar++)
        FitParamsSecondFit[ipar] += polFit_ReducedPoints->GetParameter(ipar);
    }

    //polFit_ReducedPoints->Write();
    delete polFit_ReducedPoints;
    delete gr_ReducedLnLik;
  }

  delete polFit;
  delete gr_LnLik;
}

void getMinimum(vector< vector< vector<double> > > LnLikArray, vector< vector<double> > SF, vector<std::string> name, vector<double> norm, double Lumi, TFile* outFile, double Var[], double FitMin, double FitMax, bool getIndivSampleMin, std::string pseudoTitle = ""){

  if(pseudoTitle != "") Lumi = 1;
  std::cout << "\n Considered a total of : " << LnLikArray.size() << " samples" << std::endl;
  vector<double> summedEntries; summedEntries.clear();
  for(int iFile = 0; iFile < LnLikArray.size(); iFile++){
    std::cout << " Considering : " << name[iFile] << std::endl;
    std::cout << " Stored a total of : " << LnLikArray[iFile].size() << " arrays " << std::endl;
    std::cout << " Each event array has " << LnLikArray[iFile][0].size() << " entries " << std::endl;

    if(pseudoTitle != "") norm[iFile] = 1;
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

    double SampleEntries[100];
    for(int i = 0; i < LnLikArray[iFile][0].size(); i++)
      SampleEntries[i] = summedSampleEntries[i];

    if(getIndivSampleMin){
      TGraph* gr_sampleSum = new TGraph(LnLikArray[iFile][0].size(), Var, SampleEntries);
      TF1* polFit_sampleSum = new TF1(("polFit_SummedSampleGraph_"+name[iFile]).c_str(),"pol2",FitMin, FitMax);
      gr_sampleSum->Fit(polFit_sampleSum,"Q","",polFit_sampleSum->GetXmin(), polFit_sampleSum->GetXmax());
      std::cout << "\n Minimum for " << polFit_sampleSum->GetName() << " is : " << polFit_sampleSum->GetMinimumX() << " +- " << polFit_sampleSum->GetX(polFit_sampleSum->GetMinimum()+0.5, polFit_sampleSum->GetMinimumX(),0.2) - polFit_sampleSum->GetX(polFit_sampleSum->GetMinimum()+0.5, -0.2, polFit_sampleSum->GetMinimumX()) << endl;

      TDirectory* grdir = outFile->GetDirectory("graphs");   //Check whether directory already exists ..
      if(!grdir) grdir = outFile->mkdir("graphs");          // .. and otherwise create it!
      grdir->cd();
      gr_sampleSum->SetName(("FittedGraph_"+name[iFile]).c_str());
      gr_sampleSum->SetTitle(("Graph containing the summed events for "+name[iFile]).c_str());
      gr_sampleSum->Write();
    }
  }
  
  double Entries[100];
  for(int i = 0; i < LnLikArray[0][0].size(); i++)
    Entries[i] = summedEntries[i];

  TDirectory* grdir = outFile->GetDirectory("graphs");   //Check whether directory already exists ..
  if(!grdir) grdir = outFile->mkdir("graphs");          // .. and otherwise create it!
  grdir->cd();
  TGraph* gr_totalSum = new TGraph(LnLikArray[0][0].size(), Var, Entries);
  TF1* polFit_totalSum = new TF1(("polFit_SummedTotalGraph"+pseudoTitle).c_str(),"pol2",FitMin, FitMax);
  gr_totalSum->Fit(polFit_totalSum,"Q","",polFit_totalSum->GetXmin(), polFit_totalSum->GetXmax());
  std::cout << "\n Minimum for " << polFit_totalSum->GetName() << " is : " << polFit_totalSum->GetMinimumX() << " +- " << polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+0.5, polFit_totalSum->GetMinimumX(),0.2) - polFit_totalSum->GetX(polFit_totalSum->GetMinimum()+0.5, -0.2, polFit_totalSum->GetMinimumX()) << endl;

  gr_totalSum->SetName(("FittedGraph_AllSamples"+pseudoTitle).c_str());
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
  // 2-..) Weight files which will be considered

  double LikCut = 100;
  if( argc >= 2)
    LikCut = atoi(argv[1]);
  std::cout << " - Applied likelihood cut value is : " << LikCut << std::endl;

  vector<string> inputFiles;
  if( argc >= 3){
    for(int iFile = 2; iFile < argc; iFile++){
      inputFiles.push_back(string(argv[iFile]).c_str());
      std::cout << " - Stored file name is : " << inputFiles[inputFiles.size()-1] << std::endl;
    }
  }

  TFile* outputFile = new TFile(("Events_MC/OutFile_"+timestamp()+".root").c_str(),"RECREATE");   //So what if also Data is added?
  outputFile->cd();

  //Considered values and corresponding XS-values
  const int NrConfigs = 9; 
  std::stringstream ssNrConfigs; ssNrConfigs << NrConfigs; std::string sNrConfigs = ssNrConfigs.str();
  double Var[NrConfigs]     = {-0.2,     -0.15,   -0.1,    -0.05,   0.0,     0.05,    0.1,     0.15,    0.2    };
  double MGXSCut[NrConfigs] = {0.947244, 1.13624, 1.36448, 1.63952, 1.96892, 2.36027, 2.82111, 3.35903, 3.98157};

  double FitMin = -0.15, FitMax = 0.15;
  double LnLik[NrConfigs] = {0.0};

  //vector< vector<TH1F> > indivLnLik;
  vector< vector< vector<double> > > indivLnLik;
  vector< std::string > sampleName; sampleName.clear();
  vector< double > normFactor;     normFactor.clear();
  vector< vector<double> > scaleFactor;

  //Store the TH1Fs in a map such that they can easily be looped over!
  std::map<string,TH1F*> histo1D;  

  //** Loop over the different input files and read out the necessary info **//
  double Luminosity = 19646.8;
  for(int iWeightFile = 0; iWeightFile < inputFiles.size(); iWeightFile++){

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
    int SMConfig = 99, consEvts = 0;
    double MaxLik = 0, MinLik = 9999;

    vector<double> FitParametersFirstFit, FitParametersSecondFit;
    FitParametersFirstFit.clear(); FitParametersSecondFit.clear();
    bool doFits = false, allPointsInFit = false;       //Before changing this to true the issue of passing on the fit parameters should be resolved!!
  
    int NrBins = 8; 
    double halfBinWidth = (Var[NrConfigs-1]- Var[0])/((double) NrBins*2.0);
    float xLow = Var[0] - halfBinWidth, xHigh = Var[NrConfigs-1] + halfBinWidth; 
    histo1D["SMLikValue"] = new TH1F("SMLikelihoodValue","Distribution of likelihood value at gR = 0.0",500,30,90);
    histo1D["SMLikValue_AfterCut"] = new TH1F("SMLikelihoodValue_AfterCut","Distribution of likelihood value at gR = 0.0 (after cut of ..)",500,30,90);
    histo1D["MCScaleFactor"] = new TH1F("MCScaleFactor","Scale factor for MC sample", 100,0,2);
    histo1D["Luminosity"] = new TH1F("Luminosity","Luminosity used", 100,15000,23000);
    histo1D["NormFactor"] = new TH1F("NormFactor","Norm factor for MC sample", 100,0,0.00001);

    while( std::getline(ifs,line)){
      std::istringstream iss(line);

      if( iss >> evt >> config >> tf >> weight >> weightUnc >> MCScaleFactor ){

        if(config == 1 && (consEvts+1) % 2000 == 0) std::cout << " Looking at event : " << consEvts+1 << flush<<"\r";
        stringstream ssEvt; ssEvt << evt; string sEvt = ssEvt.str();

        //--- Initialize the event-per-event variables! ---//          
        if(config == 1){
          MaxLik = 0; MinLik = 9999; //Reset this variable at the beginning of each event!
          SMConfig = 99;

          if(consEvts == 0){
            //Loop over all histo1Ds and add the sampleName!
            if(histo1D.size() > 0){
              for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
                TH1F *temp = it->second;
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

          //Plot the SM likelihood value:
          histo1D["SMLikValue"]->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);
          histo1D["MCScaleFactor"]->Fill(MCScaleFactor);
          histo1D["Luminosity"]->Fill(Luminosity);
          histo1D["NormFactor"]->Fill(NormFactor);

          //Need to make sure event is rejected for all normalisations otherwise counter is wrong, therefore nrNorms-1 is uses which corresponds to acceptance norm!
          //  --> But then last bin is not yet filled for all norms (so using average will be difficult ...)
          if( LnLik[SMConfig] > LikCut){ if(LikCut == 100) std::cout << " ******** \n ERROR: Should not reject any event when cut-value = 100 !!! \n ******** \n " << std::endl; continue; }
          consEvts++;   //Count the number of full events!

          vector<double> indivEvtLnLik;
          indivEvtLnLik.clear();

          for(int i=0; i< NrConfigs; i++)
            indivEvtLnLik.push_back(LnLik[i]);
          indivSampleLnLik.push_back(indivEvtLnLik);
          sampleSF.push_back(MCScaleFactor);

          histo1D["SMLikValue_AfterCut"]->Fill(LnLik[SMConfig], MCScaleFactor*Luminosity*NormFactor);

          //-- Send the array containing the ln(weights) to the predefined function to define the TGraph, fit this, detect the deviation points and fit again! --//
          if(doFits){
            calculateFit(LnLik, sEvt, consEvts, allPointsInFit, NrConfigs, FitParametersFirstFit, FitParametersSecondFit, Var, FitMin, FitMax);
          }
        }
      }
      else{
        sampleNameSet = true;
        std::cout << " In ELSE -> Line = " << line << std::endl;
        sSampleName = line.substr(0, line.find(" "));
        sNormFactor = line.substr(sSampleName.length(), line.length());
        NormFactor = atof(sNormFactor.c_str());
        std::cout << " * What is sample name : " << sSampleName << std::endl;
        sampleName.push_back(sSampleName);
        std::cout << " * What is norm-factor : " << NormFactor << std::endl;    //To notice: 1/normFactor = lumi of sample!
        normFactor.push_back(NormFactor);
      }
    }
    ifs.close();

    //Push back all the vectors containing event-per-event information!
    indivLnLik.push_back(indivSampleLnLik);
    scaleFactor.push_back(sampleSF);

    std::cout << " Studied a total of " << consEvts << " events !" << std::endl;

    if(histo1D.size() > 0){
      TDirectory* th1dir = outputFile->GetDirectory(("1D_histo_"+sampleName[iWeightFile]).c_str());   //Check whether directory already exists ..
      if(!th1dir) th1dir = outputFile->mkdir(("1D_histo_"+sampleName[iWeightFile]).c_str());          // .. and otherwise create it!
      th1dir->cd();
      for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
      }
    }

  }//End of looping over the different weight files

  int nrRandomSamples = 10;
  vector< vector< vector<double> > > pseudoDataLnLik;
  vector< vector<double> > scaleFactorPD; 
  for(int iSample = 0; iSample < indivLnLik.size(); iSample++){

    vector< vector<double> > pseudoSampleLnLik; pseudoSampleLnLik.clear();
    vector<double> pseudoSF;
    for(int iEvt = 0; iEvt < indivLnLik[iSample].size(); iEvt++){
      double random = ((double) rand() / (RAND_MAX));
      if(random <= Luminosity*normFactor[iSample]){
        pseudoSampleLnLik.push_back(indivLnLik[iSample][iEvt]);
        pseudoSF.push_back(scaleFactor[iSample][iEvt]);
      }
    }
    pseudoDataLnLik.push_back(pseudoSampleLnLik);
    scaleFactorPD.push_back(pseudoSF);
    std::cout << " Stored for dataset " << sampleName[iSample] << " a total of " << pseudoDataLnLik[iSample].size() << std::endl;
  }
  getMinimum(pseudoDataLnLik, scaleFactorPD, sampleName, normFactor, Luminosity, outputFile, Var, FitMin, FitMax, false, "_PseudoSample");

  bool getMinForIndivSamples = false;
  getMinimum(indivLnLik, scaleFactor, sampleName, normFactor, Luminosity, outputFile, Var, FitMin, FitMax, getMinForIndivSamples);

  outputFile->Close();
  cout << "\n It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << "s to run the program \n" << endl;
  return 0;
  }
