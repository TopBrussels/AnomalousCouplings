//#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
#include "TTree.h"
#include "TBranch.h"
#include "TNtuple.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TDirectory.h"
#include "TClonesArray.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include <sys/stat.h>  //Needed for mkdir option

#include "TLorentzVector.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!

#include "AnomalousCouplings/PersonalClasses/interface/AnomCoupLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/BTagStudy.h"
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"
#include "AnomalousCouplings/PersonalClasses/interface/ExtraEvtSelCuts.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

using namespace std;
using namespace reweight;   //Need this in order to use LumiReWeighting!

//Convert integer to string!
template <typename T>
string NumberToString ( T Number ){
    stringstream ss;
    ss << Number;
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

//Check if a file exists or not!
bool file_exist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
}

int main (int argc, char *argv[])
{
  clock_t start = clock();
  cout << "************************************************************" << endl;
  cout << "  Beginning of the program for analyzing the Light Trees !  " << endl;
  cout << "************************************************************" << endl;

  //Get the timestamp to add this in the output root file!
  std::string time = timestamp();

  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  //--------------------------------//
  //  Get NrEvts from command line  //
  //--------------------------------//
  int nrEvtsComLine = -1;
  if(argc >= 2 && string(argv[1]) != "-1")
    nrEvtsComLine = atoi(argv[1]);

  //------------------------------//
  //  Get the systematic from CL  //
  //------------------------------//
  //1) Systematic applied in nTuple!
  std::string nTupleSyst = "Nominal";
  if(argc >= 3) nTupleSyst = string(argv[2]); 
  std::cout << " - Will be using nTupleSyst " << nTupleSyst << " ! " << std::endl;

  //2) Systematic applied locally
  std::string systematic = "Nominal";
  if(argc >= 4) systematic = string(argv[3]);
  std::cout << " - Will be using systematic " << systematic << " ! " << std::endl;
  
  if(systematic != "Nominal" && nTupleSyst != "Nominal"){ std::cout << " ERROR: Cannot apply a local systematic on a nTuple systematic ....!! " << std::endl; exit(-1);}

  //------------------------//
  //  Verbosity for output  //
  //  -->Used in bTagStudy  //
  //------------------------//
  int verbosity                 = 1;
  //1 Main Info
  //2 mcParticlesMatchin Info
  //3 
  //4 Info for each event
  //5 Debug

  //---------------------------//
  //  Run specific parts only  //
  //---------------------------//
  bool getLHCOOutput = true;
  bool getEvtNrInfo = false;
  bool splitLeptonChargeLHCO = false;
  bool getCorrectAndWrongLHCO = true; 
  bool savePDF = false;
  bool bTagChoiceMade = true;  
  bool getMassFits = false;
  bool bTagPlotsMade = true;         //Decide whether the bTag distributions for the efficiencies should be created or read from file
  bool createTFTree = false;         //Decide whether the TFTree should be made or just skipped during this run! 
  bool fillTFAfterCuts = false;      //Decide whether the TFTree is filled before or after the additional event selection cuts (b-tag, chi-sq and Mw-Mt)
  bool onlyTTbar = false;
  bool onlyMuonChannel = true;       //Set this to true in order to reject all electron channel events!
  bool runLocally = false;           //Differentiate between local and m-machines running for the location of the input samples

  //In case the matching or scaling systematics are considered, only the ttbar sample should be considered
  if(nTupleSyst == "matchingdown_Nominal" || nTupleSyst == "matchingup_Nominal" || nTupleSyst == "scaledown_Nominal" || nTupleSyst == "scaleup_Nominal") onlyTTbar = true;
  if(systematic != "Nominal"){ getLHCOOutput = false; getEvtNrInfo = true;}

  //Give correct name to TF light ntuples!
  std::string tfFill = "";
  if(fillTFAfterCuts) tfFill = "_AfterExtraCuts";

  //Mention in name whether only ttbar sample is considered
  std::string onlyTTbarTitle = "";
  if(onlyTTbar) onlyTTbarTitle = "_OnlyTTbar";

  //Safety check: Can only create the mass plots if onlyTTbar samples are considered!
  if(!onlyTTbar && getMassFits){std::cout << " Will not perform the mass fits since all samples are considered!" << std::endl; getMassFits = false;}
  if(getMassFits) onlyTTbarTitle = "_MadeMassFits";

  //In case the bTagPlots still need to be made, only ttbarSemiLept sample should be considered!
  if(!bTagPlotsMade) onlyTTbar = true;
  if(!bTagPlotsMade) onlyTTbarTitle = "_MadeBTagPlots";

  //Set the correct bTagPlots ROOT file!
  std::string bTagPlotsOutput = "PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_WithTTbar_OnlyTTbarSemiLept_12Bins.root";
  //std::string bTagPlotsOutput = "PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_NoTTbar_OnlyTTbarSemiLept_12Bins.root";
  //std::string bTagPlotsOutput = "PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_NoTTbar_AllDataSets_12Bins.root";

  int stop = 0;
  if( file_exist(bTagPlotsOutput.c_str()) ){
    if(!bTagPlotsMade){
      std::cout << " \n According to the bTagPlotsMade boolean the plots are not yet made but the output file has been found ... " << std::endl;
      std::cout << "  --> Maybe the output file needs to be updated: " << bTagPlotsOutput << std::endl;
      std::cout << "  ==> Will stop the script until this issue is fixed ! \n" << std::endl;
      stop = 1;
    }
    else
      std::cout << " - Will apply the b-tag efficiencies from file : " << bTagPlotsOutput << std::endl;
  }
  else{
    if(bTagPlotsMade){
      std::cout << "\n According to the bTagPlotsMade boolean the plots should be created but the output file cannot be found ... " << std::endl;
      std::cout << " --> Maybe the wrong output file has been given : " << bTagPlotsOutput << std::endl;
      std::cout << " ==> Will stop the script until this issue is fixed ! \n" << std::endl;
      stop = 1;
    }
    else
      std::cout << " - Will create the b-tag efficiencies and store them in : " << bTagPlotsOutput << std::endl;
  }
  if(stop == 1) return 0;

  //-- Specific b-tag stuff!
  int ChosenBTag = 5;  //-->Corresponds to two T b-tags and no light veto!
  int NrBTags = 9;
  string bTitle[9] = {"LooseTags","LooseTagsLVeto","MediumTags","MediumTagsMVeto","MediumTagsLVeto","TightTags","TightTagsTVeto","TightTagsMVeto","TightTagsLVeto"};
  string whichBTags = "AllBTags";
  if(bTagChoiceMade){ NrBTags = 1; bTitle[0] = bTitle[ChosenBTag]; whichBTags = bTitle[ChosenBTag];}

  //-------------------------//
  //   Which systematics ?   //
  //-------------------------//
  std::string doLeptonSFShift = "Nominal";   //Other options are "Minus" and "Plus" (Capital letters are needed!)
  std::string doLumiWeightShift = "Nominal"; //Other options are "Minus" and "Plus"
  if(systematic != "Nominal" && systematic != "bTagMinus" && systematic != "bTagPlus" && systematic != "misTagMinus" && systematic != "misTagPlus" && systematic != "MuonSFMinus" && systematic != "MuonSFPlus" && systematic != "ElecSFMinus" && systematic != "ElecSFPlus"){
    std::cout << "   *****  Given systematic is not allowed! " << std::endl;
    std::cout << "   *****  Possibilities are : Nominal, bTagMinus, bTagPlus, MuonSFMinus, MuonSFPlus, ElecSFMinus and ElecSFPlus " << std::endl;
    exit(-1);
  }

  //b-tag systematics
  int syst_btag = 0; //+1 for SF up, -1 for SF down
  if(systematic == "bTagMinus")     syst_btag = -1;
  else if(systematic == "bTagPlus") syst_btag = 1;
  int syst_mistag = 0; //+1 for SF up, -1 for SF down
  if(systematic == "misTagMinus")     syst_mistag = -1;
  else if(systematic == "misTagPlus") syst_mistag = 1;
        
  //lepton SF systematics
  string syst_muonSF = "Nominal";
  if(systematic == "MuonSFMinus")     syst_muonSF = "Minus";
  else if(systematic == "MuonSFPlus") syst_muonSF = "Plus";
  string syst_elecSF = "Nominal";
  if(systematic == "ElecSFMinus")     syst_elecSF = "Minus";
  else if(systematic == "ElecSFPlus") syst_elecSF = "Plus";

  //----------------------------//
  // Input & output information //
  //----------------------------//
  //ROOT file for storing plots
  string pathPNG = "PlotsMacro_Light";
  TFile* outputFile = new TFile((pathPNG+"/AnomCoup_Analysis_"+nTupleSyst+""+onlyTTbarTitle+"_"+systematic+"_"+whichBTags+"_"+time+".root").c_str(),"RECREATE");

  //Which datasets should be considered
  vector<string> inputFiles;
  vector<Dataset*> datasets;
  std::string inputFileDir = "/user/aolbrech/PBS_ScriptRunning/";
  if(onlyTTbar){
    if(nTupleSyst == "Nominal" && systematic == "Nominal"){
      if(!runLocally) inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_TTbarJets_SemiLept_"+nTupleSyst+".root").c_str());
      else            inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_SemiLept_Nominal.root");
    }
    else{    //This corresponds to matching & scaling which is all TTbar (not only SemiLept)!
      if(!runLocally) inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_TTbarJets_"+nTupleSyst+".root").c_str());
      else            inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_Nominal.root");
    }
  }
  else{ 
    if(!runLocally){
      if(nTupleSyst == "Nominal" && systematic == "Nominal"){
        inputFiles.push_back((inputFileDir+"DATA/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root").c_str());                                                  //Winter14_v8
        //inputFiles.push_back("/user/aolbrech/PBS_ScriptRunning/Results/RESULTS_AnomCoup_13112015_121239/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root");  //Winter14_v5
      }
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_ZJets_1jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_ZJets_2jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_ZJets_3jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_ZJets_4jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_WJets_1jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_WJets_2jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_WJets_3jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_WJets_4jets_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_sChannel_t_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_sChannel_tbar_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_tWChannel_t_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_tWChannel_tbar_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_tChannel_t_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_SingleTop_tChannel_tbar_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_TTbarJets_SemiLept_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_TTbarJets_FullLept_"+nTupleSyst+".root").c_str());
      inputFiles.push_back((inputFileDir+""+nTupleSyst+"/AnomCoupLight_TTbarJets_FullHadr_"+nTupleSyst+".root").c_str());
    }
    else{
      inputFiles.push_back("LightTree/AnomCoupLight_Data_Mu_Merged_22Jan2013_Nominal.root");    //Winter14_v8
      inputFiles.push_back("LightTree/AnomCoupLight_ZJets_1jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_ZJets_2jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_ZJets_3jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_ZJets_4jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_WJets_1jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_WJets_2jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_WJets_3jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_WJets_4jets_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_sChannel_t_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_sChannel_tbar_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_tWChannel_t_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_tWChannel_tbar_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_tChannel_t_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_SingleTop_tChannel_tbar_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_SemiLept_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_FullLept_Nominal.root");
      inputFiles.push_back("LightTree/AnomCoupLight_TTbarJets_FullHadr_Nominal.root");
    }
  }
  if(verbosity > 0) std::cout << " - All ROOT files loaded " << std::endl;
	
  //-------------------------//
  // Set analysis luminosity //
  //-------------------------//
  float Luminosity = 19646.840; //IsoMu24_eta2p1 trigger for semiMu case!  
  if(verbosity > 0) std::cout << " - Analysis luminosity is set to : " << Luminosity << std::endl;

  //-------------------//
  // Load all datasets //
  //-------------------//
  for(unsigned int iDataSet=0; iDataSet<inputFiles.size(); iDataSet++){
    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");

    TTree* inConfigTree = (TTree*) inputFile->Get("configTreeLightFile");
    TBranch* d_br = (TBranch*) inConfigTree->GetBranch("Dataset");
    TClonesArray* tc_dataset = new TClonesArray("Dataset",0);
    d_br->SetAddress(&tc_dataset);
    inConfigTree->GetEvent(0);
    Dataset* dataSet = (Dataset*) tc_dataset->At(0);
    //XSection = dataSet->Xsection();

    int color = 0;
    if(dataSet->Name().find("QCD") == 0)             color = kYellow;
    if(dataSet->Name().find("TT") == 0){             color = kRed+1; dataSet->SetTitle("t#bar{t}");}        //Giving the different ttbar samples the same title will merge them in the MSPlot!
    if(dataSet->Name().find("TTbarJets_Other") == 0) color = kRed-7;
    if(dataSet->Name().find("ST") == 0 || 
       dataSet->Name().find("SingleTop") == 0){      color = kMagenta; dataSet->SetTitle("Single top"); }
    if(dataSet->Name().find("WJets") == 0){          color = kGreen-3; dataSet->SetTitle("W#rightarrowl#nu");}      
    if(dataSet->Name().find("ZJets") == 0){          color = kAzure-2; dataSet->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");}
    if(dataSet->Name().find("Data_Mu") == 0){        color = kBlack; dataSet->SetTitle("Data"); dataSet->SetName("Data_Mu_Merged"); dataSet->SetEquivalentLuminosity(Luminosity);}

    Dataset* tmpDS = new Dataset(dataSet->Name(), dataSet->Title(), dataSet->DoIt(), color, dataSet->LineStyle(), dataSet->LineWidth(), dataSet->NormFactor(), dataSet->Xsection());
    tmpDS->SetEquivalentLuminosity( dataSet->EquivalentLumi() );
    datasets.push_back( tmpDS );
  }
  if(verbosity > 0) std::cout << " - Datasets initialized " << std::endl;

  //--------------------------//
  //  Initialize MSPlots      //
  //  --> Needs dataset info  //
  //--------------------------//
  map<string,MultiSamplePlot*> MSPlot;
  string leptFlavs[2]={"_mu","_el"};

  //-- Histograms which should be made separately for the two possible decay channels
  for(int ii = 0; ii < 2; ii++){
    string leptFlav = leptFlavs[ii];
    if(onlyMuonChannel && leptFlav == "_el") continue;

    //-- Histograms which should contain the considered btag name!
    for(int ibTag = 0; ibTag < NrBTags; ibTag++){
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# light jets");
      MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# selected jets");
      MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# b-tagged jets");
      MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav] = new MultiSamplePlot(datasets, "nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+ leptFlav,10, -0.5, 9.5, "# light jets");
    }
      
    MSPlot["HadrTopMass"+leptFlav] = new MultiSamplePlot(datasets,"HadrTopMass"+leptFlav, 175, 30, 320,"Top mass (GeV)");
    MSPlot["HadrTopMass_AllCuts"+leptFlav] = new MultiSamplePlot(datasets,"HadrTopMass_AllCuts"+leptFlav, 175, 30, 320,"Top mass (GeV)");
    MSPlot["HadrWMass"+leptFlav] = new MultiSamplePlot(datasets,"HadrWMass"+leptFlav, 90, 10, 180,"W mass (GeV)");
    MSPlot["HadrWMass_AllCuts"+leptFlav] = new MultiSamplePlot(datasets,"HadrWMass_AllCuts"+leptFlav, 90, 10, 180,"W mass (GeV)");
    MSPlot["MlbMqqbChiSq"+leptFlav] = new MultiSamplePlot(datasets,"MlbMqqbChiSq"+leptFlav, 50, 0, 35,"Mlb-Mqqb chi-sq");
    MSPlot["MlbMqqbChiSq_OnlyChiSq"+leptFlav] = new MultiSamplePlot(datasets,"MlbMqqbChiSq_OnlyChiSq"+leptFlav, 50, 0, 35,"Mlb-Mqqb chi-sq");
    MSPlot["MlbMqqbChiSq_AllCuts"+leptFlav] = new MultiSamplePlot(datasets,"MlbMqqbChiSq_AllCuts"+leptFlav, 50, 0, 35,"Mlb-Mqqb chi-sq");

    MSPlot["nSelectedJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nSelectedJets_BeforeBTag"+leptFlav,10, -0.5, 9.5, "# selected jets");
    MSPlot["partFlav_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets,"partFlav_BeforeBTag"+leptFlav, 5, -0.5, 4.5, "Parton flavour of selected jets");

    MSPlot["nBJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nBJets_BeforeBTag"+leptFlav, 10, -0.5, 9.5, "# b-jets");
    MSPlot["nBJets_BTag"+leptFlav]       = new MultiSamplePlot(datasets, "nBJets_BTag"+leptFlav,       10, -0.5, 9.5, "# b-jets");
    MSPlot["nBJets_BTagSF"+leptFlav]     = new MultiSamplePlot(datasets, "nBJets_BTagSF"+leptFlav,     10, -0.5, 9.5, "# b-jets");
    MSPlot["nBJets_AllCuts"+leptFlav]    = new MultiSamplePlot(datasets, "nBJets_AllCuts"+leptFlav,    10, -0.5, 9.5, "# b-jets");
    MSPlot["nLightJets_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nLightJets_BeforeBTag"+leptFlav, 10, -0.5, 9.5, "# light jets");
    MSPlot["nLightJets_BTag"+leptFlav]       = new MultiSamplePlot(datasets, "nLightJets_BTag"+leptFlav,       10, -0.5, 9.5, "# light jets");
    MSPlot["nLightJets_BTagSF"+leptFlav]     = new MultiSamplePlot(datasets, "nLightJets_BTagSF"+leptFlav,     10, -0.5, 9.5, "# light jets");
    MSPlot["nLightJets_AllCuts"+leptFlav]    = new MultiSamplePlot(datasets, "nLightJets_AllCuts"+leptFlav,    10, -0.5, 9.5, "# light jets");

    MSPlot["JetsPt_partFlavGluon_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavLight_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavLight_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with light parton flavour");
    MSPlot["JetsPt_partFlavB_BeforeBTag"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavB_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with b parton flavour");
    MSPlot["JetsPt_partFlavC_BeforeBTag"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavC_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with c parton flavour");
    MSPlot["JetsPt_partFlavUndef_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetsPt_partFlavUndef_BeforeBTag"+leptFlav, 120, 0, 350, "Jet p_T for all jets with undefined parton flavour");

    MSPlot["JetsPt_partFlavGluon_BTag"+leptFlav]       = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BTag"+leptFlav,       120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavGluon_BTagSF"+leptFlav]     = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_BTagSF"+leptFlav,     120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["JetsPt_partFlavGluon_AllCuts"+leptFlav]    = new MultiSamplePlot(datasets, "JetsPt_partFlavGluon_AllCuts"+leptFlav,    120, 0, 350, "Jet p_T for all jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_BTag"+leptFlav]    = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_BTag"+leptFlav,    120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_BTagSF"+leptFlav]  = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_BTagSF"+leptFlav,  120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");
    MSPlot["BJetPt_partFlavGluon_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "BJetPt_partFlavGluon_AllCuts"+leptFlav, 120, 0, 350, "Jet p_T for b-tagged jets with gluon parton flavour");

    MSPlot["LeptonPt_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforePU"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_PU_NoLeptonSF"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BeforeBTag"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_BTagSF"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_AllCuts"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonPt_AllCuts_noMuSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonPt_AllCuts_noMuSF"+leptFlav, 60, 20, 150,"Lepton p_{T} (GeV)");
    MSPlot["LeptonMass_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonMass_AllCuts"+leptFlav, 100, 0, 2,"Lepton Mass (GeV)");
    MSPlot["LeptonEta"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonEta_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta_BTagSF"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonEta_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonEta_AllCuts"+leptFlav, 50, -2.6, 2.6,"Lepton #eta");
    MSPlot["LeptonCharge"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge"+leptFlav, 10, -2, 2,"Lepton charge");
    MSPlot["LeptonCharge_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge_BTagSF"+leptFlav, 10, -2, 2,"Lepton charge");
    MSPlot["LeptonCharge_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "LeptonCharge_AllCuts"+leptFlav, 10, -2, 2,"Lepton charge");

    MSPlot["JetPt_LeadingJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforePU"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_PU_NoLeptonSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BeforeBTag"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_BTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_AllCuts"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_LeadingJet_AllCuts_noBTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_LeadingJet_AllCuts_noBTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");

    MSPlot["JetPt_SecondJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet_BeforePU"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_SecondJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet_PU_NoLeptonSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_SecondJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet_BeforeBTag"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_SecondJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_SecondJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet_AllCuts"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_SecondJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_SecondJet_BTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");

    MSPlot["JetPt_ThirdJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet_BeforePU"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_ThirdJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet_PU_NoLeptonSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_ThirdJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet_BeforeBTag"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_ThirdJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_ThirdJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet_AllCuts"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_ThirdJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_ThirdJet_BTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");

    MSPlot["JetPt_FourthJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet_BeforePU"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FourthJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet_PU_NoLeptonSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FourthJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet_BeforeBTag"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FourthJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FourthJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet_AllCuts"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FourthJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FourthJet_BTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");

    MSPlot["JetPt_FifthJet_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet_BeforePU"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FifthJet_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet_PU_NoLeptonSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FifthJet_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet_BeforeBTag"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FifthJet"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FifthJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet_AllCuts"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");
    MSPlot["JetPt_FifthJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetPt_FifthJet_BTagSF"+leptFlav, 85, 20, 270,"Jet p_{T} (GeV)");

    MSPlot["JetMass_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetMass_LeadingJet_AllCuts"+leptFlav, 75, 0, 50,"Jet p_{T} (GeV)");
    MSPlot["JetEta_LeadingJet"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet"+leptFlav, 50, -2.6, 2.6,"Jet #eta");
    MSPlot["JetEta_LeadingJet_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet_BTagSF"+leptFlav, 50, -2.6, 2.6,"Jet #eta");
    MSPlot["JetEta_LeadingJet_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "JetEta_LeadingJet_AllCuts"+leptFlav, 50, -2.6, 2.6,"Jet #eta");

    MSPlot["nPV_BeforePU"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforePU"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_PU_NoLeptonSF"+leptFlav] = new MultiSamplePlot(datasets, "nPV_PU_NoLeptonSF"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BeforeBTag"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV"+leptFlav] = new MultiSamplePlot(datasets, "nPV"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_BTagSF"+leptFlav] = new MultiSamplePlot(datasets, "nPV_BTagSF"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_AllCuts"+leptFlav] = new MultiSamplePlot(datasets, "nPV_AllCuts"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_AllCuts_noPU"+leptFlav] = new MultiSamplePlot(datasets, "nPV_AllCuts_noPU"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");
    MSPlot["nPV_AllCuts_noLumi"+leptFlav] = new MultiSamplePlot(datasets, "nPV_AllCuts_noLumi"+leptFlav, 36, 0.5, 35.5,"Number of primary vertices");

    MSPlot["CSVDistr_BeforeBTag"+leptFlav] = new MultiSamplePlot(datasets, "CSVDistribution_BeforeBTag"+leptFlav, 100, -1, 2, "CSV discriminant");
    MSPlot["CSVDistr_BJets"+leptFlav] = new MultiSamplePlot(datasets, "CSVDistr_BJets"+leptFlav,50,-1,2,"CSV discriminant");
    MSPlot["CSVDistr_LightJets"+leptFlav] = new MultiSamplePlot(datasets, "CSVDistr_LightJets"+leptFlav,50,-1,2,"CSV discriminant");
  }

  //-----------------------//
  // Initialize histograms //
  //-----------------------//
  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;  

  //--------------------------------//
  // Lumi reweighting and lepton SF //
  //--------------------------------//
  LumiReWeighting LumiWeights;
  LumiReWeighting LumiWeightsUp;
  LumiReWeighting LumiWeightsDown;
  LumiWeights = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/nominal.root", "pileup","pileup");
  LumiWeightsUp = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_up.root","pileup","pileup");
  LumiWeightsDown = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_Mu/sys_down.root", "pileup", "pileup");
  cout << " - LumiReWeighting instantiated ... " << endl;
  
  // initialize lepton SF (ROOT files taken from: IsoMu24_eta2p1)
  LeptonTools* leptonTools = new LeptonTools(false);
  //leptonTools->readMuonSF("PersonalClasses/Calibrations/LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "PersonalClasses/Calibrations/LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root"); //Order: ID/Iso/Trig files!
  leptonTools->readMuonSF("PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_Run2012ReReco_53X.root","PersonalClasses/Calibrations/LeptonSF/MuonEfficiencies_ISO_Run_2012ReReco_53X.root","PersonalClasses/Calibrations/LeptonSF/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root");
  leptonTools->readElectronSF();

  //------------------------//
  // B-tagging scale factor //
  //------------------------//
  //Method 1a of https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods will be used!

  //BTagWeightTools *bTagTool = new BTagWeightTools("PersonalClasses/Calibrations/BTagSF/SFb-pt_NOttbar_payload_EPS13.txt","CSVT");   //Standard EPS2013 is used!
  BTagWeightTools *bTagTool = new BTagWeightTools("PersonalClasses/Calibrations/BTagSF/SFb-pt_WITHttbar_payload_EPS13.txt","CSVT");   //Standard EPS2013 is used!

  //During a first run the plots need to be created!
  if(!bTagPlotsMade){
    bTagTool->InitializeMCEfficiencyHistos(12,30.,360.,2);  //How to get these histo-values?
  }
  else{
    bTagTool->ReadMCEfficiencyHistos(bTagPlotsOutput); //"PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_noTTbar.root");
  }

  //-----------------------//
  // Load personal classes //
  //-----------------------//
  float Mlb  = 107.7945, S_Mlb  = 32.4255;
  float Mqqb = 175.0311, S_Mqqb = 17.0589;
  float MW   = 83.6161,  S_MW   = 10.2171;
  ExtraEvtSelCuts extraEvtSelCuts(Mqqb, S_Mqqb, MW, S_MW, bTagChoiceMade, 3, 2);
  BTagStudy bTagStudy(verbosity, datasets, bTagChoiceMade, ChosenBTag, Mlb, S_Mlb, Mqqb, S_Mqqb);

  //--------------------------------------//
  // Loop on datasets for actual analysis //
  //--------------------------------------//
  for(unsigned int iDataSet = 0; iDataSet < inputFiles.size(); iDataSet++){

    TFile* inputFile = new TFile(inputFiles[iDataSet].c_str(),"READ");  

    TTree* inLightTree = (TTree*) inputFile->Get("LightTree");
    TBranch* light_br = (TBranch*) inLightTree->GetBranch("TheAnomCoupLight");
    AnomCoupLight* light = 0;
    light_br->SetAddress(&light);

    //Number of events that will be used in the "loop on events"
    int nrEvts = 0;
    if(nrEvtsComLine == -1) nrEvts = inLightTree->GetEntries();
    else                    nrEvts = nrEvtsComLine;
  
    Dataset* dataSet = datasets[iDataSet];//(Dataset*) tc_dataset->At(0);
    string dsName = dataSet->Name();
    string dsTitle = dataSet->Title();
    if(verbosity > 0) std::cout << "\n   *** Looking at dataset "<< dsName << " (" << iDataSet+1 << "/" << inputFiles.size() << ") with " << nrEvts << " events (total = " << inLightTree->GetEntries() << ") !" << std::endl;

    //Set some booleans in order to know whether the considered sample is Data or TTbarSemiLept!
    bool isData = false, isTTbarSemiLept = false;
    if( dsName.find("Data") == 0 || dsName.find("data") == 0 || dsName.find("DATA") == 0 )
      isData = true;
    if( dsName.find("TTbarJets_SemiLept") == 0 )
      isTTbarSemiLept = true;

    dsName += "_"+nTupleSyst;

    //**** Create the 1D and 2D histograms for each dataset specific! ****//
    histo1D["BTagWeight_"+dsName] = new TH1F(("BTagWeight_"+dsName).c_str(),("b-tag reweighting distribution for "+dsTitle).c_str(),200,0,2);
    histo1D["lumiWeight_"+dsName] = new TH1F(("lumiWeight_"+dsName).c_str(),("lumi reweighting distribution for "+dsTitle).c_str(),200,0,2);
    histo1D["leptonSF_"+dsName] = new TH1F(("leptonSF_"+dsName).c_str(), ("lepton scale factor for "+dsTitle).c_str(),200,0,2);

    histo2D["bTagWeight_vs_BJetPt_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_"+dsName).c_str(),("btag weight versus b-jet PT for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavB_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavB_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = b-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavC_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavC_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = c-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_BJetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_BJetPt_partFlavLight_"+dsName).c_str(),("btag weight versus b-jet PT (partFlav = light jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_"+dsName).c_str(),("btag weight versus light jet PT for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavB_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavB_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = b-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavC_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavC_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = c-jet) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_LightJetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_LightJetPt_partFlavLight_"+dsName).c_str(),("btag weight versus light jet PT (partFlav = light) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo2D["bTagWeight_vs_JetPt_partFlavLight_"+dsName] = new TH2F(("bTagWeight_vs_JetPt_partFlavLight_"+dsName).c_str(),("btag weight versus jet PT (partFlav = light) for "+dsTitle).c_str(),800,0,800,150,0,1);
    histo1D["partonFlav_BJets_"+dsName] = new TH1F(("partonFlav_BJets_"+dsName).c_str(),("parton flavour for b-jets ("+dsTitle+")").c_str(),101,-0.5,100.5);
    histo1D["partonFlav_LightJets_"+dsName] = new TH1F(("partonFlav_LightJets_"+dsName).c_str(),("parton flavour for light jets ("+dsTitle+")").c_str(),101,-0.5,100.5);
    
    histo1D["nPartons_PartFlav5_"+dsName] = new TH1F(("nPartons_PartFlav5_"+dsName).c_str(), ("Number of partons with parton-flavour 5 ("+dsTitle+")").c_str(),6, -0.5, 5.5);
    
    histo2D["MlbMqqbChiSq_vs_nPv_"+dsName] = new TH2F(("MlbMqqbChiSq_vs_nPv_"+dsName).c_str(), ("Mlb-Mqqb #chi^{2} versus nr of primary vertices for "+dsTitle).c_str(), 40, 0, 40, 50,0, 10); 
    histo2D["HadrTopMass_vs_nPv_"+dsName] = new TH2F(("HadrTopMass_vs_nPv_"+dsName).c_str(), ("Hadronic top mass versus nr of primary vertices for "+dsTitle).c_str(), 40, 0, 40, 175, 30, 320);
    histo2D["HadrWMass_vs_nPv_"+dsName] = new TH2F(("HadrWMass_vs_nPv_"+dsName).c_str(), ("Hadronic w-boson mass versus nr of primary vertices for "+dsTitle).c_str(), 40, 0, 40, 120, 15, 175);
    histo2D["MlbMqqbChiSq_vs_LeptPt_"+dsName] = new TH2F(("MlbMqqbChiSq_vs_LeptPt_"+dsName).c_str(), ("Mlb-Mqqb #chi^{2} versus lepton pT for "+dsTitle).c_str(), 75, 0, 200, 50, 0, 10);
    histo2D["HadrTopMass_vs_LeptPt_"+dsName] = new TH2F(("HadrTopMass_vs_LeptPt_"+dsName).c_str(), ("Hadronic top mass versus lepton pt for "+dsTitle).c_str(), 75, 0, 200, 175, 30, 320);
    histo2D["HadrWMass_vs_LeptPt_"+dsName] = new TH2F(("HadrWMass_vs_LeptPt_"+dsName).c_str(), ("Hadronic w-boson mass versus lepton pt for "+dsTitle).c_str(), 75, 0, 200, 120, 15, 175);
    histo2D["MlbMqqbChiSq_vs_JetPt_"+dsName] = new TH2F(("MlbMqqbChiSq_vs_JetPt_"+dsName).c_str(), ("Mlb-Mqqb #chi^{2} versus jet pT for "+dsTitle).c_str(), 150, 0, 350, 50, 0, 10);
    histo2D["HadrTopMass_vs_JetPt_"+dsName] = new TH2F(("HadrTopMass_vs_JetPt_"+dsName).c_str(), ("Hadronic top mass versus jet pT for "+dsTitle).c_str(), 150, 0, 350, 175, 30, 320);
    histo2D["HadrWMass_vs_JetPt_"+dsName] = new TH2F(("HadrWMass_vs_JetPt_"+dsName).c_str(), ("Hadronic w-boson mass versus jet pT for "+dsTitle).c_str(), 150, 0, 350, 120, 15, 175);
    
    histo2D["MlbMqqbChiSq_vs_HadrTopMass_"+dsName] = new TH2F(("MlbMqqbChiSq_vs_HadrTopMass_"+dsName).c_str(), ("Mlb-Mqqb #chi^{2} versus hadronic top mass for "+dsTitle).c_str(), 100, 150, 200, 50, 0, 10);
    histo2D["MlbMqqbChiSq_vs_HadrWMass_"+dsName] = new TH2F(("MlbMqqbChiSq_vs_HadrWMass_"+dsName).c_str(), ("Mlb-Mqqb #chi^{2} versus hadronic w-boson mass for "+dsTitle).c_str(), 120, 40, 140, 50, 0, 10);
    histo2D["HadrTopMass_vs_HadrWMass_"+dsName] = new TH2F(("HadrTopMass_vs_HadrWMass_"+dsName).c_str(),("Hadronic top mass versus hadronic w-boson mass for "+dsTitle).c_str(), 120, 15, 175, 175, 30, 320);

    //Only want these to be filled for TTbarJets_SemiMu!  (what about FullLept and FullHadr ??)
    int BothB = 0;
    double tagEffPartFlavB = 0, tagEffPartFlavNotB = 0, tagEffPartFlavC = 0, tagEffPartFlavL = 0;
    int nrPartFlavBJets = 0, nrPartFlavNotBJets = 0, nrPartFlavCJets = 0, nrPartFlavLJets = 0;
    int lightIsMatched = 0, evtHasThreeLights = 0;
    int firstIsMatched = 0;
    int firstLightIsFirst = 0, firstLightIsSecond = 0, firstLightIsElse = 0;
    int lightIsFirst = 0, lightIsSecond = 0, lightIsThird = 0;
    int firstLightMatchedThreeExist = 0, firstLightIsThird = 0;
    int firstJetChosen = 0, secondJetChosen = 0, thirdJetChosen = 0;
    int firstJetChosenAndCorrect = 0, secondJetChosenAndCorrect = 0, thirdJetChosenAndCorrect = 0;

    int chiSqLowerThan5 = 0, chiSqLowerThan10 = 0, chiSqLowerThan50 = 0;
    int chiSqLowerThan5Correct = 0, chiSqLowerThan10Correct = 0, chiSqLowerThan50Correct = 0;
    double SFcount = 0.;

    if(isTTbarSemiLept){
      histo1D["PartFlavBJets_BTag"+dsName] = new TH1F(("PartFlavourBJets_BTag"+dsName).c_str(), ("Parton flavour of b-jets for "+dsTitle+" events after 2T CSV tags").c_str(), 10, -0.5, 9.5);
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(1,"Both b");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(2,"b and light (u)");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(3,"b and light (d)");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(4,"b and light (s)");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(5,"b and c");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(6,"b and gluon");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(7,"b and undef");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(8,"Both light (udscg)");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(9,"Both undef");
      histo1D["PartFlavBJets_BTag"+dsName]->GetXaxis()->SetBinLabel(10,"else ..");

      histo1D["PartFlavLightJets_BTag"+dsName] = new TH1F(("PartFlavourLightJets_BTag"+dsName).c_str(), ("Parton flavour of light jets for "+dsTitle+" events after 2T CSV tags").c_str(), 12, -0.5, 11.5);
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(1,"Both light (uds)");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(2,"Both gluon");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(3,"Light and gluon");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(4,"Light and c");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(5,"Gluon and c");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(6,"Light and b");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(7,"Gluon and b");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(8,"One undef");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(9,"Both b");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(10,"Both c");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(11,"b and c");
      histo1D["PartFlavLightJets_BTag"+dsName]->GetXaxis()->SetBinLabel(12,"else ..");

      histo2D["BTaggedJets_BTagEffvsPt"]   = new TH2F("BTaggedJets_BTagEffvsPt",  "B-tag efficiency versus pT for b-tagged jets",                        120,0,400,75,0,1);
      histo2D["LightJets_BTagEffvsPt"]     = new TH2F("LightJets_BTagEffvsPt",    "B-tag efficiency versus pT for light jets",                           120,0,400,75,0,1);
      histo2D["PartFlavB_BTagEffvsPt"]     = new TH2F("PartFlavB_BTagEffvsPt",    "B-tag efficiency versus pT for jets with parton flavour of b",        120,0,400,75,0,1);
      histo2D["PartFlavC_BTagEffvsPt"]     = new TH2F("PartFlavC_BTagEffvsPt",    "B-tag efficiency versus pT for jets with parton flavour of c",        120,0,400,75,0,1);
      histo2D["PartFlavLight_BTagEffvsPt"] = new TH2F("PartFlavLight_BTagEffvsPt","B-tag efficiency versus pT for jets with parton flavour of light",    120,0,400,75,0,1);
      histo2D["PartFlavGluon_BTagEffvsPt"] = new TH2F("PartFlavGluon_BTagEffvsPt","B-tag efficiency versus pT for jets with parton flavour of gluon",    120,0,400,75,0,1);
      histo2D["PartFlavUndef_BTagEffvsPt"] = new TH2F("PartFlavUndef_BTagEffvsPt","B-tag efficiency versus pT for jets with parton flavour of undefined",120,0,400,75,0,1);

      histo1D["BTagWeight_SelJetsOnly"] = new TH1F("BTagWeight_SelJetsOnly","B-tag weight for the four selected jets only",200,0,2);
      histo1D["BTagWeight_NonSelJetsOnly"] = new TH1F("BTagWeight_NonSelJetsOnly","B-tag weight for the non selected jets only",200,0,2);
      histo1D["BTagWeight_AllJets"] = new TH1F("BTagWeight_AllJets","B-tag weight for all jets in the event",200,0,2);

      histo1D["CSVDiscr_LightJets_PartFlavNob"] = new TH1F("CSVDiscr_LightJets_PartFlavNob", ("CSV discriminant of light jets for "+dsTitle+" events with parton flavour not 5").c_str(), 75, -1, 2);
      histo1D["CSVDiscr_LightJets_PartFlavb"]   = new TH1F("CSVDiscr_LightJets_PartFlavb",   ("CSV discriminant of light jets for "+dsTitle+" events with parton flavour 5").c_str(),     75, -1, 2);

      histo1D["CSVDiscr_CorrectBJets"]     = new TH1F("CSVDiscr_CorrectBJets",    "CSV discriminant for the two correct b-jets in the event",    100,-0.005,1.005);
      histo1D["CSVDiscr_CorrectLightJets"] = new TH1F("CSVDiscr_CorrectLightJets","CSV discriminant for the two correct light jets in the event",100,-0.005,1.005);
      histo1D["BTagEff_partFlavB"]    = new TH1F("BTagEff_partFlavB",   "b-tag efficiency for jets with parton flavour of b",     100, 0, 0.6);
      histo1D["BTagEff_partFlavC"]    = new TH1F("BTagEff_partFlavC",   "b-tag efficiency for jets with parton flavour of c",     100, 0, 0.1);
      histo1D["BTagEff_partFlavL"]    = new TH1F("BTagEff_partFlavL",   "b-tag efficiency for jets with parton flavour of udsg",  100, 0, 0.01);
      histo1D["BTagEff_partFlavNotB"] = new TH1F("BTagEff_partFlavNotB","b-tag efficiency for jets with parton flavour of udscg", 100, 0, 0.1);

      histo1D["Mlb_CorrTT_"+dsName]        = new TH1F(("Mlb_CorrTT_"+dsName).c_str(),       ("mass_{l,b} for the actual "+dsTitle+" events").c_str(),     200,  0, 250);
      histo1D["HadrMTop_CorrTT_"+dsName]   = new TH1F(("HadrMTop_CorrTT_"+dsName).c_str(),  ("Hadronic m_{t} for the actual "+dsTitle+" events").c_str(), 200, 50, 350);
      histo1D["HadrMW_CorrTT_"+dsName]     = new TH1F(("HadrMW_CorrTT_"+dsName).c_str(),    ("Hadronic m_{W} for the actual "+dsTitle+" events").c_str(), 200,  0, 250);
      histo2D["MW_vs_MTop_CorrTT_"+dsName] = new TH2F(("MW_vs_MTop_CorrTT_"+dsName).c_str(),("Hadronic m_{W} versus m_{t} for the actual "+dsTitle+" jet combinations").c_str(),200,50,350,200,0,250);

      histo1D["LightJetIndex_Correct"] = new TH1F("LightJetIndex_Correct","Light jet index of the correct light jets", 10, -0.5, 10.5);
      histo1D["LightJetIndex_Wrong"] = new TH1F("LightJetIndex_Wrong",    "Light jet index of the wrong light jets",   10, -0.5, 10.5);
      histo1D["ChosenJetCombi"] = new TH1F("ChosenJetCombi","Number of times the different jet-combinations are chosen", 3, -0.5, 2.5);
      histo1D["CorrectJetCombi"] = new TH1F("CorrectJetCombi","Number of times the different jet-combinations are correct", 3, -0.5, 2.5);

      histo1D["CorrectLightJetIndex"] = new TH1F("CorrectLightJetIndex","Correct jet index of the light jets", 10, -0.5, 9.5);
      histo1D["ChosenLightJetIndex"]  = new TH1F("ChosenLightJetIndex", "Chosen jet index of the light jets",  10, -0.5, 9.5);
      histo1D["LightJetsCorrComb"]    = new TH1F("LightJetsCorrComb",   "Correct jet combination for the light jets", 6, -0.5, 5.5);
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(1,"0 & 1 correct");
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(2,"0 & 2 correct");
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(3,"1 & 2 correct");
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(4,"Only 0 correct");
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(5,"Only 1 correct");
      histo1D["LightJetsCorrComb"]->GetXaxis()->SetBinLabel(6,"Only 2 correct");
      histo1D["ThirdJetInfo"] = new TH1F("ThirdJetInfo", "Events with more than two light jets",4, -0.5, 3.5);
      histo1D["ThirdJetInfo"]->GetXaxis()->SetBinLabel(1,"All events");
      histo1D["ThirdJetInfo"]->GetXaxis()->SetBinLabel(2,"Light jet not in leading 2");
      histo1D["ThirdJetInfo"]->GetXaxis()->SetBinLabel(3,"Light jet is 3th one");
      histo1D["ThirdJetInfo"]->GetXaxis()->SetBinLabel(4,"Chosen jet is 3th one");

      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]        = new TH1F(("Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(),  ("mass_{l,b} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
        histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]   = new TH1F(("HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(), ("Hadronic m_{t} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 50, 350);
        histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]     = new TH1F(("HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(),   ("Hadronic m_{W} for the actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 0, 250);
        histo2D["MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]] = new TH2F(("MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]).c_str(), ("Hadronic m_{W} versus m_{t} for actual "+dsTitle+" events -- "+bTitle[ibTag]).c_str(), 200, 50, 350, 200, 0,250);
      }
    }
    
    for(int ii = 0; ii < 2; ii++){
      string leptFlav = leptFlavs[ii];
      histo1D["CosTheta_Gen_SelEvts_"+dsName+""+leptFlav] = new TH1F(("CosTheta_Gen_SelEvts_"+dsName+""+leptFlav).c_str(),("Cos #theta^{*}_{gen} distribution for the selected "+dsTitle+" events").c_str(),200,-1,1);
    }
    histo1D["CosTheta_Gen_SelEvts_"+dsName] = new TH1F(("CosTheta_Gen_SelEvts_"+dsName).c_str(),("Cos #theta^{*}_{gen} distribution for the selected "+dsTitle+" events").c_str(),200,-1,1);

    TFLight* tfLight_mu = 0;
    //TFLight* tfLight_el = 0;
    TTree* TFLightTree;
    TFile* TFLightFile = 0;
    if(isTTbarSemiLept && createTFTree){
      std::cout << " --> Going into this loop for dataset : " << dsName  << std::endl;
      //Initialize LightTuple (TFTree) specific stuff:
      TFLightFile = new TFile(("TFTree/TFLight_"+dsName+"_"+systematic+""+tfFill+".root").c_str(),"RECREATE");  
      TFLightTree = new TTree("TFLightTree",("Tree containing the TFLight information for "+dsName+" at "+timestamp()).c_str());
      TFLightTree->Branch("TheTFLight_muCh","TFLight",&tfLight_mu);
      //TFLightTree->Branch("TheTFLight_elCh","TFLight",&tfLight_el);
    }

    //-----------------------//
    // Load personal classes //
    //-----------------------//
    bTagStudy.InitializeDataSet(dsName);
    LHCOOutput lhcoOutput(verbosity, getLHCOOutput, splitLeptonChargeLHCO, getCorrectAndWrongLHCO);
    extraEvtSelCuts.Initialize(bTitle[0], dsName);
    lhcoOutput.Initialize("Reco", onlyMuonChannel, dsName);

    ofstream EvtNrMatching;
    if(bTagChoiceMade && (getLHCOOutput || getEvtNrInfo)){
      EvtNrMatching.open(("MadWeightInput/AnalyzerOutput/MWEventNrMatching_"+dsName+"_"+systematic+".txt").c_str());
      EvtNrMatching << "  Event Nr     Extra cuts     Decay channel      MW Number     TTbar splitting file      MW Nr (split)    MC scale factor " << endl;
      EvtNrMatching << "   * Lumi = " << Luminosity << endl;
      EvtNrMatching << "   * NormFactor = " << datasets[iDataSet]->NormFactor() << endl;
    }

    // --------------------------- //
    //  Start looping over events  //
    // --------------------------- //    
    int nSelectedMu = 0, nSelectedEl = 0;
    for(unsigned int iEvt = 0; iEvt < nrEvts; iEvt++){
      inLightTree->GetEvent(iEvt);
      if(iEvt%5000 == 0)
	std::cout<<"    Processing the "<<iEvt<<"th event ("<<((double)iEvt/(double)nrEvts)*100<<"%)"<<" -> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\r";

      //*** Start with corrections ***//
      // Beam scraping and PU reweighting
      double lumiWeight = 1;  
      if(!isData){   
        if(doLumiWeightShift == "Nominal")    lumiWeight = LumiWeights.ITweight( (int)light->nTruePU() );
        else if(doLumiWeightShift == "Plus")  lumiWeight = LumiWeightsUp.ITweight( (int)light->nTruePU() );
        else if(doLumiWeightShift == "Minus") lumiWeight = LumiWeightsDown.ITweight( (int)light->nTruePU() );
      }

      int eventId = light->eventID();
      int runId = light->runID();
      int nTruePU = light->nTruePU();
      int nPrimVertices = light->nPV();
      float scaleFactor = light->scaleFactor();

      vector<float> bTagCSV = light->CSVbTag();
      vector<TLorentzVector> selJets = light->selectedJets();
      vector<int> partFlavour = light->selectedJetsPartonFlavour();
      vector<int> partFlavAbs; partFlavAbs.clear();
      for(int ii = 0; ii < partFlavour.size(); ii++) partFlavAbs.push_back(abs(partFlavour[ii]));
      TLorentzVector selLepton = light->selectedLepton();
      TLorentzVector MET = light->met();
      int decayCh = light->decayChannel();  //0 = semiMu and 1 = semiEl
      std::string leptChannel;
      if(decayCh == 0) leptChannel = "_mu";
      else if(decayCh == 1) leptChannel = "_el";
      float leptCharge = light->leptonCharge();
      vector<int> corrJetCombi = light->correctJetCombi();    //0 = LeptB, 1 = HadrB, 2 = Quark1 & 3 = Quark2
      float genCosTheta = light->genCosTh();

      if(onlyMuonChannel && decayCh == 1) continue;

      //Distributions before PU reweighting is taken into account!
      MSPlot["nPV_BeforePU"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_BeforePU"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_SecondJet_BeforePU"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_ThirdJet_BeforePU"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_FourthJet_BeforePU"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      if(selJets.size() > 4) MSPlot["JetPt_FifthJet_BeforePU"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_BeforePU"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      scaleFactor = scaleFactor * lumiWeight;
      histo1D["lumiWeight_"+dsName]->Fill(lumiWeight);

      MSPlot["nPV_PU_NoLeptonSF"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_SecondJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_ThirdJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_FourthJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      if(selJets.size() > 4) MSPlot["JetPt_FifthJet_PU_NoLeptonSF"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_PU_NoLeptonSF"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      //Get the correct lepton scalefactor!
      double leptonSF = 1;
      if( !isData ){
        if(decayCh == 0 ) leptonSF = leptonTools->getMuonSF(selLepton.Eta(), selLepton.Pt(), syst_muonSF);
        else              leptonSF = leptonTools->getElectronSF(selLepton.Eta(), selLepton.Pt(), syst_elecSF);
      }
      scaleFactor = scaleFactor * leptonSF;
      histo1D["leptonSF_"+dsName]->Fill(leptonSF); 

      //Fill the b-tag histo's in case they do not yet exist
      if(!bTagPlotsMade){
        bTagTool->FillMCEfficiencyHistos(selJets, partFlavour, bTagCSV);
        continue;  //In case the histo's still need to be filled, do not continue with the analysis!
      }

      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      //ooOOooOOoo      Reading out nTuples done           ooOOooOOoo
      //ooOOooOOoo-----------------------------------------ooOOooOOoo
      //ooOOooOOoo      Start of actual analysis           ooOOooOOoo
      //ooooooooOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOooooooooooooOOOOO
      bTagStudy.CalculateJets(selJets, bTagCSV, corrJetCombi, selLepton, datasets[iDataSet], Luminosity*scaleFactor);
 
      //MSPlots with number of jets information before requiring at least two b-jets and at least 2 light jets!
      MSPlot["nSelectedJets_BeforeBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);

      MSPlot["nPV_BeforeBTag"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_LeadingJet_BeforeBTag"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_SecondJet_BeforeBTag"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_ThirdJet_BeforeBTag"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["JetPt_FourthJet_BeforeBTag"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      if(selJets.size() > 4) MSPlot["JetPt_FifthJet_BeforeBTag"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
      MSPlot["LeptonPt_BeforeBTag"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);

      //Store the expected Mlb, Mqqb, Mt and MW distributions
      if(isTTbarSemiLept){
        if(corrJetCombi[0] != 9999 && corrJetCombi[1] != 9999 && corrJetCombi[2] != 9999 && corrJetCombi[3] != 9999){
          histo1D["Mlb_CorrTT_"+dsName]->Fill( (selLepton+selJets[corrJetCombi[0]]).M());
          histo1D["HadrMTop_CorrTT_"+dsName]->Fill( (selJets[corrJetCombi[1]]+selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M() );
          histo1D["HadrMW_CorrTT_"+dsName]->Fill( (selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M());
          histo2D["MW_vs_MTop_CorrTT_"+dsName]->Fill( (selJets[corrJetCombi[1]]+selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M(), (selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M() );
        }
      }
    
      //Store the CSV discriminant for the actual light and b-jets such that this can be plotted together!
      if(isTTbarSemiLept){
        if(corrJetCombi[0] != 9999) histo1D["CSVDiscr_CorrectBJets"]->Fill(bTagCSV[corrJetCombi[0]]);
        if(corrJetCombi[1] != 9999) histo1D["CSVDiscr_CorrectBJets"]->Fill(bTagCSV[corrJetCombi[1]]);
        if(corrJetCombi[2] != 9999) histo1D["CSVDiscr_CorrectLightJets"]->Fill(bTagCSV[corrJetCombi[2]]);
        if(corrJetCombi[3] != 9999) histo1D["CSVDiscr_CorrectLightJets"]->Fill(bTagCSV[corrJetCombi[3]]);
      }
 
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        vector<int> selJetCombi = bTagStudy.getIndices(ibTag);

        for(int iJet = 0; iJet < selJets.size(); iJet++)
          MSPlot["CSVDistr_BeforeBTag"+leptChannel]->Fill(bTagCSV[iJet], datasets[iDataSet], true, Luminosity*scaleFactor);

        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_BeforeBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);

        MSPlot["nBJets_BeforeBTag"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_BeforeBTag"+leptChannel]->Fill(bTagStudy.getNrLightJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(!isData){
            if(partFlavour[ii] == 1 || partFlavour[ii] == 2 || partFlavour[ii] == 3 ){
              MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(1, datasets[iDataSet], true, Luminosity*scaleFactor);
              MSPlot["JetsPt_partFlavLight_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
            }
            else if(partFlavour[ii] == 4 ){
              MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(2, datasets[iDataSet], true, Luminosity*scaleFactor);
              MSPlot["JetsPt_partFlavC_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
            }
            else if(partFlavour[ii] == 5 ){
              MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(3, datasets[iDataSet], true, Luminosity*scaleFactor);
              MSPlot["JetsPt_partFlavB_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
            }
            else if(partFlavour[ii] == 21 ){
              MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(4, datasets[iDataSet], true, Luminosity*scaleFactor);
              MSPlot["JetsPt_partFlavGluon_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
            }
            else{
              MSPlot["partFlav_BeforeBTag"+leptChannel]->Fill(0, datasets[iDataSet], true, Luminosity*scaleFactor);
              MSPlot["JetsPt_partFlavUndef_BeforeBTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
            }
          }
        }
      
        if(isTTbarSemiLept && !fillTFAfterCuts && createTFTree){  //Combine both muon and electron events!
          tfLight_mu = new TFLight();
            
          double fullScaleFactor = scaleFactor;   
          tfLight_mu->setFullScaleFactor(fullScaleFactor);
          tfLight_mu->setSelectedJets(selJets);
          tfLight_mu->setSelectedLepton(selLepton);
          tfLight_mu->setCorrectJetCombi(corrJetCombi);
	
          //Store the information needed for the TF (but only has value when dataset is ttbar)
          tfLight_mu->setGenVectorLight1( light->genVectorLight1() );
          tfLight_mu->setGenVectorLight2( light->genVectorLight2() );
          tfLight_mu->setGenVectorHadrB(  light->genVectorHadrB()  );
          tfLight_mu->setGenVectorLeptB(  light->genVectorLeptB()  );
          tfLight_mu->setGenVectorLepton( light->genVectorLepton() );

          TFLightTree->Fill();
          delete tfLight_mu;
          //----  End of Tree file filling (for TF's after evtSel)  ----//
        }

        //*****************//  
        // Apply the b-tag //
        //*****************//  
        if( bTagStudy.getNrBTaggedJets(ibTag) < 2 || bTagStudy.getNrLightJets(ibTag) < 2 ) continue;
        int BLeptIndex = bTagStudy.getBLeptIndex(ibTag);
        int BHadrIndex = bTagStudy.getBHadrIndex(ibTag);
        int Light1Index = bTagStudy.getLight1Index(ibTag);
        int Light2Index = bTagStudy.getLight2Index(ibTag);

        //--- Check whether the event is correctly reconstructed  ---// 
        //---  (jetCombi is initialized to 9999 for all dataSets) ---//
        int CWUIndex = 999;
        if(corrJetCombi[0] != 9999 && corrJetCombi[1] != 9999 && corrJetCombi[2] != 9999 && corrJetCombi[3] != 9999){
          if( selJetCombi[0] == corrJetCombi[0] && selJetCombi[1] == corrJetCombi[1]  &&
             (selJetCombi[2] == corrJetCombi[2] || selJetCombi[2] == corrJetCombi[3]) &&
             (selJetCombi[3] == corrJetCombi[2] || selJetCombi[3] == corrJetCombi[3]) )
            CWUIndex = 0;
          else
            CWUIndex = 1; 
        }
        else
          CWUIndex = 2;

        //How often is the third light jet one of the correct ones ...
        if(isTTbarSemiLept){
          histo1D["CorrectLightJetIndex"]->Fill(corrJetCombi[2]);
          histo1D["CorrectLightJetIndex"]->Fill(corrJetCombi[3]);
          histo1D["ChosenLightJetIndex"]->Fill(Light1Index);
          histo1D["ChosenLightJetIndex"]->Fill(Light2Index);

          vector<int> lightJet = bTagStudy.getLightJets(ibTag);
          if(                            (lightJet[0]==corrJetCombi[2] && lightJet[1]==corrJetCombi[3]) || (lightJet[0]==corrJetCombi[3] && lightJet[1]==corrJetCombi[2])) histo1D["LightJetsCorrComb"]->Fill(0);
          else if(lightJet.size() > 2 && (lightJet[0]==corrJetCombi[2] && lightJet[2]==corrJetCombi[3]) || (lightJet[0]==corrJetCombi[3] && lightJet[2]==corrJetCombi[2])) histo1D["LightJetsCorrComb"]->Fill(1);
          else if(lightJet.size() > 2 && (lightJet[1]==corrJetCombi[2] && lightJet[2]==corrJetCombi[3]) || (lightJet[1]==corrJetCombi[3] && lightJet[2]==corrJetCombi[2])) histo1D["LightJetsCorrComb"]->Fill(2);
          else if(lightJet[0] == corrJetCombi[2] || lightJet[0] == corrJetCombi[3]) histo1D["LightJetsCorrComb"]->Fill(3);
          else if(lightJet[1] == corrJetCombi[2] || lightJet[1] == corrJetCombi[3]) histo1D["LightJetsCorrComb"]->Fill(4);
          else if(lightJet.size() > 2 && (lightJet[2] == corrJetCombi[2] || lightJet[2] == corrJetCombi[3]) ) histo1D["LightJetsCorrComb"]->Fill(5);

          //Get some percentages to understand the 3rd light jet improvement
          if( corrJetCombi[2] != 9999 && corrJetCombi[3] != 9999){  //Look only at events with matched light jets

            if( (Light1Index == lightJet[0] && Light2Index == lightJet[1]) || (Light1Index == lightJet[1] && Light2Index == lightJet[0]) ) histo1D["ChosenJetCombi"]->Fill(0);
            if( (Light1Index == lightJet[1] && Light2Index == lightJet[2]) || (Light1Index == lightJet[2] && Light2Index == lightJet[1]) ) histo1D["ChosenJetCombi"]->Fill(1);
            if( (Light1Index == lightJet[0] && Light2Index == lightJet[2]) || (Light1Index == lightJet[2] && Light2Index == lightJet[0]) ) histo1D["ChosenJetCombi"]->Fill(2);

            if( (corrJetCombi[2] == lightJet[0] && corrJetCombi[3] == lightJet[1]) || (corrJetCombi[2] == lightJet[1] && corrJetCombi[3] == lightJet[0]) ) histo1D["CorrectJetCombi"]->Fill(0);
            if( (corrJetCombi[2] == lightJet[1] && corrJetCombi[3] == lightJet[2]) || (corrJetCombi[2] == lightJet[2] && corrJetCombi[3] == lightJet[1]) ) histo1D["CorrectJetCombi"]->Fill(1);
            if( (corrJetCombi[2] == lightJet[0] && corrJetCombi[3] == lightJet[2]) || (corrJetCombi[2] == lightJet[2] && corrJetCombi[3] == lightJet[0]) ) histo1D["CorrectJetCombi"]->Fill(2);

            if(lightJet.size() > 2 && ( lightJet[0] == Light1Index || lightJet[0] == Light2Index)){
              firstJetChosen++;
              if(lightJet[0] == corrJetCombi[2] || lightJet[0] == corrJetCombi[3]) firstJetChosenAndCorrect++;
            }
            if(lightJet.size() > 2 && ( lightJet[1] == Light1Index || lightJet[1] == Light2Index)){
              secondJetChosen++;
              if(lightJet[1] == corrJetCombi[2] || lightJet[1] == corrJetCombi[3]) secondJetChosenAndCorrect++;
            }
            if(lightJet.size() > 2 && (lightJet[2] == Light1Index || lightJet[2] == Light2Index)){
              thirdJetChosen++;
              if(lightJet[2] == corrJetCombi[2] || lightJet[2] == corrJetCombi[3]) thirdJetChosenAndCorrect++;
            }
            
            for(int iLight = 0; iLight < lightJet.size(); iLight++){
              if(lightJet[iLight] == corrJetCombi[2] || lightJet[iLight] == corrJetCombi[3] ) histo1D["LightJetIndex_Correct"]->Fill(iLight);
              else                                                                            histo1D["LightJetIndex_Wrong"]->Fill(iLight);
            }
    
            lightIsMatched++;
            if(corrJetCombi[2] == lightJet[0] || corrJetCombi[3] == lightJet[0]) lightIsFirst++;        //How often is the light jet the leading jet
            if(corrJetCombi[2] == lightJet[1] || corrJetCombi[3] == lightJet[1]) lightIsSecond++;       //How often is the light jet the second jet
            if(lightJet.size() > 2){
              evtHasThreeLights++;
              if(corrJetCombi[2] == lightJet[2] || corrJetCombi[3] == lightJet[2]) lightIsThird++;      //How often is the light the third jet
            }
          }
          if( corrJetCombi[2] != 9999){
            firstIsMatched++;
            if(corrJetCombi[2] == lightJet[0])      firstLightIsFirst++;
            else if(corrJetCombi[2] == lightJet[1]) firstLightIsSecond++;
            else{
              firstLightIsElse++;
              if(lightJet.size() > 2){
                firstLightMatchedThreeExist++;
                if(corrJetCombi[2] == lightJet[2]) firstLightIsThird++;
              }
            }
            
          }

          if(lightJet.size() > 2){
            histo1D["ThirdJetInfo"]->Fill(0);
            if( (corrJetCombi[2] != lightJet[0] & corrJetCombi[2] != lightJet[1] && corrJetCombi[2] != 9999) || (corrJetCombi[3] != lightJet[0] && corrJetCombi[3] != lightJet[1] && corrJetCombi[3] != 9999) ) histo1D["ThirdJetInfo"]->Fill(1);
            if( corrJetCombi[2] == lightJet[2] || corrJetCombi[3] == lightJet[2]) histo1D["ThirdJetInfo"]->Fill(2);
            if(Light1Index == lightJet[2] || Light2Index == lightJet[2]) histo1D["ThirdJetInfo"]->Fill(3);
          }
        }

        if(partFlavour[BLeptIndex] == 21) MSPlot["BJetPt_partFlavGluon_BTag"+leptChannel]->Fill( selJets[BLeptIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(partFlavour[BHadrIndex] == 21) MSPlot["BJetPt_partFlavGluon_BTag"+leptChannel]->Fill( selJets[BHadrIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_BTag"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        }

        MSPlot["nBJets_BTag"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_BTag"+leptChannel]->Fill(bTagStudy.getNrLightJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        //vector<TLorentzVector> selJetsAfterBTag;
        //vector<int> partFlavAfterBTag;
        //vector<float> bTagCSVAfterBTag;
        //selJetsAfterBTag.push_back(selJets[BLeptIndex]); partFlavAfterBTag.push_back(partFlavour[BLeptIndex]); bTagCSVAfterBTag.push_back(bTagCSV[BLeptIndex]);
        //selJetsAfterBTag.push_back(selJets[BHadrIndex]); partFlavAfterBTag.push_back(partFlavour[BHadrIndex]); bTagCSVAfterBTag.push_back(bTagCSV[BHadrIndex]);
        //selJetsAfterBTag.push_back(selJets[Light1Index]); partFlavAfterBTag.push_back(partFlavour[Light1Index]); bTagCSVAfterBTag.push_back(bTagCSV[Light1Index]);
        //selJetsAfterBTag.push_back(selJets[Light2Index]); partFlavAfterBTag.push_back(partFlavour[Light2Index]); bTagCSVAfterBTag.push_back(bTagCSV[Light2Index]);

//        if(decayCh == 0) nSelectedMu += 1;
//        else if(decayCh == 1) nSelectedEl += 1;

        MSPlot["nPV"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_LeadingJet"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_SecondJet"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_ThirdJet"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_FourthJet"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(selJets.size() > 4) MSPlot["JetPt_FifthJet"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonPt"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetEta_LeadingJet"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonEta"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonCharge"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

        //Get the bTag scaleFactor!
        double BTagWeight = 1;
        double BTagWeight2b = 1;
        if(bTagPlotsMade && !isData){
          BTagWeight = bTagTool->getMCEventWeight(selJets, partFlavour, bTagCSV, syst_btag, syst_mistag);

 /*         //Take into account the two b-tags (following 1b this time ...)
          vector<float> bTagSFs, bTagEffs;
          bTagSFs.clear(); bTagEffs.clear();
          double w0TagsMC = 1;
          double w0TagsData = 1;
          //std::cout << " w0Tags is = " << endl;
          for(int iJet = 0; iJet < selJets.size(); iJet++){
            bTagSFs.push_back(bTagTool->getSF(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet], "CSVT", syst_btag, syst_mistag));
            bTagEffs.push_back(bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]));
            w0TagsMC = w0TagsMC * (1-bTagEffs[iJet]);
            w0TagsData = w0TagsData * (1-bTagEffs[iJet]*bTagSFs[iJet]);
            //std::cout << "  " << iJet << ") Jet with flavour " << partFlavAfterBTag[iJet] << ": 1 - " << bTagEffs[iJet] << " * " << bTagSFs[iJet] << " = " <<  1-bTagEffs[iJet]*bTagSFs[iJet] << " ==> w0 = " << w0TagsData << std::endl;
          }

          //Now calculate the event weight w(>= 2 btags):
          double w1TagIndivMC[4] = {0, 0, 0, 0};
          double w1TagIndivData[4] = {0, 0, 0, 0};
          double w1TagMC = 0;
          double w1TagData = 0;
          //std::cout << " w1Tag = " << endl;
          for(int i = 0; i< selJets.size(); i++){
            w1TagIndivMC[i] = bTagEffs[i];
            w1TagIndivData[i] = bTagSFs[i]*bTagEffs[i];
            //std::cout << "  " << i << ") Jet with flavour " << partFlavAfterBTag[i] << ": " << bTagSFs[i] << " * " << bTagEffs[i] << " = " << w1TagIndivData[i] << std::endl;
            for(int j = 0; j < selJets.size(); j++){
              if(j != i) w1TagIndivMC[i] = w1TagIndivMC[i] * (1-bTagEffs[j]);
              if(j != i) w1TagIndivData[i] = w1TagIndivData[i] * (1-bTagSFs[j]*bTagEffs[j]);
              //if(j != i) std::cout << "   " << j << ") Jet with flavour " << partFlavAfterBTag[j] << ": 1- " << bTagSFs[j] << " * " << bTagEffs[j] << " = " << 1-bTagSFs[j]*bTagEffs[j] << " ==> w1 = " << w1TagIndivData[i] << std::endl;
            }
            w1TagMC = w1TagMC + w1TagIndivMC[i];
            w1TagData = w1TagData + w1TagIndivData[i];
            //std::cout << "  ==> w1 = " << w1TagData << std::endl;
          }
          BTagWeight2b = 1.0-w0TagsData-w1TagData;
          //BTagWeight2b = (1-w0TagsData-w1TagData)/(1-w0TagsMC-w1TagData);
          float DataBTag = 1-(1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3]) - ( bTagEffs[0]*bTagSFs[0]*((1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3])) + bTagEffs[1]*bTagSFs[1]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3])) + bTagEffs[2]*bTagSFs[2]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])) + bTagEffs[3]*bTagSFs[3]*((1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])));
          float MCBTag = 1-(1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3]) - ( bTagEffs[0]*((1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3])) + bTagEffs[1]*((1-bTagEffs[0])*(1-bTagEffs[2])*(1-bTagEffs[3])) + bTagEffs[2]*((1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])) + bTagEffs[3]*((1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])));
          float DataBTag1 = 1-(1-bTagEffs[0]*bTagSFs[0])*(1-bTagEffs[1]*bTagSFs[1])*(1-bTagEffs[2]*bTagSFs[2])*(1-bTagEffs[3]*bTagSFs[3]);
          float MCBTag1 = 1-(1-bTagEffs[0])*(1-bTagEffs[1])*(1-bTagEffs[2])*(1-bTagEffs[3]);
          //std::cout << " Full formula (2btag) 'Data' is " << DataBTag << " and 'MC' is " << MCBTag << " ==> Ratio is : " << DataBTag/MCBTag << std::endl;
          //std::cout << " Full formula (1btag) 'Data' is " << DataBTag1 << " and 'MC' is " << MCBTag1 << " ==> Ratio is : " << DataBTag1/MCBTag1 << std::endl; 
          BTagWeight2b = DataBTag/MCBTag;
         */
          if(isTTbarSemiLept){

            for(int iJet = 0; iJet < selJets.size(); iJet++){
              if(partFlavAbs[iJet] == 5){ 
                histo1D["BTagEff_partFlavB"]->Fill(bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet])); 
                tagEffPartFlavB += bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]); 
                nrPartFlavBJets++;
              }
              else{
                histo1D["BTagEff_partFlavNotB"]->Fill(bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]));
                tagEffPartFlavNotB += bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]); 
                nrPartFlavNotBJets++;
                if(partFlavAbs[iJet] == 4){
                  histo1D["BTagEff_partFlavC"]->Fill(bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]));
                  tagEffPartFlavC += bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]); 
                  nrPartFlavCJets++;
                }
                else{
                  histo1D["BTagEff_partFlavL"]->Fill(bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]));
                  tagEffPartFlavL += bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]); 
                  nrPartFlavLJets++;
                }
              } 
            }
          }
        }
        histo1D["BTagWeight_"+dsName]->Fill(BTagWeight);

        int nrPartFlavBs = 0;
        for(int i = 0; i < selJets.size(); i++){
          if(partFlavAbs[i] == 5) nrPartFlavBs++;
        }
        histo1D["nPartons_PartFlav5_"+dsName]->Fill(nrPartFlavBs);

        if(isTTbarSemiLept && bTagPlotsMade){
          if      (partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 5){                                                                     histo1D["PartFlavBJets_BTag"+dsName]->Fill(0); BothB++;}
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 1) || (partFlavAbs[BLeptIndex] == 1 && partFlavAbs[BHadrIndex] == 5))   histo1D["PartFlavBJets_BTag"+dsName]->Fill(1);
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 2) || (partFlavAbs[BLeptIndex] == 2 && partFlavAbs[BHadrIndex] == 5))   histo1D["PartFlavBJets_BTag"+dsName]->Fill(2);
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavour[BHadrIndex] == 3) || (partFlavAbs[BLeptIndex] == 3 && partFlavAbs[BHadrIndex] == 5))   histo1D["PartFlavBJets_BTag"+dsName]->Fill(3);
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 4) || (partFlavAbs[BLeptIndex] == 4 && partFlavAbs[BHadrIndex] == 5))   histo1D["PartFlavBJets_BTag"+dsName]->Fill(4);
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 21) || (partFlavAbs[BLeptIndex] == 21 && partFlavAbs[BHadrIndex] == 5)) histo1D["PartFlavBJets_BTag"+dsName]->Fill(5);
          else if((partFlavAbs[BLeptIndex] == 5 && partFlavAbs[BHadrIndex] == 0) || (partFlavAbs[BLeptIndex] == 0 && partFlavAbs[BHadrIndex] == 5))   histo1D["PartFlavBJets_BTag"+dsName]->Fill(6);
          else if((partFlavAbs[BLeptIndex] > 0 && partFlavAbs[BLeptIndex] <= 21) && (partFlavAbs[BHadrIndex] > 0 && partFlavAbs[BHadrIndex] <= 21 ))  histo1D["PartFlavBJets_BTag"+dsName]->Fill(7);
          else if(partFlavAbs[BLeptIndex] == 0 && partFlavAbs[BHadrIndex] == 0)                                                                       histo1D["PartFlavBJets_BTag"+dsName]->Fill(8);
          else{            histo1D["PartFlavBJets_BTag"+dsName]->Fill(9);}

          if      (partFlavAbs[Light1Index] == 0 || partFlavAbs[Light2Index] == 0)                                                                        histo1D["PartFlavLightJets_BTag"+dsName]->Fill(7);
          else if (partFlavAbs[Light1Index] < 4 && partFlavAbs[Light2Index] < 4)                                                                          histo1D["PartFlavLightJets_BTag"+dsName]->Fill(0);
          else if (partFlavAbs[Light1Index] == 21 && partFlavAbs[Light2Index] == 21)                                                                      histo1D["PartFlavLightJets_BTag"+dsName]->Fill(1);
          else if((partFlavAbs[Light1Index] == 21 && partFlavAbs[Light2Index] < 4) || (partFlavAbs[Light1Index] < 4 && partFlavAbs[Light2Index] == 21))   histo1D["PartFlavLightJets_BTag"+dsName]->Fill(2);
          else if((partFlavAbs[Light1Index] == 4 && partFlavour[Light2Index] < 4) || (partFlavAbs[Light1Index] < 4 && partFlavAbs[Light2Index] == 4))     histo1D["PartFlavLightJets_BTag"+dsName]->Fill(3);
          else if((partFlavAbs[Light1Index] == 4 && partFlavour[Light2Index] == 21) || (partFlavAbs[Light1Index] == 21 && partFlavAbs[Light2Index] == 4)) histo1D["PartFlavLightJets_BTag"+dsName]->Fill(4);
          else if((partFlavAbs[Light1Index] == 5 && partFlavAbs[Light2Index] < 4) || (partFlavAbs[Light1Index] < 4 && partFlavAbs[Light2Index] == 5))     histo1D["PartFlavLightJets_BTag"+dsName]->Fill(5);
          else if((partFlavAbs[Light1Index] == 5 && partFlavAbs[Light2Index] == 21) || (partFlavAbs[Light1Index] == 21 && partFlavAbs[Light2Index] == 5)) histo1D["PartFlavLightJets_BTag"+dsName]->Fill(6);
          else if(partFlavAbs[Light1Index] == 5 && partFlavAbs[Light2Index] == 5)                                                                         histo1D["PartFlavLightJets_BTag"+dsName]->Fill(8);
          else if(partFlavAbs[Light1Index] == 4 && partFlavAbs[Light2Index] == 4)                                                                         histo1D["PartFlavLightJets_BTag"+dsName]->Fill(9);
          else if((partFlavAbs[Light1Index] == 5 && partFlavAbs[Light2Index] == 4) || (partFlavAbs[Light1Index] == 4 && partFlavAbs[Light2Index] == 5))   histo1D["PartFlavLightJets_BTag"+dsName]->Fill(10);
          else{              histo1D["PartFlavLightJets_BTag"+dsName]->Fill(11);}

          //Draw the CSV discriminant for some type of light jets:
          if(partFlavAbs[Light1Index] != 5 && partFlavAbs[Light2Index] != 5){histo1D["CSVDiscr_LightJets_PartFlavNob"]->Fill(bTagCSV[Light1Index]); histo1D["CSVDiscr_LightJets_PartFlavNob"]->Fill(bTagCSV[Light2Index]);}
          else                                                              {histo1D["CSVDiscr_LightJets_PartFlavb"]->Fill(bTagCSV[Light1Index]); histo1D["CSVDiscr_LightJets_PartFlavb"]->Fill(bTagCSV[Light2Index]);}

          for(int iJet = 0; iJet < selJets.size(); iJet++){
            if(iJet == BLeptIndex || iJet == BHadrIndex)         histo2D["BTaggedJets_BTagEffvsPt"]->Fill(  selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(iJet == Light1Index || iJet == Light2Index)       histo2D["LightJets_BTagEffvsPt"]->Fill(    selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(partFlavAbs[iJet] == 5)                           histo2D["PartFlavB_BTagEffvsPt"]->Fill(    selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(partFlavAbs[iJet] == 4)                           histo2D["PartFlavC_BTagEffvsPt"]->Fill(    selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(partFlavAbs[iJet] == 21)                          histo2D["PartFlavGluon_BTagEffvsPt"]->Fill(selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(partFlavAbs[iJet] != 0 && partFlavAbs[iJet] < 4 ) histo2D["PartFlavLight_BTagEffvsPt"]->Fill(selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );
            if(partFlavAbs[iJet] == 0)                           histo2D["PartFlavUndef_BTagEffvsPt"]->Fill(selJets[iJet].Pt(), bTagTool->getTagEff(selJets[iJet].Pt(), selJets[iJet].Eta(), partFlavour[iJet]) );

          }

          //Check how much the non-selected events contribute to the weight
          float BTagWeightMCSelJets = 1,    BTagWeightDataSelJets = 1;
          float BTagWeightMCNonSelJets = 1, BTagWeightDataNonSelJets = 1;
          float BTagWeightMCAllJets = 1,    BTagWeightDataAllJets = 1;
          /*for(int i = 0; i < selJets.size(); i++){

            //std::cout << i << ") In analyzer: tagEff = " << bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]) << " and SF = " << bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag) << " ... with CSV value of " << bTagCSV[i] << std::endl;
            if(bTagCSV[i] >= 0.898){
              BTagWeightMCAllJets =  BTagWeightMCAllJets*bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]);
              BTagWeightDataAllJets = BTagWeightDataAllJets*bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag);
            }
            else{
              BTagWeightMCAllJets = BTagWeightMCAllJets*(1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]));
              BTagWeightDataAllJets = BTagWeightDataAllJets*(1.0-(bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag)));
            }
            //std::cout << "  ==> BTagWeightMCAllJets = " << BTagWeightMCAllJets << " and BTagWeightDataAllJets = " << BTagWeightDataAllJets << std::endl;

            if(i == BLeptIndex || i == BHadrIndex || i == Light1Index || i == Light2Index){
              if(bTagCSV[i] >= 0.898){
                BTagWeightMCSelJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]); 
                BTagWeightDataSelJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag); 
                //BTagWeightMCAllJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]);
                //BTagWeightDataAllJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag);
              }
              else{
                BTagWeightMCSelJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])); 
                BTagWeightDataSelJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag)); 
                //BTagWeightMCAllJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]));
                //BTagWeightDataAllJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag));
              }
            }
            else{
              if(bTagCSV[i] >= 0.898){
                BTagWeightMCNonSelJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]); 
                BTagWeightDataNonSelJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag); 
                //BTagWeightMCAllJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]);
                //BTagWeightDataAllJets *= bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag);
              }
              else{
                BTagWeightMCNonSelJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])); 
                BTagWeightDataNonSelJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag)); 
                //BTagWeightMCAllJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i]));
                //BTagWeightDataAllJets *= (1.0-bTagTool->getTagEff(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i])*bTagTool->getSF(selJets[i].Pt(), selJets[i].Eta(), partFlavour[i], "CSVT", syst_btag, syst_mistag));
              }
            }
          }*/
          histo1D["BTagWeight_SelJetsOnly"]->Fill(BTagWeightDataSelJets/BTagWeightMCSelJets);
          histo1D["BTagWeight_NonSelJetsOnly"]->Fill(BTagWeightDataNonSelJets/BTagWeightMCNonSelJets);
          histo1D["BTagWeight_AllJets"]->Fill(BTagWeightDataAllJets/BTagWeightMCAllJets);
        }

        if(BHadrIndex == 0 || BLeptIndex == 0){
          histo2D["bTagWeight_vs_BJetPt_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          if(partFlavAbs[0] == 5)      histo2D["bTagWeight_vs_BJetPt_partFlavB_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else if(partFlavAbs[0] == 4) histo2D["bTagWeight_vs_BJetPt_partFlavC_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else                              histo2D["bTagWeight_vs_BJetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
        }
        if(Light1Index == 0 || Light2Index == 0){
          histo2D["bTagWeight_vs_LightJetPt_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          if(partFlavAbs[0] == 5)      histo2D["bTagWeight_vs_LightJetPt_partFlavB_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else if(partFlavAbs[0] == 4) histo2D["bTagWeight_vs_LightJetPt_partFlavC_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
          else                         histo2D["bTagWeight_vs_LightJetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);
        }
        if(partFlavAbs[0] != 4 && partFlavAbs[0] != 5) histo2D["bTagWeight_vs_JetPt_partFlavLight_"+dsName]->Fill(selJets[0].Pt(),BTagWeight);

        histo1D["partonFlav_BJets_"+dsName]->Fill(partFlavAbs[BLeptIndex]);
        histo1D["partonFlav_BJets_"+dsName]->Fill(partFlavAbs[BHadrIndex]);
        histo1D["partonFlav_LightJets_"+dsName]->Fill(partFlavAbs[Light1Index]);
        histo1D["partonFlav_LightJets_"+dsName]->Fill(partFlavAbs[Light2Index]);
        scaleFactor = scaleFactor * BTagWeight;

        if(partFlavour[BLeptIndex] == 21) MSPlot["BJetPt_partFlavGluon_BTagSF"+leptChannel]->Fill( selJets[BLeptIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(partFlavour[BHadrIndex] == 21) MSPlot["BJetPt_partFlavGluon_BTagSF"+leptChannel]->Fill( selJets[BHadrIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        for(int ii = 0; ii < selJets.size(); ii++){
          if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_BTagSF"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        }
        MSPlot["nBJets_BTagSF"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_BTagSF"+leptChannel]->Fill(bTagStudy.getNrLightJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nPV_BTagSF"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_LeadingJet_BTagSF"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_SecondJet_BTagSF"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_ThirdJet_BTagSF"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetPt_FourthJet_BTagSF"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(selJets.size() > 4) MSPlot["JetPt_FifthJet_BTagSF"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonPt_BTagSF"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["JetEta_LeadingJet_BTagSF"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonEta_BTagSF"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["LeptonCharge_BTagSF"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

        MSPlot["CSVDistr_BJets"+leptChannel]->Fill(bTagCSV[BHadrIndex], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_BJets"+leptChannel]->Fill(bTagCSV[BLeptIndex], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_LightJets"+leptChannel]->Fill(bTagCSV[Light1Index], datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["CSVDistr_LightJets"+leptChannel]->Fill(bTagCSV[Light2Index], datasets[iDataSet], true, Luminosity*scaleFactor);

        //Check whether the Mlb, Mqqb, Mt and MW depends a lot on the considered b-tag and on the fact whether it is applied!
        if(isTTbarSemiLept){
          if(corrJetCombi[0] != 9999 && corrJetCombi[1] != 9999 && corrJetCombi[2] != 9999 && corrJetCombi[3] != 9999){
            histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selLepton+selJets[corrJetCombi[0]]).M());
            histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selJets[corrJetCombi[1]]+selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M() );
            histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill( (selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M());
            histo2D["MW_vs_MTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fill((selJets[corrJetCombi[1]]+selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M(),(selJets[corrJetCombi[2]]+selJets[corrJetCombi[3]]).M());
          }
        }

        //Identical MSPlots with number of jets information after requiring at least two b-jets and at least 2 light jets!
        MSPlot["nSelectedJets_"+bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( selJets.size(),                    datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nBTaggedJets_"+ bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["nLightJets_"+   bTitle[ibTag]+"_AfterBTag"+leptChannel]->Fill( bTagStudy.getNrLightJets(ibTag),   datasets[iDataSet], true, Luminosity*scaleFactor);      

        //Check the reconstruction fraction for different chi-sq values
        if(isTTbarSemiLept && corrJetCombi[2] != 9999 && corrJetCombi[3] != 9999){
          if(bTagStudy.getMlbMqqbChiSq(ibTag) < 5){ 
            chiSqLowerThan5++;
            if( (Light1Index == corrJetCombi[2] && Light2Index == corrJetCombi[3]) || (Light1Index == corrJetCombi[3] && Light2Index == corrJetCombi[2]) ) chiSqLowerThan5Correct++;

          }
          if(bTagStudy.getMlbMqqbChiSq(ibTag) < 10){
            chiSqLowerThan10++;
            if( (Light1Index == corrJetCombi[2] && Light2Index == corrJetCombi[3]) || (Light1Index == corrJetCombi[3] && Light2Index == corrJetCombi[2]) ) chiSqLowerThan10Correct++;
          }
          else{
            chiSqLowerThan50++;
            if( (Light1Index == corrJetCombi[2] && Light2Index == corrJetCombi[3]) || (Light1Index == corrJetCombi[3] && Light2Index == corrJetCombi[2]) ) chiSqLowerThan50Correct++;
          }
        }

        bool CutsSurvived = extraEvtSelCuts.KeepEvent(selLepton, selJets, selJetCombi, bTagStudy.getMlbMqqbChiSq(ibTag), CWUIndex, decayCh);
        MSPlot["MlbMqqbChiSq"+leptChannel]->Fill( bTagStudy.getMlbMqqbChiSq(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        if(bTagStudy.getMlbMqqbChiSq(ibTag) < 10) MSPlot["MlbMqqbChiSq_OnlyChiSq"+leptChannel]->Fill( bTagStudy.getMlbMqqbChiSq(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["HadrTopMass"+leptChannel]->Fill( (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor);
        MSPlot["HadrWMass"+leptChannel]->Fill( (selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor);

        //2D plots in order to understand whether there is a correlation between the cut-values and nPV, lepton and jet Pt!
        histo2D["MlbMqqbChiSq_vs_nPv_"+dsName]->Fill(nPrimVertices, bTagStudy.getMlbMqqbChiSq(ibTag));
        histo2D["HadrTopMass_vs_nPv_"+dsName]->Fill(nPrimVertices, (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
        histo2D["HadrWMass_vs_nPv_"+dsName]->Fill(nPrimVertices, (selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
        histo2D["MlbMqqbChiSq_vs_LeptPt_"+dsName]->Fill(selLepton.Pt(), bTagStudy.getMlbMqqbChiSq(ibTag));
        histo2D["HadrTopMass_vs_LeptPt_"+dsName]->Fill(selLepton.Pt(), (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
        histo2D["HadrWMass_vs_LeptPt_"+dsName]->Fill(selLepton.Pt(), (selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
        for(int ii = 0; ii < 4; ii++){
          histo2D["MlbMqqbChiSq_vs_JetPt_"+dsName]->Fill(selJets[selJetCombi[ii]].Pt(), bTagStudy.getMlbMqqbChiSq(ibTag));
          histo2D["HadrTopMass_vs_JetPt_"+dsName]->Fill(selJets[selJetCombi[ii]].Pt(), (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
          histo2D["HadrWMass_vs_JetPt_"+dsName]->Fill(selJets[selJetCombi[ii]].Pt(), (selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M());
        }
    
        //Also understand the correlation between the three cut variables
        histo2D["MlbMqqbChiSq_vs_HadrTopMass_"+dsName]->Fill((selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), bTagStudy.getMlbMqqbChiSq(ibTag) );
        histo2D["MlbMqqbChiSq_vs_HadrWMass_"+dsName]->Fill((selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), bTagStudy.getMlbMqqbChiSq(ibTag) );
        histo2D["HadrTopMass_vs_HadrWMass_"+dsName]->Fill((selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M() );

        //Now get the MSPlots after this additional cuts
        if( CutsSurvived ){
          if(decayCh == 0) nSelectedMu += 1;
          else if(decayCh == 1) nSelectedEl += 1;

          MSPlot["MlbMqqbChiSq_AllCuts"+leptChannel]->Fill( bTagStudy.getMlbMqqbChiSq(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["HadrTopMass_AllCuts"+leptChannel]->Fill( (selJets[selJetCombi[1]]+selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["HadrWMass_AllCuts"+leptChannel]->Fill( (selJets[selJetCombi[2]]+selJets[selJetCombi[3]]).M(), datasets[iDataSet], true, Luminosity*scaleFactor);

          MSPlot["nLightJets_AllCuts"+leptChannel]->Fill(bTagStudy.getNrLightJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["nBJets_AllCuts"+leptChannel]->Fill(bTagStudy.getNrBTaggedJets(ibTag), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["nPV_AllCuts"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["nPV_AllCuts_noPU"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, Luminosity*scaleFactor/lumiWeight);
          if(!isData)MSPlot["nPV_AllCuts_noLumi"+leptChannel]->Fill( nPrimVertices, datasets[iDataSet], true, scaleFactor/(datasets[iDataSet]->NormFactor()));
          SFcount += scaleFactor;
          MSPlot["JetPt_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetPt_LeadingJet_AllCuts_noBTagSF"+leptChannel]->Fill( selJets[0].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor/BTagWeight);
          MSPlot["JetPt_SecondJet_AllCuts"+leptChannel]->Fill( selJets[1].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetPt_ThirdJet_AllCuts"+leptChannel]->Fill( selJets[2].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetPt_FourthJet_AllCuts"+leptChannel]->Fill( selJets[3].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          if(selJets.size() > 4) MSPlot["JetPt_FifthJet_AllCuts"+leptChannel]->Fill( selJets[4].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetMass_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].M(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonPt_AllCuts"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonPt_AllCuts_noMuSF"+leptChannel]->Fill( selLepton.Pt(), datasets[iDataSet], true, Luminosity*scaleFactor/leptonSF);
          MSPlot["LeptonMass_AllCuts"+leptChannel]->Fill( selLepton.M(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["JetEta_LeadingJet_AllCuts"+leptChannel]->Fill( selJets[0].Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonEta_AllCuts"+leptChannel]->Fill( selLepton.Eta(), datasets[iDataSet], true, Luminosity*scaleFactor);
          MSPlot["LeptonCharge_AllCuts"+leptChannel]->Fill( leptCharge, datasets[iDataSet], true, Luminosity*scaleFactor);

          if(partFlavour[BLeptIndex] == 21) MSPlot["BJetPt_partFlavGluon_AllCuts"+leptChannel]->Fill( selJets[BLeptIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          if(partFlavour[BHadrIndex] == 21) MSPlot["BJetPt_partFlavGluon_AllCuts"+leptChannel]->Fill( selJets[BHadrIndex].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          for(int ii = 0; ii < selJets.size(); ii++){
            if(partFlavour[ii] == 21) MSPlot["JetsPt_partFlavGluon_AllCuts"+leptChannel]->Fill(selJets[ii].Pt(), datasets[iDataSet], true, Luminosity*scaleFactor);
          }

          if(isTTbarSemiLept && fillTFAfterCuts && createTFTree){  //Combine muon and electron channel!
            tfLight_mu = new TFLight();
            
            double fullScaleFactor = scaleFactor;   
            tfLight_mu->setFullScaleFactor(fullScaleFactor);
            tfLight_mu->setSelectedJets(selJets);
            tfLight_mu->setSelectedLepton(selLepton);
            tfLight_mu->setCorrectJetCombi(corrJetCombi);
	
            //Store the information needed for the TF (but only has value when dataset is ttbar)
            tfLight_mu->setGenVectorLight1( light->genVectorLight1() );
            tfLight_mu->setGenVectorLight2( light->genVectorLight2() );
            tfLight_mu->setGenVectorHadrB(  light->genVectorHadrB()  );
            tfLight_mu->setGenVectorLeptB(  light->genVectorLeptB()  );
            tfLight_mu->setGenVectorLepton( light->genVectorLepton() );

            TFLightTree->Fill();
            delete tfLight_mu;
            //----  End of Tree file filling (for TF's after evtSel)  ----//
          }
          /*else if(isTTbarSemiLept && decayCh == 1 && createTFTree){
           tfLight_el = new TFLight();
            
            double fullScaleFactor = scaleFactor;   
            tfLight_el->setFullScaleFactor(fullScaleFactor);
            tfLight_el->setSelectedJets(selJets);
            tfLight_el->setSelectedLepton(selLepton);
            tfLight_el->setDecayChannel(decayCh);
            tfLight_el->setLeptonCharge(leptCharge);
            tfLight_el->setCorrectJetCombi(corrJetCombi);
            tfLight_el->setMET(MET);
	
            //Store the information needed for the TF (but only has value when dataset is ttbar)
            tfLight_el->setGenVectorLight1( light->genVectorLight1() );
            tfLight_el->setGenVectorLight2( light->genVectorLight2() );
            tfLight_el->setGenVectorHadrB(  light->genVectorHadrB()  );
            tfLight_el->setGenVectorLeptB(  light->genVectorLeptB()  );
            tfLight_el->setGenVectorLepton( light->genVectorLepton() );

            TFLightTree->Fill();
            delete tfLight_el;
            //----  End of Tree file filling (for TF's after evtSel)  ----//
          }*/
        }
 
        //Write out the LHCO output!
        if( (getLHCOOutput || getEvtNrInfo) && bTagChoiceMade){

          //Keep track of the original event number and the madweight numbers:
          EvtNrMatching << "\n    " << iEvt;
          if( iEvt < 10) EvtNrMatching << " "; if( iEvt < 100) EvtNrMatching << " "; if( iEvt < 1000) EvtNrMatching << " "; if( iEvt < 10000) EvtNrMatching << " "; if( iEvt < 100000) EvtNrMatching << " ";
          EvtNrMatching << "        " << CutsSurvived << "     ";

          //if(CutsSurvived) histo1D["CosTheta_SelEvts"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
          // --> Cannot add this since the neutrino is not completely reconstructed and thus the leptonic top is not known ...
          if(CutsSurvived){
            lhcoOutput.StoreRecoInfo(selLepton, selJets, selJetCombi, decayCh, leptCharge, EvtNrMatching, CWUIndex);
            EvtNrMatching << "            " << scaleFactor;
          }
          
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 0) histo1D["CosTheta_Gen_SelEvts_"+dsName+"_mu"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0 && decayCh == 1) histo1D["CosTheta_Gen_SelEvts_"+dsName+"_el"]->Fill(genCosTheta);
          if(CutsSurvived && genCosTheta != 2.0) histo1D["CosTheta_Gen_SelEvts_"+dsName]->Fill(genCosTheta);
        }
        else if((getLHCOOutput || getEvtNrInfo) && !bTagChoiceMade){
          cout << " ERROR : Not possible to write the lhco file when all b-tag options are still being considered!!  ==> Setting boolean to FALSE !" << endl;
          getLHCOOutput = false;
        }
      }

    }//End of loop on events
    std::cout<<"    Processed all "<<nrEvts<<" events  --> # selected: "<<nSelectedMu<<" (mu+jets) "<<nSelectedEl<<" (e+jets)"<< flush<<"\n";

    //------------------------------//
    //  Calculate the ChiSq values  //
    //------------------------------//
    if(isTTbarSemiLept && getMassFits){
      std::cout << " -----------------     Fitting the different mass distributions    --------------------- \n" << endl;   
      histo1D["Mlb_CorrTT_"+dsName]->Fit("gaus","Q");     
      histo1D["HadrMTop_CorrTT_"+dsName]->Fit("gaus","Q");
      histo1D["HadrMW_CorrTT_"+dsName]->Fit("gaus","Q");
      histo1D["Mlb_CorrTT_"+dsName]->Fit("gaus","Q","", histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                        histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      histo1D["HadrMTop_CorrTT_"+dsName]->Fit("gaus","Q","", histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                             histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      histo1D["HadrMW_CorrTT_"+dsName]->Fit("gaus","Q","",histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2),
                                                          histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2));
      
      //Write out the mass values!
      std::cout << "   ** Mlb  = " << histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)      << " +- " << histo1D["Mlb_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** Mqqb = " << histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["HadrMTop_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << endl;
      std::cout << "   ** MW   = " << histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["HadrMW_CorrTT_"+dsName]->GetFunction("gaus")->GetParameter(2) << std::endl;
    
      //Output for the different b-tags
      for(int ibTag = 0; ibTag < NrBTags; ibTag++){
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q"); histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q");
        histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q","", 
                                                      histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                                      histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                          histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                          histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
        histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->Fit("gaus","Q","",
                                               histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)-histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2),
                                               histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)+histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2));
      
        //Write out the mass values!
        cout << "\n   ** Mlb  -- " << bTitle[ibTag] << " = " << histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)    << " +- " << histo1D["Mlb_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** Mqqb -- " << bTitle[ibTag] << " = " << histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1) << " +- " << histo1D["HadrMTop_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << endl;
        cout << "   ** MW   -- " << bTitle[ibTag] << " = " << histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(1)   << " +- " << histo1D["HadrMW_CorrTT_"+dsName+"_"+bTitle[ibTag]]->GetFunction("gaus")->GetParameter(2) << std::endl;
      }
    }

    //Some output on the b-tag efficiencies
    if(isTTbarSemiLept){
      cout << " Mean value obtained for the b-tag efficiency for jets with parton flavour b is : " << (double)(tagEffPartFlavB/nrPartFlavBJets) << " ("<< tagEffPartFlavB << "/" << nrPartFlavBJets << ")" << endl;
      cout << " Mean value obtained for the b-tag efficiency for jets with parton flavour udscg is : " << (double)(tagEffPartFlavNotB/nrPartFlavNotBJets) << " ("<< tagEffPartFlavNotB << "/" << nrPartFlavNotBJets << ")" << endl;
      cout << " Mean value obtained for the b-tag efficiency for jets with parton flavour c is : " << (double)(tagEffPartFlavC/nrPartFlavCJets) << " ("<< tagEffPartFlavC << "/" << nrPartFlavCJets << ")" << endl;
      cout << " Mean value obtained for the b-tag efficiency for jets with parton flavour udsg is : " << (double)(tagEffPartFlavL/nrPartFlavLJets) << " ("<< tagEffPartFlavL << "/" << nrPartFlavLJets << ")" << endl;

      cout << "\n ******* Some information about the light jet assignment ******* " << endl;
      cout << "   * Number of events with the light jets matched : " << lightIsMatched << endl;
      cout << "    - How often is the light jet the leading jet : " << lightIsFirst << endl;
      cout << "    - How often is the light jet the second jet : " << lightIsSecond << endl;
      cout << "   * How often are more than 2 light jets present in the event: " << evtHasThreeLights << endl;
      cout << "    - How often is the light jet the third jet : " << lightIsThird << endl;
  
      cout << "\n ******* First light jet separate ******* " << endl;
      cout << "   * Number of events with first light jet matched : " << firstIsMatched << endl;
      cout << "     - How often is the first light jet the leading one : " << firstLightIsFirst << endl;
      cout << "     - How often is the first light jet the second  one : " << firstLightIsSecond << endl;
      cout << "     - How often is the first light jet another one     : " << firstLightIsElse << endl;
      cout << "      * Number of events with a third light jet : " << firstLightMatchedThreeExist << endl;
      cout << "        - How often is this third light jet the correct one : " << firstLightIsThird << endl;

      cout << "\n ******** Some percentages ******** " << endl;
      cout << "   * How often is when the first jet is chosen this also correct  : " << (float)(firstJetChosenAndCorrect)/(float)(firstJetChosen) << endl;
      cout << "   * How often is when the second jet is chosen this also correct : " << (float)(secondJetChosenAndCorrect)/(float)secondJetChosen << endl;
      cout << "   * How often is when the third jet is chosen this also correct  : " << (float)(thirdJetChosenAndCorrect)/(float)thirdJetChosen << endl;

      cout << "\n ******** Reconstruction efficiency for different chi-sq bins ********** " << endl;
      cout << "   * Fraction of correct light jets for chisq < 5  : " << (float)(chiSqLowerThan5Correct)/(float)(chiSqLowerThan5) << " (Nr events is " << chiSqLowerThan5 << ") " << endl;
      cout << "   * Fraction of correct light jets for chisq < 10 : " << (float)(chiSqLowerThan10Correct)/(float)(chiSqLowerThan10) << " (Nr events is " << chiSqLowerThan10 << ") " << endl;
      cout << "   * Fraction of correct light jets for chisq < 50 : " << (float)(chiSqLowerThan50Correct)/(float)(chiSqLowerThan50) << " (Nr events is " << chiSqLowerThan50 << ") " << endl;
    }
    cout << "    --> Overall event count using the SF (no lumi) is : " << SFcount << endl;

    //--- Get output from bTagStudy class ---//
    if(isTTbarSemiLept) bTagStudy.ReturnBTagTable();
    bTagStudy.CreateHistograms(outputFile, savePDF, pathPNG, iDataSet);  //Security is added inside class such that MSPlots are only written when all datasets are considered!

    //--- Get output from LHCOOutput class ---//
    if(getLHCOOutput && bTagChoiceMade) lhcoOutput.WriteLHCOPlots(outputFile);

    //--- Get output form the ExtraEvtSel class ---//
    extraEvtSelCuts.StoreCutInfluence(outputFile);

    //---- Close the EventNrMatching output file for the considered dataset  ----//
    if(bTagChoiceMade && (getLHCOOutput || getEvtNrInfo)) EvtNrMatching.close();

    //---- Store the dataset specific 1D and 2D histograms ----//
    if(histo1D.size() > 0){
      TDirectory* th1dir = outputFile->GetDirectory("1D_histograms");   //Check whether directory already exists ..
      if(!th1dir) th1dir = outputFile->mkdir("1D_histograms");          // .. and otherwise create it!
      th1dir->cd();
      TDirectory* ds1Dir = th1dir->mkdir(dsName.c_str());
      ds1Dir->cd();
      for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
      }
    }

    if(histo2D.size() > 0){
      TDirectory* th2dir = outputFile->GetDirectory("2D_histograms_graphs");
      if(!th2dir) th2dir = outputFile->mkdir("2D_histograms_graphs");
      th2dir->cd();
      TDirectory* ds2Dir = th2dir->mkdir(dsName.c_str());
      ds2Dir->cd();
      for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
        TH2F *temp = it->second;
        temp->Write();
      }
    }
    histo1D.clear();
    histo2D.clear();

    //------ Store the TF tree ------//
    if(isTTbarSemiLept && createTFTree){
      TFLightFile->cd();
      TTree* configTreeTFLightFile = new TTree("configTreeTFLightFile","configuration Tree in TFLight File");
      TClonesArray* tcdatasettflightfile = new TClonesArray("Dataset",1);
      configTreeTFLightFile->Branch("Dataset","TClonesArray",&tcdatasettflightfile);
      new ((*tcdatasettflightfile)[0]) Dataset(*datasets[iDataSet]);
      //TClonesArray* tcanaenvlightfile = new TClonesArray("AnalysisEnvironment",1);   --> Needed?
      //configTreeLightFile->Branch("AnaEnv","TClonesArray",&tcanaenvlightfile);
      //new ((*tcanaenvlightfile)[0]) AnalysisEnvironment(anaEnv);

      configTreeTFLightFile->Fill();
      configTreeTFLightFile->Write();
      TFLightTree->Write();
      TFLightFile->Close();
      //----  End of storing Tree  ----//
    }
    delete TFLightFile;  //Will need to delete this for each dataset since it is initialized for each of them!

    inputFile->Close();
    delete inputFile;

    if(isTTbarSemiLept) std::cout << "\n Number of events with both b-jets actually parton flavour = b is : " << BothB << std::endl;
  }//End of loop on datasets

  //Store the b-tag histograms:
  if(!bTagPlotsMade)
    bTagTool->WriteMCEfficiencyHistos(bTagPlotsOutput); //"PersonalClasses/Calibrations/BTagSF/BTagWeightPlots_CSVT_noTTbar.root");

  /////////////////////////
  // Write out the plots //
  /////////////////////////
  outputFile -> cd();
  mkdir((pathPNG+"/MSPlots").c_str(),0777);

  if(verbosity > 0) std::cout << "\n - Making the plots in the LightAnomCoupAnalyzer file: " << outputFile->GetName() << std::endl;
  if(MSPlot.size() > 0 && !onlyTTbar){
    TDirectory* msdir = outputFile->mkdir("MSPlots");
    msdir->cd(); 
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
      MultiSamplePlot *temp = it->second;
      string name = it->first;
      temp->Draw(name, 1, false, false, false, 1);
      temp->Write(outputFile, name, savePDF, (pathPNG+"/MSPlots/").c_str(), "pdf");
    }
  }

  delete leptonTools;
  delete outputFile;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
