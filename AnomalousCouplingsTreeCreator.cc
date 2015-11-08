#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "PersonalClasses/Style.C"                                                 //CHECK if this works!
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TLorentzVector.h"

//Specific code for anomalous couplings analysis:
#include "AnomalousCouplings/PersonalClasses/interface/LHCOOutput.h"
#include "AnomalousCouplings/PersonalClasses/interface/AnomCoupLight.h"
#include "AnomalousCouplings/PersonalClasses/interface/KinematicFunctions.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

int main (int argc, char *argv[]){

  TFile *fout = new TFile ("PlotsMacro/AnomCouplings.root", "RECREATE");  
  clock_t start = clock();
  
  cout << "*************************************************************" << endl;
  cout << " Beginning of the program for creating the AnomCoupl Trees ! " << endl;
  cout << "*************************************************************" << endl;
  
  //SetStyle if needed
  //setTDRStyle(); 
  setMyStyle();

  /////////////////////
  // Configuration
  /////////////////////

  //Nr events
  int NrEvts = -1;
  if(argc >= 2 && string(argv[1]) != "-1"){
    NrEvts = atoi(argv[1]);
    cout << " ** Will only look at " << NrEvts << " events ! " << endl;
  }

  //xml file
  string xmlFileName ="PersonalClasses/config/myAnomCouplConfig.xml";
  if (argc >= 3)                            //--> Needed for the PBS script!!
    xmlFileName = string(argv[2]);
  const char *xmlfile = xmlFileName.c_str();
  cout << " ** Used config file: " << xmlfile << endl;

  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);

  //////////////////////////
  // Verbosity for output //
  //////////////////////////
  int verbosity                 = 1;
  //0 muet
  //1 Main Info
  //2 mcParticlesMatchin Info
  //3 
  //4 Info for each event
  //5 Debug
 
  ///////////////////////////////
  //  Run specific parts only  //
  ///////////////////////////////
  bool getLHCOOutput = false;
  bool getEventInfo = true;
  bool saveAsPDF = false;
  bool splitLeptonCharge = false;
  bool getCorrectWrong = false;

  /////////////////////////////
  //   Which systematics ?   //
  /////////////////////////////
  std::string systematic = "Nominal";
  if(argc >= 4)
    systematic = string(argv[3]);
  std::cout << " ** Will be using systematic : " << systematic << std::endl;
  if(systematic != "Nominal" && systematic != "JESPlus" && systematic != "JESMinus" && systematic != "JERPlus" && systematic != "JERMinis" && systematic != "bTagPlus" && systematic != "bTagMinus"){
    std::cout << " *****  Given systematic is not allowed! " << std::endl;
    std::cout << " *****  Possibilities are : Nominal, JESPlus, JESMinus, JERPlus, JERMinus, bTagPlus, bTagMinus " << std::endl;
    exit(-1);
  }
 
  ////////////////////////////////////
  /// AnalysisEnvironment  
  ////////////////////////////////////
  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
 
  /////////////////////
  // Load Datasets
  /////////////////////
  TTreeLoader treeLoader;
  cout << " - Load datasets ..." << endl;
  vector < Dataset* > datasets;
  vector < Dataset* > datasetsMu;
  vector < Dataset* > datasetsEl;
  vector < Dataset* > datasetsPlot;
  treeLoader.LoadDatasets (datasets, xmlfile);
  for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  //Set luminosity correct!
  bool foundMu = false,  foundEl = false;    //Boolean used to differentiate between muon-data and electron-date in case the two channels are ran at the same time!
  float Luminosity = oldLuminosity, LuminosityMu = oldLuminosity, LuminosityEl = oldLuminosity;
  for (unsigned int d = 0; d < datasets.size (); d++) {
    
    if(Luminosity > datasets[d]->EquivalentLumi() ){
      Luminosity = datasets[d]->EquivalentLumi();
      LuminosityMu = datasets[d]->EquivalentLumi(); LuminosityEl = datasets[d]->EquivalentLumi(); //This way EventSelectionTable info is set correctly !
    }
    string dataSetName = datasets[d]->Name();
    
    if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
      LuminosityMu = datasets[d]->EquivalentLumi();
      foundMu=true;
    }  
    if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
      LuminosityEl = datasets[d]->EquivalentLumi();
      foundEl=true;  
    }  
    
    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // when only using one channel
      Luminosity = datasets[d]->EquivalentLumi();
    }   

    if(dataSetName.find("QCD") == 0)                                      datasets[d]->SetColor(kYellow);
    if(dataSetName.find("TT") == 0)                                       datasets[d]->SetColor(kRed+1);
    if(dataSetName.find("TTbarJets_Other") == 0)                          datasets[d]->SetColor(kRed-7);
    if(dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") == 0) datasets[d]->SetColor(kMagenta);
    if(dataSetName.find("WJets") == 0){                                   datasets[d]->SetColor(kGreen-3); datasets[d]->SetTitle("W#rightarrowl#nu");}      
    if(dataSetName.find("ZJets") == 0){                                   datasets[d]->SetColor(kAzure-2); datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");}
  }
  
  if(!foundMu && !foundEl && Luminosity != oldLuminosity) cout << "   --> Changed analysis environment luminosity to "<< Luminosity << endl;
  else {
    if(LuminosityMu != oldLuminosity) cout << "   --> Muon PD: changed analysis environment luminosity to "<< LuminosityMu << endl;
    if(LuminosityEl != oldLuminosity) cout << "   --> Electron PD: changed analysis environment luminosity to "<< LuminosityEl << endl;
  }
  
  //Now that it has been identified whether the two channels are ran at the same time, create a separate datasets vector not containing the other DATA sample!
  if (foundMu) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) ) 
	datasetsMu.push_back(datasets[d]);
    }
  }
  if (foundEl) {
    for (unsigned int d = 0; d < datasets.size (); d++) {
      string dataSetName = datasets[d]->Name();
      if ( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) ) 
	datasetsEl.push_back(datasets[d]);
    }
  }

  ////////////////////////////////////
  /// Normal Plots (TH1F* and TH2F*)
  ////////////////////////////////////
  //All histograms can be defined as pre-programmed maps which makes definitions and looping easier
  map<string,TH1F*> histo1D;     
  map<string,TH2F*> histo2D;  
  
  histo1D["Pt_genEvtMuon"] = new TH1F("Pt_genEvtMuon","Pt distribution for genEvt muon",100,0,150);
  histo1D["Mass_genEvtMuon"] = new TH1F("Mass_genEvtMuon","Mass distribution for genEvt muon",100,0,0.2);
  histo1D["Pt_genEvtElec"] = new TH1F("Pt_genEvtElec","Pt distribution for genEvt electron",100,0,150);
  histo1D["Mass_genEvtElec"] = new TH1F("Mass_genEvtElec","Mass distribution for genEvt electron",100,0,0.001);

  histo1D["StCosTheta"] = new TH1F("StCosTheta","StCosTheta",200,-1,1);
  histo1D["CorrectBLeptCSVDiscr"] = new TH1F("CorrectBLeptCSVDiscr","CorrectBLeptCSVDiscr",400,-2.5,1.5);
  histo1D["CorrectBHadrCSVDiscr"] = new TH1F("CorrectBHadrCSVDiscr","CorrectBHadrCSVDiscr",400,-2.5,1.5);
  histo1D["CorrectQuark1CSVDiscr"] = new TH1F("CorrectQuark1CSVDiscr","CorrectQuark1CSVDiscr",400,-2.5,1.5);
  histo1D["CorrectQuark2CSVDiscr"] = new TH1F("CorrectQuark2CSVDiscr","CorrectQuark2CSVDiscr",400,-2.5,1.5);

  histo1D["Quark1JetNumber"] = new TH1F("Quark1JetNumber","Quark1JetNumber",12,-1.5,10.5);
  histo1D["Quark2JetNumber"] = new TH1F("Quark2JetNumber","Quark2JetNumber",12,-1.5,10.5);
  histo1D["nSelectedEvents"] = new TH1F("nSelectedEvents","nSelectedEvents",8, 0.5, 8.5);
  std::string BinLabel[8] = {"#mu+ selected","#mu- selected","#mu selected","gen #mu identified","e+ selected","e- selected","e selected","gen e identified"};
  for(int iBin = 1; iBin <= 8; iBin++)
    histo1D["nSelectedEvents"]->GetXaxis()->SetBinLabel(iBin,BinLabel[iBin-1].c_str());

  ////////////////////////////////////
  /// MultiSamplePlot                   --> Should be removed here since localgrid running is not compatible with MSPlots!
  ////////////////////////////////////
  map<string,MultiSamplePlot*> MSPlot;
  /*
  MSPlot["InitJets_pT_jet1_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet1_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet2_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet2_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet3_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet3_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");
  MSPlot["InitJets_pT_jet4_beforeEvtSel"] = new MultiSamplePlot(datasets, "InitJets_pT_jet4_beforeEvtSel", 60, -150, 650, "p_{T} (GeV)");

  MSPlot["InitJets_METPt"] = new MultiSamplePlot(datasets, "InitJets_METPt", 60,0,300,"p_{T} (GeV) ");
  MSPlot["InitJets_Pt_jet1"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet1",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet2"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet2",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet3"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet3",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_Pt_jet4"] = new MultiSamplePlot(datasets, "InitJets_Pt_jet4",60,0,400,"p_{T} (GeV)");
  MSPlot["InitJets_METPt_METTypeOneCorrected"] = new MultiSamplePlot(datasets, "InitJets_METPt_METTypeOneCorrected", 60,0,300, "p_{T} (GeV)");
  MSPlot["InitJets_METPt_JerSmearingApplied"] = new MultiSamplePlot(datasets, "InitJets_METPt_JerSmearingApplied", 60,0,300, "p_{T} (GeV)");

  MSPlot["nEventsAfterCutsSemiMu"] = new MultiSamplePlot(datasets, "nEventsAfterCutsSemiMu",15, -0.5, 14.5, "#events after each cut");

  //string leptFlavs[3]={"_other","_mu","_el"};
  string leptFlavs[2]={"_mu","_el"};
  for(int ii = 0; ii < 2; ii++){
    string leptFlav = leptFlavs[ii];
    MSPlot["Selected_Events_pT_jet1"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet1"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_jet2"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet2"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_jet3"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet3"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_jet4"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_jet4"+leptFlav, 60, -150, 650, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_4leadingjets"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_4leadingjets"+leptFlav,60, 0, 600, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_alljets"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_alljets"+leptFlav, 60, 0, 600, "p_{T} (GeV)");
    MSPlot["Selected_Events_pT_lepton"+leptFlav] = new MultiSamplePlot(datasets, "Selected_Events_pT_lepton"+leptFlav,150,-100,350,"p_{t} (GeV)");
  }

  MSPlot["lumiWeights"] = new MultiSamplePlot(datasets, "lumiWeights", 200,-1,2, "lumiWeight obtained from nTruePU");
  MSPlot["nTruePU"] = new MultiSamplePlot(datasets, "nTruePU",100, 0, 75,"number of true PU interactions");
  */  

  ////////////////////////////////////
  /// Selection table
  ////////////////////////////////////
  vector<string> CutsSelecTableSemiMu;
  CutsSelecTableSemiMu.push_back(string("Preselected"));
  CutsSelecTableSemiMu.push_back(string("Trigged"));
  CutsSelecTableSemiMu.push_back(string("Good PV"));
  CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  vector<string> CutsSelecTableSemiEl;
  CutsSelecTableSemiEl.push_back(string("Preselected"));
  CutsSelecTableSemiEl.push_back(string("Trigged"));
  CutsSelecTableSemiEl.push_back(string("Good PV"));
  CutsSelecTableSemiEl.push_back(string("1 selected electron"));
  CutsSelecTableSemiEl.push_back(string("Veto muon"));
  CutsSelecTableSemiEl.push_back(string("Veto 2nd electron from Z-decay"));
  CutsSelecTableSemiEl.push_back(string("Conversion veto"));
  
  char LabelNJets[100];
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-3);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-2);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets-1);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));
  sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  CutsSelecTableSemiMu.push_back(string(LabelNJets));
  CutsSelecTableSemiEl.push_back(string(LabelNJets));

  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  selecTableSemiMu.SetLuminosity(LuminosityMu);
  SelectionTable selecTableSemiEl(CutsSelecTableSemiEl, datasets);
  selecTableSemiEl.SetLuminosity(LuminosityEl);

  if (verbose > 0) cout << " - SelectionTable instantiated ..." << endl;

  ////////////////////////
  // PileUp Reweighting //   --> Need this for the event selection table!
  ////////////////////////
  LumiReWeighting LumiWeights;
  LumiWeights = LumiReWeighting("PersonalClasses/Calibrations/PUReweighting/pileup_MC_Summer12_S10.root", "PersonalClasses/Calibrations/PUReweighting/pileup_2012Data53X_UpToRun208686_pixelcorrv2_Mu/nominal.root", "pileup", "pileup");
  cout << " - LumiReWeighting instantiated ... " << endl;

  ////////////////////////////////////
  //	Loop on datasets
  ////////////////////////////////////
  if (verbose > 0) cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  
  for (unsigned int d = 0; d < datasets.size (); d++) {
    if (verbose > 1) cout << "   * Dataset " << d << ": " << datasets[d]->Name () << " with " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    int iFile = -1;
    string previousFilename = "", dataSetName = datasets[d]->Name();
    
    int nSelectedMuPos = 0, nSelectedMuNeg = 0;
    int nSelectedElPos = 0, nSelectedElNeg = 0;
    int nGenMu = 0, nGenEl = 0;
    
    //open files and load
    treeLoader.LoadDataset (datasets[d], anaEnv);
    bool isData = false;
    if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
      isData = true;
    
    /////////////////////////////////////
    /// Initialize JEC factors            --> Updated on 3/11/2015 (Recommendations from: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC)
    /////////////////////////////////////
    vector<JetCorrectorParameters> vCorrParam;
    if( isData ){// Data!
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_DATA_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_DATA_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_DATA_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_DATA_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else{
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_MC_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_MC_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("PersonalClasses/Calibrations/JECFiles/Summer13_V4_MC_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("PersonalClasses/Calibrations/JECFiles/Summer13_V4_MC_Uncertainty_AK5PFchs.txt");

    //True means also redo L1        
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
   
    /////////////////////
    //  Used classes   //
    /////////////////////
    LHCOOutput lhcoOutput(verbose, getLHCOOutput, splitLeptonCharge, getCorrectWrong); 
    if(dataSetName.find("TTbarJets") == 0) lhcoOutput.Initialize("Gen");
    KinematicFunctions kinFunctions;  //Variable accessible in KinematicFunctions using kinFunctions.CosTheta(TLorentzVector *Top, TLorentzVector *WLept, TLorentzVector *lepton)
    AnomCoupLight* anomCoupLight = 0;

    //Initialize LightTuple (AnomCoupTree) specific stuff:
    TFile* LightFile = new TFile(("LightTree/AnomCoupLight_"+dataSetName+"_"+systematic+".root").c_str(),"RECREATE");  
    TTree* LightTree = new TTree("LightTree","Tree containing the AnomCoup information");
    LightTree->Branch("TheAnomCoupLight","AnomCoupLight",&anomCoupLight);

    /////////////////////////////////////////
    //  LHCO Output files + GeneratorInfo  //
    /////////////////////////////////////////
    ofstream EventInfoFile;
    if(getEventInfo == true){	
      EventInfoFile.open("MadWeightInput/AnalyzerOutput/EventNumberInformation.lhco");
      EventInfoFile << " Event Number      Lepton Type       Event selection       selectedChannelNumber " << endl;
    }
    unsigned int NumberCorrectEvents = 0; //Counts the number of semi-leptonic events

    ////////////////////////////////////
    //	loop on events
    ////////////////////////////////////
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;

    unsigned int NrEvtsToRunOver = NrEvts;
    if(NrEvts == -1) NrEvtsToRunOver = datasets[d]->NofEvtsToRunOver();
    std::cout << "  ---> Will process " << NrEvtsToRunOver << " events for this dataset ! " << std::endl;
    for (unsigned int ievt = 0; ievt < NrEvtsToRunOver; ievt++){
      
      if(verbosity > 3) std::cout << " Looking at event : " << ievt << std::endl;    
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      vector < TRootGenJet* > genjets;
      vector<int> jetCombi(4,9999);;   //Define this here and just initialize to 9999 such that it also can be used for other datasets!
      
      if(ievt%1000 == 0)
	std::cout<<"    - Processing the "<<ievt<<"th event (" << ((double)ievt/(double)NrEvtsToRunOver)*100  << "%)" << " +> # selected: " << nSelectedMuPos+nSelectedMuNeg << " (mu+jets) " << nSelectedElPos+nSelectedElNeg << " (e+jets)" << flush<<"\r";
      
      ////////////////
      // LOAD EVENT //
      ////////////////
      //TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);  //Use uncorrected jets ...
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
      
      if( !isData ) {
        genjets = treeLoader.LoadGenJet(ievt,false);
        sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      vector<TRootMCParticle*> mcParticles;      
      if(dataSetName.find("TTbarJets") == 0){	
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);  
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }

      // check with genEvent which ttbar channel it is
      if(dataSetName.find("TTbarJets") == 0)  {
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	if( genEvt->isSemiLeptonic(TRootGenEvent::kMuon) ) {
          histo1D["Mass_genEvtMuon"]->Fill(genEvt->lepton().M());
          histo1D["Pt_genEvtMuon"]->Fill(genEvt->lepton().Pt());
          histo1D["nSelectedEvents"]->Fill(4);
          nGenMu++;
	}
	else if( genEvt->isSemiLeptonic(TRootGenEvent::kElec) ) {
          histo1D["Mass_genEvtElec"]->Fill(genEvt->lepton().M());
          histo1D["Pt_genEvtElec"]->Fill(genEvt->lepton().Pt());
          histo1D["nSelectedEvents"]->Fill(8);
          nGenEl++;
	}
      }
      
      /////////////////////////////////
      // DETERMINE EVENT SCALEFACTOR //
      /////////////////////////////////
      float scaleFactor = 1.;
      
      // Load the GenEvent and calculate the branching ratio correction
      /*if(dataSetName.find("TTbarJets") == 0){
        std::cout << " In this if loop for WJets ? --> " << dataSetName << std::endl;
	TRootGenEvent* genEvt = treeLoader.LoadGenEvent(ievt,false);
	if( genEvt->isSemiLeptonic() )      scaleFactor *= (0.108*9.)*(0.676*1.5);   --> Added this in given XS!
	else if( genEvt->isFullHadronic() ) scaleFactor *= (0.676*1.5)*(0.676*1.5);
	else if( genEvt->isFullLeptonic() ) scaleFactor *= (0.108*9.)*(0.108*9.); 
      }*/

      //------------------------//
      // Start with corrections //
      //------------------------//

      ////////////////////////////////////////
      //  Beam scraping and PU reweighting     --> Will be moved to nTuple analyzer!
      ////////////////////////////////////////
      double lumiWeight = 1;  
      if(!isData)   
        lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
      
      //////////////////////////////////////
      // Apply Jet Corrections on-the-fly //   
      //////////////////////////////////////
      jetTools->unCorrectMETTypeOne(init_jets, mets[0], isData);
      jetTools->correctJets(init_jets, event->kt6PFJets_rho(), isData);      //Why is this kt6PFJets rho information needed?
      jetTools->correctMETTypeOne(init_jets, mets[0], isData);
//      MSPlot["InitJets_METPt_METTypeOneCorrected"]->Fill(mets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor*lumiWeight);

      ///////////////////
      // TRIGGER SETUP //
      ///////////////////
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if(previousFilename != currentFilename){
	previousFilename = currentFilename;
	iFile++;
	cout<<"    - File changed!!! => iFile = "<<iFile<<endl;
      }
      
      int currentRun = event->runId();
      if(previousRun != currentRun) {
        previousRun = currentRun;
	
	//semi-mu
	if(dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
	  
	  if( event->runId() >= 190456 && event->runId() <= 190738 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
	  else if( event->runId() >= 190782 && event->runId() <= 193621)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
	  else if(event->runId() >= 193834  && event->runId() <= 196531 )
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
	  else if( event->runId() >= 198022  && event->runId() <= 199608)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
	  else if( event->runId() >= 199698 && event->runId() <= 209151)
	    itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
	  else
	    cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
	  if( itriggerSemiMu == 9999 ){
	    cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
	    exit(-1);
	  }
        }
        else{ 
	  itriggerSemiMu = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun); // Summer12 DR53X
	}
	
	// semi-electron
        if(dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0 ) {
	  
	  if( event->runId() <= 190738 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
	  else if( event->runId() >= 190782 && event->runId() <= 191411 )
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);
	  else if( event->runId() >= 191695 && event->runId() <= 196531)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
	  else if( event->runId() >= 198049 && event->runId() <= 208686)
	    itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
	  //else if( event->runId() > 208686)
	  //     itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
	  else { 
	    cout << "Unknown run for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(1);
	  }
	  if( itriggerSemiEl == 9999 ){	    
	    cout << "itriggerSemiEl == 9999 for SemiEl HLTpath selection: " << event->runId() << endl;
	    exit(-1);
	  }
        }
        else	  
	  itriggerSemiEl = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun); // Summer12 DR53X 
      }
      
      if (itriggerSemiMu == 9999 && itriggerSemiEl == 9999) {
	cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT IN RUN " << event->runId() << endl;
	exit(1);
      }
      
      ///////////////////////////////////////////////////////////////////////////////////////
      //   JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA         //
      //     - JES Corrections: The nominal corrections are already applied at PAT level   //
      //       so these tools should only be used for studies of the effect of systematics //
      ///////////////////////////////////////////////////////////////////////////////////////
      //JER Smearing:
      if(!isData) {
        if( systematic == "JERPlus")  jetTools->correctJetJER(init_jets, genjets, mets[0], "plus", false); //false means don't use old numbers but newer ones (~ 8TeV recommendation of 2014!) 
        if( systematic == "JERMinus") jetTools->correctJetJER(init_jets, genjets, mets[0], "minus", false);
        if( systematic == "Nominal")  jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
	if( systematic == "JESPlus")  jetTools->correctJetJESUnc(init_jets, mets[0], "plus", 1);  //last integer (1) = nSigma
	if( systematic == "JESMinus") jetTools->correctJetJESUnc(init_jets, mets[0], "minus", 1);
      }

      ////////////////////////////////////////////////////////
      // Access particle information before event selection //
      // Write this information to LHCO Output for MW       //
      ////////////////////////////////////////////////////////
      if(dataSetName.find("TTbarJets") == 0){
        lhcoOutput.StoreGenInfo(mcParticles);
        if(lhcoOutput.GenEventContentCorrect())
          histo1D["StCosTheta"]->Fill(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
      }

      //Accessing information and store in EventInfoFile    
      std::string leptonTypeString[5] = {"muPlus","muMinus","elPlus","elMinus","notIdentified"};      //Not possible to define this in header file of LHCOOutput ...
      
      if( getEventInfo == true){
        EventInfoFile << "     " << ievt << "         ";
        if( lhcoOutput.GenEventContentCorrect() ){
	  NumberCorrectEvents++;
          EventInfoFile << "     " << leptonTypeString[lhcoOutput.getLeptonType()] << "    ";
	}//End of Gen event with correct content    
        else    
          EventInfoFile << "        0      ";  //Output when event content is wrong!
      }
      
      /////////////////////
      // EVENT SELECTION //
      /////////////////////
      
      //Declare selection instance    
      Selection selection(init_jets, init_muons, init_electrons, mets, event->kt6PFJets_rho());
      selection.setJetCuts(30,2.5,0.01,1.,0.98,0.3,0.1);     
      selection.setMuonCuts(26,2.1,0.12,0.2,0.3,1,0.5,5,0); 
      selection.setElectronCuts(32,2.5,0.1,0.02,0.5,0.3,0); 
      selection.setLooseMuonCuts(10,2.5,0.2);
      selection.setLooseElectronCuts(20,2.5,0.15,0.);

      bool triggedSemiMu = false;
      bool triggedSemiEl = false;

      if( ! (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) )
        triggedSemiMu = treeLoader.EventTrigged (itriggerSemiMu);
      if( ! (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) )
        triggedSemiEl = treeLoader.EventTrigged (itriggerSemiEl);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);

      vector<TRootJet*> selectedJets;
      vector<TRootMuon*> selectedMuons;
      vector<TRootElectron*> selectedElectrons;
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10,2.5,0.2);
      vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectrons(20,2.5,0.15);
      vector<TRootElectron*> vetoElectronsSemiEl = selection.GetSelectedLooseElectrons(20,2.5,0.15);

      selectedJets = selection.GetSelectedJets(true);
      selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
      selectedElectrons = selection.GetSelectedElectrons(selectedJets);
      
      //////////////////////
      // Event selection  //
      //////////////////////
      bool eventselectedSemiMu = false, eventselectedSemiEl = false;
      enum DecayChannel_t {semiMu, semiEl};
      DecayChannel_t decayChannel;

      // semi-mu selection
      selecTableSemiMu.Fill(d,0,scaleFactor*lumiWeight);              
      if (triggedSemiMu) {
	selecTableSemiMu.Fill(d,1,scaleFactor*lumiWeight);
	if (isGoodPV) {
	  selecTableSemiMu.Fill(d,2,scaleFactor*lumiWeight);
	  if (selectedMuons.size() == 1) {
	    selecTableSemiMu.Fill(d,3,scaleFactor*lumiWeight);
	    if( vetoMuons.size() == 1 ) {
	      selecTableSemiMu.Fill(d,4,scaleFactor*lumiWeight);
	      if (vetoElectronsSemiMu.size() == 0) {
		selecTableSemiMu.Fill(d,5,scaleFactor*lumiWeight);
		if (selectedJets.size() >= 1) {
		  selecTableSemiMu.Fill(d,6,scaleFactor*lumiWeight);
		  if (selectedJets.size() >= 2) {
		    selecTableSemiMu.Fill(d,7,scaleFactor*lumiWeight);
		    if (selectedJets.size() >= 3) {
		      selecTableSemiMu.Fill(d,8,scaleFactor*lumiWeight);
		      if (selectedJets.size() >= 4) {
			selecTableSemiMu.Fill(d,9,scaleFactor*lumiWeight);
		 	eventselectedSemiMu = true;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      // semi-elec selection
      selecTableSemiEl.Fill(d,0,scaleFactor*lumiWeight);
      if( triggedSemiEl) {
	selecTableSemiEl.Fill(d,1,scaleFactor*lumiWeight);
	if (isGoodPV ) {
	  selecTableSemiEl.Fill(d,2,scaleFactor*lumiWeight);
	  if( selectedElectrons.size() == 1 ) {
	    selecTableSemiEl.Fill(d,3,scaleFactor*lumiWeight);
	    if( vetoMuons.size() == 0 ) {
	      selecTableSemiEl.Fill(d,4,scaleFactor*lumiWeight);
	      if (vetoElectronsSemiEl.size() == 1) {
		selecTableSemiEl.Fill(d,5,scaleFactor*lumiWeight);
		selecTableSemiEl.Fill(d,6,scaleFactor*lumiWeight);
		if( selectedJets.size()>=1 ) {
		  selecTableSemiEl.Fill(d,7,scaleFactor*lumiWeight);
		  if( selectedJets.size()>=2 ) {
		    selecTableSemiEl.Fill(d,8,scaleFactor*lumiWeight);
		    if( selectedJets.size()>=3 ) {
		      selecTableSemiEl.Fill(d,9,scaleFactor*lumiWeight);
		      if( selectedJets.size()>=4 ) {
			selecTableSemiEl.Fill(d,10,scaleFactor*lumiWeight);
			eventselectedSemiEl=true;			
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      
      //--- Only continue with events passing one of the two event selections!  ---//
      if (!eventselectedSemiMu && !eventselectedSemiEl){
        if(getEventInfo == true) EventInfoFile << "        failed       "<<endl;
        continue;
      }
      
      if(getEventInfo == true) EventInfoFile << "        passed            ";  
      //--- do some data-mc ---//
      if (!foundMu && !foundEl) datasetsPlot = datasets;    //Data!
      string leptonFlav="_other";
      TLorentzVector* selectedLepton;
      float LeptonRecoCharge;
      if (eventselectedSemiMu){
        selectedLepton = (TLorentzVector*)selectedMuons[0]; LeptonRecoCharge = selectedMuons[0]->charge();
        decayChannel = semiMu;
        histo1D["nSelectedEvents"]->Fill(3);
        if(LeptonRecoCharge > 0)     { nSelectedMuPos++; histo1D["nSelectedEvents"]->Fill(1);}
        else if(LeptonRecoCharge < 0){ nSelectedMuNeg++; histo1D["nSelectedEvents"]->Fill(2);}
	datasetsPlot = datasetsMu; Luminosity = LuminosityMu;
        leptonFlav="_mu";
      }
      else if (eventselectedSemiEl){
        selectedLepton = (TLorentzVector*)selectedElectrons[0]; LeptonRecoCharge = selectedElectrons[0]->charge();
        decayChannel = semiEl;
        histo1D["nSelectedEvents"]->Fill(7);
        if(LeptonRecoCharge > 0)     { nSelectedElPos++; histo1D["nSelectedEvents"]->Fill(5);}
        else if(LeptonRecoCharge < 0){ nSelectedElNeg++; histo1D["nSelectedEvents"]->Fill(6);}
	datasetsPlot = datasetsEl; Luminosity = LuminosityEl;
        leptonFlav="_el";
      }
      
      ////////////////////////////////////////////////////////////////////
      //   Use genEvent information to get the correct event topology   //
      ////////////////////////////////////////////////////////////////////
      TLorentzVector genLight1, genLight2, genHadrB, genLeptB, genLepton;
      if(dataSetName.find("TTbarJets_SemiLept") == 0){
        vector<TRootMCParticle> mcParticlesMatching;      	
        vector< pair<unsigned int, unsigned int> > JetPartonPair, ISRJetPartonPair; // First one is jet number, second one is mcParticle number
	
	//First index is the jet number, the second one the mcParticle number:
	pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_;
	leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int, unsigned int>(9999,9999);
	vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
	bool muPlusFromTop = false, muMinusFromTop = false, elPlusFromTop = false, elMinusFromTop = false;
	if(verbosity>1) cout << " Looking at mcParticlesMatching " << endl;

	for(unsigned int i=0; i<mcParticles.size(); i++){
	  if(mcParticles[i]->status() != 3) continue;
	  
	  //Muon identification:
	  if(mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6){
            if(muMinusFromTop) cerr<<"muMinusFromTop already true"<<endl; muMinusFromTop = true;
          }
	  if(mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6){
            if(muPlusFromTop) cerr<<"muPlusFromTop already true"<<endl; muPlusFromTop = true;
          }
	  
	  //Electron identification:
	  if(mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6){
            if(elMinusFromTop) cerr<<"elMinusFromTop already true"<<endl; elMinusFromTop = true;
          }
	  if(mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6){
            if(elPlusFromTop) cerr<<"elPlusFromTop already true"<<endl; elPlusFromTop = true;
          }
	  
	  if( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 ){
	    mcParticlesTLV.push_back(*mcParticles[i]); mcParticlesMatching.push_back(*mcParticles[i]);
	  }
	}

        //--- Control checks ---//
	if(muPlusFromTop && muMinusFromTop) cerr<<"muPlusFromTop and muMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;
	if(elPlusFromTop && elMinusFromTop) cerr<<"elPlusFromTop and elMinusFromTop are both true ?!\nCheck if you are using the right sample..."<<endl;

        //--- MATCHING ---//	
	for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);  //Need selectedJets for radiation stuff!
	
	JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
	if(matching.getNumberOfAvailableCombinations() != 1) cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<"  This should be equal to 1 !!!"<<endl;
	for(unsigned int i=0; i<mcParticlesTLV.size(); i++){
	  int matchedJetNumber = matching.getMatchForParton(i, 0);
	  if(matchedJetNumber != -1) JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
	}
	
	for(unsigned int itJPPair=0; itJPPair<JetPartonPair.size(); itJPPair++){
	  unsigned int mcPartNr = JetPartonPair[itJPPair].second;
	  
	  if(verbosity > 3){ 
	    std::cout <<" Jet number " << JetPartonPair[itJPPair].first << " can be matched with the mcParticle number " << JetPartonPair[itJPPair].second << " which is of the type : " << mcParticlesMatching[mcPartNr].type() << std::endl;
	  }

          //--- Recover the light jets ---//	  
	  if( fabs(mcParticlesMatching[mcPartNr].type()) < 6 ){
	    if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == -24 && mcParticlesMatching[mcPartNr].grannyType() == -6 ) || 
                ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == 24  && mcParticlesMatching[mcPartNr].grannyType() == 6 )  ){
	      if(hadronicWJet1_.first == 9999) 
		hadronicWJet1_ = JetPartonPair[itJPPair];
	      else if(hadronicWJet2_.first == 9999) 
		hadronicWJet2_ = JetPartonPair[itJPPair];
	      else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
	    }
	  }

          //--- Recover the b-jets ---//
	  if( fabs(mcParticlesMatching[mcPartNr].type()) == 5 ){
	    if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == -6 ) || 
                ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == 6  ) )
	      hadronicBJet_ = JetPartonPair[itJPPair];
	    
	    else if( ( (muPlusFromTop  || elPlusFromTop)  && mcParticlesMatching[mcPartNr].motherType() == 6 ) || 
                     ( (muMinusFromTop || elMinusFromTop) && mcParticlesMatching[mcPartNr].motherType() == -6) )
	      leptonicBJet_ = JetPartonPair[itJPPair];
	  }
	}
	
        //--- Initialize jetCombi vector which has the jet number information ---//
	jetCombi[0] = leptonicBJet_.first;
	jetCombi[1] = hadronicBJet_.first;
	jetCombi[2] = hadronicWJet1_.first;
	jetCombi[3] = hadronicWJet2_.first;
	if(verbose > 3) cout<<" Index of BHadronic: "<<jetCombi[1]<<" , Index of BLeptonic: "<<jetCombi[0]<<" , Index of quark1: "<<jetCombi[2]<<" & Index of quark2: "<<jetCombi[3]<<endl;
	
        //--- Make CSV discriminant and jet number plots ---//	
	if(jetCombi[0] != 9999) histo1D["CorrectBLeptCSVDiscr"]->Fill(selectedJets[jetCombi[0]]->btag_combinedSecondaryVertexBJetTags());
	else                    histo1D["CorrectBLeptCSVDiscr"]->Fill(-2);
	if(jetCombi[1] != 9999) histo1D["CorrectBHadrCSVDiscr"]->Fill(selectedJets[jetCombi[1]]->btag_combinedSecondaryVertexBJetTags());
	else                    histo1D["CorrectBHadrCSVDiscr"]->Fill(-2);
	if(jetCombi[2] != 9999){histo1D["CorrectQuark1CSVDiscr"]->Fill(selectedJets[jetCombi[2]]->btag_combinedSecondaryVertexBJetTags()); histo1D["Quark1JetNumber"]->Fill(jetCombi[2]);}
	else                   {histo1D["CorrectQuark1CSVDiscr"]->Fill(-2);                                                                histo1D["Quark1JetNumber"]->Fill(-1);}
	if(jetCombi[3] != 9999){histo1D["CorrectQuark2CSVDiscr"]->Fill(selectedJets[jetCombi[3]]->btag_combinedSecondaryVertexBJetTags()); histo1D["Quark2JetNumber"]->Fill(jetCombi[3]);}
	else                   {histo1D["CorrectQuark2CSVDiscr"]->Fill(-2);                                                                histo1D["Quark2JetNumber"]->Fill(-1);}
 	
	//Working on generator level (i.e. parton level):  
	if(jetCombi[0]!=9999 && jetCombi[1]!=9999 && jetCombi[2]!=9999 && jetCombi[3]!=9999){
          genLight1 = (TLorentzVector) mcParticlesMatching[hadronicWJet1_.second];
          genLight2 = (TLorentzVector) mcParticlesMatching[hadronicWJet2_.second];
          genHadrB = (TLorentzVector) mcParticlesMatching[hadronicBJet_.second];
          genLeptB = (TLorentzVector) mcParticlesMatching[leptonicBJet_.second];
          genLepton = (TLorentzVector) *lhcoOutput.getGenLepton();

	}//End of matched particles reconstructed
      
      }//End of TTbarJets!

      //------- Fill the Tree file for the LightAnomCoupAnalyzer file -------//
      anomCoupLight = new AnomCoupLight();

      std::vector<float> CSVbTagValues;
      vector<TLorentzVector> SelectedJets;
      for(int iJet = 0; iJet < selectedJets.size(); iJet++){
        CSVbTagValues.push_back(selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags());
        SelectedJets.push_back(*selectedJets[iJet]);
      }

      anomCoupLight->setEventId(event->eventId());
      anomCoupLight->setRunId(event->runId());
      anomCoupLight->setLumiBlockId(event->lumiBlockId());
      anomCoupLight->setNPV(vertex.size());
      anomCoupLight->setNTruePU(event->nTruePU());
      anomCoupLight->setScaleFactor(scaleFactor);

      anomCoupLight->setSelectedJets(SelectedJets);
      anomCoupLight->setBTagCSV(CSVbTagValues);
      anomCoupLight->setSelectedLepton(*selectedLepton);
      anomCoupLight->setDecayChannel(decayChannel);
      anomCoupLight->setLeptonCharge(LeptonRecoCharge);
      anomCoupLight->setCorrectJetCombi(jetCombi);
      anomCoupLight->setMET(*mets[0]);
      if(dataSetName.find("TTbarJets") == 0 && lhcoOutput.GenEventContentCorrect()) anomCoupLight->setGenCosTheta(kinFunctions.CosTheta(lhcoOutput.getGenLeptTop(), lhcoOutput.getGenLeptW(), lhcoOutput.getGenLepton()));
      else anomCoupLight->setGenCosTheta(2.);
	
      //Store the information needed for the TF (but only has value when dataset is ttbar)
      anomCoupLight->setGenVectorLight1( genLight1 );
      anomCoupLight->setGenVectorLight2( genLight2 );
      anomCoupLight->setGenVectorHadrB( genHadrB );
      anomCoupLight->setGenVectorLeptB( genLeptB );
      anomCoupLight->setGenVectorLepton( genLepton);

      LightTree->Fill();
      delete anomCoupLight;
      //----  End of Tree file filling (LightAnomCoupAnalyzer)  ----//

    } //loop on events
  
    if(dataSetName.find("TTbarJets_SemiLept") == 0){
      cout << "\n -> " << nSelectedMuPos << " mu+ and " << nSelectedMuNeg << " mu- events are selected on " << nGenMu << " ==> " << (float)(nSelectedMuPos+nSelectedMuNeg)/(float)nGenMu << endl;
      cout <<   " -> " << nSelectedElPos << " el+ and " << nSelectedElNeg << " el- events are selected on " << nGenEl << " ==> " << (float)(nSelectedElPos+nSelectedElNeg)/(float)nGenEl << endl;
    }
    else{
      cout << "\n -> " << nSelectedMuPos << " mu+ and " << nSelectedMuNeg << " mu- events are selected on " << NrEvtsToRunOver << " ==> " << (float)(nSelectedMuPos+nSelectedMuNeg)/(float)NrEvtsToRunOver << endl;
      cout <<   " -> " << nSelectedElPos << " el+ and " << nSelectedElNeg << " el- events are selected on " << NrEvtsToRunOver << " ==> " << (float)(nSelectedElPos+nSelectedElNeg)/(float)NrEvtsToRunOver << endl;
    }
    std::cout << " " << std::endl;

    //--------------------------------//
    // Store the gen-level LHCO plots //
    //--------------------------------//
    lhcoOutput.WriteLHCOPlots(fout);

    //------ Store the LightAnomCoupAnalyzer tree ------//
    LightFile->cd();
    TTree* configTreeLightFile = new TTree("configTreeLightFile","configuration Tree in Light File");
    TClonesArray* tcdatasetlightfile = new TClonesArray("Dataset",1);
    configTreeLightFile->Branch("Dataset","TClonesArray",&tcdatasetlightfile);
    new ((*tcdatasetlightfile)[0]) Dataset(*datasets[d]);
    //TClonesArray* tcanaenvlightfile = new TClonesArray("AnalysisEnvironment",1);   --> Needed?
    //configTreeLightFile->Branch("AnaEnv","TClonesArray",&tcanaenvlightfile);
    //new ((*tcanaenvlightfile)[0]) AnalysisEnvironment(anaEnv);

    configTreeLightFile->Fill();
    configTreeLightFile->Write();
    LightTree->Write();
    LightFile->Close();
    delete LightFile;
    //----  End of storing Tree  ----//

    //////////////
    // CLEANING //
    //////////////
    if(getEventInfo == true) EventInfoFile.close();
    
    if (jecUnc) delete jecUnc;
    if (jetTools) delete jetTools;
  
    //important: free memory
    treeLoader.UnLoadDataset();    
  }  //loop on datasets
  
  //Once everything is filled ...
  if(verbose > 0) cout << "We ran over all the data ;-)" << endl;
  
  /////////////////////////
  // Write out the plots //
  /////////////////////////
  string pathPNG = "PlotsMacro";

  //Write away the MSPlots, the histo1D's and the histo2D's!!
  fout -> cd();
  mkdir((pathPNG+"/MSPlots").c_str(),0777);

  TDirectory* msdir = fout->mkdir("MSPlots");
  msdir->cd(); 
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 0, false, false, false, 1);     //string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSSignal 
    temp->Write(fout, name, saveAsPDF, (pathPNG+"/MSPlots/").c_str(), "pdf");
  }

  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  TDirectory* th2dir = fout->mkdir("2D_histograms_graphs");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){
    TH2F *temp = it->second;
    temp->Write();
  }

  fout->cd();
  configTree->Fill();
  configTree->Write();

  //Selection tables
  selecTableSemiMu.TableCalculator(false,true,true,true,true);
  string selectiontableMu = "EventSelectionResults/AnalyzerOutput/SelectionTable_BTAG_SEMIMU.tex";
  selecTableSemiMu.Write(selectiontableMu.c_str());
  selecTableSemiEl.TableCalculator(false, true, true, true, true);
  string selectiontableEl = "EventSelectionResults/AnalyzerOutput/SelectionTable_BTAG_SEMIEL.tex";
  selecTableSemiEl.Write(selectiontableEl.c_str());

  fout->Close(); 
   
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  //delete configTree;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "           hasn't crashed yet ;-)           " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}
