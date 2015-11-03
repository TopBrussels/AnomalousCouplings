#include "../interface/TFCreation.h"

TFCreation::TFCreation(int nEtaBins, std::string etaCons, bool doFits){

  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);

  ///////////////////////////////////////
  //  Declare the used fit functions!  //
  ///////////////////////////////////////
  //
  // 1) Double Gaussian --> its range depends on the jet/lepton energy range (hence, the Y-axis)
  doubleGaussianFit = new TF1("doubleGaussianFit","[5]*(1/(TMath::Sqrt(2*TMath::Pi())*(TMath::Sqrt(TMath::Power([1],2))+TMath::Sqrt(TMath::Power([4],2))*[2])))*(TMath::Exp(-TMath::Power((x-[0]),2)/(2*TMath::Power([1],2)))+[2]*TMath::Exp(-TMath::Power((x-[3]),2)/(2*TMath::Power([4],2))))");
  nParsFit_ = doubleGaussianFit->GetNpar();
  std::string parnames[6]={"a1","a2","a3","a4","a5","Ampl"};
  std::string ParName[6] = {"Mean first gaussian", "Width first gaussian","Relative constant gaussians","Mean second gaussian","Width second gaussian","Amplitude"};

  if(nParsFit_ != sizeof(parnames)/sizeof(parnames[0])) std::cout << " ERROR : Difference between number of parameters and defined array --> Also check the header file than !! " << std::endl;

  for(int ipar = 0; ipar < nParsFit_; ipar++){
    parnames_[ipar] = parnames[ipar];
    ParName_[ipar] = ParName[ipar];
    doubleGaussianFit->SetParName(ipar,parnames[ipar].c_str());
  }
 
  //2) Calorimeter Energy formula (ai = ai0 + ai1*Ep + ai2*sqrt(Ep)) --> its range depends on the part energy range (hence, the X-axis)
  caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");
  if(doFits){ caloFitFile = new TFile(("TFInformation/CaloEnergyFitFunctions"+etaCons+".root").c_str(),"RECREATE");}

  //Store the EtaBin Title and Name for the histograms!
  nEtaBins_ = nEtaBins;
  EtaBin[0] = ""; EtaTitle[0] = "";
  EtaValues[1] = 0.; EtaValues[2] = 0.375; EtaValues[3] = 0.750; EtaValues[4] = 1.450; EtaValues[5] = 2.5;
  EtaValue[1] = "0"; EtaValue[2] = "0p375"; EtaValue[3] = "0p750"; EtaValue[4] = "1p450"; EtaValue[5] = "2p5";
  for(int ii = 1; ii <= 4; ii++){
    EtaBin[ii] = "_Eta_"+EtaValue[ii]+"_"+EtaValue[ii+1];
    EtaTitle[ii] = " -- "+tostr(EtaValues[ii])+" < |#eta| #leq "+tostr(EtaValues[ii+1]);
  }
}

TFCreation::~TFCreation(){
  delete caloEnergyFit;
  delete doubleGaussianFit;
  caloFitFile->Close();
  delete caloFitFile;
}

void TFCreation::InitializeVariables(){

  //Actual TF distributions!
  histo1D["DeltaR_TFClass_Light1"] = new TH1F("DeltaR_TFClass_Light1","DeltaR_TFClass_Light1",200,0,0.4);
  histo1D["DeltaR_TFClass_Light2"] = new TH1F("DeltaR_TFClass_Light2","DeltaR_TFClass_Light2",200,0,0.4);
  histo1D["DeltaR_TFClass_HadrB"]  = new TH1F("DeltaR_TFClass_HadrB","DeltaR_TFClass_HadrB",200,0,0.4);
  histo1D["DeltaR_TFClass_LeptB"]  = new TH1F("DeltaR_TFClass_LeptB","DeltaR_TFClass_LeptB",200,0,0.4);
  histo1D["DeltaR_TFClass_Mu"]     = new TH1F("DeltaR_TFClass_Mu","DeltaR_TFClass_Mu",200,0,4);
  histo1D["DeltaR_TFClass_El"]     = new TH1F("DeltaR_TFClass_El","DeltaR_TFClass_El",200,0,4);

  histo1D["PtDistribution_Electron"] = new TH1F("PtDistr_El","Pt-distribution for electron",150,25,200);
  histo1D["EDistribution_Electron"] = new TH1F("EDistr_El","E-distribution for electron",250,25,350);
  histo1D["PtDistribution_Muon"] = new TH1F("PtDistr_Mu","Pt-distribution for muon",150,25,200);
  histo1D["EDistribution_Muon"] = new TH1F("EDistr_Mu","E-distribution for muon",250,25,350);
  histo1D["Mass_GenElec"] = new TH1F("MassDistr_GenEl","Mass distribution for generator-level electron",150,0,0.001);
  histo1D["Mass_GenMuon"] = new TH1F("MassDistr_GenMu","Mass distribution for generator-level muon",150, 0,0.2);
  histo1D["Mass_RecoElec"] = new TH1F("MassDistr_RecoEl","Mass distribution for reco-level electron",150,0,0.1);
  histo1D["Mass_RecoMuon"] = new TH1F("MassDistr_RecoMu","Mass distribution for reco-level muon",150,0,1);
  histo2D["Muon_Pt_vs_Eta"] = new TH2F("Muon_Pt_vs_Eta","Pt-value versus eta-value for muons",60,-4,4,200,0,200);
  histo2D["Muon_E_vs_Eta"] = new TH2F("Muon_E_vs_Eta","E-value versus eta-value for muons",60,-4,4,200,0,200);

  histo2D["Light_RecoEVsGenE"]         = new TH2F("Light_RecoEVsGenE",         "Energy of light quarks (reco vs gen)",                             150,    0,  300, 150,     0,   300);
  histo2D["Light_DiffEVsGenE"]         = new TH2F("Light_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for light quarks",           16,   25,  160, 200,   -75,    75);
  histo2D["Light_RecoPtVsGenPt"]       = new TH2F("Light_RecoPtVsGenPt",       "Transverse momentum of light quarks (reco vs gen)",                150,    0,  300, 150,     0,   300);
  histo2D["Light_DiffPtVsGenPt"]       = new TH2F("Light_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for light quarks",        10,   30,  125, 100,   -30,    35);
  histo2D["Light_RecoThetaVsGenTheta"] = new TH2F("Light_RecoThetaVsGenTheta", "Polar angle distribution of light quarks (reco vs gen)",            60,    0, 3.15,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenTheta"] = new TH2F("Light_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for light quarks", 10,  0.1,  3.1, 100, -0.15,  0.15);
  histo2D["Light_RecoThetaVsGenE"]     = new TH2F("Light_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for light quarks",         120,    0,  300,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenE"]     = new TH2F("Light_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for light quarks",      15,   30,  185, 150,  -0.1,   0.1);
  histo2D["Light_RecoThetaVsGenPt"]    = new TH2F("Light_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for light quarks",       120,    0,  300,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenPt"]    = new TH2F("Light_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{Tgen} for light quarks",     10,   30,  165, 150, -0.12,  0.12);
  histo2D["Light_RecoPhiVsGenPhi"]     = new TH2F("Light_RecoPhiVsGenPhi",     "Azimuthal angle distribution of light quarks (reco vs gen)",        60, -3.2,  3.2,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenPhi"]     = new TH2F("Light_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for light quarks",     10, -3.2,  3.2, 100,  -0.2,   0.2);
  histo2D["Light_RecoPhiVsGenE"]       = new TH2F("Light_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for light quarks",       150,    0,  300,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenE"]       = new TH2F("Light_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for light quarks",        15,   30,  185, 100, -0.12,  0.12);
  histo2D["Light_RecoPhiVsGenPt"]      = new TH2F("Light_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for light quarks",     150,    0,  300,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenPt"]      = new TH2F("Light_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for light quarks",      10,   30,  165, 100, -0.15,  0.15);

  histo2D["BJet_RecoEVsGenE"]         = new TH2F("BJet_RecoEVsGenE",         "Energy of b-jets (reco vs gen level)",                        150,    0,  300, 150,     0,  300);
  histo2D["BJet_DiffEVsGenE"]         = new TH2F("BJet_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for b-jets",            18,   30,  230, 350,  -100,  100);
  histo2D["BJet_RecoPtVsGenPt"]       = new TH2F("BJet_RecoPtVsGenPt",       "Transverse momentum of b-jets (reco vs gen level)",           150,    0,  300, 150,     0,  300);
  histo2D["BJet_DiffPtVsGenPt"]       = new TH2F("BJet_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for b-jets",         10,   30,  150, 100,   -35,   50);
  histo2D["BJet_RecoThetaVsGenTheta"] = new TH2F("BJet_RecoThetaVsGenTheta", "Polar angle distribution of b-jets (reco vs gen)",             60,    0, 3.15,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenTheta"] = new TH2F("BJet_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for b-jets",  10,  0.1,  3.1, 100, -0.15, 0.15);
  histo2D["BJet_RecoThetaVsGenE"]     = new TH2F("BJet_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus energy E_{gen} for b-jets",   120,    0,  300,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenE"]     = new TH2F("BJet_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for b-jets",       17,   30,  240, 150, -0.09, 0.09);
  histo2D["BJet_RecoThetaVsGenPt"]    = new TH2F("BJet_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus energy P_{T,gen} for b-jets", 120,    0,  300,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenPt"]    = new TH2F("BJet_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for b-jets",     10,   30,  160, 150,  -0.1,  0.1);
  histo2D["BJet_RecoPhiVsGenPhi"]     = new TH2F("BJet_RecoPhiVsGenPhi",     "Azimuthal angle distribution of b-jets (reco vs gen)",         60, -3.2,  3.2,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenPhi"]     = new TH2F("BJet_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for b-jets",      10, -3.2,  3.2, 100,  -0.2,  0.2);
  histo2D["BJet_RecoPhiVsGenE"]       = new TH2F("BJet_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for b-jets",        150,    0,  300,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenE"]       = new TH2F("BJet_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for b-jets",         18,   20,  250, 100,  0.15, 0.15);
  histo2D["BJet_RecoPhiVsGenPt"]      = new TH2F("BJet_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for b-jets",      150,    0,  300,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenPt"]      = new TH2F("BJet_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for b-jets",       10,   30,  160, 100, -0.15, 0.15);

  histo2D["El_RecoEVsGenE"]         = new TH2F("El_RecoEVsGenE",         "Energy of electron (reco vs gen)",                             150,    0,  200, 150,      0,   200);
  histo2D["El_DiffEVsGenE"]         = new TH2F("El_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for electron",           25,   30,  220, 100,     -6,     6);
  histo2D["El_RecoPtVsGenPt"]       = new TH2F("El_RecoPtVsGenPt",       "Transverse momentum of electron (reco vs gen)",                100,    0,  200, 100,      0,   200);
  histo2D["El_DiffPtVsGenPt"]       = new TH2F("El_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for electron",        10,   30,  105, 100,     -6,     6);
  histo2D["El_RecoThetaVsGenTheta"] = new TH2F("El_RecoThetaVsGenTheta", "Polar angle distribution of electron (reco vs gen)",            60,    0, 3.15,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenTheta"] = new TH2F("El_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for electron", 10,    0, 3.15, 100,  -0.15,  0.15);
  histo2D["El_RecoThetaVsGenE"]     = new TH2F("El_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for electron",         100,    0,  200,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenE"]     = new TH2F("El_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for electron",      20,   30,  220, 100, -0.015, 0.015);
  histo2D["El_RecoThetaVsGenPt"]    = new TH2F("El_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for electron",       100,    0,  200,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenPt"]    = new TH2F("El_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for electron",    10,   30,  130, 100,  -0.02,  0.02);
  histo2D["El_RecoPhiVsGenPhi"]     = new TH2F("El_RecoPhiVsGenPhi",     "Azimuthal angle distribution of electron (reco vs gen)", 	  60, -3.2,  3.2,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenPhi"]     = new TH2F("El_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for electron", 	  10, -3.2,  3.2,  75,  -0.15,  0.15);
  histo2D["El_RecoPhiVsGenE"]       = new TH2F("El_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for electron",       150,    0,  250,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenE"]       = new TH2F("El_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for electron",        18,   30,  200, 120, -0.008, 0.008);
  histo2D["El_RecoPhiVsGenPt"]      = new TH2F("El_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for electron",     150,    0,  250,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenPt"]      = new TH2F("El_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for electron",      10,   30,  130, 120, -0.015, 0.015);
	
  histo2D["Mu_RecoInvEVsGenInvE"]    = new TH2F("Mu_RecoInvEVsGenInvE",   "Inverse of energy of muon (reco vs gen)",                                  100,     0,  0.05, 100,      0,   0.05);
  histo2D["Mu_DiffInvEVsGenInvE"]    = new TH2F("Mu_DiffInvEVsGenInvE",   "#frac{1}{E} difference (gen-reco) versus #frac{1}{E_{gen}} for muon",       15, 0.001,  0.03, 120, -0.001,  0.001);
  histo2D["Mu_RecoInvPtVsGenInvPt"]  = new TH2F("Mu_RecoInvPtVsGenInvPt", "Inverse of transverse momentum of muon (reco vs gen)",                     100,     0,  0.05, 100,      0,   0.05);
  histo2D["Mu_DiffInvPtVsGenInvPt"]  = new TH2F("Mu_DiffInvPtVsGenInvPt", "#frac{1}{P_{T}} difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon", 10, 0.005, 0.035, 120,-0.0018, 0.0018);
  histo2D["Mu_RecoEVsGenE"]          = new TH2F("Mu_RecoEVsGenE",         "Energy of muon (reco vs gen)",                                             150,     0,   200, 150,      0,    200);
  histo2D["Mu_DiffEVsGenE"]          = new TH2F("Mu_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for muon",                           14,    20,   160, 125,     -6,      6);
  histo2D["Mu_RecoPtVsGenPt"]        = new TH2F("Mu_RecoPtVsGenPt",       "Transverse momentum of muon (reco vs gen)",                                150,     0,   200, 150,      0,    200);
  histo2D["Mu_DiffPtVsGenPt"]        = new TH2F("Mu_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for muon",                        10,    26,   150, 100,    -10,     10);
  histo2D["Mu_RecoThetaVsGenTheta"]  = new TH2F("Mu_RecoThetaVsGenTheta", "Polar angle distribution of muon (reco vs gen)",                            60,     0,  3.15,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenTheta"]  = new TH2F("Mu_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for muon",                 10,     0,  3.15, 100,  -0.15,   0.15);
  histo2D["Mu_RecoThetaVsGenInvE"]   = new TH2F("Mu_RecoThetaVsGenInvE",  "Polar angle #theta_{rec} versus #frac{1}{E_{gen}} for muon",               100,     0,  0.05,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenInvE"]   = new TH2F("Mu_DiffThetaVsGenInvE",  "#theta difference (gen-reco) versus #frac{1}{E_{gen}} for muon",            15, 0.001,  0.03, 150, -0.004,  0.004);
  histo2D["Mu_RecoThetaVsGenInvPt"]  = new TH2F("Mu_RecoThetaVsGenInvPt", "Polar angle #theta_{rec} versus #frac{1}{P_{T,gen}} for muon",             100,     0,  0.05,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenInvPt"]  = new TH2F("Mu_DiffThetaVsGenInvPt", "#theta difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon",          10, 0.005, 0.035, 150, -0.004,  0.004);
  histo2D["Mu_RecoThetaVsGenE"]      = new TH2F("Mu_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for muon",                         150,     0,   200,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenE"]      = new TH2F("Mu_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for muon",                      10,    30,   150, 100,   -0.1,    0.1);
  histo2D["Mu_RecoThetaVsGenPt"]     = new TH2F("Mu_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for muon",                       150,     0,   200,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenPt"]     = new TH2F("Mu_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for muon",                    10,    30,   150, 100,   -0.1,    0.1);
  histo2D["Mu_RecoPhiVsGenPhi"]      = new TH2F("Mu_RecoPhiVsGenPhi",     "Azimuthal angle distribution of muon (reco vs gen)", 		       60,  -3.2,   3.2,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenPhi"]      = new TH2F("Mu_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for muon", 		       10,  -3.2,   3.2,  75,   -0.2,    0.2);
  histo2D["Mu_RecoPhiVsGenInvE"]     = new TH2F("Mu_RecoPhiVsGenInvE",    "Azimuthal angle #phi_{rec} versus #frac{1}{E_{gen}} for muon",             100,     0,  0.05,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenInvE"]     = new TH2F("Mu_DiffPhiVsGenInvE",    "#phi difference (gen-reco) versus #frac{1}{E_{gen}} for muon",              15, 0.001,  0.03, 120, -0.005,  0.005);
  histo2D["Mu_RecoPhiVsGenInvPt"]    = new TH2F("Mu_RecoPhiVsGenInvPt",   "Azimuthal angle #phi_{rec} versus #frac{1}{P_{T,gen}} for muon",           100,     0,  0.05,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenInvPt"]    = new TH2F("Mu_DiffPhiVsGenInvPt",   "#phi difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon",            10, 0.005, 0.035, 120, -0.006,  0.006);
  histo2D["Mu_RecoPhiVsGenE"]        = new TH2F("Mu_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for muon",                       150,     0,   200,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenE"]        = new TH2F("Mu_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for muon",                        10,    26,   150, 100,   -0.2,    0.2);
  histo2D["Mu_RecoPhiVsGenPt"]       = new TH2F("Mu_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for muon",                     150,     0,   200,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenPt"]       = new TH2F("Mu_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for muon",                      10,    26,   150, 100,   -0.2,    0.2);

  //Store all the existing histo2D's in a new collection!
  map<string,TH2F*> histo2DCopy = histo2D;
  for(int iEta = 1; iEta <= 4; iEta++){

    for(std::map<std::string,TH2F*>::const_iterator it = histo2DCopy.begin(); it != histo2DCopy.end(); it++){
      TH2F *temp = it->second;
	
      //Need different Ymax and Ymin for Mu_DiffThetaVsGenInvPt histogram!
      double YMax = temp->GetYaxis()->GetXmax();
      double YMin = temp->GetYaxis()->GetXmin();
      int NYBins = (int)(temp->GetYaxis()->GetNbins());

      //if( string(temp->GetName()) == "Mu_DiffThetaVsGenInvPt"){ YMax = 0.015; YMin = -0.015;}
      if( string(temp->GetName()) == "Mu_DiffEVsGenE" ){
        YMax = 4.0; YMin = -3.5;
        if(iEta == 1 || iEta == 2) NYBins = NYBins*1.1;
        else if(iEta == 4){ NYBins = NYBins*0.9; }
      }
      else if(string(temp->GetName()) == "Light_DiffEVsGenE" && iEta == 4) NYBins = NYBins*0.95;

      if(iEta == 4)
	histo2D[temp->GetName()+EtaBin[iEta]] = new TH2F( (temp->GetName()+EtaBin[iEta]).c_str(), (temp->GetTitle()+EtaTitle[iEta]).c_str(), temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXmin(), temp->GetXaxis()->GetXmax(), NYBins*0.9, YMin*1.4, YMax*1.4 );
      else if(iEta == 3)
	histo2D[temp->GetName()+EtaBin[iEta]] = new TH2F( (temp->GetName()+EtaBin[iEta]).c_str(), (temp->GetTitle()+EtaTitle[iEta]).c_str(), temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXmin(), temp->GetXaxis()->GetXmax(), NYBins, YMin*1.2, YMax*1.2);
      else
        histo2D[temp->GetName()+EtaBin[iEta]] = new TH2F( (temp->GetName()+EtaBin[iEta]).c_str(), (temp->GetTitle()+EtaTitle[iEta]).c_str(), temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXmin(), temp->GetXaxis()->GetXmax(), NYBins, YMin, YMax);      
    }
  }
}

void TFCreation::FillHistograms(TLorentzVector* hadrWJet1, TLorentzVector* hadrWJet2, TLorentzVector* hadrBJet, TLorentzVector* leptBJet, TLorentzVector* lepton, TLorentzVector* selHadrWJet1, TLorentzVector* selHadrWJet2, TLorentzVector* selHadrBJet, TLorentzVector* selLeptBJet, TLorentzVector* selLepton, int enumDecayChannel){

  bool isSemiMu = false, isSemiEl = false;
  if(enumDecayChannel == 0) isSemiMu = true;
  else if(enumDecayChannel == 1) isSemiEl = true;

  //Should use Pt information in stead of E!   ... BUT ... doesn't work in MadWeight as expected ....
  // --> Both concepts are identical in the case of CaloJets, but not in the case of PF
  // --> PF uses massive objects to construct particles!

  histo1D["DeltaR_TFClass_Light1"]->Fill( hadrWJet1->DeltaR(*selHadrWJet1) );
  histo1D["DeltaR_TFClass_Light2"]->Fill( hadrWJet2->DeltaR(*selHadrWJet2) );
  histo1D["DeltaR_TFClass_HadrB"]->Fill( hadrBJet->DeltaR(*selHadrBJet) );
  histo1D["DeltaR_TFClass_LeptB"]->Fill( leptBJet->DeltaR(*selLeptBJet) );
  if(isSemiMu) histo1D["DeltaR_TFClass_Mu"]->Fill( lepton->DeltaR(*selLepton) );
  if(isSemiEl) histo1D["DeltaR_TFClass_El"]->Fill( lepton->DeltaR(*selLepton) );

  if(isSemiEl) histo1D["PtDistribution_Electron"]->Fill(lepton->Pt());
  if(isSemiEl) histo1D["EDistribution_Electron"]->Fill(lepton->E());
  if(isSemiMu) histo1D["PtDistribution_Muon"]->Fill(lepton->Pt());
  if(isSemiMu) histo1D["EDistribution_Muon"]->Fill(lepton->E());
  if(isSemiEl){histo1D["Mass_RecoElec"]->Fill(selLepton->M()); histo1D["Mass_GenElec"]->Fill(lepton->M());}
  if(isSemiMu){histo1D["Mass_RecoMuon"]->Fill(selLepton->M()); histo1D["Mass_GenMuon"]->Fill(lepton->M());}

  if(isSemiMu){histo2D["Muon_Pt_vs_Eta"]->Fill(selLepton->Eta(), selLepton->Pt()); histo2D["Muon_E_vs_Eta"]->Fill(selLepton->Eta(), selLepton->E());}

  for(int iEtaCase = 0; iEtaCase < 2; iEtaCase++){
    // Select the correct eta-bin which should be used!
    int whichEtaBin = 0; 
    int useEtaBinWJet1 = 99, useEtaBinWJet2 = 99, useEtaBinHadrB = 99, useEtaBinLeptB = 99, useEtaBinLepton = 99;
    if(iEtaCase == 1){
      for(int iEta = 1; iEta <= 4; iEta++){
        if(abs(selHadrWJet1->Eta()) <= EtaValues[iEta+1] && abs(selHadrWJet1->Eta()) > EtaValues[iEta]) useEtaBinWJet1 = iEta;
        if(abs(selHadrWJet2->Eta()) <= EtaValues[iEta+1] && abs(selHadrWJet2->Eta()) > EtaValues[iEta]) useEtaBinWJet2 = iEta;
        if(abs(selHadrBJet->Eta())  <= EtaValues[iEta+1] && abs(selHadrBJet->Eta()) > EtaValues[iEta])  useEtaBinHadrB = iEta;
        if(abs(selLeptBJet->Eta())  <= EtaValues[iEta+1] && abs(selLeptBJet->Eta()) > EtaValues[iEta])  useEtaBinLeptB = iEta;
        if(abs(selLepton->Eta())    <= EtaValues[iEta+1] && abs(selLepton->Eta()) > EtaValues[iEta])    useEtaBinLepton = iEta;
      }
    }
    //Fill histograms for first light jet in case eta-splitting is desired!
    if(iEtaCase == 1) whichEtaBin = useEtaBinWJet1;
    histo2D["Light_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrWJet1->E(),    selHadrWJet1->E()     );
    histo2D["Light_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrWJet1->Pt(),   selHadrWJet1->Pt()     );
    histo2D["Light_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrWJet1->Theta(),selHadrWJet1->Theta() );
    histo2D["Light_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrWJet1->E(),    selHadrWJet1->Theta() );
    histo2D["Light_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrWJet1->Pt(),    selHadrWJet1->Theta() );
    histo2D["Light_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrWJet1->Phi(),  selHadrWJet1->Phi()   );
    histo2D["Light_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrWJet1->E(),    selHadrWJet1->Phi()   );
    histo2D["Light_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrWJet1->Pt(),    selHadrWJet1->Phi()   );   

    histo2D["Light_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrWJet1->E(),    hadrWJet1->E()     - selHadrWJet1->E()     );
    histo2D["Light_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrWJet1->Pt(),   hadrWJet1->Pt()    - selHadrWJet1->Pt()     );
    histo2D["Light_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrWJet1->Theta(),hadrWJet1->Theta() - selHadrWJet1->Theta() );
    histo2D["Light_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrWJet1->E(),    hadrWJet1->Theta() - selHadrWJet1->Theta() );
    histo2D["Light_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrWJet1->Pt(),    hadrWJet1->Theta() - selHadrWJet1->Theta() );
    histo2D["Light_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrWJet1->Phi(),  hadrWJet1->DeltaPhi(*selHadrWJet1)   );
    histo2D["Light_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrWJet1->E(),    hadrWJet1->DeltaPhi(*selHadrWJet1)   );
    histo2D["Light_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrWJet1->Pt(),    hadrWJet1->DeltaPhi(*selHadrWJet1)   );

    //Fill histograms for second light jet!
    if(iEtaCase == 1) whichEtaBin = useEtaBinWJet2;
    histo2D["Light_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrWJet2->E(),    selHadrWJet2->E()     );
    histo2D["Light_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrWJet2->Pt(),   selHadrWJet2->Pt()     );
    histo2D["Light_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrWJet2->Theta(),selHadrWJet2->Theta() );
    histo2D["Light_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrWJet2->E(),    selHadrWJet2->Theta() );
    histo2D["Light_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrWJet2->Pt(),    selHadrWJet2->Theta() );
    histo2D["Light_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrWJet2->Phi(),  selHadrWJet2->Phi()   );
    histo2D["Light_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrWJet2->E(),    selHadrWJet2->Phi()   );
    histo2D["Light_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrWJet2->Pt(),    selHadrWJet2->Phi()   );

    histo2D["Light_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrWJet2->E(),    hadrWJet2->E()     - selHadrWJet2->E()     );
    histo2D["Light_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrWJet2->Pt(),   hadrWJet2->Pt()    - selHadrWJet2->Pt()    );
    histo2D["Light_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrWJet2->Theta(),hadrWJet2->Theta() - selHadrWJet2->Theta() );
    histo2D["Light_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrWJet2->E(),    hadrWJet2->Theta() - selHadrWJet2->Theta() );
    histo2D["Light_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrWJet2->Pt(),   hadrWJet2->Theta() - selHadrWJet2->Theta() );
    histo2D["Light_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrWJet2->Phi(),  hadrWJet2->DeltaPhi(*selHadrWJet2)   );
    histo2D["Light_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrWJet2->E(),    hadrWJet2->DeltaPhi(*selHadrWJet2)   );
    histo2D["Light_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrWJet2->Pt(),   hadrWJet2->DeltaPhi(*selHadrWJet2)   );

    //Fill histograms for hadronic b-jet
    if(iEtaCase == 1) whichEtaBin = useEtaBinHadrB;
    histo2D["BJet_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrBJet->E(),    selHadrBJet->E()     );
    histo2D["BJet_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrBJet->Pt(),   selHadrBJet->Pt()    );
    histo2D["BJet_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrBJet->Theta(),selHadrBJet->Theta() );
    histo2D["BJet_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrBJet->E(),    selHadrBJet->Theta() );
    histo2D["BJet_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrBJet->Pt(),   selHadrBJet->Theta() );
    histo2D["BJet_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrBJet->Phi(),  selHadrBJet->Phi());
    histo2D["BJet_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrBJet->E(),    selHadrBJet->Phi());
    histo2D["BJet_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrBJet->Pt(),   selHadrBJet->Phi());

    histo2D["BJet_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        hadrBJet->E(),    hadrBJet->E()     - selHadrBJet->E()     );
    histo2D["BJet_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      hadrBJet->Pt(),   hadrBJet->Pt()    - selHadrBJet->Pt()    );
    histo2D["BJet_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(hadrBJet->Theta(),hadrBJet->Theta() - selHadrBJet->Theta() );
    histo2D["BJet_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    hadrBJet->E(),    hadrBJet->Theta() - selHadrBJet->Theta() );
    histo2D["BJet_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   hadrBJet->Pt(),   hadrBJet->Theta() - selHadrBJet->Theta() );
    histo2D["BJet_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    hadrBJet->Phi(),  hadrBJet->DeltaPhi(*selHadrBJet)   );
    histo2D["BJet_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      hadrBJet->E(),    hadrBJet->DeltaPhi(*selHadrBJet)   );
    histo2D["BJet_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     hadrBJet->Pt(),   hadrBJet->DeltaPhi(*selHadrBJet)   );

    //Fill histograms for leptonic b-jet
    if(iEtaCase == 1) whichEtaBin = useEtaBinLeptB;
    histo2D["BJet_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        leptBJet->E(),    selLeptBJet->E()     );
    histo2D["BJet_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      leptBJet->Pt(),   selLeptBJet->Pt()    );
    histo2D["BJet_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(leptBJet->Theta(),selLeptBJet->Theta() );
    histo2D["BJet_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    leptBJet->E(),    selLeptBJet->Theta() );
    histo2D["BJet_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   leptBJet->Pt(),   selLeptBJet->Theta() );
    histo2D["BJet_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    leptBJet->Phi(),  selLeptBJet->Phi()   );
    histo2D["BJet_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      leptBJet->E(),    selLeptBJet->Phi()   );
    histo2D["BJet_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     leptBJet->Pt(),   selLeptBJet->Phi()   );

    histo2D["BJet_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        leptBJet->E(),    leptBJet->E()     - selLeptBJet->E()     );
    histo2D["BJet_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      leptBJet->Pt(),   leptBJet->Pt()    - selLeptBJet->Pt()     );
    histo2D["BJet_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(leptBJet->Theta(),leptBJet->Theta() - selLeptBJet->Theta() );
    histo2D["BJet_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    leptBJet->E(),    leptBJet->Theta() - selLeptBJet->Theta() );
    histo2D["BJet_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   leptBJet->Pt(),   leptBJet->Theta() - selLeptBJet->Theta() );
    histo2D["BJet_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    leptBJet->Phi(),  leptBJet->DeltaPhi(*selLeptBJet)   );
    histo2D["BJet_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      leptBJet->E(),    leptBJet->DeltaPhi(*selLeptBJet)   );
    histo2D["BJet_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     leptBJet->Pt(),   leptBJet->DeltaPhi(*selLeptBJet)   );

    //Fill histograms for lepton!
    if(iEtaCase == 1) whichEtaBin = useEtaBinLepton;
    if(isSemiEl){
      histo2D["El_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        lepton->E(),    selLepton->E()     );
      histo2D["El_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      lepton->Pt(),   selLepton->Pt()    );
      histo2D["El_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(lepton->Theta(),selLepton->Theta() );
      histo2D["El_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    lepton->E(),    selLepton->Theta() );
      histo2D["El_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   lepton->Pt(),   selLepton->Theta() );
      histo2D["El_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      lepton->E(),    selLepton->Phi()   );
      histo2D["El_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     lepton->Pt(),   selLepton->Phi()   );
      histo2D["El_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    lepton->Phi(),  selLepton->Phi()   );
      
      histo2D["El_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        lepton->E(),    lepton->E()     - selLepton->E()     );
      histo2D["El_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      lepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
      histo2D["El_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(lepton->Theta(),lepton->Theta() - selLepton->Theta() );
      histo2D["El_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    lepton->E(),    lepton->Theta() - selLepton->Theta() );
      histo2D["El_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   lepton->Pt(),   lepton->Theta() - selLepton->Theta() );
      histo2D["El_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    lepton->Phi(),  lepton->DeltaPhi(*selLepton)   );
      histo2D["El_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      lepton->E(),    lepton->DeltaPhi(*selLepton)   );
      histo2D["El_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     lepton->Pt(),   lepton->DeltaPhi(*selLepton)   );
    }
    else if(isSemiMu){
      float InvEgenMu = 1./lepton->E();
      float InvErecMu = 1./selLepton->E();
      float InvPtgenMu = 1./lepton->Pt();
      float InvPtrecMu = 1./selLepton->Pt();
      histo2D["Mu_RecoInvEVsGenInvE"+EtaBin[whichEtaBin]]->Fill(  InvEgenMu,      InvErecMu         );
      histo2D["Mu_RecoEVsGenE"+EtaBin[whichEtaBin]]->Fill(        lepton->E(),    selLepton->E()    );
      histo2D["Mu_RecoInvPtVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(InvPtgenMu,     InvPtrecMu        );
      histo2D["Mu_RecoPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      lepton->Pt(),   selLepton->Pt()   );
      histo2D["Mu_RecoThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(lepton->Theta(),selLepton->Theta());
      histo2D["Mu_RecoThetaVsGenInvE"+EtaBin[whichEtaBin]]->Fill( InvEgenMu,      selLepton->Theta());
      histo2D["Mu_RecoThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    lepton->E(),    selLepton->Theta());
      histo2D["Mu_RecoThetaVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(InvPtgenMu,     selLepton->Theta());
      histo2D["Mu_RecoThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   lepton->Pt(),   selLepton->Theta());
      histo2D["Mu_RecoPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    lepton->Phi(),  selLepton->Phi()  );
      histo2D["Mu_RecoPhiVsGenInvE"+EtaBin[whichEtaBin]]->Fill(   InvEgenMu,      selLepton->Phi()  );
      histo2D["Mu_RecoPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      lepton->E(),    selLepton->Phi()  );
      histo2D["Mu_RecoPhiVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(  InvPtgenMu,     selLepton->Phi()  );
      histo2D["Mu_RecoPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     lepton->Pt(),   selLepton->Phi()  );

      histo2D["Mu_DiffInvEVsGenInvE"+EtaBin[whichEtaBin]]->Fill(  InvEgenMu,      InvEgenMu       - InvErecMu          );
      histo2D["Mu_DiffEVsGenE"+EtaBin[whichEtaBin]]->Fill(        lepton->E(),    lepton->E()     - selLepton->E()     );
      histo2D["Mu_DiffInvPtVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(InvPtgenMu,     InvPtgenMu      - InvPtrecMu         );
      histo2D["Mu_DiffPtVsGenPt"+EtaBin[whichEtaBin]]->Fill(      lepton->Pt(),   lepton->Pt()    - selLepton->Pt()    );
      histo2D["Mu_DiffThetaVsGenTheta"+EtaBin[whichEtaBin]]->Fill(lepton->Theta(),lepton->Theta() - selLepton->Theta() );
      histo2D["Mu_DiffThetaVsGenInvE"+EtaBin[whichEtaBin]]->Fill( InvEgenMu,      lepton->Theta() - selLepton->Theta() );
      histo2D["Mu_DiffThetaVsGenE"+EtaBin[whichEtaBin]]->Fill(    lepton->E(),    lepton->Theta() - selLepton->Theta() );
      histo2D["Mu_DiffThetaVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(InvPtgenMu,     lepton->Theta() - selLepton->Theta() );
      histo2D["Mu_DiffThetaVsGenPt"+EtaBin[whichEtaBin]]->Fill(   lepton->Pt(),   lepton->Theta() - selLepton->Theta() );
      histo2D["Mu_DiffPhiVsGenPhi"+EtaBin[whichEtaBin]]->Fill(    lepton->Phi(),  lepton->DeltaPhi(*selLepton)         );
      histo2D["Mu_DiffPhiVsGenInvE"+EtaBin[whichEtaBin]]->Fill(   InvEgenMu,      lepton->DeltaPhi(*selLepton)         );
      histo2D["Mu_DiffPhiVsGenE"+EtaBin[whichEtaBin]]->Fill(      lepton->E(),    lepton->DeltaPhi(*selLepton)         );
      histo2D["Mu_DiffPhiVsGenInvPt"+EtaBin[whichEtaBin]]->Fill(  InvPtgenMu,     lepton->DeltaPhi(*selLepton)         );
      histo2D["Mu_DiffPhiVsGenPt"+EtaBin[whichEtaBin]]->Fill(     lepton->Pt(),   lepton->DeltaPhi(*selLepton)         );
    }
  }
}

void TFCreation::CalculateTFFromFile(string fitHistoName, bool useStartValues, int histoNr, bool useROOTClass, bool useStartArray, float startValues[], bool changeFitRange, TFile* file, int whichEtaBin, TFile* readFile){
 
  //Select the correct Eta-bin histogram!
  if(nEtaBins_ == 1) EtaBin[1] = EtaBin[0];
  TH2F* fitHisto = (TH2F*) readFile->Get( ("2D_histograms_graphs/"+fitHistoName+""+EtaBin[whichEtaBin]).c_str() );

  TDirectory* th2dir;
  if(file->GetDirectory("2D_histograms_graphs") == 0) th2dir = file->mkdir("2D_histograms_graphs");
  else{                                               th2dir = file->GetDirectory("2D_histograms_graphs");}
  //Save the 2D histogram used for the fit!
  th2dir->cd(); fitHisto->Write(); file->cd();             
    
  TDirectory* histoFitDir = file->mkdir(fitHisto->GetName());
  histoFitDir->cd();

  caloEnergyFit->SetRange( fitHisto->GetXaxis()->GetXmin(), fitHisto->GetXaxis()->GetXmax() );

  //Initialize the start values if asked
  startValuesArray = startValues;
  if(useStartValues) SetStartValuesDoubleGaussian(histoNr, useStartArray, string(fitHisto->GetName()));         //Can only be done after that doubleGaussianFit is initialized!

  //Choose the correct fit method:
  TObjArray aSlices;
  if(useROOTClass){
    std::cout << " ERROR : Not able to use ROOT Class any more (not adapted to the use of TGraphs ...)" << std::endl;
  }
  else
    FitSliceClassCode(fitHisto, changeFitRange, whichEtaBin);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //   Now histogram with all parameters needs to be fitted with Calorimeter Energy formula   //
  //////////////////////////////////////////////////////////////////////////////////////////////
  std::string fitHistName = string(fitHisto->GetName());
  unsigned int fitHistSize = fitHistName.size();

  for( int ipar = 0; ipar < nParsFit_; ipar++ ){
    
    //Set the name and function formula of this parameter fit correctly!
    caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");          //Only expect the calorimeter behavior for the sigma's of the gaussians! 
    if(ipar == 2){
      if(fitHistName.find("Light") < fitHistSize) caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x"); 
      if(fitHistName.find("BJet") < fitHistSize)  caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x");
    }
    for(int ii = 0; ii < caloEnergyFit->GetNpar(); ii++) caloEnergyFit->SetParName(ii, ( parnames_[ipar]+tostr(ii)).c_str() );

    FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[0];                            
    FitMax_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetN()-1]; //-1 is because bins go from 0 to N (Keep overflow when it is well determined!!)

    //Exclude starting bins for the second gaussian determination
    // --> Have to make sure this is really the second gaussian describing the tails!!
    if(fitHistName.find("Light_DiffEVsGenE") < fitHistSize){
      if(fitHistName.find("Eta_0_") <= fitHistSize && (ipar == 3||ipar == 4) )     FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[1];
      else if(fitHistName.find("Eta_0p375_") <= fitHistSize && (ipar == 3||ipar == 4) ) FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[2];
      else if(fitHistName.find("Eta_0p75") <= fitHistSize && (ipar == 3||ipar == 4) )   FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[1];
    }
    if(fitHistName.find("BJet_DiffEVsGenE") < fitHistSize){
      if( fitHistName.find("Eta_0") <= fitHistSize && (ipar == 3||ipar == 4) ) FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[2];
    }

    //In case of muons the last parameters correspond to the narrow gaussian (not for the non-binned hist)!!
    if(fitHistName.find("Mu_DiffEVsGenE") < fitHistSize){
      if(fitHistName.find("Eta") > fitHistSize && (ipar == 3 || ipar ==4) ) FitMin_[nParsFit_*whichEtaBin+ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[2];
      else if(fitHistName.find("Eta_0p75") < fitHistSize && (ipar == 3 || ipar == 4) ) FitMax_[ipar] = grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetX()[grE_ParamFit[nParsFit_*whichEtaBin+ipar]->GetN()-2];
    } 

    //for(int ii = 0; ii < 3; ii++) caloEnergyFit->SetParName(ii, ( parnames_[ipar]+tostr(ii)).c_str() ); //Name here since different for each doubleGaussian parameter!
    caloEnergyFit->SetName( (string(fitHisto->GetName())+"_"+parnames_[ipar]+"_Fit").c_str() );
    caloEnergyFit->SetTitle((ParName_[ipar]+" for "+string(fitHisto->GetName())).c_str());

    grE_ParamFit[nParsFit_*whichEtaBin+ipar]->Fit(caloEnergyFit, "Q", "",FitMin_[nParsFit_*whichEtaBin+ipar], FitMax_[nParsFit_*whichEtaBin+ipar]);
    caloEnergyFit->SetRange(0,250);
    AllCaloEnergyFits[nParsFit_*whichEtaBin+ipar] = *caloEnergyFit;       //caloEnergyFit is a pointer, but each member of the array should point to the corresponding value of the TF1!
    grE_ParamFit[nParsFit_*whichEtaBin+ipar]->Write();
//    caloEnergyFit->SaveAs(("TFInformation/Plots/"+string(caloEnergyFit->GetName())+".pdf").c_str());
    caloFitFile->cd(); caloEnergyFit->Write(); file->cd(); histoFitDir->cd(); 
  }
  //PlotDlbGaus(fitHisto,file);
							  
}

void TFCreation::FitSliceClassCode(TH2F* histoFit, bool ChangeFitRange, int etaBin){
  //------------------------------------------------------------------------------------------//
  // Main difference with the Root class FitSlicesY() is the plotting of histograms !        
  // In the Root class the distribution of each histogram is not given!
  // --> Useful to use the own code when needing control histograms
  //
  // Other difference between the two codes have been removed!
  // Originally the treatment of the overflow bin was different, but is now made similar!
  // Create one histogram for each function parameter -> 5 histograms for each 2D plot

  TH1D* h_chi2 = new TH1D( (string(histoFit->GetName())+"_chi2").c_str(), (string(histoFit->GetName())+": #chi^{2} distribution for "+string(doubleGaussianFit->GetExpFormula())).c_str(), histoFit->GetXaxis()->GetNbins(), histoFit->GetXaxis()->GetXmin(), histoFit->GetXaxis()->GetXmax() );

  //Loop on all bins in X, generate a projection along Y and fit each bin separately!
  int cut = 0; // require a minimum number of bins in the slice to be filled --> Should this ever be larger than 0 ??
  int nbins = histoFit->GetXaxis()->GetNbins();
  int nActiveBins = 0;
  vector<double> xValue, yValue, xError, yError;
  std::string histoName = string(histoFit->GetName());
  int lastBin = nbins+1;
  
  //Combine some bins!
  int binStart[10] = {50};
  int binEnd[10] = {50};
  int nCombBins = 0;
  if(     histoName.find("BJet_DiffEVsGenE") <= histoName.size() ){ 
    if(histoName.find("Eta") > histoName.size() ){            binStart[0] = 14; binEnd[0] = 15; binStart[1] = 16; binEnd[1] = 17; binStart[2] = 18; binEnd[2] = 19;}
    else if(histoName.find("Eta_0_")    <= histoName.size()){ binStart[0] = 12; binEnd[0] = 14;}
    else if(histoName.find("Eta_0p375") <= histoName.size()){ binStart[0] = 14; binEnd[0] = 15; binStart[1] = 16; binEnd[1] = 18;}
    else if(histoName.find("1p45") <= histoName.size()     ){ binStart[0] = 13; binEnd[0] = 14; binStart[1] = 15; binEnd[1] = 16; binStart[2] = 17; binEnd[2] = 18;}
  }
  else if(histoName.find("Light_DiffEVsGenE") <= histoName.size()){
    if(histoName.find("Eta") > histoName.size() ){            binStart[0] = 1; binEnd[0] = 2;}
    if(histoName.find("Eta_0_")         <= histoName.size()){ binStart[0] = 1; binEnd[0] = 2; binStart[1] = 13; binEnd[1] = 14; binStart[2] = 15; binEnd[2] = 17; lastBin = 15;}
    else if(histoName.find("Eta_0p375") <= histoName.size()){ binStart[0] = 1; binEnd[0] = 2; binStart[1] = 13; binEnd[1] = 14; binStart[2] = 15; binEnd[2] = 17; lastBin = 15;}
    else if(histoName.find("Eta_0p75")  <= histoName.size()){ binStart[0] = 2; binEnd[0] = 3; binStart[1] = 15; binEnd[1] = 17; lastBin = 15;}
    else if(histoName.find("Eta_1p45")  <= histoName.size()){ binStart[0] = 6; binEnd[0] = 8; binStart[1] = 9;  binEnd[1] = 10;}
  }
  else if(histoName.find("El_DiffEVsGenE") <= histoName.size() ){
    //if(histoName.find("Eta") > histoName.size() ){ binStart[0] = 15; binEnd[0] = 16; binStart[1] = 17; binEnd[1] = 19; binStart[2] = 20; binEnd[2] = 25;}
    if(histoName.find("Eta_0_") <= histoName.size() ){ binStart[0] = 13; binEnd[0] = 14; binStart[1] = 15; binEnd[1] = 17;}
  }
  else if(histoName.find("Mu_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){   binStart[0] = 13; binEnd[0] = 15;}
    else if(histoName.find("Eta_0_") <= histoName.size() ){   binStart[0] = 8; binEnd[0] = 10;}
    else if(histoName.find("Eta_0p375") <= histoName.size() ){binStart[0] = 9; binEnd[0] = 10; binStart[1] = 11; binEnd[1] = 12;}
    else if(histoName.find("Eta_0p75") <= histoName.size() ){ binStart[0] = 2; binEnd[0] = 3;  binStart[1] = 10; binEnd[1] = 11; binStart[2] = 12; binEnd[2] = 14;}
    else if(histoName.find("Eta_1p45") <= histoName.size() ){ binStart[0] = 5; binEnd[0] = 6;  binStart[1] = 11; binEnd[1] = 12; binStart[2] = 13; binEnd[2] = 15;}
  }

  //Skip some bins!
  int binToSkip[10] = {50};
  int nSkippedBins = 0;
  if( histoName.find("BJet_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){ binToSkip[0] = 1; lastBin = 18;}
    if(histoName.find("Eta_0_")    <= histoName.size() ){ binToSkip[0] = 15; binToSkip[1] = 16; binToSkip[2] = 17; binToSkip[3] = 18; binToSkip[4] = 19; lastBin = 12;}
    else if(histoName.find("Eta_0p75")  <= histoName.size() )  binToSkip[0] = 1;
    else if(histoName.find("Eta_1p45")  <= histoName.size() ){ binToSkip[0] = 1;  binToSkip[1] = 2; binToSkip[2] = 3; binToSkip[3] = 4; binToSkip[4] = 5;}
  }
  else if(histoName.find("Light_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta_0p75") <= histoName.size() ){ binToSkip[0] = 1;}
    else if(histoName.find("Eta_1p45") <= histoName.size() ){binToSkip[0] = 1; binToSkip[1] = 2; binToSkip[2] = 3; binToSkip[3] = 4; binToSkip[4] = 5;}
  }
  else if(histoName == "El_DiffEVsGenE_Eta_0_0p375"){binToSkip[0] = 18; binToSkip[1] = 19; binToSkip[2] = 20; binToSkip[3] = 21; binToSkip[4] = 22; binToSkip[5] = 23; binToSkip[6] = 24; binToSkip[7] = 25; binToSkip[8] = 26; lastBin = 15;}
  else if(histoName.find("Mu_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta_0_") <= histoName.size() ){        binToSkip[0] = 11; binToSkip[1] = 12; binToSkip[2] = 13; binToSkip[3] = 14; binToSkip[4] = 15; lastBin = 8;}
    else if(histoName.find("Eta_0p375") <= histoName.size() ){binToSkip[0] = 1; binToSkip[1] = 13; binToSkip[2] = 14; binToSkip[3] = 15; lastBin = 11;}
    else if(histoName.find("Eta_0p75") <= histoName.size() ){ binToSkip[0] = 1; }
    else if(histoName.find("Eta_1p45") <= histoName.size() ){ binToSkip[0] = 1; binToSkip[1] = 2; binToSkip[2] = 3; binToSkip[3] = 4; lastBin = 13;}
  }
  
  //Vector that will contain the string information to add in the LaTeX file containing all individual fitted histograms!
  indivFitHistos.clear();
  indivFitHistos.push_back("\\newpage \n\\subsection{Individual fit histograms} \n");

  for(int bin=1;bin <= nbins+1;bin ++) {
    string projection_title = "sliceYbin"+tostr(bin)+"_"+string(histoFit->GetName());

    TH1D *hp = 0;
    double binLow, binHigh;
    if(bin == binStart[nCombBins]){
      for(int ii = 1; ii <= binEnd[nCombBins] - binStart[nCombBins]; ii ++) projection_title = "sliceYbin"+tostr(bin)+"And"+tostr(bin+ii)+"_"+string(histoFit->GetName());
      binLow = histoFit->GetXaxis()->GetBinLowEdge(bin);
      binHigh = histoFit->GetXaxis()->GetBinUpEdge(binEnd[nCombBins]);
      hp = histoFit->ProjectionY(projection_title.c_str(),bin,binEnd[nCombBins],"e");
    }
    else if(bin > binStart[nCombBins-1] && bin <= binEnd[nCombBins-1] && nCombBins != 0){
      //std::cout << " Excluded histogram for bin : " << bin << std::endl;
    }
    else if(bin == binToSkip[nSkippedBins]){
      nSkippedBins++;
      //std::cout << " Skipping bin : " << bin  << std::endl;
    }
    else{
      hp = histoFit->ProjectionY(projection_title.c_str(),bin,bin,"e");
      binLow = histoFit->GetXaxis()->GetBinLowEdge(bin);
      binHigh = histoFit->GetXaxis()->GetBinUpEdge(bin);
      if(histoName == "Mu_DiffEVsGenE" && bin == 15) hp->SetAxisRange(-8,8,"X");
    }

    //Histogram doesn't have any memory space ...
    if(hp == 0){ continue;}
    if( float(hp->GetEntries()) <= 0){ delete hp; continue;} //|| float(hp->GetEntries()) < cut) {delete hp; continue;}

    doubleGaussianFit->SetName((projection_title+"Fitted").c_str());
    hp->SetTitle(("Y-axis projection of #DeltaE (gen-reco) for E_{gen} #in ["+tostr(binLow)+","+tostr(binHigh)+"]").c_str());
    hp->GetXaxis()->SetTitle("#DeltaE (gen-reco)"); hp->GetYaxis()->SetTitle("Normalized entries"); hp->GetYaxis()->SetTitleOffset(1.3);
    double ActualFitRange[2] = {(double) (histoFit->GetYaxis())->GetXmin(), (double) (histoFit->GetYaxis())->GetXmax()};
    if(ChangeFitRange == true){
      std::vector<double> fitBin = SetFitRange(histoName, bin, ActualFitRange);
      ActualFitRange[0] = fitBin[0]; ActualFitRange[1] = fitBin[1];
    }

    doubleGaussianFit->SetParLimits(0, -10,  10);
    doubleGaussianFit->SetParLimits(1, 0.0,  15);
    doubleGaussianFit->SetParLimits(2, 0.0, 0.2);
    doubleGaussianFit->SetParLimits(3, -20,  20);
    doubleGaussianFit->SetParLimits(4,   0,  30);
    if( histoName.find("BJet_DiffEVsGenE") <= histoName.size() || histoName.find("Light_DiffEVsGenE") <= histoName.size()){
      if(bin > 4) {
        doubleGaussianFit->SetParLimits(1, 0.0,  25);
        doubleGaussianFit->SetParLimits(2, 0.0, 1.0);
        doubleGaussianFit->SetParLimits(3, -50,  50);
        doubleGaussianFit->SetParLimits(4,   0,  80);
      }
    }
    hp->Fit(doubleGaussianFit,"Q","",ActualFitRange[0],ActualFitRange[1]);

    int npfits = doubleGaussianFit->GetNumberFitPoints();              //WHAT IS THIS .... ???
    if(npfits > nParsFit_ && npfits >= cut) {
      //Fill the TGraph for each parameter with the obtained Fit parameter and its uncertainty
      //--> Each bin in this histogram represents a bin range in x-axis of considered 2D histogram!
      for(int ipar=0; ipar<nParsFit_; ipar++ ){
        if( !( (histoName == "Light_DiffPtVsGenPt" || histoName.find("Light_DiffPtVsGenPt_Eta_0") <= histoName.size() ) && bin == 2) &&
            !( histoName == "Mu_DiffInvPtVsGenInvPt_Eta_1p45_2p5" && ( bin == 11 || bin == 10) ) &&
            !( histoName == "Mu_DiffInvPtVsGenInvPt" && bin == 1 ) &&                            
            !( histoName.find("Light_DiffPtVsGenPt_Eta_1p45") <= histoName.size() && bin == 7)   &&
            !( histoName.find("Mu_DiffInvPtVsGenInvPt_Eta_0p75") <= histoName.size() && bin == 11) ){ // &&

          if(ipar == 0) nActiveBins += 1;
 
          if( bin == binStart[nCombBins] ){
            xValue.push_back( (histoFit->GetXaxis()->GetBinCenter(bin)+histoFit->GetXaxis()->GetBinCenter(bin+1))/2);
            float binCenter = 0;
            for(int ii = 0; ii <= binEnd[nCombBins] - binStart[nCombBins]; ii++) binCenter += histoFit->GetXaxis()->GetBinWidth(bin+ii);
            binCenter = binCenter/(binEnd[nCombBins] - binStart[nCombBins]+1);         //Need this 1/n factor since the width is from bin low edge to bin high edge!
            xError.push_back( binCenter); 
            yValue.push_back(doubleGaussianFit->GetParameter(ipar));
            yError.push_back(doubleGaussianFit->GetParError(ipar));
            if(ipar == nParsFit_-1) nCombBins++;
          }
          else{
            xValue.push_back(histoFit->GetXaxis()->GetBinCenter(bin+1/2));
            xError.push_back(histoFit->GetXaxis()->GetBinWidth(bin)/2 );
            yValue.push_back(doubleGaussianFit->GetParameter(ipar));
            yError.push_back(doubleGaussianFit->GetParError(ipar));
          }
	}

        //Store the TGraph once all the bins have been considered!
        if(bin == lastBin){

          //Select the information for each parameter separately and store it in a new vector!
          vector<double> xValue_, yValue_, xError_, yError_;
          xValue_.clear(); yValue_.clear(); yError_.clear();
          for(int iAct = 0; iAct < nActiveBins; iAct++){
            xValue_.push_back(xValue[ipar+nParsFit_*iAct]);
            //xError_.push_back(xError[ipar+nParsFit_*iAct]);    //Problem since fit really sees it as an error ...
            yValue_.push_back(yValue[ipar+nParsFit_*iAct]);
            yError_.push_back(yError[ipar+nParsFit_*iAct]);
          }

          grE_ParamFit[nParsFit_*etaBin+ipar] = new TGraphErrors(nActiveBins, &xValue_[0], &yValue_[0], 0, &yError_[0]);
          grE_ParamFit[nParsFit_*etaBin+ipar]->GetXaxis()->SetTitle(histoFit->GetXaxis()->GetTitle());
          grE_ParamFit[nParsFit_*etaBin+ipar]->SetTitle((ParName_[ipar]+" for "+histoName).c_str());
          grE_ParamFit[nParsFit_*etaBin+ipar]->SetName((string(histoName+"_"+parnames_[ipar])).c_str());
          grE_ParamFit[nParsFit_*etaBin+ipar]->Draw("AP");
          grE_ParamFit[nParsFit_*etaBin+ipar]->SetMarkerStyle(1);
        }
      }
      //Save hchi2 histogram!
      h_chi2->Fill(histoFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetChisquare()/(npfits-nParsFit_));
    
//      if( bin == nbins/2 ){
//	if(abs(doubleGaussianFit->GetParameter(1)) < abs(doubleGaussianFit->GetParameter(4))){ NarrowGaus[0] = 0; NarrowGaus[1] = 1; NarrowGaus[2] = 2; WideGaus[0] = 3; WideGaus[1] = 4; WideGaus[2] = 5;}
//	else{  NarrowGaus[0] = 3; NarrowGaus[1] = 4; NarrowGaus[2] = 5; WideGaus[0] = 0; WideGaus[1] = 1; WideGaus[2] = 2;}
//      }
        
    }
    TCanvas* canv_IndivFit = new TCanvas("canv_IndivFit","canv_IndivFit");
    canv_IndivFit->SetName(hp->GetName()); canv_IndivFit->SetTitle(hp->GetTitle());
    canv_IndivFit->cd();
    hp->Draw();
    int n = 200;
    double x[200], yDbl[200], ySngl[200];
    for(int ii = 0; ii < n; ii ++){
      x[ii] = -100+ii;
      yDbl[ii] = doubleGaussianFit->GetParameter(5)*(1/(TMath::Sqrt(2*TMath::Pi())*(TMath::Sqrt(TMath::Power(doubleGaussianFit->GetParameter(1),2))+TMath::Sqrt(TMath::Power(doubleGaussianFit->GetParameter(4),2))*doubleGaussianFit->GetParameter(2))))*(TMath::Exp(-TMath::Power((x[ii]-doubleGaussianFit->GetParameter(0)),2)/(2*TMath::Power(doubleGaussianFit->GetParameter(1),2)))+doubleGaussianFit->GetParameter(2)*TMath::Exp(-TMath::Power((x[ii]-doubleGaussianFit->GetParameter(3)),2)/(2*TMath::Power(doubleGaussianFit->GetParameter(4),2))));
      ySngl[ii] = doubleGaussianFit->GetParameter(5)*(1/(TMath::Sqrt(2*TMath::Pi())*(TMath::Sqrt(TMath::Power(doubleGaussianFit->GetParameter(1),2)))))*(TMath::Exp(-TMath::Power((x[ii]-doubleGaussianFit->GetParameter(0)),2)/(2*TMath::Power(doubleGaussianFit->GetParameter(1),2))));;
    }
/*    if( histoName.find("Light_DiffEVsGenE") <= histoName.size() ){ // || histoName.find("Light_DiffEVsGenE") <= histoName.size() ){
      std::cout << " ****  Looking at bin : " << bin << " **** " << endl;
      if(histoName.find("Eta") > histoName.size() ){
        std::cout << " Parameters of the double Gaussian fit are : " << endl;
        std::cout << "  * Narrow (mean & width)  : " << doubleGaussianFit->GetParameter(0) << " , " << doubleGaussianFit->GetParameter(1) << endl;
        std::cout << "  * Relative constant      : " << doubleGaussianFit->GetParameter(2) << endl;
        std::cout << "  * Wide (mean & width)    : " << doubleGaussianFit->GetParameter(3) << " , " << doubleGaussianFit->GetParameter(4) << endl;
        std::cout << "  * Amplitude              : " << doubleGaussianFit->GetParameter(5) << endl;
      }
    }*/
    TGraph* gr_TFDistDblG = new TGraph(n,x,yDbl);
    gr_TFDistDblG->Draw("C");
    gr_TFDistDblG->SetLineColor(6); gr_TFDistDblG->SetLineStyle(2);
    TGraph* gr_TFDistSnglG = new TGraph(n,x,ySngl);
    gr_TFDistSnglG->Draw("C");
    gr_TFDistSnglG->SetLineColor(3); gr_TFDistSnglG->SetLineStyle(10);
    TLegend* leg_Canv = new TLegend(0.1, 0.7, 0.4, 0.9);
    leg_Canv->AddEntry(hp,"Actual ttbar entries fitted with Dbl gaus","lpe");
    leg_Canv->AddEntry(gr_TFDistDblG,"Overal Dbl gauss function","l");
    leg_Canv->AddEntry(gr_TFDistSnglG,"Narrow gaus of Dbl gaus","l");
    leg_Canv->Draw();
    canv_IndivFit->Write();
    if(histoName.find("BJet_DiffEVsGenE") <= histoName.size() || histoName.find("Light_DiffEVsGenE") <= histoName.size() ){
      canv_IndivFit->SaveAs(("TFInformation/Plots/"+string(canv_IndivFit->GetName())+".pdf").c_str());
      indivFitHistos.push_back(("\\includegraphics[width = 0.50 \\textwidth]{Plots/"+string(canv_IndivFit->GetName())+"} \n").c_str());
    }
    //hp->Write();
    delete hp;
  }//loop over bins!
  h_chi2->Write();
  delete h_chi2;
}

std::vector<double> TFCreation::SetFitRange(std::string histoName, unsigned int iBin, double startRange[2]){
   
  std::vector<double> BinnedFitRange;
  BinnedFitRange.clear();
  double FitRangeBinNeg = startRange[0], FitRangeBinPos = startRange[1];

  if(histoName.find("BJet_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FitRangeNeg[9] = {-43, -43, -42, -52, -55, -55, -60, -64, -66}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[9] = { 20,  20,  27,  33,  42,  45,  57,  62,  68}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    if(histoName.find("Eta_0_") <= histoName.size()){
      double FitRangeNeg[11] = {-13, -20, -26, -30, -30, -30, -30, -32, -34, -34, -35}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[11] = {  7,  15,  15,  20,  23,  25,  25,  28,  30,  32,  32}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p375") <= histoName.size()){
      double FitRangeNeg[7] = {-16, -20, -30, -32, -33, -35, -35}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[7] = {  4,  11,  21,  28,  30,  32,  35}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p75") <= histoName.size()){
      double FitRangeNeg[6] = {-10, -16, -30, -37, -40, -40}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[6] = { 10,   7,  13,  24,  36,  42}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1p45") <= histoName.size()){
      double FitRangeNeg[10] = {-10, -10, -10, -10, -10, -32, -40, -45, -45, -47}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[10] = { 10,  10,  10,  10,  10,  18,  30,  40,  47,  48}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
  } //End of BJet_DiffE histo

  if(histoName.find("El_DiffEVsGenE") <= histoName.size() ){
    if( histoName.find("Eta") > histoName.size() ){
      double FitRangeNeg[16] = { -3,  -4, -4,  -4,  -4, -4.5, -4.5,  -4, -4.5, -4.5,  -5, -5, -5.5, -6, -6, -6}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[16] = {2.5, 2.5,  3, 3.5, 3.5,    4,  4.5, 4.5,    5,  5.5, 5.5,  6,    6,  6,  6,  6}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else{
      double FitRangeNeg[15] = {-1.3, -1.5, -2.4, -2.5, -2.5, -2.7, -2.7, -2.7, -3.0, -2.5, -2.5, -2.5, -3.0, -3.0, -3.0}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[15] = { 0.9,  1.7,  1.5,  1.6,  1.8,  2.3,  2.7,  2.5,  2.5,  3.0,  3.0,  3.0,  3.5,  3.5,  3.7}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
  } //End of El_DiffE GenE histo

  if(histoName.find("Light_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FitRangeNeg[7] = {-49, -50, -45, -50, -53, -61, -70}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[7] = {  9,  11,  17,  22,  28,  35,  41}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size()){
      double FitRangeNeg[8] = {-16, -10, -21, -23, -25, -30, -32, -35}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[8] = {  6,   6,  14,  20,  21,  26,  33,  35}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p375") <= histoName.size()){
      double FitRangeNeg[8] = {-15, -13, -17, -25, -28, -32, -32, -34}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[8] = {  4,   5,  11,  18,  24,  28,  32,  33}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p75") <= histoName.size()){
      double FitRangeNeg[9] = {-13, -21, -20, -21, -30, -33, -36, -40, -40}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[9] = {  4,   6,   5,   9,  17,  25,  34,  40,  43}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1p45") <= histoName.size() ){
      double FitRangeNeg[17] = {-10, -10, -10, -10, -10, -25, -25, -25, -20, -20, -20, -17, -17, -19, -20, -18, -18}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[17] = {  4,   4,   4,   4,   4,  10,  10,  10,  18,  12,  22,  26,  32,  35,  35,  35,  32}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
  } //End of Light_DiffE GenE histo

  if(histoName.find("Mu_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size()){
      double FitRangeNeg[8] = {-0.5, -1.0, -1.8, -2.3, -3.2, -4.0, -4.5, -5.5}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[8] = { 0.7,  1.0,  1.3,  1.7,  2.2,  2.8,  3.6,  4.2}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() ){
      double FitRangeNeg[5] = {-1.0, -1.5, -2.0, -2.4, -3.0}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[5] = { 1.0,  1.3,  2.2,  3.0,  3.3}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p375") <= histoName.size() ){
      double FitRangeNeg[5] = {-1.0, -1.3, -1.6, -2.2, -3.4}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[5] = { 1.0,  1.4,  2.0,  2.7,  3.2}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0p75") <= histoName.size() ){
      double FitRangeNeg[6] = {-1.0, -2.5, -1.6, -3.0, -3.6, -4.0}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[6] = { 1.0,  2.5,  1.6,  3.2,  3.6,  4.0}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1p45") <= histoName.size() ){
      double FitRangeNeg[12] = {-0.5, -0.5, -0.5, -0.5, -3.5, -3.5, -4.5, -4.5, -6.0, -6.5, -6.5, -7.5}; if(iBin <= sizeof(FitRangeNeg)/sizeof(FitRangeNeg[0])) FitRangeBinNeg = FitRangeNeg[iBin-1];
      double FitRangePos[12] = { 0.5,  0.5,  0.5,  0.5,  4.0,  3.0,  5.0,  6.0,  7.0,  7.0,  8.0,  8.0}; if(iBin <= sizeof(FitRangePos)/sizeof(FitRangePos[0])) FitRangeBinPos = FitRangePos[iBin-1];
    }
  }

  BinnedFitRange.push_back(FitRangeBinNeg); BinnedFitRange.push_back(FitRangeBinPos);

  return  BinnedFitRange;
}
void TFCreation::SetStartValuesDoubleGaussian(int whichHisto, bool useStartArray, std::string histoName){

  if(useStartArray == true){
    for(int ipar = 0; ipar < nParsFit_; ipar++){
      doubleGaussianFit->SetParameter(ipar, startValuesArray[ipar]);
      if(histoName.find("_Eta_") <= histoName.size() && ipar == 2 ){ doubleGaussianFit->SetParameter(ipar, startValuesArray[ipar]/4.); }
    }
  }
  else{
    //if(whichHisto==1 || whichHisto==4 || whichHisto == 7){ // for E transfer function of JETS (and elec -- added as test ...)
    //  float StartValues[] = {-8,18,0,0,8.6};          //First three values are for the first broad gaussian (central, sigma and constant value respectively)
    //                                                       //Second three values are the same for the second narrow gaussian
    //  for(int ipar = 0; ipar < nParsFit_; ipar++)
    //	doubleGaussianFit->SetParameter(ipar,StartValues[ipar]);
    //}
    if(whichHisto==1 || whichHisto==4 || whichHisto == 7){ // for Pt transfer function of JETS (and elec -- added as test ...)
      float StartValues[] = {-8,18,0,0,8.6};          //First three values are for the first broad gaussian (central, sigma and constant value respectively)
                                                           //Second three values are the same for the second narrow gaussian
      for(int ipar = 0; ipar < nParsFit_; ipar++)
	doubleGaussianFit->SetParameter(ipar,StartValues[ipar]);
    }
    else if (whichHisto==0 || whichHisto==2 || whichHisto==3 || whichHisto==5 || whichHisto == 6 || whichHisto == 8) { //for theta and phi transfer functions of JETS (and elec)
      float StartValues[] = {0,0.038,0,0.004,0.011};
      for(int ipar = 0; ipar < nParsFit_; ipar++)
	doubleGaussianFit->SetParameter(ipar, StartValues[ipar]);
    }
    //else if (whichHisto==10){ //for 1/E transfer function of muons
    //  float StartValues[] = {-0.0008,0.001,0,-0.0001,0.0001};
    //  for(int ipar = 0; ipar < nParsFit_; ipar++)
    //	doubleGaussianFit->SetParameter(ipar, StartValues[ipar]);
    //}
    else if (whichHisto==10){ //for 1/Pt transfer function of muons
      float StartValues[] = {-0.0008,0.001,0,-0.0001,0.0001};
      for(int ipar = 0; ipar < nParsFit_; ipar++)
	doubleGaussianFit->SetParameter(ipar, StartValues[ipar]);
    }
    else if (whichHisto==9 || whichHisto==11) { //for theta, phi transfer function of muons
      float StartValues[] = {0.0,0.01,0,0,0.001};
      for(int ipar = 0; ipar < nParsFit_; ipar++)
	doubleGaussianFit->SetParameter(ipar, StartValues[ipar]);
    }
  }
} 

void TFCreation::WriteTF(ostream &myTFTable, ostream &myTransferCard, ostream &myTF, ostream &myLaTeX, std::string kinVar, std::string partName, int TFDep){ 

  //Simple text files which can be used to copy the fit formula (with parameters) to a ROOT file
  //ofstream myFitFormula;
  //myFitFormula.open("TFInformation/TF_user.dat");

  std::string pVar[2] = {"p(0)","pt(p)"};
  std::string pexpVar[2] = {"pexp(0)","pt(pexp)"};

  //Is Pt or E dependent considered
  int whichDep = TFDep;

  string TFDependency[3]  =     {"","*dsqrt("+pVar[whichDep]+")",   "*"+pVar[whichDep]   };
  string TFDependencyConst[6] = {"","*"+pVar[whichDep],             "*"+pVar[whichDep]+"**2","*"+pVar[whichDep]+"**3","*"+pVar[whichDep]+"**4","*"+pVar[whichDep]+"**5"};
  string WidthDependency[3] =   {"","*dsqrt("+pexpVar[whichDep]+")","*"+pexpVar[whichDep]};
  //if(partName == "muon" && whichDep == 1){
  //  TFDependencyWidth[1] = "*dsqrt(1d0/"+pVar[whichDep]+")"; TFDependencyWidth[2] = "*1d0/"+pVar[whichDep]+")";
  //  TFDependency[1] = "*1d0/"+pVar[whichDep]; TFDependency[2] = "*1d0/"+pVar[whichDep]+"**2"; TFDependency[3] = "*1d0/"+pVar[whichDep]+"**3"; TFDependency[4] = "*1d0/"+pVar[whichDep]+"**4)";
  //  WidthDependency[1] = "*dsqrt(1d0/"+pexpVar[whichDep]+")"; WidthDependency[2] = "*1d0/"+pexpVar[whichDep]+")";
  //}

  ostream *TransferCard = &myTransferCard;
  ostream *TF = &myTF;
  ostream *LaTeX = &myLaTeX;

  *LaTeX << "\n\\newpage \n\\section{"<< partName << " for " << kinVar <<"}" << endl;

  int dummyCounter = 0;
  vector<std::string> WidthText;
  for(int iEta = 1; iEta <= nEtaBins_; iEta++){

    if(nEtaBins_ != 1){
      //Need this eta separation both for the TF itself and for the width!
      if(iEta == 1)     { WidthText.push_back("\n      ENDIF \n    </tf> \n    <width>"); WidthText.push_back("\n      IF( ABS(eta(pexp)) .LE. 0.375) THEN "); *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 2){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 0.375 .AND. ABS(eta(pexp)) .LE. 0.75) THEN "); *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 3){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 0.75 .AND. ABS(eta(pexp)) .LE. 1.45) THEN ");  *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 4){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 1.45 .AND. ABS(eta(pexp)) .LE. 2.5) THEN ");   *TF << WidthText[WidthText.size()-1];}
      
      if(iEta != 1) *LaTeX << " \\caption{Fit distributions for " << EtaValues[iEta-1] << " $<$ $\\vert \\eta \\vert$ $<$ " << EtaValues[iEta] << "}\n \\end{figure}\n \\newpage \n";
      *LaTeX << " \\subsection{Eta between "<< EtaValues[iEta] << " and " << EtaValues[iEta+1] << "}\n";
    }
    else WidthText.push_back("\n    </tf> \n    <width>");
        
    *LaTeX << " \\begin{figure}[h!b] " << endl;
    std::string provFormula = "", provMin = "", provMax = "";
    for(int ipar = 0; ipar < nParsFit_-1; ipar++){
      int NrConsideredCaloPars = AllCaloEnergyFits[iEta*nParsFit_+ipar].GetNpar();

      TCanvas *canv = new TCanvas(("canv_"+partName+""+EtaBin[iEta]).c_str(),"canvas");
      canv->cd();
      grE_ParamFit[iEta*nParsFit_+ipar]->Draw("AP");
      grE_ParamFit[iEta*nParsFit_+ipar]->GetXaxis()->SetTitle(("Generator parton energy (GeV) -- Fit between "+tostr(FitMin_[iEta*nParsFit_+ipar])+" and "+tostr(FitMax_[iEta*nParsFit_+ipar])).c_str());
      AllCaloEnergyFits[iEta*nParsFit_+ipar].SetRange(FitMin_[iEta*nParsFit_+ipar], FitMax_[iEta*nParsFit_+ipar]);
      AllCaloEnergyFits[iEta*nParsFit_+ipar].Draw("same");
      canv->SaveAs(("TFInformation/Plots/"+string(AllCaloEnergyFits[iEta*nParsFit_+ipar].GetName())+".pdf").c_str());
      delete canv;
      *LaTeX << "  \\includegraphics[width = 0.45 \\textwidth]{Plots/" << AllCaloEnergyFits[iEta*nParsFit_+ipar].GetName() << ".pdf}";
      if(ipar == 1 || ipar == 3) *LaTeX << " \\\\ " << endl;
      else                       *LaTeX << " " << endl;

      float FitMinValue = 0, FitMaxValue = 0;
      for(int icalpar = 0; icalpar < NrConsideredCaloPars; icalpar++){
	dummyCounter++;
	if(icalpar == 0) myTFTable<<ParName_[ipar]<<" & $a_{" <<ipar <<icalpar <<"}$ = "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<<" $\\pm$ "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParError(icalpar);
	else             myTFTable<<                " & $a_{" <<ipar <<icalpar <<"}$ = "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<<" $\\pm$ "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParError(icalpar);

	*TransferCard<< dummyCounter << "     " << AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<< "     # " << ParName_[ipar] << endl;

	if(icalpar==0){
//          *TF << "\n        prov"<<ipar+1<<"=(#"<<dummyCounter;
          provFormula += "\n          prov"+tostr(ipar+1)+"=(#"+tostr(dummyCounter);
          if(ipar == 1 || ipar == 4){
            WidthText.push_back("\n        prov"+tostr(ipar+1)+"=(#"+tostr(dummyCounter));
          }
        }

	if(icalpar != 0 && ipar == 2){ provFormula += "+#"+tostr(dummyCounter)+TFDependencyConst[icalpar];}
	else if(icalpar != 0){         provFormula += "+#"+tostr(dummyCounter)+TFDependency[icalpar];     
          if(ipar == 1 || ipar == 4) WidthText.push_back("+#"+tostr(dummyCounter)+WidthDependency[icalpar]);
        }
        if( icalpar == NrConsideredCaloPars-1){ provFormula += " )"; if(ipar == 1 || ipar == 4) WidthText.push_back(" )");}

        //Now add the cut-off!
        if(     NrConsideredCaloPars == 3 && icalpar == 1){
          FitMinValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*sqrt(FitMin_[iEta*nParsFit_+ipar]);
          FitMaxValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*sqrt(FitMax_[iEta*nParsFit_+ipar]);
        }
        else if(NrConsideredCaloPars == 3 && icalpar == 2){
          FitMinValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*FitMin_[iEta*nParsFit_+ipar];
          FitMaxValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*FitMax_[iEta*nParsFit_+ipar];
        }
        else{
          FitMinValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*pow(FitMin_[iEta*nParsFit_+ipar],icalpar);
          FitMaxValue += AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)*pow(FitMax_[iEta*nParsFit_+ipar],icalpar);
        }

        if(icalpar == NrConsideredCaloPars-1){
          provMin += "\n          prov"+tostr(ipar+1)+" = "+tostr(AllCaloEnergyFits[iEta*nParsFit_+ipar].Eval(FitMin_[iEta*nParsFit_+ipar]));
          provMax += "\n          prov"+tostr(ipar+1)+" = "+tostr(AllCaloEnergyFits[iEta*nParsFit_+ipar].Eval(FitMax_[iEta*nParsFit_+ipar],icalpar));
//          *TF << "\n        IF( p(0) .LE. " << tostr(FitMin_[iEta*nParsFit_+ipar]) << ") THEN " << endl;
//          *TF << "          prov"<< ipar+1 << "=" << tostr(AllCaloEnergyFits[iEta*nParsFit_+ipar].Eval(FitMin_[iEta*nParsFit_+ipar])) << "\n        ENDIF " << endl;
//          *TF << "\n        IF( p(0) .GE. " << tostr(FitMax_[iEta*nParsFit_+ipar]) << ") THEN " << endl;
//          *TF << "          prov"<<ipar+1 << " = " << tostr(AllCaloEnergyFits[iEta*nParsFit_+ipar].Eval(FitMax_[iEta*nParsFit_+ipar],icalpar)) << "\n        ENDIF " << endl;
        }
      }
      myTFTable << "\\\\" << endl;
      //delete grE_ParamFit[ipar];   //--> Will probably not work for more eta-bins ... !!
    }
    *TF << "\n        IF( p(0) .LE. " << tostr(FitMin_[iEta*nParsFit_]) << ") THEN ";
    *TF << provMin;
    *TF << "\n        ELSE IF( p(0) .GE. " << tostr(FitMax_[iEta*nParsFit_]) << ") THEN ";
    *TF << provMax;
    *TF << "\n        ELSE";
    *TF << "  " << provFormula << "\n        END IF " << endl;

    myTFTable << " \\hline" << endl;
    *TransferCard << " " << endl;  //Need a white line between the different eta-blocks!
        
    if(kinVar == "PT" || kinVar == "E"){
      if(whichDep == 0 || (whichDep == 1 && partName != "muon") ){
        *TF << "\n        tf=(exp(-("+pVar[whichDep]+"-"+pexpVar[whichDep]+"-prov1)**2/2d0/prov2**2))                  !first gaussian\n";
        *TF <<   "        tf=tf+prov3*(exp(-("+pVar[whichDep]+"-"+pexpVar[whichDep]+"-prov4)**2/2d0/prov5**2))         !second gaussian\n";
        *TF <<   "        tf=tf*((1d0/dsqrt(2d0*pi))/(dsqrt(prov2*prov2)+prov3*dsqrt(prov5*prov5)))   !normalisation";
      }
      else{
        *TF << "\n\n        tf=(exp(-(1d0/"+pVar[whichDep]+"-1d0/"+pexpVar[whichDep]+"-prov1)**2/2d0/prov2**2))          !first gaussian\n";
        *TF <<     "        tf=tf+prov3*(exp(-(1d0/"+pVar[whichDep]+"-1d0/"+pexpVar[whichDep]+"-prov4)**2/2d0/prov5**2)) !second gaussian\n";
        *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(dsqrt(prov2*prov2)+prov3*dsqrt(prov5*prov5)))   !normalisation";
      }
    }
    else if(kinVar == "THETA"){
      *TF << "\n\n        tf=(exp(-(theta(p)-theta(pexp)-prov1)**2/2d0/prov2**2))          !first gaussian\n";
      *TF <<     "        tf=tf+prov3*(exp(-(theta(p)-theta(pexp)-prov4)**2/2d0/prov5**2)) !second gaussian\n";
      *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(dsqrt(prov2*prov2)+prov3*dsqrt(prov5*prov5)))   !normalisation";
    }
    else if(kinVar == "PHI"){
      *TF << "\n\n        tf=(exp(-(phi(p)-phi(pexp)-prov1)**2/2d0/prov2**2))           !first gaussian\n";
      *TF <<     "        tf=tf+prov3*(exp(-(phi(p)-phi(pexp)-prov4)**2/2d0/prov5**2))  !second gaussian\n";
      *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(dsqrt(prov2*prov2)+prov3*dsqrt(prov5*prov5)))               !normalisation";
    }
    
    WidthText.push_back("\n \n        width = max(prov2, prov5) ");
    if(nEtaBins_ != 1 && iEta == 4){
      WidthText.push_back("\n      ENDIF \n    </width> \n  </variable>");
      for(unsigned int ii = 0; ii < WidthText.size(); ii++) *TF << WidthText[ii];
    }
    else if(nEtaBins_ == 1){
      WidthText.push_back("\n    </width> \n  </variable>");
      for(unsigned int ii = 0; ii < WidthText.size(); ii++) *TF << WidthText[ii];
    }
  }

  //Will need to decide what to do when eta-bins are considered --> Probably best to make a subsection structure then and show each eta-bins separately!
  if(nEtaBins_ == 1) *LaTeX << " \\caption{Fit distributions} " << endl;
  else               *LaTeX << " \\caption{Fit distributions for " << EtaValues[nEtaBins_] << " $<$ $\\vert \\eta \\vert$ $<$ " << EtaValues[nEtaBins_+1] << "}\n";
  *LaTeX << "\\end{figure} " << endl;

  for(int ii = 0; ii < indivFitHistos.size(); ii++)
    *LaTeX << indivFitHistos[ii];
  *LaTeX << " " << endl;

}

void TFCreation::WritePlots(TFile* outfile){
  outfile->cd();

  TDirectory* th1dir = outfile->mkdir("1D_histograms");
  th1dir->cd();
  for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
  }
  TDirectory* th2dir = outfile->mkdir("2D_histograms_graphs");
  th2dir->cd();
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
    TH2F *temp = it->second;
    temp->Write();
  }
  outfile->cd(); //Step out from 2D_histograms_graphs directory!
}
