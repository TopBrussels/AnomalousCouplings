#include "../interface/TFCreation.h"

TFCreation::TFCreation(int nEtaBins){
  ///////////////////////////////////////
  //  Declare the used fit functions!  //
  ///////////////////////////////////////
  //
  // 1) Double Gaussian --> its range depends on the jet/lepton energy range (hence, the Y-axis)
  doubleGaussianFit = new TF1("doubleGaussianFit","(1/(TMath::Sqrt(2*TMath::Pi())*([1]+[2]*[4])))*(TMath::Exp(-TMath::Power((x-[0]),2)/(2*TMath::Power([1],2)))+[2]*TMath::Exp(-TMath::Power((x-[3]),2)/(2*TMath::Power([4],2))))");
  nParsFit_ = doubleGaussianFit->GetNpar();
  std::string parnames[5]={"a1","a2","a3","a4","a5"};
  std::string ParName[5] = {"Mean broad gaussian", "Width broad gaussian","Relative Constant gaussians","Mean narrow gaussian","Width narrow gaussian"};

  if(nParsFit_ != sizeof(parnames)/sizeof(parnames[0])) std::cout << " ERROR : Difference between number of parameters and defined array --> Also check the header file than !! " << std::endl;

  for(int ipar = 0; ipar < nParsFit_; ipar++){
    parnames_[ipar] = parnames[ipar];
    ParName_[ipar] = ParName[ipar];
    doubleGaussianFit->SetParName(ipar,parnames[ipar].c_str());
  }
 
  //2) Calorimeter Energy formula (ai = ai0 + ai1*Ep + ai2*sqrt(Ep)) --> its range depends on the part energy range (hence, the X-axis)
  caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");
  
  //Store the EtaBin Title and Name for the histograms!
  nEtaBins_ = nEtaBins;
  EtaBin[0] = ""; EtaTitle[0] = "";
  if(nEtaBins_ == 4){         
    EtaValues [1] = 0.; EtaValues[2] = 0.375; EtaValues[3] = 0.750; EtaValues[4] = 1.450; EtaValues[5] = 2.5;
    for(int ii = 1; ii <= nEtaBins_; ii++){
      EtaBin[ii] = "_Eta_"+tostr(EtaValues[ii])+"_"+tostr(EtaValues[ii+1]);
      EtaTitle[ii] = " -- "+tostr(EtaValues[ii])+" < |#eta| #leq "+tostr(EtaValues[ii+1]);
    }
  }
  else if (nEtaBins_ != 1) std::cout << " Wrong choice for NrEtaBins --> Can only be 1 or 4 !! " << std::endl;
}

TFCreation::~TFCreation(){
  delete caloEnergyFit;
  delete doubleGaussianFit;
}

void TFCreation::InitializeVariables(){
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
  histo2D["Light_DiffEVsGenE"]         = new TH2F("Light_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for light quarks",           10,   30,  125, 100,   -30,    35);
  histo2D["Light_RecoPtVsGenPt"]       = new TH2F("Light_RecoPtVsGenPt",       "Transverse momentum of light quarks (reco vs gen)",                150,    0,  300, 150,     0,   300);
  histo2D["Light_DiffPtVsGenPt"]       = new TH2F("Light_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for light quarks",        10,   30,  125, 100,   -30,    35);
  histo2D["Light_RecoThetaVsGenTheta"] = new TH2F("Light_RecoThetaVsGenTheta", "Polar angle distribution of light quarks (reco vs gen)",            60,    0, 3.15,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenTheta"] = new TH2F("Light_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for light quarks", 10,  0.1,  3.1, 100, -0.15,  0.15);
  histo2D["Light_RecoThetaVsGenE"]     = new TH2F("Light_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for light quarks",         120,    0,  300,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenE"]     = new TH2F("Light_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for light quarks",      10,   30,  165, 150, -0.12,  0.12);
  histo2D["Light_RecoThetaVsGenPt"]    = new TH2F("Light_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for light quarks",       120,    0,  300,  60,     0,  3.15);
  histo2D["Light_DiffThetaVsGenPt"]    = new TH2F("Light_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{Tgen} for light quarks",     10,   30,  165, 150, -0.12,  0.12);
  histo2D["Light_RecoPhiVsGenPhi"]     = new TH2F("Light_RecoPhiVsGenPhi",     "Azimuthal angle distribution of light quarks (reco vs gen)",        60, -3.2,  3.2,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenPhi"]     = new TH2F("Light_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for light quarks",     10, -3.2,  3.2, 100,  -0.2,   0.2);
  histo2D["Light_RecoPhiVsGenE"]       = new TH2F("Light_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for light quarks",       150,    0,  300,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenE"]       = new TH2F("Light_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for light quarks",        10,   30,  165, 100, -0.15,  0.15);
  histo2D["Light_RecoPhiVsGenPt"]      = new TH2F("Light_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for light quarks",     150,    0,  300,  60,  -3.2,   3.2);
  histo2D["Light_DiffPhiVsGenPt"]      = new TH2F("Light_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for light quarks",      10,   30,  165, 100, -0.15,  0.15);

  histo2D["BJet_RecoEVsGenE"]         = new TH2F("BJet_RecoEVsGenE",         "Energy of b-jets (reco vs gen level)",                        150,    0,  300, 150,     0,  300);
  histo2D["BJet_DiffEVsGenE"]         = new TH2F("BJet_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for b-jets",            10,   30,  150, 100,   -35,   50);
  histo2D["BJet_RecoPtVsGenPt"]       = new TH2F("BJet_RecoPtVsGenPt",       "Transverse momentum of b-jets (reco vs gen level)",           150,    0,  300, 150,     0,  300);
  histo2D["BJet_DiffPtVsGenPt"]       = new TH2F("BJet_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for b-jets",         10,   30,  150, 100,   -35,   50);
  histo2D["BJet_RecoThetaVsGenTheta"] = new TH2F("BJet_RecoThetaVsGenTheta", "Polar angle distribution of b-jets (reco vs gen)",             60,    0, 3.15,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenTheta"] = new TH2F("BJet_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for b-jets",  10,  0.1,  3.1, 100, -0.15, 0.15);
  histo2D["BJet_RecoThetaVsGenE"]     = new TH2F("BJet_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus energy E_{gen} for b-jets",   120,    0,  300,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenE"]     = new TH2F("BJet_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for b-jets",       10,   30,  160, 150,  -0.1,  0.1);
  histo2D["BJet_RecoThetaVsGenPt"]    = new TH2F("BJet_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus energy P_{T,gen} for b-jets", 120,    0,  300,  60,     0, 3.15);
  histo2D["BJet_DiffThetaVsGenPt"]    = new TH2F("BJet_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for b-jets",     10,   30,  160, 150,  -0.1,  0.1);
  histo2D["BJet_RecoPhiVsGenPhi"]     = new TH2F("BJet_RecoPhiVsGenPhi",     "Azimuthal angle distribution of b-jets (reco vs gen)",         60, -3.2,  3.2,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenPhi"]     = new TH2F("BJet_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for b-jets",      10, -3.2,  3.2, 100,  -0.2,  0.2);
  histo2D["BJet_RecoPhiVsGenE"]       = new TH2F("BJet_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for b-jets",        150,    0,  300,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenE"]       = new TH2F("BJet_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for b-jets",         10,   30,  160, 100,  0.15, 0.15);
  histo2D["BJet_RecoPhiVsGenPt"]      = new TH2F("BJet_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for b-jets",      150,    0,  300,  60,  -3.2,  3.2);
  histo2D["BJet_DiffPhiVsGenPt"]      = new TH2F("BJet_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for b-jets",       10,   30,  160, 100, -0.15, 0.15);

  histo2D["El_RecoEVsGenE"]         = new TH2F("El_RecoEVsGenE",         "Energy of electron (reco vs gen)",                             100,    0,  200, 100,      0,   200);
  histo2D["El_DiffEVsGenE"]         = new TH2F("El_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for electron",           10,   30,  105, 100,     -6,     6);
  histo2D["El_RecoPtVsGenPt"]       = new TH2F("El_RecoPtVsGenPt",        "Transverse momentum of electron (reco vs gen)",                100,    0,  200, 100,      0,   200);
  histo2D["El_DiffPtVsGenPt"]       = new TH2F("El_DiffPtVsGenPt",        "Pt difference (gen-reco) versus P_{T,gen} for electron",        10,   30,  105, 100,     -6,     6);
  histo2D["El_RecoThetaVsGenTheta"] = new TH2F("El_RecoThetaVsGenTheta", "Polar angle distribution of electron (reco vs gen)",            60,    0, 3.15,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenTheta"] = new TH2F("El_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for electron", 10,    0, 3.15, 100,  -0.15,  0.15);
  histo2D["El_RecoThetaVsGenE"]     = new TH2F("El_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for electron",         100,    0,  200,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenE"]     = new TH2F("El_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for electron",      10,   30,  130, 100,  -0.02,  0.02);
  histo2D["El_RecoThetaVsGenPt"]    = new TH2F("El_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for electron",       100,    0,  200,  60,      0,  3.15);
  histo2D["El_DiffThetaVsGenPt"]    = new TH2F("El_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for electron",    10,   30,  130, 100,  -0.02,  0.02);
  histo2D["El_RecoPhiVsGenPhi"]     = new TH2F("El_RecoPhiVsGenPhi",     "Azimuthal angle distribution of electron (reco vs gen)", 	  60, -3.2,  3.2,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenPhi"]     = new TH2F("El_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for electron", 	  10, -3.2,  3.2,  75,  -0.15,  0.15);
  histo2D["El_RecoPhiVsGenE"]       = new TH2F("El_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for electron",       150,    0,  250,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenE"]       = new TH2F("El_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for electron",        10,   30,  130, 120, -0.015, 0.015);
  histo2D["El_RecoPhiVsGenPt"]      = new TH2F("El_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for electron",     150,    0,  250,  60,   -3.2,   3.2);
  histo2D["El_DiffPhiVsGenPt"]      = new TH2F("El_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for electron",      10,   30,  130, 120, -0.015, 0.015);
	
  histo2D["Mu_RecoInvEVsGenInvE"]    = new TH2F("Mu_RecoInvEVsGenInvE",   "Inverse of energy of muon (reco vs gen)",                                  100,     0,  0.05, 100,      0,   0.05);
  histo2D["Mu_DiffInvEVsGenInvE"]    = new TH2F("Mu_DiffInvEVsGenInvE",   "#frac{1}{E} difference (gen-reco) versus #frac{1}{E_{gen}} for muon",       10, 0.005, 0.035, 120,-0.0018, 0.0018);
  histo2D["Mu_RecoInvPtVsGenInvPt"]  = new TH2F("Mu_RecoInvPtVsGenInvPt", "Inverse of transverse momentum of muon (reco vs gen)",                     100,     0,  0.05, 100,      0,   0.05);
  histo2D["Mu_DiffInvPtVsGenInvPt"]  = new TH2F("Mu_DiffInvPtVsGenInvPt", "#frac{1}{P_{T}} difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon", 10, 0.005, 0.035, 120,-0.0018, 0.0018);
  histo2D["Mu_RecoEVsGenE"]          = new TH2F("Mu_RecoEVsGenE",         "Energy of muon (reco vs gen)",                                             150,     0,   200, 150,      0,    200);
  histo2D["Mu_DiffEVsGenE"]          = new TH2F("Mu_DiffEVsGenE",         "E difference (gen-reco) versus E_{gen} for muon",                           10,    26,   150, 100,    -10,     10);
  histo2D["Mu_RecoPtVsGenPt"]        = new TH2F("Mu_RecoPtVsGenPt",       "Transverse momentum of muon (reco vs gen)",                                150,     0,   200, 150,      0,    200);
  histo2D["Mu_DiffPtVsGenPt"]        = new TH2F("Mu_DiffPtVsGenPt",       "Pt difference (gen-reco) versus P_{T,gen} for muon",                        10,    26,   150, 100,    -10,     10);
  histo2D["Mu_RecoThetaVsGenTheta"]  = new TH2F("Mu_RecoThetaVsGenTheta", "Polar angle distribution of muon (reco vs gen)",                            60,     0,  3.15,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenTheta"]  = new TH2F("Mu_DiffThetaVsGenTheta", "#theta difference (gen-reco) versus #theta_{gen} for muon",                 10,     0,  3.15, 100,  -0.15,   0.15);
  histo2D["Mu_RecoThetaVsGenInvE"]   = new TH2F("Mu_RecoThetaVsGenInvE",  "Polar angle #theta_{rec} versus #frac{1}{E_{gen}} for muon",               100,     0,  0.05,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenInvE"]   = new TH2F("Mu_DiffThetaVsGenInvE",  "#theta difference (gen-reco) versus #frac{1}{E_{gen}} for muon",            10, 0.005, 0.035, 150, -0.004,  0.004);
  histo2D["Mu_RecoThetaVsGenInvPt"]  = new TH2F("Mu_RecoThetaVsGenInvPt", "Polar angle #theta_{rec} versus #frac{1}{P_{T,gen}} for muon",             100,     0,  0.05,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenInvPt"]  = new TH2F("Mu_DiffThetaVsGenInvPt", "#theta difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon",          10, 0.005, 0.035, 150, -0.004,  0.004);
  histo2D["Mu_RecoThetaVsGenE"]      = new TH2F("Mu_RecoThetaVsGenE",     "Polar angle #theta_{rec} versus E_{gen} for muon",                         150,     0,   200,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenE"]      = new TH2F("Mu_DiffThetaVsGenE",     "#theta difference (gen-reco) versus E_{gen} for muon",                      10,    30,   150, 100,   -0.1,    0.1);
  histo2D["Mu_RecoThetaVsGenPt"]     = new TH2F("Mu_RecoThetaVsGenPt",    "Polar angle #theta_{rec} versus P_{T,gen} for muon",                       150,     0,   200,  60,      0,   3.15);
  histo2D["Mu_DiffThetaVsGenPt"]     = new TH2F("Mu_DiffThetaVsGenPt",    "#theta difference (gen-reco) versus P_{T,gen} for muon",                    10,    30,   150, 100,   -0.1,    0.1);
  histo2D["Mu_RecoPhiVsGenPhi"]      = new TH2F("Mu_RecoPhiVsGenPhi",     "Azimuthal angle distribution of muon (reco vs gen)", 		       60,  -3.2,   3.2,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenPhi"]      = new TH2F("Mu_DiffPhiVsGenPhi",     "#phi difference (gen-reco) versus #phi_{gen} for muon", 		       10,  -3.2,   3.2,  75,   -0.2,    0.2);
  histo2D["Mu_RecoPhiVsGenInvE"]     = new TH2F("Mu_RecoPhiVsGenInvE",    "Azimuthal angle #phi_{rec} versus #frac{1}{E_{gen}} for muon",             100,     0,  0.05,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenInvE"]     = new TH2F("Mu_DiffPhiVsGenInvE",    "#phi difference (gen-reco) versus #frac{1}{E_{gen}} for muon",              10, 0.005, 0.035, 120, -0.006,  0.006);
  histo2D["Mu_RecoPhiVsGenInvPt"]    = new TH2F("Mu_RecoPhiVsGenInvPt",   "Azimuthal angle #phi_{rec} versus #frac{1}{P_{T,gen}} for muon",           100,     0,  0.05,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenInvPt"]    = new TH2F("Mu_DiffPhiVsGenInvPt",   "#phi difference (gen-reco) versus #frac{1}{P_{T,gen}} for muon",            10, 0.005, 0.035, 120, -0.006,  0.006);
  histo2D["Mu_RecoPhiVsGenE"]        = new TH2F("Mu_RecoPhiVsGenE",       "Azimuthal angle #phi_{rec} versus E_{gen} for muon",                       150,     0,   200,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenE"]        = new TH2F("Mu_DiffPhiVsGenE",       "#phi difference (gen-reco) versus E_{gen} for muon",                        10,    26,   150, 100,   -0.2,    0.2);
  histo2D["Mu_RecoPhiVsGenPt"]       = new TH2F("Mu_RecoPhiVsGenPt",      "Azimuthal angle #phi_{rec} versus P_{T,gen} for muon",                     150,     0,   200,  60,   -3.2,    3.2);
  histo2D["Mu_DiffPhiVsGenPt"]       = new TH2F("Mu_DiffPhiVsGenPt",      "#phi difference (gen-reco) versus P_{T,gen} for muon",                      10,    26,   150, 100,   -0.2,    0.2);

  //Initialize the different eta bins! 
  if(nEtaBins_ == 4){         
    //Store all the existing histo2D's in a new collection!
    map<string,TH2F*> histo2DCopy = histo2D;
    for(int iEta = 1; iEta <= nEtaBins_; iEta++){

      for(std::map<std::string,TH2F*>::const_iterator it = histo2DCopy.begin(); it != histo2DCopy.end(); it++){
	TH2F *temp = it->second;
	
	//Need different Ymax and Ymin for Mu_DiffThetaVsGenInvPt histogram!
	double YMax = temp->GetYaxis()->GetXmax();
	double YMin = temp->GetYaxis()->GetXmin();
	if( string(temp->GetName()) == "Mu_DiffThetaVsGenInvPt"){ YMax = 0.015; YMin = -0.015;}

	if(iEta != nEtaBins_)  //Last eta-bins should have fewer bins to take into account the lower statistics!
	  histo2D[temp->GetName()+EtaBin[iEta]] = new TH2F( (temp->GetName()+EtaBin[iEta]).c_str(), (temp->GetTitle()+EtaTitle[iEta]).c_str(), temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXmin(), temp->GetXaxis()->GetXmax(), (int)(temp->GetYaxis()->GetNbins()*0.75), YMin, YMax);
	else
	  histo2D[temp->GetName()+EtaBin[iEta]] = new TH2F( (temp->GetName()+EtaBin[iEta]).c_str(), (temp->GetTitle()+EtaTitle[iEta]).c_str(), temp->GetXaxis()->GetNbins(), temp->GetXaxis()->GetXmin(), temp->GetXaxis()->GetXmax(), (int)(temp->GetYaxis()->GetNbins()), YMin*1.2, YMax*1.2 );                     
                    
      }
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

  // Select the correct eta-bin which should be used!
  int whichEtaBin = 0; 
  int useEtaBinWJet1 = 99, useEtaBinWJet2 = 99, useEtaBinHadrB = 99, useEtaBinLeptB = 99, useEtaBinLepton = 99;
  if(nEtaBins_ != 1){
    for(int iEta = 1; iEta <= nEtaBins_; iEta++){
      if(abs(selHadrWJet1->Eta()) <= EtaValues[iEta+1] && abs(selHadrWJet1->Eta()) > EtaValues[iEta]) useEtaBinWJet1 = iEta;
      if(abs(selHadrWJet2->Eta()) <= EtaValues[iEta+1] && abs(selHadrWJet2->Eta()) > EtaValues[iEta]) useEtaBinWJet2 = iEta;
      if(abs(selHadrBJet->Eta())  <= EtaValues[iEta+1] && abs(selHadrBJet->Eta()) > EtaValues[iEta])  useEtaBinHadrB = iEta;
      if(abs(selLeptBJet->Eta())  <= EtaValues[iEta+1] && abs(selLeptBJet->Eta()) > EtaValues[iEta])  useEtaBinLeptB = iEta;
      if(abs(selLepton->Eta())    <= EtaValues[iEta+1] && abs(selLepton->Eta()) > EtaValues[iEta])    useEtaBinLepton = iEta;
    }
  }

  //Fill histograms for first light jet in case eta-splitting is desired!
  if(nEtaBins_ != 1) whichEtaBin = useEtaBinWJet1;
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
  if(nEtaBins_ != 1) whichEtaBin = useEtaBinWJet2;
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
  if(nEtaBins_ != 1) whichEtaBin = useEtaBinHadrB;
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
  if(nEtaBins_ != 1) whichEtaBin = useEtaBinLeptB;
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
  if(nEtaBins_ != 1) whichEtaBin = useEtaBinLepton;
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

void TFCreation::CalculateTFFromFile(string fitHistoName, bool useStartValues, int histoNr, bool useROOTClass, bool useStartArray, float startValues[], bool changeFitRange, TFile* file, int whichEtaBin, TFile* readFile){
 
  //Select the correct Eta-bin histogram!
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
  hlist = new TH1D*[nParsFit_+1];
  TObjArray aSlices;
  if(useROOTClass){
    fitHisto->FitSlicesY(doubleGaussianFit, 0, -1, 0, "", &aSlices);
    for(int ipar = 0; ipar <= nParsFit_; ipar++) hlist[ipar] = (TH1D*) aSlices[ipar];
  }
  else
    FitSliceClassCode(fitHisto, changeFitRange);

  //////////////////////////////////////////////////////////////////////////////////////////////
  //   Now histogram with all parameters needs to be fitted with Calorimeter Energy formula   //
  //////////////////////////////////////////////////////////////////////////////////////////////
  for( int ipar = 0; ipar < nParsFit_; ipar++ ){
    if(ipar == 0 || ipar == 2 || ipar == 3 || ipar == 5){
      caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");    //Quartic function as fit!
      for(int ii = 0; ii < 5; ii++) caloEnergyFit->SetParName(ii, ( parnames_[ipar]+tostr(ii)).c_str() );
    }
    else{
      caloEnergyFit = new TF1("caloEnergyFit", "[0]+[1]*sqrt(x)+[2]*x");          //Only expect the calorimeter behavior for the sigma's of the gaussians! 
      for(int ii = 0; ii < 3; ii++) caloEnergyFit->SetParName(ii, ( parnames_[ipar]+tostr(ii)).c_str() );
    }

    double FitMax = fitHisto->GetXaxis()->GetXmax();
    double FitMin = fitHisto->GetXaxis()->GetXmin();
    if( string(fitHisto->GetName()) == "Mu_DiffInvPtVsGenInvPt_Eta_1.45_2.5"){
      for(int ii = 1; ii < 11; ii++){  //Only go until bin 10 since overflow bin is always excluded from fit (and otherwise edge of this bin is taken as max!)
	if( hlist[ipar]->GetBinContent(ii) == 0.) FitMax = hlist[ipar]->GetXaxis()->GetBinLowEdge(ii);
      }
    }
		
    for(int ii = 0; ii < 3; ii++) caloEnergyFit->SetParName(ii, ( parnames_[ipar]+tostr(ii)).c_str() ); //Name here since different for each doubleGaussian parameter!
    caloEnergyFit->SetName( (string(fitHisto->GetName())+"_"+parnames_[ipar]+"_Fit").c_str() );
    hlist[ipar]->SetName( (string(fitHisto->GetName())+"_"+parnames_[ipar]+"_PointsAndFit").c_str() );

    hlist[ipar]->Fit(caloEnergyFit, "Q","",FitMin, FitMax);
    AllCaloEnergyFits[nParsFit_*whichEtaBin+ipar] = *caloEnergyFit;       //caloEnergyFit is a pointer, but each member of the array should point to the corresponding value of the TF1!
    hlist[ipar]->Write();                    
  }
  hlist[nParsFit_]->Write();
  //PlotDlbGaus(fitHisto,file);
							  
  delete [] hlist;
}

void TFCreation::FitSliceClassCode(TH2F* histoFit, bool ChangeFitRange){
  //------------------------------------------------------------------------------------------//
  // Main difference with the Root class FitSlicesY() is the plotting of histograms !        
  // In the Root class the distribution of each hlist histogram is not given!
  // --> Useful to use the own code when needing control histograms
  //
  // Other difference between the two codes have been removed!
  // Originally the treatment of the overflow bin was different, but is now made similar!
  // Create one histogram for each function parameter -> 5 histograms for each 2D plot

  for(int ipar=0 ; ipar < nParsFit_; ipar++){

    float hlistMax = histoFit->GetXaxis()->GetXmax() + ((histoFit->GetXaxis()->GetXmax()-histoFit->GetXaxis()->GetXmin())/histoFit->GetXaxis()->GetNbins());	
    hlist[ipar] = new TH1D( (string(histoFit->GetName())+"_"+parnames_[ipar]).c_str(), (string(histoFit->GetName())+" : Fitted value of "+parnames_[ipar]).c_str(), histoFit->GetXaxis()->GetNbins()+1, histoFit->GetXaxis()->GetXmin(), hlistMax);
    hlist[ipar]->GetXaxis()->SetTitle(histoFit->GetXaxis()->GetTitle());
  }
  hlist[nParsFit_] = new TH1D( (string(histoFit->GetName())+"_chi2").c_str(), (string(histoFit->GetName())+": #chi^{2} distribution for "+string(doubleGaussianFit->GetExpFormula())).c_str(), histoFit->GetXaxis()->GetNbins(), histoFit->GetXaxis()->GetXmin(), histoFit->GetXaxis()->GetXmax() );

  //Loop on all bins in X, generate a projection along Y and fit each bin separately!
  int cut = 0; // require a minimum number of bins in the slice to be filled --> Should this ever be larger than 0 ??
  int nbins = histoFit->GetXaxis()->GetNbins();

  for(int bin=1;bin <= nbins+1;bin ++) {
    string projection_title = string(histoFit->GetName())+"_sliceYbin"+tostr(bin);

    TH1D *hp;
    if(string(histoFit->GetName()).find("Mu_DiffInvPtVsGenInvPt_Eta_1.45") <= string(histoFit->GetName()).size() && bin == 8)
      hp = histoFit->ProjectionY(projection_title.c_str(),8,9,"e");            
    else if(string(histoFit->GetName()).find("Mu_DiffInvPtVsGenInvPt_Eta_1.45") <= string(histoFit->GetName()).size() && bin == 9)    //Temporary fix ... Need to figure out how to enlarge 1 bin!
      hp = histoFit->ProjectionY(projection_title.c_str(),8,9,"e");
    else if(string(histoFit->GetName()).find("El_DiffPtVsGenPt_Eta_1.45") <= string(histoFit->GetName()).size() && bin == 8)
      hp = histoFit->ProjectionY(projection_title.c_str(),8,9,"e");            
    else if(string(histoFit->GetName()).find("El_DiffPtVsGenPt_Eta_1.45") <= string(histoFit->GetName()).size() && bin == 9)
      hp = histoFit->ProjectionY(projection_title.c_str(),8,9,"e");
    else
      hp = histoFit->ProjectionY(projection_title.c_str(),bin,bin,"e");

    //Histogram doesn't have any memory space ...
    if(hp == 0) continue;
    if( float(hp->GetEntries()) <= 0){ delete hp; continue;} //|| float(hp->GetEntries()) < cut) {delete hp; continue;}
    std::string histoName = string(histoFit->GetName());

    doubleGaussianFit->SetName((projection_title+"Fitted").c_str());
    double ActualFitRange[2];
    if(ChangeFitRange == false){                                    ActualFitRange[0] = (double) (histoFit->GetYaxis())->GetXmin(); ActualFitRange[1] = (double) (histoFit->GetYaxis())->GetXmax(); }
    else{ std::vector<double> fitBin = SetFitRange(histoName, bin); ActualFitRange[0] = fitBin[0] ;                                 ActualFitRange[1] = fitBin[1];                                  }

    //Do the actual fit:
    hp->Fit(doubleGaussianFit,"Q","",ActualFitRange[0],ActualFitRange[1]);

    int npfits = doubleGaussianFit->GetNumberFitPoints();              //WHAT IS THIS .... ???
    if(npfits > nParsFit_ && npfits >= cut) {

      //Fill the hlist histogram for each parameter with the obtained Fit parameter and its uncertainty
      //--> Each bin in this histogram represents a bin range in x-axis of considered 2D histogram!
      for(int ipar=0; ipar<nParsFit_; ipar++ ){
        if( !( (histoName == "Light_DiffPtVsGenPt" || histoName.find("Light_DiffPtVsGenPt_Eta_0") <= histoName.size() ) && bin == 2) &&
            !( histoName == "Mu_DiffInvPtVsGenInvPt_Eta_1.45_2.5" && ( bin == 11 || bin == 10) ) &&
            !( histoName == "Mu_DiffInvPtVsGenInvPt" && bin == 1 ) &&                            
            !( histoName.find("Light_DiffPtVsGenPt_Eta_1.45") <= histoName.size() && bin == 7)   &&
            !( histoName.find("Mu_DiffInvPtVsGenInvPt_Eta_0.75") <= histoName.size() && bin == 11) ){

	  hlist[ipar]->Fill(histoFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetParameter(ipar));
	  hlist[ipar]->SetBinError( (int) (bin+1/2) ,doubleGaussianFit->GetParError(ipar)); //WHY +1/2 .... (Is bin size always equal to 1 .. )?
	}
      }
      //Save hchi2 histogram as extra hlist!
      hlist[nParsFit_]->Fill(histoFit->GetXaxis()->GetBinCenter(bin+1/2),doubleGaussianFit->GetChisquare()/(npfits-nParsFit_));
    
//      if( bin == nbins/2 ){
//	if(abs(doubleGaussianFit->GetParameter(1)) < abs(doubleGaussianFit->GetParameter(4))){ NarrowGaus[0] = 0; NarrowGaus[1] = 1; NarrowGaus[2] = 2; WideGaus[0] = 3; WideGaus[1] = 4; WideGaus[2] = 5;}
//	else{  NarrowGaus[0] = 3; NarrowGaus[1] = 4; NarrowGaus[2] = 5; WideGaus[0] = 0; WideGaus[1] = 1; WideGaus[2] = 2;}
//      }
        
    }
    hp->Write();
    delete hp;
  }//loop over bins!
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

void TFCreation::WriteTF(ostream &myTFTable, ostream &myTransferCard, ostream &myTF, std::string kinVar, std::string partName){ 

  std::string pVar[2] = {"p(0)","pt(p)"};
  std::string pexpVar[2] = {"pexp(0)","pt(pexp)"};

  //Is Pt or E dependent considered
  int whichDep = 0;

  string TFDependencyWidth[3]  = {"","*dsqrt("+pVar[whichDep]+")","*"+pVar[whichDep]+")"};
  string TFDependency[5] = {"","*"+pVar[whichDep],"*"+pVar[whichDep]+"**2","*"+pVar[whichDep]+"**3","*"+pVar[whichDep]+"**4)"};
  string WidthDependency[3] = {"","*dsqrt("+pexpVar[whichDep]+")","*"+pexpVar[whichDep]+")"};
  if(partName == "muon"){
    TFDependencyWidth[1] = "*dsqrt(1d0/"+pVar[whichDep]+")"; TFDependencyWidth[2] = "*1d0/"+pVar[whichDep]+")";
    TFDependency[1] = "*1d0/"+pVar[whichDep]; TFDependency[2] = "*1d0/"+pVar[whichDep]+"**2"; TFDependency[3] = "*1d0/"+pVar[whichDep]+"**3"; TFDependency[4] = "*1d0/"+pVar[whichDep]+"**4)";
    WidthDependency[1] = "*dsqrt(1d0/"+pexpVar[whichDep]+")"; WidthDependency[2] = "*1d0/"+pexpVar[whichDep]+")";
  }

  ostream *TransferCard = &myTransferCard;
  ostream *TF = &myTF;
  int dummyCounter = 0;
  vector<std::string> WidthText;
  for(int iEta = 1; iEta <= nEtaBins_; iEta++){

    if(nEtaBins_ != 1){
      //
      //Need this eta separation both for the TF itself and for the width!
      if(iEta == 1)     { WidthText.push_back("\n      ENDIF \n    </tf> \n    <width>"); WidthText.push_back("\n      IF( ABS(eta(pexp)) .LE. 0.375) THEN "); *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 2){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 0.375 .AND. ABS(eta(pexp)) .LE. 0.75) THEN "); *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 3){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 0.75 .AND. ABS(eta(pexp)) .LE. 1.45) THEN ");  *TF << WidthText[WidthText.size()-1];}
      else if(iEta == 4){ WidthText.push_back("\n      ENDIF \n \n      IF( ABS(eta(pexp)) .GT. 1.45 .AND. ABS(eta(pexp)) .LE. 2.5) THEN ");   *TF << WidthText[WidthText.size()-1];}
    }
    else WidthText.push_back("\n    </tf> \n    <width>");
        
    for(int ipar = 0; ipar < nParsFit_; ipar++){
      int NrConsideredCaloPars;
      if(ipar == 0 || ipar == 2 || ipar == 3) NrConsideredCaloPars = 5;
      else NrConsideredCaloPars = 3;

      for(int icalpar = 0; icalpar < NrConsideredCaloPars; icalpar++){
	dummyCounter++;
	if(icalpar == 0) myTFTable<<ParName_[ipar]<<" & $a_{" <<ipar <<icalpar <<"}$ = "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<<"$\\pm$"<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParError(icalpar);
	else              myTFTable<<                 " & $a_{" <<ipar <<icalpar <<"}$ = "<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<<"$\\pm$"<<AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParError(icalpar);

	*TransferCard<< dummyCounter << "     " << AllCaloEnergyFits[iEta*nParsFit_+ipar].GetParameter(icalpar)<< "     # " << ParName_[ipar] << endl;

	if(icalpar==0){
          *TF << "\n        prov"<<ipar+1<<"=(#"<<dummyCounter; 
          if(ipar == 1 || ipar == 4){
            WidthText.push_back("\n        prov"+tostr(ipar+1)+"=(#"+tostr(dummyCounter));
          }
        }
	if(icalpar != 0 && (ipar == 0||ipar == 2||ipar == 3)) *TF << "+#"<<dummyCounter<<TFDependency[icalpar];
	if(icalpar != 0 && (ipar == 1||ipar == 4)){           *TF <<"+#"<<dummyCounter<<TFDependencyWidth[icalpar];
          WidthText.push_back("+#"+tostr(dummyCounter)+WidthDependency[icalpar]);
        }
	if(icalpar==NrConsideredCaloPars-1 && ipar == 2) *TF << "\n        prov"<<ipar+1<<"=max(0,prov"<<ipar+1<<")";
      }
      myTFTable << "\\\\" << endl;
    }

    myTFTable << " \\hline" << endl;
    *TransferCard << " " << endl;  //Need a white line between the different eta-blocks!
        
    if(kinVar == "PT"){
      if(partName != "muon"){
        *TF << "\n\n        tf=(exp(-("+pVar[whichDep]+"-"+pexpVar[whichDep]+"-prov1)**2/2d0/prov2**2))          !first gaussian\n";
        *TF <<     "        tf=tf+prov3*(exp(-("+pVar[whichDep]+"-"+pexpVar[whichDep]+"-prov4)**2/2d0/prov5**2)) !second gaussian\n";
        *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2+prov3*prov5))                                      !normalisation";
      }
      else{
        *TF << "\n\n        tf=(exp(-(1d0/"+pVar[whichDep]+"-1d0/"+pexpVar[whichDep]+"-prov1)**2/2d0/prov2**2))          !first gaussian\n";
        *TF <<     "        tf=tf+prov3*(exp(-(1d0/"+pVar[whichDep]+"-1d0/"+pexpVar[whichDep]+"-prov4)**2/2d0/prov5**2)) !second gaussian\n";
        *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2+prov3*prov5))                                              !normalisation";
      }
    }
    else if(kinVar == "THETA"){
      *TF << "\n\n        tf=(exp(-(theta(p)-theta(pexp)-prov1)**2/2d0/prov2**2))          !first gaussian\n";
      *TF <<     "        tf=tf+prov3*(exp(-(theta(p)-theta(pexp)-prov4)**2/2d0/prov5**2)) !second gaussian\n";
      *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2+prov3*prov5))                  !normalisation";
    }
    else if(kinVar == "PHI"){
      *TF << "\n\n        tf=(exp(-(phi(p)-phi(pexp)-prov1)**2/2d0/prov2**2))           !first gaussian\n";
      *TF <<     "        tf=tf+prov3*(exp(-(phi(p)-phi(pexp)-prov4)**2/2d0/prov5**2))  !second gaussian\n";
      *TF <<     "        tf=tf*((1d0/dsqrt(2d0*pi))/(prov2+prov3*prov5))               !normalisation";
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
}

void TFCreation::PlotDlbGaus(TH2F* fitHisto, TFile* plotsFile){
/*
  const int EPars = 15;
  float EGenValues[EPars] = {10,15,20,30,40,55,70,85,100,115,130,145,160,180,200};

  TCanvas* canvasSame = new TCanvas((string(fitHisto->GetName())+"_StackCanvas_WideAndNarrowGaussian").c_str(),"Stacked canvas for wide and narrow gaussian (only fitted E values)");
  TLegend* sameLegend = new TLegend(0.55,0.7,0.95,0.9);
  canvasSame->Divide(2,3);
  int ptCounter=0;

  if( string(fitHisto->GetName()).find("VsGenInvPt") <= string(fitHisto->GetName()).size()){
    float InvEGenValues[EPars] = {0.1,0.0667, 0.05, 0.033, 0.025, 0.01818, 0.014, 0.0118, 0.01, 0.00869, 0.00769, 0.00689, 0.00625, 0.00556, 0.005};
    for(int ii = 0; ii < EPars; ii++) EGenValues[ii] = InvEGenValues[ii];
  }

  for(int iGenE = 0; iGenE < EPars; iGenE++){
    TH1F* DblGausPlot = new TH1F("DblGausPlot","Double Gaussian distribution using the fit parameters",200,(fitHisto->GetYaxis()->GetXmin())*2,(fitHisto->GetYaxis()->GetXmax())*2);
    DblGausPlot->SetTitle( (string(DblGausPlot->GetTitle())+" (E of parton = "+tostr(EGenValues[iGenE])+")").c_str());
    DblGausPlot->SetName( (string(fitHisto->GetName())+"_DblGausPlot_GenE"+tostr(EGenValues[iGenE])).c_str());
        
    float CaloParGenE[nParsFit_]={0,0,0,0,0,0};
    for(int ipar = 0; ipar < nParsFit_; ipar++){
      if(ipar == 0 || ipar == 2 || ipar == 3){
	for(int icalo = 0; icalo < 5; icalo++) CaloParGenE[ipar] += AllCaloEnergyFits[ipar].GetParameter(icalo)*pow(EGenValues[iGenE],icalo);                
      }
      else{
	for(int icalo = 0; icalo < 3; icalo++) CaloParGenE[ipar] += AllCaloEnergyFits[ipar].GetParameter(icalo)*pow(EGenValues[iGenE],(double) (icalo/2.));
      }
    }
        
    if(CaloParGenE[2] < 0) CaloParGenE[2] = 0;
    if(CaloParGenE[5] < 0) CaloParGenE[5] = 0;    
    float Sqrt2Pi = 2.506628275;

    for(int iBin = 0; iBin <= 200; iBin++)
      DblGausPlot->SetBinContent(iBin,(1/Sqrt2Pi)*(1/(CaloParGenE[1]*CaloParGenE[2] + CaloParGenE[4]*CaloParGenE[5]))*CaloParGenE[2]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[0]),2)/(2*pow(CaloParGenE[1],2))))+CaloParGenE[5]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[3]),2)/(2*pow(CaloParGenE[4],2)))));
    DblGausPlot->Write();

    if(iGenE > 3 && iGenE < 10){
      ptCounter++;    //Counter needed for the divide option of the canvas!

      TH1F* GausNarrow = new TH1F("GausNarrow","Gaussian distribution of the narrow fit",200,(fitHisto->GetYaxis()->GetXmin())*2,(fitHisto->GetYaxis()->GetXmax())*2);
      GausNarrow->SetTitle( (string(GausNarrow->GetTitle())+" (E of parton = "+tostr(EGenValues[iGenE])+")").c_str());
      GausNarrow->SetName( (string(GausNarrow->GetName())+"_NarrowGausPlot_GenE"+tostr(EGenValues[iGenE])).c_str());
      TH1F* GausWide = new TH1F("GausWide","Gaussian distribution of the second fit",200,(fitHisto->GetYaxis()->GetXmin())*2,(fitHisto->GetYaxis()->GetXmax())*2);
      GausWide->SetTitle( (string(GausWide->GetTitle())+" (E of parton = "+tostr(EGenValues[iGenE])+")").c_str());
      GausWide->SetName( (string(GausWide->GetName())+"_WideGausPlot_GenE"+tostr(EGenValues[iGenE])).c_str());
      //TH1F* GausSum = new TH1F("GausSum","Gaussian distribution of the sum",200,(fitHisto->GetYaxis()->GetXmin())*2,(fitHisto->GetYaxis()->GetXmax())*2);
 
      for(int iBin = 0; iBin <= 200; iBin++){
	GausNarrow->SetBinContent(iBin,(1/Sqrt2Pi)*(1/(CaloParGenE[1]*CaloParGenE[2] + CaloParGenE[4]*CaloParGenE[5]))*CaloParGenE[NarrowGaus[2]]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[NarrowGaus[0]]),2)/(2*pow(CaloParGenE[NarrowGaus[1]],2)))));
	GausWide->SetBinContent(iBin,(1/Sqrt2Pi)*(1/(CaloParGenE[1]*CaloParGenE[2] + CaloParGenE[4]*CaloParGenE[5]))*CaloParGenE[WideGaus[2]]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[WideGaus[0]]),2)/(2*pow(CaloParGenE[WideGaus[1]],2)))));
	//GausSum->SetBinContent(iBin,(1/Sqrt2Pi)*(1/(CaloParGenE[1]*CaloParGenE[2] + CaloParGenE[4]*CaloParGenE[5]))*(CaloParGenE[2]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[0]),2)/(2*pow(CaloParGenE[1],2))))+CaloParGenE[5]*(exp(-pow((DblGausPlot->GetXaxis()->GetBinCenter(iBin)-CaloParGenE[3]),2)/(2*pow(CaloParGenE[4],2))))));
      }
      double GausMax = GausNarrow->GetMaximum();
      if(GausWide->GetMaximum() > GausNarrow->GetMaximum()){ GausMax = GausWide->GetMaximum();}// std::cout << "Wrong maximum for histogram " << fitHisto->GetName() << " ! " << std::endl;}

      //gStyle->SetOptStat(0);
      canvasSame->cd(ptCounter);
      GausNarrow->SetLineColor(kRed);
      //GausNarrow->SetMaximum(GausSum->GetMaximum());
      GausNarrow->Draw();
      GausWide->SetLineColor(kGreen);
      GausWide->Draw("same");
      //GausSum->SetLineColor(kOrange);
      //GausSum->Draw("same");

      sameLegend->Clear();
      sameLegend->SetHeader(fitHisto->GetName());
      sameLegend->AddEntry(GausNarrow, ("Narrow gaussian (mean a_{"+tostr(NarrowGaus[0])+"}, sigma a_{"+tostr(NarrowGaus[1])+"} and amplitude a_{"+tostr(NarrowGaus[2])+"})").c_str(), "l");
      sameLegend->AddEntry(GausWide,   ("Wide gaussian   (mean a_{"+tostr(WideGaus[0])+"}, sigma a_{"+tostr(WideGaus[1])+"} and amplitude a_{"+tostr(WideGaus[2])+"})").c_str(), "l");
      //sameLegend->AddEntry(GausSum, "Sum of both gaussians ","l");
      sameLegend->Draw();

      if(iGenE == 9) canvasSame->Write();            
    }
  }
*/
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

std::vector<double> TFCreation::SetFitRange(std::string histoName, int iBin){
    
  std::vector<double> BinnedFitRange;
  BinnedFitRange.clear();
  double FitRangeBinNeg = 0., FitRangeBinPos = 0.;

  if(histoName.find("BJet_DiffPhiVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0") <= histoName.size() ){            //Same fit ranges for first three eta-bins and full distribution!
      double FullFitRangeNeg[11] = {-0.1, -0.12, -0.12, -0.1, -0.1, -0.08, -0.08, -0.08, -0.05, -0.05, -0.05}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.1,  0.12,  0.12,  0.1,  0.1,  0.08,  0.08,  0.08,  0.05,  0.05,  0.05}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size()){
      double FullFitRangeNeg[11] = {-0.15, -0.15, -0.12, -0.12, -0.1, -0.1, -0.08, -0.08, -0.06, -0.06, -0.06}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.15,  0.15,  0.12,  0.12,  0.1,  0.1,  0.08,  0.08,  0.06,  0.06,  0.06}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of BJet_DiffPhi GenE histo

  if(histoName.find("BJet_DiffPhiVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0") <= histoName.size() ){            //Same fit ranges for first three eta-bins and full distribution!
      double FullFitRangeNeg[11] = {-0.12, -0.12, -0.12, -0.1, -0.1, -0.08, -0.08, -0.08, -0.05, -0.05, -0.05}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.12,  0.12,  0.12,  0.1,  0.1,  0.08,  0.08,  0.08,  0.05,  0.05,  0.05}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size()){
      double FullFitRangeNeg[11] = {-0.15, -0.15, -0.12, -0.12, -0.1, -0.1, -0.08, -0.08, -0.06, -0.06, -0.06}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.15,  0.15,  0.12,  0.12,  0.1,  0.1,  0.08,  0.08,  0.06,  0.06,  0.06}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of BJet_DiffPhi GenPt histo

  if(histoName.find("BJet_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-18, -18, -20, -22, -22, -25, -25, -28, -28, -28, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 10,  20,  25,  30,  30,  30,  30,  30,  30,  30,  40}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -18, -20, -22, -22, -25, -25, -28, -28, -28, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];    //Difference for first bin!
      double FullFitRangePos[11] = { 10,  20,  25,  30,  30,  30,  30,  30,  30,  30,  40}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-20, -20, -20, -25, -25, -25, -30, -30, -35, -35, -35}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 10,  20,  25,  25,  30,  20,  20,  20,  20,  20,  20}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of BJet_DiffE histo
  if(histoName.find("BJet_DiffPtVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-18, -18, -20, -22, -22, -25, -25, -28, -28, -28, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 10,  20,  25,  30,  30,  30,  30,  30,  30,  30,  40}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -18, -20, -22, -22, -25, -25, -28, -28, -28, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];    //Difference for first bin!
      double FullFitRangePos[11] = { 10,  20,  25,  30,  30,  30,  30,  30,  30,  30,  40}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-20, -20, -20, -25, -25, -25, -30, -30, -35, -35, -35}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 10,  20,  25,  25,  30,  20,  20,  20,  20,  20,  20}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of BJet_DiffPt GenPt histo

  if(histoName.find("BJet_DiffThetaVsGenE") <= histoName.size() ){ FitRangeBinNeg = -0.1; FitRangeBinPos =  0.1; } //End of BJet_DiffTheta GenE histo

  if(histoName.find("BJet_DiffThetaVsGenPt") <= histoName.size() ){ FitRangeBinNeg = -0.1; FitRangeBinPos =  0.1; } //End of BJet_DiffTheta GenPt histo

  if(histoName.find("El_DiffPhiVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){ FitRangeBinNeg = -0.012; FitRangeBinPos =  0.012; }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.01, -0.01, -0.01}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.01,  0.01,  0.01}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){ FitRangeBinNeg = -0.015; FitRangeBinPos = 0.015; }
  } //End of El_DiffPhi GenE histo

  if(histoName.find("El_DiffPhiVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){ FitRangeBinNeg = -0.012; FitRangeBinPos =  0.012; }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.012, -0.01, -0.01, -0.01}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.012,  0.01,  0.01,  0.01}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){ FitRangeBinNeg = -0.015; FitRangeBinPos = 0.015; }
  } //End of El_DiffPhi GenPt histo

  if(histoName.find("El_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta_1") > histoName.size() ){ FitRangeBinNeg = -4; FitRangeBinPos = 5; }
    else if( histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-6, -6, -6, -6, -6, -6, -6.5, -6.5, -6.5, -6.5, -6.5}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 4,  5,  5,  6,  6,  6,  6.5,  6.5,  6.5,  6.5,  6.5}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of El_DiffE GenE histo

  if(histoName.find("El_DiffPtVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta_1") > histoName.size() ){ FitRangeBinNeg = -4; FitRangeBinPos = 5; }
    else if( histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-6, -6, -6, -6, -6, -6, -6.5, -6.5, -6.5, -6.5, -6.5}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 4,  5,  5,  6,  6,  6,  6.5,  6.5,  6.5,  6.5,  6.5}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of El_DiffE GenPt histo

  if(histoName.find("El_DiffThetaVsGenE") <= histoName.size() ){
    if(histoName.find("Eta_1") > histoName.size() ){ FitRangeBinNeg = -0.018; FitRangeBinPos = 0.018; }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.012, -0.012, -0.012, -0.01, -0.01, -0.007, -0.007, -0.007, -0.005, -0.005, -0.005}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.012,  0.012,  0.012,  0.01,  0.01,  0.007,  0.007,  0.007,  0.005,  0.005,  0.005}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of El_DiffTheta GenE histo

  if(histoName.find("El_DiffThetaVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta_1") > histoName.size() ){ FitRangeBinNeg = -0.018; FitRangeBinPos = 0.018; }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.012, -0.012, -0.012, -0.01, -0.01, -0.007, -0.007, -0.007, -0.005, -0.005, -0.005}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.012,  0.012,  0.012,  0.01,  0.01,  0.007,  0.007,  0.007,  0.005,  0.005,  0.005}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of El_DiffTheta GenPt histo

  if(histoName.find("Light_DiffPhiVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){ 
      double FullFitRangeNeg[11] = {-0.14, -0.14, -0.1, -0.1, -0.1, -0.1, -0.08, -0.08, -0.08, -0.08, -0.08}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.14,  0.14,  0.1,  0.1,  0.1,  0.1,  0.08,  0.08,  0.08,  0.08,  0.08};  FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.12, -0.12, -0.1, -0.1, -0.1, -0.1, -0.08, -0.08, -0.08, -0.08, -0.08}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.12,  0.12,  0.1,  0.1,  0.1,  0.1,  0.08,  0.08,  0.08,  0.08,  0.08}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.15, -0.15, -0.15, -0.12, -0.12, -0.1, -0.08, -0.08, -0.07, -0.07, -0.07}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.15,  0.15,  0.15,  0.12,  0.12,  0.1,  0.08,  0.08,  0.07,  0.07,  0.07}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Light_DiffPhi GenE histo

  if(histoName.find("Light_DiffPhiVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() || histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){ 
      double FullFitRangeNeg[11] = {-0.14, -0.14, -0.1, -0.1, -0.1, -0.1, -0.08, -0.08, -0.08, -0.08, -0.08}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.14,  0.14,  0.1,  0.1,  0.1,  0.1,  0.08,  0.08,  0.08,  0.08,  0.08};  FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.12, -0.12, -0.1, -0.1, -0.1, -0.1, -0.08, -0.08, -0.08, -0.08, -0.08}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.12,  0.12,  0.1,  0.1,  0.1,  0.1,  0.08,  0.08,  0.08,  0.08,  0.08}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.15, -0.15, -0.15, -0.12, -0.12, -0.1, -0.08, -0.08, -0.07, -0.07, -0.07}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.15,  0.15,  0.15,  0.12,  0.12,  0.1,  0.08,  0.08,  0.07,  0.07,  0.07}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Light_DiffPhi GenPt histo

  if(histoName.find("Light_DiffEVsGenE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -18, -20, -20, -22, -22, -25, -28, -28, -28, -28}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  12,  18,  22,  25,  25,  28,  35,  35,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -20, -20, -20, -22, -22, -22, -25, -25, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  14,  18,  20,  24,  24,  25,  30,  30,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-16, -18, -20, -20, -22, -22, -22, -30, -30, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  12,  18,  22,  27,  27,  27,  35,  35,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];      //Check influence of fit (-35 was used in previous code ...)
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -20, -20, -20, -22, -22, -22, -25, -25, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  14,  18,  20,  24,  24,  30,  30,  30,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Light_DiffE GenE histo

  if(histoName.find("Light_DiffPtVsGenPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -18, -20, -20, -22, -22, -25, -28, -28, -28, -28}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  12,  18,  22,  25,  25,  28,  35,  35,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -20, -20, -20, -22, -22, -22, -25, -25, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  14,  18,  20,  24,  24,  25,  30,  30,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-16, -18, -20, -20, -22, -22, -22, -30, -30, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  12,  18,  22,  27,  27,  27,  35,  35,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];      //Check influence of fit (-35 was used in previous code ...)
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-15, -20, -20, -20, -22, -22, -22, -25, -25, -30, -30}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = {  7,  14,  18,  20,  24,  24,  30,  30,  30,  35,  35}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Light_DiffPt GenPt histo

  if(histoName.find("Light_DiffThetaVsGenE") <= histoName.size() ){ FitRangeBinNeg = -0.12; FitRangeBinPos = 0.12; } //End of Light_DiffTheta GenPt histo

  if(histoName.find("Light_DiffThetaVsGenPt") <= histoName.size() ){ FitRangeBinNeg = -0.12; FitRangeBinPos = 0.12; } //End of Light_DiffTheta GenPt histo

  if(histoName.find("Mu_DiffPhiVsGenInvE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0025, -0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0025,  0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.003, -0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.003,  0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.002, -0.0025, -0.003, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.002,  0.0025,  0.003,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.005, -0.005, -0.006}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.005,  0.005,  0.006}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0026, -0.003, -0.0039, -0.0039, -0.0039, -0.0039, -0.0039, -0.0045, -0.0045, -0.0045, -0.0045}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0026,  0.003,  0.0039,  0.0039,  0.0039,  0.0039,  0.0039,  0.0045,  0.0045,  0.0045,  0.0045}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Mu_DiffPhi GenInvE histo

  if(histoName.find("Mu_DiffPhiVsGenInvPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0025, -0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0025,  0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.003, -0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.003,  0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.002, -0.0025, -0.003, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005, -0.005}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.002,  0.0025,  0.003,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005,  0.005}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.003, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.004, -0.005, -0.005, -0.006}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.003,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.004,  0.005,  0.005,  0.006}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0026, -0.003, -0.0039, -0.0039, -0.0039, -0.0039, -0.0039, -0.0045, -0.0045, -0.0045, -0.0045}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0026,  0.003,  0.0039,  0.0039,  0.0039,  0.0039,  0.0039,  0.0045,  0.0045,  0.0045,  0.0045}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Mu_DiffPhi GenInvPt histo

  if(histoName.find("Mu_DiffInvEVsGenInvE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size()){ //0.00045 & 0.00025 working 
      double FullFitRangeNeg[11] = {-0.0005, -0.0012, -0.0012, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0004,  0.001,   0.001,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015 }; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.75") <= histoName.size() ){ FitRangeBinNeg = -0.0018; FitRangeBinPos = 0.0018; }
    else if(histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.001, -0.0012, -0.0012, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.005,  0.008,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015 }; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0015, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017}; FitRangeBinNeg = -0.0022; //FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.001,   0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015}; FitRangeBinPos = 0.0022; //FullFitRangePos[iBin-1];
    }        
  } //End of Mu_DiffInvE GenInvE histo

  if(histoName.find("Mu_DiffInvPtVsGenInvPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size()){ //0.00045 & 0.00025 working 
      double FullFitRangeNeg[11] = {-0.0005, -0.0012, -0.0012, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.0004,  0.001,   0.001,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015 }; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.75") <= histoName.size() ){ FitRangeBinNeg = -0.0018; FitRangeBinPos = 0.0018; }
    else if(histoName.find("Eta_0.375") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.001, -0.0012, -0.0012, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015, -0.0015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.005,  0.008,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015,   0.0015 }; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.0015, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017, -0.0017}; FitRangeBinNeg = -0.0022; //FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.001,   0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015,  0.0015}; FitRangeBinPos = 0.0022; //FullFitRangePos[iBin-1];
    }        
  } //End of Mu_DiffInvPt GenInvPt histo

  if(histoName.find("Mu_DiffThetaVsGenInvE") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-0.002, -0.0025, -0.0025, -0.003, -0.003, -0.003, -0.0035, -0.0035, -0.0035, -0.0035, -0.0035}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.002,  0.0025,  0.0025,  0.003,  0.003,  0.003,  0.0035,  0.0035,  0.0035,  0.0035,  0.0035}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() || histoName.find("1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.008, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.008,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.008, -0.01, -0.01, -0.01, -0.012, -0.012, -0.012, -0.015, -0.015, -0.015, -0.015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.008,  0.01,  0.01,  0.01,  0.012,  0.012,  0.012,  0.015,  0.015,  0.015,  0.015}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Mu_DiffTheta GenInvE histo

  if(histoName.find("Mu_DiffThetaVsGenInvPt") <= histoName.size() ){
    if(histoName.find("Eta") > histoName.size() ){
      double FullFitRangeNeg[11] = {-0.002, -0.0025, -0.0025, -0.003, -0.003, -0.003, -0.0035, -0.0035, -0.0035, -0.0035, -0.0035}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.002,  0.0025,  0.0025,  0.003,  0.003,  0.003,  0.0035,  0.0035,  0.0035,  0.0035,  0.0035}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0_") <= histoName.size() || histoName.find("Eta_0.375") <= histoName.size() || histoName.find("1.45") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.008, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01, -0.01}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.008,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
    else if(histoName.find("Eta_0.75") <= histoName.size() ){
      double FullFitRangeNeg[11] = {-0.008, -0.01, -0.01, -0.01, -0.012, -0.012, -0.012, -0.015, -0.015, -0.015, -0.015}; FitRangeBinNeg = FullFitRangeNeg[iBin-1];
      double FullFitRangePos[11] = { 0.008,  0.01,  0.01,  0.01,  0.012,  0.012,  0.012,  0.015,  0.015,  0.015,  0.015}; FitRangeBinPos = FullFitRangePos[iBin-1];
    }
  } //End of Mu_DiffTheta GenInvPt histo

  BinnedFitRange.push_back(FitRangeBinNeg); BinnedFitRange.push_back(FitRangeBinPos);

  return  BinnedFitRange;
}

void TFCreation::CalculateTF(bool drawHistos, bool doFits, bool useROOTClass, bool useStartValues){
  TFile* file = new TFile("TFInformation/PlotsForTransferFunctions.root","RECREATE");
  file->cd();

  if(drawHistos == true) WritePlots(file);
  if(doFits == true){

    ///////////////////////////////////////////
    //  Choose the correct histogram to fit  //
    ///////////////////////////////////////////
    TH2F* histoForFit = new TH2F();
    for (unsigned int f=0; f<24;f++) {				
				
      switch(f){
      case 0:
	histoForFit=histo2D["BJet_DiffPhiVsGenE"]; break;
      case 1:
	histoForFit=histo2D["BJet_DiffEVsGenE"]; break;
      case 2:
	histoForFit=histo2D["BJet_DiffThetaVsGenE"]; break;
      case 3:
	histoForFit=histo2D["El_DiffPhiVsGenE"]; break;
      case 4:
	histoForFit=histo2D["El_DiffEVsGenE"]; break;
      case 5:
	histoForFit=histo2D["El_DiffThetaVsGenE"]; break;
      case 6:
	histoForFit=histo2D["Light_DiffPhiVsGenE"]; break;
      case 7:
	histoForFit=histo2D["Light_DiffEVsGenE"]; break;
      case 8:
	histoForFit=histo2D["Light_DiffThetaVsGenE"]; break;
      case 9:
	histoForFit=histo2D["Mu_DiffPhiVsGenInvE"]; break;
      case 10:
	histoForFit=histo2D["Mu_DiffInvEVsGenInvE"]; break;
      case 11:
	histoForFit=histo2D["Mu_DiffThetaVsGenInvE"]; break;
      case 12:
	histoForFit=histo2D["BJet_DiffPhiVsGenPt"]; break;
      case 13:
	histoForFit=histo2D["BJet_DiffPtVsGenPt"]; break;
      case 14:
	histoForFit=histo2D["BJet_DiffThetaVsGenPt"]; break;
      case 15:
	histoForFit=histo2D["El_DiffPhiVsGenPt"]; break;
      case 16:
	histoForFit=histo2D["El_DiffPtVsGenPt"]; break;
      case 17:
	histoForFit=histo2D["El_DiffThetaVsGenPt"]; break;
      case 18:
	histoForFit=histo2D["Light_DiffPhiVsGenPt"]; break;
      case 19:
	histoForFit=histo2D["Light_DiffPtVsGenPt"]; break;
      case 20:
	histoForFit=histo2D["Light_DiffThetaVsGenPt"]; break;
      case 21:
	histoForFit=histo2D["Mu_DiffPhiVsGenInvPt"]; break;
      case 22:
	histoForFit=histo2D["Mu_DiffInvPtVsGenInvPt"]; break;
      case 23:
	histoForFit=histo2D["Mu_DiffThetaVsGenInvPt"]; break;
      }					

      TDirectory* histoFitDir = file->mkdir(histoForFit->GetName());
      histoFitDir->cd();
 
      hlist = new TH1D*[nParsFit_+1];
      if(useStartValues)
	SetStartValuesDoubleGaussian(f, false, string(histoForFit->GetName()));   //false means that normal start values are being used!

      TObjArray aSlices;
      if(useROOTClass){
	//Fit using the FitSliceY function of TF1!
	histoForFit->FitSlicesY(doubleGaussianFit, 0, -1, 0, "", &aSlices);
	for(int ipar = 0; ipar <= nParsFit_; ipar++)
	  hlist[ipar] = (TH1D*) aSlices[ipar];
      }
      else
	FitSliceClassCode(histoForFit,false);

      //////////////////////////////////////////////////////////////////////////////////////////////
      //   Now histogram with all parameters needs to be fitted with Calorimeter Energy formula   //
      //////////////////////////////////////////////////////////////////////////////////////////////
      caloEnergyFit->SetRange( histoForFit->GetXaxis()->GetXmin(), histoForFit->GetXaxis()->GetXmax() );
      for( int ipar = 0; ipar < nParsFit_; ipar++ ){

	//give names to the parameters		
	for(int jj = 0; jj < 3; jj++) caloEnergyFit->SetParName(jj, ( parnames_[ipar]+tostr(jj)).c_str() );
	caloEnergyFit->SetName( (string(histoForFit->GetName())+"_"+parnames_[ipar]+"_Fit").c_str() );
	hlist[ipar]->SetName( (string(histoForFit->GetName())+"_"+parnames_[ipar]+"_PointsAndFit").c_str() );

	hlist[ipar]->Fit(caloEnergyFit);
	hlist[ipar]->Write();            
      }
      hlist[nParsFit_]->Write();
						
    }//Loop over f						
    delete histoForFit;
    delete [] hlist;
  }                               //Boolean doFits = true
  file->Close();
  delete file;
}
