#include <iomanip>
#include "../interface/LHCOOutput.h"

#include <fstream>
#include <sstream>

LHCOOutput::LHCOOutput(int verbosity, bool writeOutput, bool splitLeptCharge, bool splitCorrectWrong){
  verbose_ = verbosity;
  writeOutput_ = writeOutput;
  splitLeptCharge_ = splitLeptCharge;     //Specifiy whether the output should be splitted into mu+ and mu- or just kept as muon!
  splitCorrectWrong_ = splitCorrectWrong;
  GenOrReco_ = "";
  CorrectGenEvtContent = false;
  leptonType = notFound;
  if(verbose_ > 3) std::cout << " In constructor of LHCOOutput! " << std::endl;
}


LHCOOutput::~LHCOOutput(){

  //Close the output files:
  for(int ii = 0; ii < 4; ii++){
    if(GenOrReco_ == "Gen" && writeOutput_) GenOutFile[ii].close();
    if(GenOrReco_ == "Reco" && writeOutput_) RecoOutFile[ii].close();
  }
  if(GenOrReco_ == "Reco" && writeOutput_){ WrongRecoMuPosFile.close(); CorrectRecoMuPosFile.close();UnmatchedRecoMuPosFile.close();}

}

void LHCOOutput::Initialize(string GenOrReco){
  if(verbose_ > 3) std::cout << " In Initialize function of LHCOOutput ... " << std::endl;    

  GenOrReco_ = GenOrReco; //Will be updated if the sample is TTbarJets and then the destructor will actually do something (for all other samples nothing will be executed in destructor)
  NumberNegativeElectrons = 0; NumberNegativeMuons = 0; NumberPositiveElectrons = 0; NumberPositiveMuons = 0;
  WrongEvtCounter = 0;
  NumberNegRecoEl = 0;    NumberNegRecoMu = 0;  NumberPosRecoEl = 0;  NumberPosRecoMu = 0; 
  NrPosRecoMuCorrect = 0; NrPosRecoMuWrong = 0; NrPosRecoMuUnmatched = 0; 

  if(GenOrReco == "Gen" && writeOutput_ == true){
    if(splitLeptCharge_){
      GenOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveMuon.lhco");
      GenOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeMuon.lhco");
      GenOutFile[2].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_PositiveElectron.lhco");
      GenOutFile[3].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_NegativeElectron.lhco");
    }
    else{
      GenOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_Muon.lhco");
      GenOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarLHCO_Electron.lhco");
    }
  }
  else if(GenOrReco == "Reco" && writeOutput_ == true){
    if(splitLeptCharge_){
      RecoOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveMuon.lhco");
      RecoOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeMuon.lhco");
      RecoOutFile[2].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_PositiveElec.lhco");
      RecoOutFile[3].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_NegativeElec.lhco");
    }
    else{
      RecoOutFile[0].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_Muon.lhco");
      RecoOutFile[1].open("MadWeightInput/AnalyzerOutput/TTbarSemiLepton_Reco_Electron.lhco");
    }

    if(splitCorrectWrong_){
      WrongRecoMuPosFile.open("MadWeightInput/AnalyzerOutput/WrongRecoEvents_PositiveMuon.lhco");
      UnmatchedRecoMuPosFile.open("MadWeightInput/AnalyzerOutput/UnmatchedRecoEvents_PositiveMuon.lhco");
      CorrectRecoMuPosFile.open("MadWeightInput/AnalyzerOutput/CorrectRecoEvents_PositiveMuon.lhco");
    }
  }

  histo1D[(GenOrReco+"_TopMassLept").c_str()] = new TH1F((GenOrReco+"_TopMassLept").c_str(),("Leptonic top-mass distribution for "+GenOrReco+" events").c_str(),250,130,210);
  histo1D[(GenOrReco+"_TopMassHadr").c_str()] = new TH1F((GenOrReco+"_TopMassHadr").c_str(),("Hadronic top-mass distribution for "+GenOrReco+" events").c_str(),250,130,210);

  histo1D[("Pt_"+GenOrReco+"_LightJets").c_str()] = new TH1F(("Pt_"+GenOrReco+"_LightJets").c_str(),("Pt-distribution for light jets ("+GenOrReco+")").c_str(),500,0,250);
  histo1D[("Pt_"+GenOrReco+"_BJets").c_str()]     = new TH1F(("Pt_"+GenOrReco+"_BJets").c_str(),    ("Pt-distribution for b-jets ("+GenOrReco+")").c_str(),500,0,250);
  histo1D[("Pt_"+GenOrReco+"_Muon").c_str()]      = new TH1F(("Pt_"+GenOrReco+"_Muon").c_str(),     ("Pt-distribution for muons ("+GenOrReco+")").c_str(),500,0,250);
  histo1D[("Pt_"+GenOrReco+"_Electron").c_str()]  = new TH1F(("Pt_"+GenOrReco+"_Electron").c_str(), ("Pt-distribution for electrons ("+GenOrReco+")").c_str(),500,0,250);
  histo1D[("Pt_"+GenOrReco+"_MET").c_str()]       = new TH1F(("Pt_"+GenOrReco+"_MET").c_str(),      ("Pt-distribution for neutrinos ("+GenOrReco+")").c_str(),500,0,250);

  histo1D[("Mass_"+GenOrReco+"_BJets").c_str()]    = new TH1F(("Mass_"+GenOrReco+"_BJets").c_str(),   ("Mass distribution for b-jets ("+GenOrReco+")").c_str(),150,-0.5,5.5);
  histo1D[("Mass_"+GenOrReco+"_Muon").c_str()]     = new TH1F(("Mass_"+GenOrReco+"_Muon").c_str(),    ("Mass distribution for muons ("+GenOrReco+")").c_str(),150,-0.1,0.2);
  histo1D[("Mass_"+GenOrReco+"_Electron").c_str()] = new TH1F(("Mass_"+GenOrReco+"_Electron").c_str(),("Mass distribution for electrons ("+GenOrReco+")").c_str(),150,-0.0001,0.001);

  histo1D[("Eta_"+GenOrReco+"_LightJets").c_str()] = new TH1F(("Eta_"+GenOrReco+"_LightJets").c_str(),("Eta-distribution for light jets ("+GenOrReco+")").c_str(),250,-5,5);
  histo1D[("Eta_"+GenOrReco+"_BJets").c_str()]     = new TH1F(("Eta_"+GenOrReco+"_BJets").c_str(),    ("Eta-distribution for b-jets ("+GenOrReco+")").c_str(),250,-5,5);
  histo1D[("Eta_"+GenOrReco+"_Muon").c_str()]      = new TH1F(("Eta_"+GenOrReco+"_Muon").c_str(),     ("Eta-distribution for muons ("+GenOrReco+")").c_str(),250,-5,5);
  histo1D[("Eta_"+GenOrReco+"_Electron").c_str()]  = new TH1F(("Eta_"+GenOrReco+"_Electron").c_str(), ("Eta-distribution for electrons ("+GenOrReco+")").c_str(),250,-5,5);
  histo1D[("Eta_"+GenOrReco+"_MET").c_str()]       = new TH1F(("Eta_"+GenOrReco+"_MET").c_str(),      ("Eta-distribution for neutrinos ("+GenOrReco+")").c_str(),250,-5,5);

  histo1D[("deltaR_"+GenOrReco+"_LightJets").c_str()]       = new TH1F(("deltaR_"+GenOrReco+"_LightJets").c_str(),      ("deltaR distribution between both light jets ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_BJets").c_str()]           = new TH1F(("deltaR_"+GenOrReco+"_BJets").c_str(),          ("deltaR distribution between both b-jets ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_LightWithB").c_str()]      = new TH1F(("deltaR_"+GenOrReco+"_LightWithB").c_str(),     ("deltaR distribution between light and b-jet ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_MuonWithJet").c_str()]     = new TH1F(("deltaR_"+GenOrReco+"_MuonWithJet").c_str(),    ("deltaR distribution between muon and light jet ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_MuonWithB").c_str()]       = new TH1F(("deltaR_"+GenOrReco+"_MuonWithB").c_str(),      ("deltaR distribution between muon and b-jet ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_ElectronWithJet").c_str()] = new TH1F(("deltaR_"+GenOrReco+"_ElectronWithJet").c_str(),("deltaR distribution between electron and jet ("+GenOrReco+")").c_str(),250,0,10);
  histo1D[("deltaR_"+GenOrReco+"_ElectronWithB").c_str()]   = new TH1F(("deltaR_"+GenOrReco+"_ElectronWithB").c_str(),  ("deltaR distribution between electron and b-jet ("+GenOrReco+")").c_str(),250,0,10);
}

void LHCOOutput::StoreGenInfo(vector<TRootMCParticle*> mcParticles){

  //Initialize EventContent counter used for selecting correct particle content!
  int EventContent[5]; //0:top; 1:b; 2: u,c,d,s; 3:W; 4:mu + neutrino
  for(int ll = 0;ll<5;ll++){EventContent[ll]=0;}
 
  //Loop over all the mcParticles
  for(unsigned int i=0; i<mcParticles.size(); i++){
    if( mcParticles[i]->status() != 3) continue;
  
    int partType=mcParticles[i]->type(); if(verbose_>4)cout<<"-->Type of mcParticle : "<<partType<<endl;
    //if( fabs(partType) == 4) std::cout << " Mass of c-quark : " << mcParticles[i]->M() << " with mother : " << mcParticles[i]->motherType() << std::endl; --> Almost lways 1.5
	
    if(fabs(partType)<7 || fabs(partType)==24 || (fabs(partType)<=14 && fabs(partType)>=11) ){ //Considering only the semileptonic particles
      int motherType=mcParticles[i]->motherType(); 
      int grannyType=mcParticles[i]->grannyType();
      if(verbose_ > 5)cout<<"Mother type of particle : "<<motherType<<", and granny type : "<<grannyType<<endl;
	  
      if(partType == 6){      Top    = (TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbose_>4) cout<<"*Particle found: Top"<<endl;    }      
      else if(partType == -6){TopBar = (TRootMCParticle*) mcParticles[i]; EventContent[0]++; if(verbose_>4) cout<<"*Particle found: AntiTop"<<endl;}	  	  

      else if(fabs(partType) == 5 && fabs(motherType) == 6){
	EventContent[1]++;
	if(partType == 5){      Bottom =    (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: Bottom"<<endl;}
	else if(partType == -5){BottomBar = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: AntiBottom"<<endl;}
      }//End of bottom particle identification  
	  
      else if(fabs(partType) == 24 && fabs(motherType) == 6){//Check correct definition!!!
	EventContent[3]++;
	if(partType == 24){      WPlus =  (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: WPlus"<<endl;}        
	else if(partType == -24){WMinus = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: WMinus"<<endl;}
      }//End of WBoson identification
	  
      else if(fabs(partType) <=4 && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[2]++;
	if(partType > 0){     Light =    (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: Light"<<endl;}
	else if(partType < 0){LightBar = (TRootMCParticle*) mcParticles[i]; if(verbose_>4) cout<<"*Particle found: AntiLight"<<endl;}
      }//End of light particle identification

      else if((fabs(partType) == 13 || fabs(partType) == 11 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[4]++;
	string LeptonType="";
	if(fabs(partType) == 13){      if(verbose_>4) LeptonType = "*Particle found: Muon";}
	else if(fabs(partType) == 11){ if(verbose_>4) LeptonType = "*Particle found: Electron";}
	Lepton = (TRootMCParticle*) mcParticles[i]; if(verbose_ > 4) cout<<LeptonType<<endl;
      }//End of lepton identification
	  
      else if((fabs(partType) == 14 || fabs(partType) == 12 ) && fabs(motherType) == 24 && fabs(grannyType) == 6){
	EventContent[4]++;
	string neutrinoType="";
	if(fabs(partType) == 14){      if(verbose_>4) neutrinoType = "*Particle found: Muon-neutrino";}
	else if(fabs(partType) == 12){ if(verbose_>4) neutrinoType = "*Particle found: Electron-neutrino";}
	NeutrinoMC = (TRootMCParticle*) mcParticles[i]; if(verbose_ > 4) cout<<neutrinoType<<endl;
      }//End of neutrino identification		
	  
    }//End of looking at semi-leptonic particles inside event ==> Semileptonic event is completely created now!	
  }//End of loop over mcParticles inside one particular event

  //std::cout << " Event content is : " << std::endl;
  //std::cout << "  * Number of top quarks is       : " << EventContent[0] << std::endl;
  //std::cout << "  * Number of b-quarks is         : " << EventContent[1] << std::endl;
  //std::cout << "  * Number of light quarks is     : " << EventContent[2] << std::endl;
  //std::cout << "  * Number of W-bosons is         : " << EventContent[3] << std::endl;
  //std::cout << "  * Number of lepton/neutrinos is : " << EventContent[4] << std::endl;

  //---  Consider only events with correct event content (b b q q l vl)  ---//
  if(EventContent[0]==2 && EventContent[1]==2 && EventContent[2]==2 && EventContent[3]==2 && EventContent[4]==2){
    CorrectGenEvtContent = true;
    vector<TLorentzVector*> LHCOVector(6);
    vector<int> MadGraphId(6,4);
    vector<float> MGBtag(6,0.0);
    MGBtag[0] = 2.0; MGBtag[3] = 2.0;     //b-jets always on position 0 and 3!
      
    if(verbose_>3){
      cout << " Event with correct event content found "       << endl;
      cout << " Mass of bottom quark    : " << Bottom->M()     << endl;
      cout << " Mass of light quark     : " << Light->M()      << endl;
      cout << " Mass of LightBar quark  : " << LightBar->M()   << endl;
      cout << " Mass of BottomBar quark : " << BottomBar->M()  << endl;
      cout << " Mass of lepton          : " << Lepton->M()     << endl;
      cout << " Mass of neutrino        : " << NeutrinoMC->M() << endl;
    }

    //Save the lepton and neutrino as TLorentzVectors
    GenLepton = (TLorentzVector*) Lepton;
    GenNeutrino = (TLorentzVector*) NeutrinoMC;

    //--- Plot the generator-level distributions for Pt, eta and dR ---//
    //---   --> Will be used for MG cuts comparison                 ---//
    histo1D["Pt_Gen_LightJets"]->Fill(Light->Pt()); histo1D["Pt_Gen_LightJets"]->Fill(LightBar->Pt());
    histo1D["Pt_Gen_BJets"]->Fill(Bottom->Pt());    histo1D["Pt_Gen_BJets"]->Fill(BottomBar->Pt()); 
    histo1D["Mass_Gen_BJets"]->Fill(Bottom->M());   histo1D["Mass_Gen_BJets"]->Fill(BottomBar->M());
    if(fabs(Lepton->type()) == 13){ histo1D["Pt_Gen_Muon"]->Fill(Lepton->Pt());     histo1D["Mass_Gen_Muon"]->Fill(Lepton->M());    } 
    if(fabs(Lepton->type()) == 11){ histo1D["Pt_Gen_Electron"]->Fill(Lepton->Pt()); histo1D["Mass_Gen_Electron"]->Fill(Lepton->M());}
    histo1D["Pt_Gen_MET"]->Fill(NeutrinoMC->Pt());

    histo1D["Eta_Gen_LightJets"]->Fill(Light->Eta()); histo1D["Eta_Gen_LightJets"]->Fill(LightBar->Eta());
    histo1D["Eta_Gen_BJets"]->Fill(Bottom->Eta());    histo1D["Eta_Gen_BJets"]->Fill(BottomBar->Eta());
    if(fabs(Lepton->type()) == 13) histo1D["Eta_Gen_Muon"]->Fill(Lepton->Eta());
    if(fabs(Lepton->type()) == 11) histo1D["Eta_Gen_Electron"]->Fill(Lepton->Eta());
    histo1D["Eta_Gen_MET"]->Fill(NeutrinoMC->Eta());

    histo1D["deltaR_Gen_LightJets"]->Fill(Light->DeltaR(*LightBar));
    histo1D["deltaR_Gen_BJets"]->Fill(Bottom->DeltaR(*BottomBar));
    histo1D["deltaR_Gen_LightWithB"]->Fill(Light->DeltaR(*Bottom)); histo1D["deltaR_Gen_LightWithB"]->Fill(Light->DeltaR(*BottomBar)); 
    histo1D["deltaR_Gen_LightWithB"]->Fill(LightBar->DeltaR(*Bottom)); histo1D["deltaR_Gen_LightWithB"]->Fill(LightBar->DeltaR(*BottomBar));
    if(fabs(Lepton->type()) == 13){
      histo1D["deltaR_Gen_MuonWithJet"]->Fill(Lepton->DeltaR(*Light)); histo1D["deltaR_Gen_MuonWithJet"]->Fill(Lepton->DeltaR(*LightBar));
      histo1D["deltaR_Gen_MuonWithB"]->Fill(Lepton->DeltaR(*Bottom)); histo1D["deltaR_Gen_MuonWithB"]->Fill(Lepton->DeltaR(*BottomBar));
    }
    if(fabs(Lepton->type()) == 11){
      histo1D["deltaR_Gen_ElectronWithJet"]->Fill(Lepton->DeltaR(*Light)); histo1D["deltaR_Gen_ElectronWithJet"]->Fill(Lepton->DeltaR(*LightBar));
      histo1D["deltaR_Gen_ElectronWithB"]->Fill(Lepton->DeltaR(*Bottom)); histo1D["deltaR_Gen_ElectronWithB"]->Fill(Lepton->DeltaR(*BottomBar));
    }
	
    //Create the lhco file for pp > t t~:
    if(Lepton->type() == 13 || Lepton->type() == 11){ //Negative lepton, hence t~ > b~ W-, W- > e/mu- ve/vm
      LHCOVector[0] = (TLorentzVector*) Bottom;
      LHCOVector[1] = (TLorentzVector*) Light;
      LHCOVector[2] = (TLorentzVector*) LightBar;
      LHCOVector[3] = (TLorentzVector*) BottomBar;
      LHCOVector[4] = (TLorentzVector*) Lepton;
      LHCOVector[5] = (TLorentzVector*) NeutrinoMC;
      MadGraphId[5] = 6;                  //MadGraph Id of MET = 6
      if(Lepton->type() == 11){           //Looking at negative electron events (index 3 for LHCO file)
  	MadGraphId[4] = 1;                //MadGraph Id of e = 1
	if(writeOutput_){
	  NumberNegativeElectrons++;
          leptonType = elMinus;      //Enum information	  --> Used anywhere??
	  if(splitLeptCharge_) LHCOEventOutput(3, GenOutFile[3], NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
          else                 LHCOEventOutput(3, GenOutFile[1], NumberNegativeElectrons+NumberPositiveElectrons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Negative electron
      else if(Lepton->type() == 13){       //Looking at negative muon events (index 1 for LHCO file)
	MadGraphId[4] = 2; //MadGraphId of mu = 2
	if(writeOutput_){
	  NumberNegativeMuons++;
          leptonType = muMinus;      //Enum information
	  if(splitLeptCharge_) LHCOEventOutput(1, GenOutFile[1], NumberNegativeMuons,LHCOVector,MadGraphId, MGBtag);
          else                 LHCOEventOutput(1, GenOutFile[0], NumberNegativeMuons+NumberPositiveMuons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Negative muon

      //--- Store information needed for cos theta* ---//
      GenLeptonicW = (TLorentzVector*) WMinus;
      GenLeptonicTop = (TLorentzVector*) TopBar;
      //GenHadronicW = (TLorentzVector*) WPlus;
      //GenHadronicTop = (TLorentzVector*) Top;

    }//Negative lepton
    else if(Lepton->type() == -13 || Lepton->type() == -11){ //Positive lepton, hence t > b W+, W+ > e/mu+ ve/vm
      LHCOVector[0] = (TLorentzVector*) Bottom; 
      LHCOVector[1] = (TLorentzVector*) Lepton;
      LHCOVector[2] = (TLorentzVector*) NeutrinoMC;
      MadGraphId[2] = 6;          //MET always on position 2
      LHCOVector[3] = (TLorentzVector*) BottomBar;
      LHCOVector[4] = (TLorentzVector*) Light;
      LHCOVector[5] = (TLorentzVector*) LightBar;
      if(Lepton->type() == -11){            //Looking at positive electron events (index 2 for LHCO file)
	MadGraphId[1] = 1;                  //MadGraphId of electron = 1
	if(writeOutput_){
	  NumberPositiveElectrons++;
          leptonType = elPlus;      //Enum information
	  if(splitLeptCharge_) LHCOEventOutput(2, GenOutFile[2], NumberPositiveElectrons,LHCOVector,MadGraphId, MGBtag);
          else                 LHCOEventOutput(2, GenOutFile[1], NumberPositiveElectrons+NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Positive electron
      else if(Lepton->type() == -13){             //Looking at positive muon events (index 0 for LHCO file)
	MadGraphId[1] = 2;                        //MadGraphId of muon = 2
	if(writeOutput_){
	  NumberPositiveMuons++;
          leptonType = muPlus;      //Enum information
	  if(splitLeptCharge_) LHCOEventOutput(0, GenOutFile[0], NumberPositiveMuons,LHCOVector,MadGraphId, MGBtag);
          else                 LHCOEventOutput(0, GenOutFile[0], NumberPositiveMuons+NumberNegativeElectrons,LHCOVector,MadGraphId, MGBtag);
	}
      }//Positive muon
	  
      //--- Store information needed for cos theta* ---//
      GenLeptonicW = (TLorentzVector*) WPlus;
      GenLeptonicTop = (TLorentzVector*) Top;
      //GenHadronicW = (TLorentzVector*) WMinus;
      //GenHadronicTop = (TLorentzVector*) TopBar;

    }//Positive lepton
  }//Correct event content found
  else{
    CorrectGenEvtContent = false;
    WrongEvtCounter++;

    //Output in case wrong event content needs to be double-checked
    if(verbose_>4){
      cout << " Number of top quarks      : " << EventContent[0] << endl;
      cout << " Number of bottom quarks   : " << EventContent[1] << endl;
      cout << " Number of light quarks    : " << EventContent[2] << endl;
      cout << " Number of W-bosons        : " << EventContent[3] << endl;
      cout << " Number of lepton/neutrino : " << EventContent[4] << endl;
    }
  }			    
}//End of class StoreGenInfo

void LHCOOutput::StoreRecoInfo(TLorentzVector lepton, vector<TLorentzVector> Jets, int bLeptIndex, int bHadrIndex, int light1Index, int light2Index, int decayChannel, float leptonCharge, vector<int> jetCombi){

  //--- Reconstruct neutrino partially ---//
  //---   (Neutrino Pt and M needed)   ---//
  if(verbose_ > 3) std::cout << " Inside StoreRecoInfo class " << std::endl;
  TLorentzVector Neutrino;
  float NeutrinoPx = -(lepton+Jets[bLeptIndex]+Jets[bHadrIndex]+Jets[light1Index]+Jets[light2Index]).Px();
  float NeutrinoPy = -(lepton+Jets[bLeptIndex]+Jets[bHadrIndex]+Jets[light1Index]+Jets[light2Index]).Py();
  Neutrino.SetPxPyPzE(NeutrinoPx, NeutrinoPy, 0.0, sqrt(NeutrinoPx*NeutrinoPx + NeutrinoPy*NeutrinoPy));
  Neutrino.SetPxPyPzE(NeutrinoPx, NeutrinoPy, 0.0, Neutrino.Pt());                        //Reset the Neutrino Energy to get the correct precision
  if(verbose_ > 3) std::cout << " Mass value for the neutrino : " << Neutrino.M() << " \n" << std::endl;
      
  //--- Initialize needed vectors ---//
  vector<TLorentzVector*> LHCORecoVector(6);
  vector<int> MadGraphRecoId(6,4);
  vector<float> MGRecoBtagId(6,0.0);
  MGRecoBtagId[0] = 1; MGRecoBtagId[3] = 1;   //Event selection is already performed and b-jets are always at the same position!!

  //--- Check whether the event is correctly reconstructed  ---// 
  //---  (jetCombi is initialized to 9999 for all dataSets) ---//
  bool jetCombiFound = false;
  bool EventCorrectlyMatched = false;
  if(jetCombi[0] != 9999 && jetCombi[1] != 9999 && jetCombi[2] != 9999 && jetCombi[3] != 9999){
    jetCombiFound = true;
    if( bLeptIndex == jetCombi[0] && bHadrIndex == jetCombi[1]    &&
       (light1Index == jetCombi[2] || light1Index == jetCombi[3]) &&
       (light2Index == jetCombi[3] || light2Index == jetCombi[3]) )
      EventCorrectlyMatched = true;
  }
   
  //--- Plot the generator-level distributions for Pt, eta and dR ---//
  //---   --> Will be used for MG cuts comparison                 ---//
  histo1D["Pt_Reco_LightJets"]->Fill(Jets[light1Index].Pt()); histo1D["Pt_Reco_LightJets"]->Fill(Jets[light2Index].Pt());
  histo1D["Pt_Reco_BJets"]->Fill(Jets[bHadrIndex].Pt());      histo1D["Pt_Reco_BJets"]->Fill(Jets[bLeptIndex].Pt()); 
  if(decayChannel == 0) histo1D["Pt_Reco_Muon"]->Fill(lepton.Pt());
  if(decayChannel == 1) histo1D["Pt_Reco_Electron"]->Fill(lepton.Pt());
  histo1D["Pt_Reco_MET"]->Fill(Neutrino.Pt());

  histo1D["Eta_Reco_LightJets"]->Fill(Jets[light1Index].Eta()); histo1D["Eta_Reco_LightJets"]->Fill(Jets[light2Index].Eta());
  histo1D["Eta_Reco_BJets"]->Fill(Jets[bHadrIndex].Eta());      histo1D["Eta_Reco_BJets"]->Fill(Jets[bLeptIndex].Eta());
  if(decayChannel == 0) histo1D["Eta_Reco_Muon"]->Fill(lepton.Eta());
  if(decayChannel == 1) histo1D["Eta_Reco_Electron"]->Fill(lepton.Eta());
  histo1D["Eta_Reco_MET"]->Fill(Neutrino.Eta());

  histo1D["deltaR_Reco_LightJets"]->Fill(Jets[light1Index].DeltaR(Jets[light2Index]));
  histo1D["deltaR_Reco_BJets"]->Fill(Jets[bHadrIndex].DeltaR(Jets[bLeptIndex]));
  histo1D["deltaR_Reco_LightWithB"]->Fill(Jets[light1Index].DeltaR(Jets[bHadrIndex])); histo1D["deltaR_Reco_LightWithB"]->Fill(Jets[light1Index].DeltaR(Jets[bLeptIndex])); 
  histo1D["deltaR_Reco_LightWithB"]->Fill(Jets[light2Index].DeltaR(Jets[bHadrIndex])); histo1D["deltaR_Reco_LightWithB"]->Fill(Jets[light2Index].DeltaR(Jets[bLeptIndex]));
  if(decayChannel == 0){
    histo1D["deltaR_Reco_MuonWithJet"]->Fill(lepton.DeltaR(Jets[light1Index])); histo1D["deltaR_Reco_MuonWithJet"]->Fill(lepton.DeltaR(Jets[light2Index]));
    histo1D["deltaR_Reco_MuonWithB"]->Fill(lepton.DeltaR(Jets[bHadrIndex]));    histo1D["deltaR_Reco_MuonWithB"]->Fill(lepton.DeltaR(Jets[bLeptIndex]));
  }
  if(decayChannel == 1){
    histo1D["deltaR_Reco_ElectronWithJet"]->Fill(lepton.DeltaR(Jets[light1Index])); histo1D["deltaR_Reco_ElectronWithJet"]->Fill(lepton.DeltaR(Jets[light2Index]));
    histo1D["deltaR_Reco_ElectronWithB"]->Fill(lepton.DeltaR(Jets[bHadrIndex]));    histo1D["deltaR_Reco_ElectronWithB"]->Fill(lepton.DeltaR(Jets[bLeptIndex]));
  }

  //---  Filling of LHCO files for reco events  ---//
  if(leptonCharge < 0.0 ){ //Negative lepton events
    LHCORecoVector[0] = &Jets[bHadrIndex]; 
    LHCORecoVector[1] = &Jets[light1Index];
    LHCORecoVector[2] = &Jets[light2Index];
    LHCORecoVector[3] = &Jets[bLeptIndex];
    LHCORecoVector[4] = &lepton;
    LHCORecoVector[5] = &Neutrino;

    MadGraphRecoId[5] = 6;
    if(decayChannel == 1){//Negative electron
      MadGraphRecoId[4] = 1;
      NumberNegRecoEl++;  
      if(splitLeptCharge_) LHCOEventOutput(3,RecoOutFile[3], NumberNegRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      else                 LHCOEventOutput(3,RecoOutFile[1], NumberNegRecoEl+NumberPosRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
    }
    if(decayChannel == 0){//Negative muon
      MadGraphRecoId[4] = 2;
      NumberNegRecoMu++;
      if(splitLeptCharge_) LHCOEventOutput(1, RecoOutFile[1], NumberNegRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      else                 LHCOEventOutput(1, RecoOutFile[0], NumberNegRecoMu+NumberPosRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
    }
  }//End of negative lepton
  else if(leptonCharge > 0.0 ){ //Positive lepton events
    LHCORecoVector[0] = &Jets[bLeptIndex];
    LHCORecoVector[1] = &lepton;
    LHCORecoVector[2] = &Neutrino;
    LHCORecoVector[3] = &Jets[bHadrIndex];
    LHCORecoVector[4] = &Jets[light1Index];
    LHCORecoVector[5] = &Jets[light2Index];

    MadGraphRecoId[2] = 6;
    if(decayChannel == 1){//Positive electron
      MadGraphRecoId[1] = 1;
      NumberPosRecoEl++;
      if(splitLeptCharge_) LHCOEventOutput(2,RecoOutFile[2], NumberPosRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);	 
      else                 LHCOEventOutput(2,RecoOutFile[1], NumberPosRecoEl+NumberNegRecoEl, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);	 
    }
    if(decayChannel == 0){//Positive muon
      MadGraphRecoId[1] = 2;
      NumberPosRecoMu++;
      if(splitLeptCharge_) LHCOEventOutput(0, RecoOutFile[0], NumberPosRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      else                 LHCOEventOutput(0, RecoOutFile[0], NumberPosRecoMu+NumberNegRecoMu, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      if(EventCorrectlyMatched == true && jetCombiFound == true && splitCorrectWrong_){
        NrPosRecoMuCorrect++; 
        LHCOEventOutput(0, CorrectRecoMuPosFile, NrPosRecoMuCorrect, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      }
      else if(EventCorrectlyMatched == false && jetCombiFound == true && splitCorrectWrong_){ 
        NrPosRecoMuWrong++;   
        LHCOEventOutput(0, WrongRecoMuPosFile, NrPosRecoMuWrong, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      }
      else if(jetCombiFound == false && splitCorrectWrong_){
        NrPosRecoMuUnmatched++;
        LHCOEventOutput(0, UnmatchedRecoMuPosFile, NrPosRecoMuUnmatched, LHCORecoVector, MadGraphRecoId, MGRecoBtagId);
      }
    }
  }//End of positive lepton

}//End of class StoreRecoInfo 

void LHCOOutput::LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag){

  if(LHCOIndex == 0 || LHCOIndex == 2)
    LeptonCharge =1;
  else if(LHCOIndex == 1 || LHCOIndex == 3)
    LeptonCharge = -1;

  if(EventNumber == 1){
    outputFile << "#</MGPGSCard> " << endl;
    outputFile << "  #  typ      eta      phi       pt   jmas  ntrk  btag   had/em  dummy  dummy " << endl;
  }

  outputFile << " 0             " << EventNumber << "        6 " << endl;  //Start of a new event

  for(int ii = 0; ii < 6; ii++){
    outputFile << "  " << setprecision(1) << ii+1;
    outputFile << "    " << setprecision(1) << MGId[ii];
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    outputFile << "  " << fixed << showpoint << setprecision(4) << vector[ii]->Eta();
    if(vector[ii]->Phi() < -3.14) outputFile << "  " << setprecision(4) << -1*(vector[ii]->Phi());
    else if(vector[ii]->Phi() < 0.0) outputFile << " " << setprecision(4) << vector[ii]->Phi()+6.28;
    else                          outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    else                      outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(MGId[ii] == 2 || MGId[ii] == 1) outputFile << "   " << setprecision(3) << LeptonCharge << "   ";       //MGId equal to 1 or 2 corresponds to muon or electron
    else                               outputFile << "    0.00";
    outputFile << " " << MGBtag[ii] << "     0.00  0.00  0.00" << endl;
  }
}

void LHCOOutput::WriteLHCOPlots(TFile* outfile){

  if(GenOrReco_ == "Gen" or GenOrReco_ == "Reco"){
    //--- Use this function to create ChiSq histograms ---//
    outfile->cd();
    if(verbose_ > 3) std::cout << " Inside WriteLHCOPlots function of LHCOOutput class ! \n Histograms will be filled in file : " << outfile->GetName() << " ********************************" << std::endl;

    TDirectory* th1dir = outfile->GetDirectory("1D_histograms_LHCOOutput");   //Check whether directory already exists ..
    if(!th1dir) th1dir = outfile->mkdir("1D_histograms_LHCOOutput");          // .. and otherwise create it!
    th1dir->cd();
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
      temp->Write();
    }
    
    TDirectory* th2dir = outfile->GetDirectory("2D_histograms_LHCOOutput");
    if(!th2dir) th2dir = outfile->mkdir("2D_histograms_LHCOOutput");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
      TH2F *temp = it->second;
      temp->Write();
    }
    outfile->cd(); 
  }
}
