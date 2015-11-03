#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

using namespace std;
using namespace TopTree;

class LHCOOutput{

public:
  LHCOOutput(int, bool, bool, bool);
  ~LHCOOutput();

  void Initialize(std::string, std::string dataSetName = "");
  void StoreGenInfo(vector<TRootMCParticle*> mcParticles);
  void StoreRecoInfo(TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int> selJetCombi, int decayChannelEnum, float leptonCharge, ofstream &EvtNrInfo, int CWUIndex); 
  void WriteLHCOPlots(TFile*);

  bool GenEventContentCorrect()    {return CorrectGenEvtContent;};
  int getLeptonType()              {return leptonType;};
  TLorentzVector* getGenLeptTop()  {return GenLeptonicTop;}; 
  TLorentzVector* getGenLeptW()    {return GenLeptonicW;};
  TLorentzVector* getGenLepton()   {return GenLepton;};
  TLorentzVector* getGenNeutrino() {return GenNeutrino;}; 
  //TLorentzVector* getGenHadrTop()  {return GenHadronicTop;}; 
  //TLorentzVector* getGenHadrW()    {return GenHadronicW;}; 

private:
  void LHCOEventOutput(float leptCharge, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag);

  TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;
  TLorentzVector *GenLeptonicTop, *GenLeptonicW, *GenLepton, *GenNeutrino;
  //TLorentzVector *GenHadronicTop, *GenHadronicW;
  unsigned int NumberNegativeElectrons, NumberNegativeMuons, NumberPositiveElectrons, NumberPositiveMuons, WrongEvtCounter;
  unsigned int NumberNegRecoEl, NumberNegRecoMu, NumberPosRecoEl, NumberPosRecoMu;
  unsigned int CWUEvtNr[3];
  bool CorrectGenEvtContent;
  ofstream MWOutFile[4], WrongGenFile, CWURecoFile[3];
  int verbose_;
  bool writeOutput_, splitLeptCharge_, splitCorrectWrong_;
  std::string GenOrReco_;

  enum LeptonType_t {muPlus, muMinus, elPlus, elMinus, notFound};
  LeptonType_t leptonType;

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
};

#endif
