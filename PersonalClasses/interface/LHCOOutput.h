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
  LHCOOutput(int, bool);
  ~LHCOOutput();

  void Initialize(std::string);
  void StoreGenInfo(vector<TRootMCParticle*> mcParticles);
  void StoreRecoInfo(TLorentzVector lepton, vector<TLorentzVector> Jets,int bLept, int bHadr, int light1, int light2, int decayChannelEnum, float leptonCharge, vector<int> jetCombi); 
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
  void LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag);

  TRootMCParticle *Top,*TopBar,*Bottom, *BottomBar,*Lepton,*NeutrinoMC,*WPlus,*WMinus,*Light,*LightBar;
  TLorentzVector *GenLeptonicTop, *GenLeptonicW, *GenLepton, *GenNeutrino;
  //TLorentzVector *GenHadronicTop, *GenHadronicW;
  unsigned int NumberNegativeElectrons, NumberNegativeMuons, NumberPositiveElectrons, NumberPositiveMuons, WrongEvtCounter;
  unsigned int NumberNegRecoEl, NumberNegRecoMu, NumberPosRecoEl, NumberPosRecoMu, NrPosRecoMuCorrect, NrPosRecoMuWrong, NrPosRecoMuUnmatched;
  bool CorrectGenEvtContent;
  ofstream GenOutFile[4], RecoOutFile[4], WrongGenFile, CorrectRecoMuPosFile, WrongRecoMuPosFile, UnmatchedRecoMuPosFile;
  int verbose_, LeptonCharge;
  bool writeOutput_;
  std::string GenOrReco_;

  enum LeptonType_t {muPlus, muMinus, elPlus, elMinus, notFound};
  LeptonType_t leptonType;

  map<string,TH1F*> histo1D;
  map<string,TH2F*> histo2D;
};

#endif
