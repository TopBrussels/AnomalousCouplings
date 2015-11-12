#ifndef AnomCoupLight_h
#define AnomCoupLight_h

#include "TObject.h"
#include "Rtypes.h"
#include "TLorentzVector.h"

using namespace std;

class AnomCoupLight : public TObject
{

  public:
    AnomCoupLight():
      TObject()
      ,eventID_(0)
      ,runID_(0)
      ,lumiBlockID_(0)
      ,nPV_(0)
      ,nTruePU_(0)
      ,scaleFactor_(0)
      ,selectedJets_()
      ,selectedJetsPartonFlavour_()
      ,CSVbTag_()
      ,selectedLepton_()
      ,decayChannel_(0)
      ,leptonCharge_(0)
      ,correctJetCombi_()
      ,met_()
      ,genCosTh_(0)
      ,genVectorLight1_()
      ,genVectorLight2_()
      ,genVectorHadrB_()
      ,genVectorLeptB_()
      ,genVectorLepton_()
      {;}
    ~AnomCoupLight() {;}

    unsigned int eventID()     const {return eventID_;}
    unsigned int runID()       const {return runID_;}
    unsigned int lumiBlockID() const {return lumiBlockID_;}
    unsigned int nPV()         const {return nPV_;}
    unsigned int nTruePU()     const {return nTruePU_;}
    float scaleFactor()        const {return scaleFactor_;}

    vector<TLorentzVector> selectedJets() const {return selectedJets_;}
    vector<int> selectedJetsPartonFlavour() const {return selectedJetsPartonFlavour_;}
    vector<float> CSVbTag() const {return CSVbTag_;}
    TLorentzVector selectedLepton() const {return selectedLepton_;}
    unsigned int decayChannel() const {return decayChannel_;}
    float leptonCharge() const {return leptonCharge_;}
    vector<int> correctJetCombi() const {return correctJetCombi_;}
    TLorentzVector met() const {return met_;}
    float genCosTh() const {return genCosTh_;}

    TLorentzVector genVectorLight1() const {return genVectorLight1_;}
    TLorentzVector genVectorLight2() const {return genVectorLight2_;}
    TLorentzVector genVectorHadrB()  const {return genVectorHadrB_;}
    TLorentzVector genVectorLeptB()  const {return genVectorLeptB_;}
    TLorentzVector genVectorLepton() const {return genVectorLepton_;}

    void setEventId(unsigned int eventID) {eventID_ = eventID;}
    void setRunId(unsigned int runID) {runID_ = runID;}
    void setLumiBlockId(unsigned int lumiBlockID) {lumiBlockID_ = lumiBlockID;}
    void setNPV(unsigned int nPV) {nPV_ = nPV;}
    void setNTruePU(unsigned int nTruePU) {nTruePU_ = nTruePU;}
    void setScaleFactor(float scaleFactor) {scaleFactor_ = scaleFactor;}

    void setSelectedJets(vector<TLorentzVector> selectedJets) {selectedJets_ = selectedJets;}
    void setSelectedJetsPartonFlavour(vector<int> selectedJetsPartonFlavour){selectedJetsPartonFlavour_ = selectedJetsPartonFlavour;}
    void setBTagCSV(vector<float> CSVbTag) {CSVbTag_ = CSVbTag;}
    void setSelectedLepton( TLorentzVector selectedLepton) {selectedLepton_ = selectedLepton;}
    void setDecayChannel(unsigned int decayChannel) {decayChannel_ = decayChannel;}
    void setLeptonCharge( float leptonCharge) {leptonCharge_ = leptonCharge;}
    void setCorrectJetCombi( vector<int> correctJetCombi) {correctJetCombi_ = correctJetCombi;}
    void setMET(TLorentzVector met) {met_ = met;}
    void setGenCosTheta(float genCosTh) {genCosTh_ = genCosTh;}

    void setGenVectorLight1(TLorentzVector genVectorLight1) {genVectorLight1_ = genVectorLight1;}
    void setGenVectorLight2(TLorentzVector genVectorLight2) {genVectorLight2_ = genVectorLight2;}
    void setGenVectorHadrB( TLorentzVector genVectorHadrB)  {genVectorHadrB_  = genVectorHadrB;}
    void setGenVectorLeptB( TLorentzVector genVectorLeptB)  {genVectorLeptB_  = genVectorLeptB;}
    void setGenVectorLepton(TLorentzVector genVectorLepton) {genVectorLepton_ = genVectorLepton;}

  protected:
    unsigned int eventID_;
    unsigned int runID_;
    unsigned int lumiBlockID_;
    unsigned int nPV_;
    unsigned int nTruePU_;
    float scaleFactor_;

    vector<TLorentzVector> selectedJets_;
    vector<int> selectedJetsPartonFlavour_;
    vector<float> CSVbTag_;
    TLorentzVector selectedLepton_;
    unsigned int decayChannel_;
    float leptonCharge_;
    vector<int> correctJetCombi_;
    TLorentzVector met_;
    float genCosTh_;

    TLorentzVector genVectorLight1_;
    TLorentzVector genVectorLight2_;
    TLorentzVector genVectorHadrB_;
    TLorentzVector genVectorLeptB_;
    TLorentzVector genVectorLepton_;

  ClassDef(AnomCoupLight,2); 
};

#endif
