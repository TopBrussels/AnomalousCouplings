#ifndef TFLight_h
#define TFLight_h

#include "TObject.h"
#include "Rtypes.h"
#include "TLorentzVector.h"

using namespace std;

class TFLight : public TObject
{

  public:
    TFLight():
      TObject()
      ,fullScaleFactor_(0)
      ,selectedJets_()
      ,selectedLepton_()
      ,correctJetCombi_()
      ,genVectorLight1_()
      ,genVectorLight2_()
      ,genVectorHadrB_()
      ,genVectorLeptB_()
      ,genVectorLepton_()
      {;}
    ~TFLight() {;}

    float fullScaleFactor()        const {return fullScaleFactor_;}
    vector<TLorentzVector> selectedJets() const {return selectedJets_;}
    TLorentzVector selectedLepton() const {return selectedLepton_;}
    vector<int> correctJetCombi() const {return correctJetCombi_;}

    TLorentzVector genVectorLight1() const {return genVectorLight1_;}
    TLorentzVector genVectorLight2() const {return genVectorLight2_;}
    TLorentzVector genVectorHadrB()  const {return genVectorHadrB_;}
    TLorentzVector genVectorLeptB()  const {return genVectorLeptB_;}
    TLorentzVector genVectorLepton() const {return genVectorLepton_;}

    void setFullScaleFactor(float fullScaleFactor) {fullScaleFactor_ = fullScaleFactor;}
    void setSelectedJets(vector<TLorentzVector> selectedJets) {selectedJets_ = selectedJets;}
    void setSelectedLepton( TLorentzVector selectedLepton) {selectedLepton_ = selectedLepton;}
    void setCorrectJetCombi( vector<int> correctJetCombi) {correctJetCombi_ = correctJetCombi;}

    void setGenVectorLight1(TLorentzVector genVectorLight1) {genVectorLight1_ = genVectorLight1;}
    void setGenVectorLight2(TLorentzVector genVectorLight2) {genVectorLight2_ = genVectorLight2;}
    void setGenVectorHadrB( TLorentzVector genVectorHadrB)  {genVectorHadrB_  = genVectorHadrB;}
    void setGenVectorLeptB( TLorentzVector genVectorLeptB)  {genVectorLeptB_  = genVectorLeptB;}
    void setGenVectorLepton(TLorentzVector genVectorLepton) {genVectorLepton_ = genVectorLepton;}

  protected:
    float fullScaleFactor_;
    vector<TLorentzVector> selectedJets_;
    TLorentzVector selectedLepton_;
    vector<int> correctJetCombi_;

    TLorentzVector genVectorLight1_;
    TLorentzVector genVectorLight2_;
    TLorentzVector genVectorHadrB_;
    TLorentzVector genVectorLeptB_;
    TLorentzVector genVectorLepton_;

  ClassDef(TFLight,2); 
};

#endif
