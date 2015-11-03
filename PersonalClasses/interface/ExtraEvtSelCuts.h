#ifndef ExtraEvtSelCuts_h
#define ExtraEvtSelCuts_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "TopTreeProducer/interface/TRootJet.h"
//#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "TopTreeAnalysisBase/Content/interface/Dataset.h"

using namespace std;
using namespace TopTree;

class ExtraEvtSelCuts{

  public:
    ExtraEvtSelCuts(float TopMassHadr, float sTopMassHadr, float WMassHadr, float sWMassHadr, bool, int chiSqCut = -1, int massSigma = -1);
    ~ExtraEvtSelCuts();

    void Initialize(std::string bTagTitle, std::string dataSetName);
    bool KeepEvent(vector<int>, TLorentzVector, vector<TLorentzVector>, vector<int>, float chiSqValue, int CWUIndex, int decayCh);
    void StoreCutInfluence(TFile*);

  private:
    float MTopHadr_, MWHadr_, sig_MTopH_, sig_MWH_;
    float ChiSqCutVal_[5], MTop_Down_[5], MTop_Up_[5], MW_Down_[5], MW_Up_[5], MassWindow_Sigmas_[5];
    int nrCuts_, chosenChiSq_, chosenMassWindow_;
    std::string histTitle_;
    bool oneBTag_;

    int OriginalEvts_[3], NrEvts_AllChosenCuts_[3];
    int NrEvts_ChiSq_[3][5], NrEvts_MT_[3][5], NrEvts_MW_[3][5], NrEvts_MComb_[3][5];

    int WrongEvts_Mu, WrongEvts_El, CorrEvts_Mu, CorrEvts_El;

    map<string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
};

#endif
