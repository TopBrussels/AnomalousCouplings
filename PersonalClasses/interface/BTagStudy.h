#ifndef BTagStudy_h
#define BTagStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"

using namespace std;
using namespace TopTree;

class BTagStudy{

  public:
    BTagStudy(int, vector<Dataset*>);
    ~BTagStudy();

    void CalculateJets(vector<TLorentzVector>, vector<float> bTagValues, vector<int> jetCombi, TLorentzVector lepton, Dataset*, float);
    void ReturnBTagTable(std::string);
    void ReturnThirdJetTable();    //Still to fill
    void CreateHistograms(TFile*, bool, std::string, int); 
 
    int getBHadrIndex(int bTagNr)  {return bHadrIndex[bTagNr];};
    int getBLeptIndex(int bTagNr)  {return bLeptIndex[bTagNr];};
    int getLight1Index4Jets(int bTagNr) {return light1Index4Jets[bTagNr];};   //Need to think how to return the light jets ... (not possible to differentiate between 4 and 5 jet case)
    int getLight2Index4Jets(int bTagNr) {return light2Index4Jets[bTagNr];};
    int getLight1Index5Jets(int bTagNr) {return light1Index5Jets[bTagNr];};  
    int getLight2Index5Jets(int bTagNr) {return light2Index5Jets[bTagNr];};

    float Mlb, Mqqb, SigmaMlb, SigmaMqqb;

    //Functions to ask everything inside the analyzer code!                          -->Is this getter still necessary?
    vector<int> getbTaggedJets(int OptionNr)   {return bTaggedJetNr[OptionNr];};
    vector<int> getNonbTaggedJets(int OptionNr){return NonbTaggedJetNr[OptionNr];};
    vector<int> getLightJets(int OptionNr)     {return LightJetNr[OptionNr];};

  private:
    void CompareJetCombi(vector<int> jetCombi, int bTagNr, int jetCase, int lightJetOne, int lightJetTwo);
    void InitializeBegin(vector<Dataset*>);
    void ResetEventArrays();
    void CalculateMlbChiSq( int bTagNr, TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int>, Dataset*, float); 
    vector<int> CalculateMqqbChiSq( int bTagNr, vector<TLorentzVector> Jets); 

    int LowestChiSqMlb[6], LowestChiSqMqqb[6];
    int bHadrIndex[6], bLeptIndex[6], light1Index4Jets[6], light2Index4Jets[6], light1Index5Jets[6], light2Index5Jets[6], light1IndexMqqb[6], light2IndexMqqb[6];
    vector<float> ChiSquaredMlb[6], ChiSquaredMqqb[6];
    int NotReconstructedEvent[2][6], CorrectlyMatched[2][3][6], atLeastOneWrongMatch[2][3][6];

    int verbose;
    float BJetWP[6], LightJetWP[6];  
    std::string OptionName[6], BTitle[6], BName[6];
    ofstream evtSelOutput[2];
    vector<int> bTaggedJetNr[6], NonbTaggedJetNr[6], LightJetNr[6];

    //4- and 5-jet case comparison
    //int EventWithTwoLightJets[6];
    //int EventWithThreeLightJets[6];
    //int EventWithTwoLightJetsAndBTagged[6];
    //int EventWithThreeLightJetsAndBTagged[6];
    //int thirdJetIsActualQuark[6][2], secondJetIsActualQuark[6][2], firstJetIsActualQuark[6][2];
    //int thirdJetIsCorrectQuark[6];
    //int thirdJetIsGoodQuark[6];    //3rd quark is one of the quarks when the two b-jets are correctly matched!

    map<std::string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
    map<string,MultiSamplePlot*> MSPlot;
    int nrDatasets_;

};

#endif
