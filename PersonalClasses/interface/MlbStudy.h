#ifndef MlbStudy_h
#define MlbStudy_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TopTreeProducer/interface/TRootJet.h"

using namespace std;
using namespace TopTree;

class MlbStudy{

  private:
    float Mlb, Mqqb, SigmaMlb, SigmaMqqb;
    int bLeptIndex, bHadrIndex;
    int LowestChiSqMlb;
    float ChiSquaredMlb[2];
  
  public:
    MlbStudy(int);
    ~MlbStudy();
    void CalculateMlbChiSq( vector<int> bJetIndex, TLorentzVector* lepton, vector<TRootJet*> Jets); 
    //void CalculateMqqbChiSq( vector<int> lightJetIndex, vector<TRootJet*> Jets); 
    void InitializeMlb();
    //void InitializeMqqb();

    int getBLept() {return bLeptIndex;};
    int getBHadr() {return bHadrIndex;}; 
    
    float getMlbFit()  {return Mlb;};
    float getMqqbFit() {return Mqqb;};

    //------------------------------------------------------------------------------------------------------------------------------------------//
/*
   void calculateChiSquared(vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, TLorentzVector* lepton, vector<TRootJet*> Jets, float MassMlb, float SigmaMlb, float MassMqqb, float SigmaMqqb);
   void getIndices(int LowestChiSqIndex);
   void calculateEfficiency(int option, vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, int NrConsideredBTagOptions, int ChiSqCutValue);
   void saveNumbers(std::string NameOfOption[6], int WhichJets, int NrOptionsConsidered, int ChosenOption, std::string StrForChiSqCutValue);
   void initializePerEvent();
   void initializeBegin(int);
   void WritePlots(TFile*);

   int chosenBHadr, chosenBLept, chosenQuark1, chosenQuark2;
   int LowestChiSq, LowestChiSq4Jets, CorrectChiSq;           //Since Mqqb Values are set artificially high, only the first two chi-sq will be the lowest one in the case of only 4 jets!
						              //Still need to separate variables since then inside the analyzer code can be chosen whether the 4- or 5-jet case will be used!
   float MlbValues[6], MqqbValues[6], ChiSquared[6];

   vector<float> h_ChiSqCorrect[6], h_ChiSqCorrectFound[6], h_ChiSqMinimum[6], h_ChiSqNotMinimum[6], h_ChiSqWrong[6];
   vector<float> h_ChiSqCorrectWhenMatched[6], h_ChiSqMinimumWhenMatched[6], h_ChiSqNotMinimumWhenMatched[6], h_ChiSqAllWhenNotMatched[6], h_ChiSqMinimumWhenCorrect[6], h_ChiSqMinimumWhenWrong[6], h_ChiSqDiffWhenWrong[6],h_ChiSqAllMatchingNotExisting[6] ;

   int NumberMatchedEvents[6], NumberNotMatchedEvents[6], CorrectOptionAvailable[6], CorrectBOptionChosen[6], WrongBOptionChosen[6], CorrectEventChosen[6],WrongEventChosen[6],CorrectMatchingNotExisting[6];
   int ThirdQuarkChosen[6], ThirdQuarkCorrectChosen[6], ThirdQuarkShouldBeChosen[6], CorrectLightJetsChosen[6], CorrectLightJetsWithThirdChosen[6];
   int FirstQuarkChosen[6], FirstQuarkCorrectChosen[6], SecondQuarkChosen[6], SecondQuarkCorrectChosen[6];
   int FirstQuarkCorrect[6], SecondQuarkCorrect[6], ThirdQuarkCorrect[6];

   int getLowestChiSqIndex()      {return LowestChiSq;};
   int getLowestChiSq4JetsIndex() {return LowestChiSq4Jets;};
   int getCorrectChiSq()          {return CorrectChiSq;};
   int getChosenBLept()           {return chosenBLept;};
   int getChosenBHadr()           {return chosenBHadr;};
   int getChosenQuark1()          {return chosenQuark1;};
   int getChosenQuark2()          {return chosenQuark2;};

   vector<float> getChiSqCorrect(int jetCase)      {return h_ChiSqCorrect[jetCase];     };
   vector<float> getChiSqCorrectFound(int jetCase) {return h_ChiSqCorrectFound[jetCase];};
   vector<float> getChiSqMinimum(int jetCase)      {return h_ChiSqMinimum[jetCase];     };
   vector<float> getChiSqNotMinimum(int jetCase)   {return h_ChiSqNotMinimum[jetCase];  };
   vector<float> getChiSqWrong(int jetCase)        {return h_ChiSqWrong[jetCase];       };

   vector<float> getChiSqCorrectWhenMatched(int jetCase)    {return h_ChiSqCorrectWhenMatched[jetCase];   };
   vector<float> getChiSqMinimumWhenMatched(int jetCase)    {return h_ChiSqMinimumWhenMatched[jetCase];   };
   vector<float> getChiSqNotMinimumWhenMatched(int jetCase) {return h_ChiSqNotMinimumWhenMatched[jetCase];};
   vector<float> getChiSqAllWhenNotMatched(int jetCase)     {return h_ChiSqAllWhenNotMatched[jetCase];    };
   vector<float> getChiSqMinimumWhenCorrect(int jetCase)    {return h_ChiSqMinimumWhenCorrect[jetCase];   };
   vector<float> getChiSqMinimumWhenWrong(int jetCase)      {return h_ChiSqMinimumWhenWrong[jetCase];     };

   vector<float> getChiSqDiffWhenWrong(int jetCase) {return h_ChiSqDiffWhenWrong[jetCase];};

   float getMlbValue(int order)     {return MlbValues[order];};
   float getMqqbValue(int order)    {return MqqbValues[order];};
   float getChiSquared(int order) {return ChiSquared[order];};

   int getNrMatchedEvents( int OptionNr) {return NumberMatchedEvents[OptionNr];};
   int getNrCorrectOptionsAvailable(int OptionNr) {return CorrectOptionAvailable[OptionNr];};
   int getNrCorrectBOptionChosen(int OptionNr) {return CorrectBOptionChosen[OptionNr];};
   int getNrWrongBOptionChosen(int OptionNr) {return WrongBOptionChosen[OptionNr];};

  private:
    map<string,TH1F*> histo1D;
    map<string,TH2F*> histo2D;
    std::string Title[3], Name[3];
*/
};

#endif
