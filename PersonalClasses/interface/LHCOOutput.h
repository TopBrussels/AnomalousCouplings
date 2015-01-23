#ifndef LHCOOutput_h
#define LHCOOutput_h

#include <iostream>
#include "TopTreeProducer/interface/TRootMCParticle.h"
#include "TLorentzVector.h"

using namespace std;
using namespace TopTree;

class LHCOOutput{

  int LeptonCharge;
public:
  void LHCOEventOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId, std::vector<float> MGBtag); //TRootMCParticle info needed?
  void LHCOEventRecoOutput(int LHCOIndex, ostream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId, std::vector<float> MGBtag); //TRootMCParticle info needed?

};

#endif