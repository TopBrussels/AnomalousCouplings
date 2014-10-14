#include <iomanip>
#include "../interface/LHCOOutput.h"

#include<fstream>
#include <sstream>

void LHCOOutput::LHCOEventOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, std::vector<TRootMCParticle*> vector, std::vector<int> MGId){

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
    else outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    else outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(ii == 1) outputFile << "   " << setprecision(3) << LeptonCharge;
    else outputFile << "    0.00";
    outputFile << "  0.00     0.00  0.00  0.00" << endl;
  }
}

void LHCOOutput::LHCOEventRecoOutput(int LHCOIndex, ofstream &outputFile, unsigned int EventNumber, std::vector<TLorentzVector*> vector, std::vector<int> MGId){

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
    else outputFile << "  " << setprecision(4) << vector[ii]->Phi();
    outputFile << "  " << setprecision(4) << vector[ii]->Pt();
    cout.setf(ios::fixed,ios::floatfield);  //Add zero to obtain the asked number of digits
    cout.precision(4);
    if(vector[ii]->M() > 0.0) outputFile << "  " << fixed << showpoint << setprecision(3) << vector[ii]->M();
    else outputFile << "  " << fixed << showpoint << setprecision(3) << " 0.00";
    if(ii == 1) outputFile << "   " << setprecision(3) << LeptonCharge;
    else outputFile << "    0.00";
    outputFile << "  0.00     0.00  0.00  0.00" << endl;
  }
}
