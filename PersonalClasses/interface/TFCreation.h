#ifndef TFCreation_h
#define TFCreation_h

#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <map>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <fstream>
#include <sstream>

using namespace std;
//using namespace TopTree;

class TFCreation{

    public:
        TFCreation(int);
        ~TFCreation();
	void InitializeVariables();
	void FillHistograms(TLorentzVector* hadrWJet1, TLorentzVector* hadrWJet2, TLorentzVector* hadrBJet, TLorentzVector* leptBJet, TLorentzVector* lepton, TLorentzVector* selHadrWJet1, TLorentzVector* selHadrWJet2, TLorentzVector* selHadrBJet, TLorentzVector* selLeptBJet, TLorentzVector* selLepton, int enumDecayChannel);
	void CalculateTF(bool, bool, bool, bool);
        void CalculateTFFromFile(string, bool, int, bool, bool, float[], bool, TFile*, int, TFile*);
	void FitSliceClassCode(TH2F*, bool);
        std::vector<double> SetFitRange(std::string, int);
	void SetStartValuesDoubleGaussian(int, bool, std::string);
	void WriteTF(ostream &output, ostream &card, ostream &TF, std::string, std::string);
        void PlotDlbGaus(TH2F*, TFile*);
	void WritePlots(TFile*);

    private:
	map<string,TH1F*> histo1D;
	map<string,TH2F*> histo2D;

	template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

	TF1 *doubleGaussianFit, *caloEnergyFit;
	TH1D **hlist;
        float* startValuesArray;
        TF1 AllCaloEnergyFits[30];
        std::string EtaBin[5];    //4 eta bins and one extra for all events!
        std::string EtaTitle[5];
        double EtaValues[6];
        int nParsFit_, nEtaBins_;
        std::string parnames_[5], ParName_[5];
        //int NarrowGaus[3];
        //int WideGaus[3];
};
#endif
//
