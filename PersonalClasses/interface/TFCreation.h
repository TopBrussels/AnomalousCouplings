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
#include "TGraphErrors.h"

#include <fstream>
#include <sstream>

using namespace std;
//using namespace TopTree;

class TFCreation{

    public:
        TFCreation(int, std::string, bool);
        ~TFCreation();
	void InitializeVariables();
	void FillHistograms(TLorentzVector* hadrWJet1, TLorentzVector* hadrWJet2, TLorentzVector* hadrBJet, TLorentzVector* leptBJet, TLorentzVector* lepton, TLorentzVector* selHadrWJet1, TLorentzVector* selHadrWJet2, TLorentzVector* selHadrBJet, TLorentzVector* selLeptBJet, TLorentzVector* selLepton, int enumDecayChannel, float scaleFactor);
        void CalculateTFFromFile(string, bool, int, bool, bool, float[], bool, TFile*, int, TFile*);
	void FitSliceClassCode(TH2F*, bool, int);
        std::vector<double> SetFitRange(std::string, unsigned int, double[]);
	void SetStartValuesDoubleGaussian(int, bool, std::string);
	void WriteTF(ostream &output, ostream &card, ostream &TF, ostream &latex, std::string, std::string, int);
	void WritePlots(TFile*);

    private:
	map<string,TH1F*> histo1D;
	map<string,TH2F*> histo2D;

	template <typename T> string tostr(const T& t) { ostringstream os; os<<t; return os.str(); }

	TF1 *doubleGaussianFit, *caloEnergyFit;
        float* startValuesArray;
        TF1 AllCaloEnergyFits[30];
        bool doFits_;

        double EtaValues[6];
        std::string EtaValue[6], EtaBin[5], EtaTitle[5];

        int nParsFit_, nEtaBins_;
        std::string parnames_[6], ParName_[6];
        TGraphErrors *grE_ParamFit[36];
        double FitMin_[36], FitMax_[36];
        //int NarrowGaus[3];
        //int WideGaus[3];
        TFile* caloFitFile;
        vector<std::string> indivFitHistos;
};
#endif
//
