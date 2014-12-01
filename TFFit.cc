//user code
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"   //Needed to load TRootMCParticle & TRootJet, which is used in TFCreation.h
#include "TFile.h"
#include "TH2.h"
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TBranch.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include <TApplication.h>

//Specific code for anomalous couplings analysis:
#include "AnomalousCouplings/PersonalClasses/interface/TFCreation.h"
#include "AnomalousCouplings/PersonalClasses/interface/TFnTuple.h"

using namespace std;
//using namespace TopTree;

int main (int argc, char **argv)
{
    TApplication theApp("App", &argc, argv); //Needed to run on local Linux!
    clock_t start = clock();
  
    cout << "***********************************************" << endl;  
    cout << " Beginning of the program for fitting the TF ! " << endl;
    cout << "*********************************************** \n" << endl;
  
    ////////////////////////////////////////////////////////////////////
    //  Choose whether created plots are used or Tree information !!  //
    ////////////////////////////////////////////////////////////////////
    bool CreateTFFromTree = true;
    bool RunFitForTF = true;
    int nEtaBins = 4;
    //Used classes
    TFCreation tfCreation(nEtaBins);

    if(CreateTFFromTree){
        //Load the TFTree information
        vector<string> inputTFRoot;
        inputTFRoot.push_back("TFInformation/TransferFunctionTree_AllEvts.root");

        for(unsigned int iDataSet = 0; iDataSet <inputTFRoot.size(); iDataSet++){
            TFile* inputTFFile = new TFile(inputTFRoot[iDataSet].c_str(),"READ");

            TTree* inputTFTree = (TTree*) inputTFFile->Get("TFTree");
            TBranch* m_br = (TBranch*) inputTFTree->GetBranch("TheTFTree");
            TFnTuple* tfNTuple = 0;
            m_br->SetAddress(&tfNTuple);

            //Set the number of selected events (for loop on events):
            int nEvent = inputTFTree->GetEntries(); 
            //int nEvent = 100000;
            std::cout << " *** Looking at dataset " << iDataSet+1 << "/" << inputTFRoot.size() << " with " << nEvent << " selected events! \n " << std::endl;

            //Initialize the TFCreation class (create all histograms):
            tfCreation.InitializeVariables(nEtaBins); //Add option of nr eta bins here!

            //Read in the TLorenztVectors:
            TLorentzVector genPart[5], recoPart[5];
            enum DecayChannel_t {isSemiMu, isSemiEl};
            DecayChannel_t decayChannel;
            for(unsigned int iEvt = 0; iEvt < nEvent; iEvt++){
                if(iEvt%10000 == 0)
                    std::cout<<"Processing the "<<iEvt<<"th event (" << ((double)iEvt/(double)inputTFTree->GetEntries())*100  << "%)" << flush<<"\r";

                inputTFTree->GetEvent(iEvt);
                genPart[0] = tfNTuple->genVectorLight1();
                genPart[1] = tfNTuple->genVectorLight2();
                genPart[2] = tfNTuple->genVectorHadrB();
                genPart[3] = tfNTuple->genVectorLeptB();
                genPart[4] = tfNTuple->genVectorLepton();

                recoPart[0] = tfNTuple->recoVectorLight1();
                recoPart[1] = tfNTuple->recoVectorLight2();
                recoPart[2] = tfNTuple->recoVectorHadrB();
                recoPart[3] = tfNTuple->recoVectorLeptB();
                recoPart[4] = tfNTuple->recoVectorLepton();

                if(genPart[4].M() <= 0.05) decayChannel = isSemiEl; //Electron channel --> decayChannel == 1
                else                       decayChannel = isSemiMu; //Muon     channel --> decayChannel == 0

                //Fill the histograms of the TFCreation class!
                tfCreation.FillHistograms( &genPart[0], &genPart[1], &genPart[2], &genPart[3], &genPart[4], &recoPart[0], &recoPart[1], &recoPart[2], &recoPart[3], &recoPart[4], decayChannel, nEtaBins);
            }//Loop on events

            TFile* fillFile = new TFile("TFInformation/PlotsForTransferFunctions_FromTree.root","RECREATE");
            std::cout << "    ----> Information writen in file : " << fillFile->GetName() << std::endl << std::endl;
            tfCreation.WritePlots(fillFile);
            fillFile->Close();
            delete fillFile;

            inputTFFile->Close();
            delete inputTFFile;
        }
    }//CreateTFFromTree = true loop

    if(RunFitForTF == true){

        std::cout << " *** Starting to perform the double Gaussian fits  " << std::endl;

        //Set which TFFile should be used
        TFile *readFile, *writeFile;
        if(CreateTFFromTree == false){
		readFile = new TFile("TFInformation/PlotsForTransferFunctions_AllEvts_UpdatedElAndMu.root","READ");
		writeFile = new TFile("TFInformation/CreatedTFFromDistributions_AllEvts_UpdatedElAndMu.root","RECREATE");
	}
        else{
		readFile = new TFile("TFInformation/PlotsForTransferFunctions_FromTree.root","READ");
		writeFile = new TFile("TFInformation/CreatedTFFromDistributions_FromTree.root","RECREATE");
	}
	//Also draw the 2D histograms!	

        //Define all histograms which need to be fitted!
        const int NrFitHistos = 12;
        int ConsHisto = 1;
        const int NrParamsDblGaus = 6;
        std::cout << " --> Will look at " << NrFitHistos << " different histograms to fit! " << std::endl;
        string HistoInfo[12][1+NrParamsDblGaus] = { "BJet_DiffPhiVsGenE",   "0.0002",  "0.022",  "8000", "0.0002",   "0.06",  "3000",
                                                    "BJet_DiffEVsGenE",         "10",    "-12", "20000",     "13",     "10", "-5000",
				                    "BJet_DiffThetaVsGenE",      "0",   "0.04",  "2000",      "0",  "0.013",  "6000",
				                    "El_DiffPhiVsGenE",          "0",  "0.006",   "600",      "0", "0.0012",  "1500",
				                    "El_DiffEVsGenE",            "0",     "-2",   "600",      "0",    "0.9",  "1500",
				                    "El_DiffThetaVsGenE",        "0",  "0.007",   "700",      "0", "0.0013",  "2500",
				                    "Light_DiffPhiVsGenE",       "0",  "0.022",  "8000", "0.0004",  "0.002",  "3000",
				                    "Light_DiffEVsGenE",         "0",      "8",  "4000",      "0",     "12",  "4000",
				                    "Light_DiffThetaVsGenE",     "0",  "-0.05",  "2000",      "0", "-0.014",  "6000",
				                    "Mu_DiffPhiVsGenInvE",       "0",  "0.008",   "700",      "0", "0.0015",   "600",
				                    "Mu_DiffInvEVsGenInvE",      "0", "0.0003",  "2000",      "0", "0.0006",   "500",
				                    "Mu_DiffThetaVsGenInvE",     "0",  "0.002",   "500",      "0", "0.0004",   "500"};

        //Set the booleans!
        bool useROOTClass = false;
        bool useStartValues = true;
        int histoNrForStartValues = NrFitHistos; //Not needed if useStartValues = false
        bool useStartArray = true;
        bool changeFitRange = true;
 
        ofstream myTFTable, myTransferCard[2], myTF[2];
        myTFTable.open("TFInformation/TransferFunctions_TABLE.txt");
        myTF[0].open("TFInformation/TF_user.dat");
        myTF[1].open("TFInformation/TF_user_etaBins.dat");
	myTransferCard[0].open("TFInformation/transfer_card_user.dat");
	myTransferCard[1].open("TFInformation/transfer_card_user_etaBins.dat");

	for(int ii = 0; ii < 2; ii++){
            myTransferCard[ii]<<"#+-----------------------------------------------------------------------+"<<endl;
	    myTransferCard[ii]<<"#|                         TRANSFER_CARD.DAT                             |"<<endl;
	    myTransferCard[ii]<<"#|                                                                       |"<<endl;
            myTransferCard[ii]<<"#|     Author: Annik Olbrechts (VUB)                                     |"<<endl;
    	    myTransferCard[ii]<<"#|             27 November 2014                                          |"<<endl;
    	    myTransferCard[ii]<<"#+-----------------------------------------------------------------------+"<<endl;
	    myTransferCard[ii]<<"#|     This file is generated automaticly by MADWEIGHT                   |"<<endl;
            myTransferCard[ii]<<"#|     card generation version: 1.0.0                                    |"<<endl;
	    myTransferCard[ii]<<"#+-----------------------------------------------------------------------+"<<endl;
	    myTransferCard[ii]<<"#|                                                                       |"<<endl;
	    myTransferCard[ii]<<"#|    To change the transfer function run ./bin/change_tf.py             |"<<endl;
	    myTransferCard[ii]<<"#|    Current parametrization :dbl_gauss_pt_jet                          |"<<endl;
	    myTransferCard[ii]<<"#|    Contains full double Gaussian for all kinematics and particles     |"<<endl;
            if(ii == 0) myTransferCard[ii]<<"#|    ** Information for all eta-bins **                                 |"<<endl;
            if(ii == 1) myTransferCard[ii]<<"#|    ** Information for the "<<nEtaBins<<" considered eta-bins separately **         |"<<endl;
	    myTransferCard[ii]<<"#+-----------------------------------------------------------------------+"<<endl;	

            myTF[ii]<<"<file>## ##################################################################"<<endl;
            myTF[ii]<<"##                                                                       ##"<<endl;
            myTF[ii]<<"##                          Matrix Element                               ##"<<endl;
            myTF[ii]<<"##                          ==============                               ##"<<endl;
            myTF[ii]<<"##                                                                       ##"<<endl;
            myTF[ii]<<"##		    Generate the transfer functions                     ##"<<endl;
            myTF[ii]<<"##	             -------------------------------                    ##"<<endl;
            myTF[ii]<<"## 	     			                                        ##"<<endl;
            myTF[ii]<<"##				                                        ##"<<endl;
            myTF[ii]<<"##    Author: Annik Olbrechts (VUB)                                      ##"<<endl;
            myTF[ii]<<"##   			                                                ##"<<endl;
            myTF[ii]<<"##    Version:     1.0.0                         		        ##"<<endl;
            myTF[ii]<<"##    Last change: 27/11/14			                        ##"<<endl;
            myTF[ii]<<"##					                                ##"<<endl;
            myTF[ii]<<"###########################################################################"<<endl;
            myTF[ii]<<"###########################################################################"<<endl;
            myTF[ii]<<"##				                                        ##"<<endl;
            myTF[ii]<<"##				                                        ##"<<endl;
            myTF[ii]<<"##    Instructions:			                                ##"<<endl;
            myTF[ii]<<"##								        ##"<<endl;
            myTF[ii]<<"##	- This program  creates transfer functions in THETA/PHI/E       ##"<<endl;
            myTF[ii]<<"##	- Those functions must be defined in f77 standard               ##"<<endl;
            myTF[ii]<<"##	- In addition to each transfer function(START_TF), you MUST give##"<<endl;
            myTF[ii]<<"##	   the typical width associated to your function (START_WIDTH)	##"<<endl;
            myTF[ii]<<"##      - If a transfer functions is not defined here it will by default ##"<<endl;
            myTF[ii]<<"##          - equals to one for neutrino/neutralino                      ##"<<endl;
            myTF[ii]<<"##          - a delta functions for the rest                             ##"<<endl;
            myTF[ii]<<"###########################################################################"<<endl;
            myTF[ii]<<"##                                                                       ##"<<endl;
            myTF[ii]<<"##   Syntax/variables:                                                   ## "<<endl;
            myTF[ii]<<"##                                                                       ##"<<endl;
            myTF[ii]<<"##  - a definition for transfer functions should define the variable tf  ##"<<endl;
            myTF[ii]<<"##    while a definition for the width shoud define the variable width   ## "<<endl;                              
            myTF[ii]<<"##	- You can use all standard f77 functions. (All variables are	##"<<endl;
            myTF[ii]<<"##		in double precision format). 	                        ##"<<endl;
            myTF[ii]<<"##	- The experimental event is  defined by the variable pexp(i)    ##"<<endl;
            myTF[ii]<<"##		i=0->3 (0->E,1->Px,2->Py,3->Pz)	                        ##"<<endl;
            myTF[ii]<<"##	- The partonic event is defined by the variable p(i)	        ##"<<endl;
            myTF[ii]<<"##		i=0->3 (0->E,1->Px,2->Py,3->Pz)			        ##"<<endl;
            myTF[ii]<<"##		sigma can not depend on those variables		        ##"<<endl;
            myTF[ii]<<"##	- You can use 10 local variables			        ##"<<endl;
            myTF[ii]<<"##		(double precision):  prov1,prov2,...,prov10	        ##"<<endl;
            myTF[ii]<<"##	- You can call specific functions on p or pexp:	                ##"<<endl;
            myTF[ii]<<"##		-pt(p)   : transverse momenta 			        ##"<<endl;
            myTF[ii]<<"##		-eta(p)  : pseudo-rapidity			        ##"<<endl;
            myTF[ii]<<"##		-rap(p)  : rapidity					##"<<endl;
            myTF[ii]<<"##		-theta(p): polar angle				        ##"<<endl;
            myTF[ii]<<"##		-phi(p)	 : azimuthal angle				##"<<endl;
            myTF[ii]<<"##	- The whole LHCO information is available.                      ##"<<endl;
            myTF[ii]<<"##              -run_number,trigger                       		##"<<endl;
            myTF[ii]<<"##		-eta_init(N),phi_init(N),pt_init(N)                     ##"<<endl;
            myTF[ii]<<"##              -j_mass(N),ntrk(N),btag(N),had_em(N)     	        ##"<<endl;
            myTF[ii]<<"##		-dummy1(N),dummy2(N)		                        ##"<<endl;
            myTF[ii]<<"##	    N is the LHCO tag(first column)	                        ##"<<endl;
            myTF[ii]<<"##		- current tag is n_lhco  				##"<<endl;
            myTF[ii]<<"##		- tag for missing ET is met_lhco			##"<<endl;
            myTF[ii]<<"##				  					##"<<endl;
            myTF[ii]<<"##	- You can incorporate parameters that will be passed through    ##"<<endl;
            myTF[ii]<<"##	        the transfert_card.dat. Those ones must have the        ##"<<endl;
            myTF[ii]<<"##		following syntax: #1,#2,#3,.. You can restart 		##"<<endl;
            myTF[ii]<<"##		the assignement for each different transfer function	##"<<endl;
            myTF[ii]<<"##	- In addition to each transfer function(tf_), you MUST give	##"<<endl;
            myTF[ii]<<"##		the typical width associated to your function (sigma_)	##"<<endl;
            myTF[ii]<<"##		This is needed for the phase space generator      	##"<<endl;
            myTF[ii]<<"##									##"<<endl;
            myTF[ii]<<"###########################################################################"<<endl;
            myTF[ii]<<"###########################################################################"<<endl;
        }

        float startValues[NrParamsDblGaus];
        for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){
            if(NrFitHistos == 1) iHisto = ConsHisto;

            //Set the correct startValues and fit the distribution
            for(int jj = 0; jj < NrParamsDblGaus; jj++) startValues[jj] = std::stof(HistoInfo[iHisto][1+jj]);
    
            for(int iEtaBin = 0; iEtaBin <= nEtaBins; iEtaBin++)
                tfCreation.CalculateTFFromFile(HistoInfo[iHisto][0], useStartValues, histoNrForStartValues, useROOTClass, useStartArray, startValues, changeFitRange, writeFile, iEtaBin, readFile);
    
            //Set the caption correct:
            string CaptionName, BlockName, PartName, KinVarName, TFName, particles, widthType, VarName;
            // -- 1) which particle
            if(HistoInfo[iHisto][0].find("Light_") == 0)    {PartName = "nonbjet";  BlockName = "TF_nonbjet_";  particles = "u,d,s,c,g"; widthType = "thin";}
            else if(HistoInfo[iHisto][0].find("BJet") == 0) {PartName = "bjet";     BlockName = "TF_bjet_";     particles = "b";         widthType = "thin";}
            else if(HistoInfo[iHisto][0].find("Mu_") == 0)  {PartName = "muon";     BlockName = "TF_muon_";     particles = "mu";        widthType = "thin";}
            else if(HistoInfo[iHisto][0].find("El_") == 0)  {PartName = "electron"; BlockName = "TF_electron_"; particles = "el";        widthType = "thin";}
            // -- 2) which kinematic variable
            if(HistoInfo[iHisto][0].find("DiffE") <= HistoInfo[iHisto][0].size())         {CaptionName = PartName+" energy";            KinVarName += "E";   VarName = "E";}
            else if(HistoInfo[iHisto][0].find("DiffTheta") <= HistoInfo[iHisto][0].size()) {CaptionName = PartName+" polar angle \\theta";            KinVarName += "THETA";VarName = "THETA";}
            else if(HistoInfo[iHisto][0].find("DiffPhi") <= HistoInfo[iHisto][0].size())   {CaptionName = PartName+" azimuthal angle \\phi";          KinVarName += "PHI";  VarName = "PHI";}
            else if(HistoInfo[iHisto][0].find("DiffInvE") <= HistoInfo[iHisto][0].size()) {CaptionName = PartName+" inverse of energy"; KinVarName += "INVE";VarName = "INVE";}
            BlockName = BlockName + KinVarName;
    
            //Write the TF's in a table and in a MadWeight card!:
            myTFTable<< endl;
            myTFTable<<" \n \\begin{table}[h!]" << endl;
            myTFTable<<"\\caption{Parameters of the transfer function for " << CaptionName << "}" << endl;
            myTFTable<<"\\label{tab::" << HistoInfo[iHisto][0] << "}" << endl;
            myTFTable<<"\\centering" << endl;
            myTFTable<<"\\begin{tabular}{c|ccc}" << endl;
            myTFTable<<"\\hline" << endl;
            myTFTable<< "Type      & $a_{i0}$ & $a_{i1}$ ($\\sqrt{E}$) & $a_{i2}$ ($E$)" << "\\\\" << endl;
            myTFTable<<"\\hline" << endl;
    
            for(int ii = 0; ii < 2; ii++){
                if(HistoInfo[iHisto][0].find("DiffPhi") <= HistoInfo[iHisto][0].size() || HistoInfo[iHisto][0].find("DiffInvPhi") <= HistoInfo[iHisto][0].size() ){   //Only write this for the first variable! 

                    myTransferCard[ii]<<"#+--------------------------------------------------------------------------------------+" <<endl;
                    myTransferCard[ii]<<"#|     Parameter for particles: "<<PartName << endl; 
                    myTransferCard[ii]<<"#|      --> Used formula: Double Gaussian fit with parameters depending on momentum" << endl;
                    myTransferCard[ii]<<"#|      --> Dependency defined as: A + B*sqrt(E) + C*E  for width of gaussians "<< endl;
		    myTransferCard[ii]<<"#|      -->                        A + B*E + C*E² + D*E³ + F*E^4"<< endl;
                    myTransferCard[ii]<<"#+--------------------------------------------------------------------------------------+" <<endl;

                    myTF[ii]<<"\n##**********************************************************************##"<<endl;
                    myTF[ii]<<"##                "<<TFName<<"                                      "<<endl;
                    myTF[ii]<<"##**********************************************************************##"<<endl;
                    myTF[ii]<<"<block name='"<<PartName<<"'>   #name can be anything"<<endl;
                    myTF[ii]<<"  <info> double gaussian with parameter depending on the energy </info>"<<endl;
                    myTF[ii]<<"  <particles> "<<particles<<" </particles>"<<endl;
                    myTF[ii]<<"  <width_type> "<<widthType<<" </width_type>"<<endl;
                    myTF[ii]<<"  #width type should be thin or large (thin is for energy accurate up to 5-10%)";
                }
                //Output for all variables!
                myTransferCard[ii]<<"BLOCK "<<BlockName << endl;
                
                myTF[ii]<<"\n  <variable name='"<<VarName<<"'>"<<endl;
                myTF[ii]<<"    <tf>";
            }
    
            tfCreation.WriteTF(myTFTable, myTransferCard[0], myTransferCard[1], myTF[0], myTF[1], nEtaBins);  
            if(HistoInfo[iHisto][0].find("DiffTheta") <= HistoInfo[iHisto][0].size() ){ myTF[0] << "\n</block> \n</file>"; myTF[1] << "\n</block> \n</file>";}
    
            myTFTable<<"\\hline" << endl;
            myTFTable<<"\\end{tabular}"<<endl;
            myTFTable<<"\\end{table} \n"<<endl;
        }             
        //Close the root file where all histograms are saved together with the output files!
        readFile->Close();
        writeFile->Close();
        myTF[0].close();
        myTF[1].close();
        myTransferCard[0].close();
        myTransferCard[1].close();
	myTFTable.close();
    
        //Delete the used pointers:
        delete readFile,writeFile;
    }//End of TF calculation when ROOT file is used!
  
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "           hasn't crashed yet ;-)           " << endl;
    cout << "********************************************" << endl;
  
    return 0;
}
