#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdlib>

/////////////////////////////////////////////////////////////
// Specify whether the stacked canvasses have to be stored //
bool storeSplittedCanvas = false; 
std::string SplittedDir = "Events_RecoTest/RecoFirstRun_50000Evts_DblGausTF/SplittedCanvasses"; 
/////////////////////////////////////////////////////////////

TFile *inFile = new TFile("Events_RecoTest/RecoFirstRun_50000Evts_DblGausTF/FitDistributions_RecoFirstRun_50000Evts_DblGausTF_50000Evts.root","READ"); 
TFile *outputFile = new TFile("Events_RecoTest/RecoFirstRun_50000Evts_DblGausTF/FitOptimizations_RecoFirstRun_50000Evts_DblGausTF_50000Evts.root","RECREATE"); 

int NrEvts = 10; 
const int xBin = 9; 
const int xFitBin = xBin*75; 
float xLow = -0.225; 
float xHigh = 0.225; 

void PaintOverflow(TH1F *h, TFile *FileToWrite, std::string dirName){     // This function draws the histogram h with an extra bin for overflows
  Int_t nx    = h->GetNbinsX()+1;
  Double_t x1 = h->GetBinLowEdge(1), bw = h->GetBinWidth(nx), x2 = h->GetBinLowEdge(nx)+bw;

  //Define a temporary histogram having an extra bin for overflows
  char newTitle[100], newName[100];
  strcpy(newTitle,h->GetTitle()); strcat(newTitle," (under- and overflow added)" );
  strcpy(newName,h->GetName());  strcat(newName, "_Flow");
  TH1F *h_tmp = new TH1F(newName, newTitle, nx, x1, x2);

  // Fill the new histogram including the extra bin for overflows
  for (Int_t i=1; i<=nx; i++)
    h_tmp->Fill(h_tmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows
  h_tmp->Fill(x1-1, h->GetBinContent(0));

  // Restore the number of entries
  h_tmp->SetEntries(h->GetEntries());

  //Set the correct path to save the file
  const char* path = 0;
  if(dirName != ""){
    TDirectory *dir = FileToWrite->GetDirectory(dirName.c_str());
    if (!dir)
      dir = FileToWrite->mkdir(dirName.c_str());
    path = dir->GetName();
  }

  FileToWrite->cd(path);
  h_tmp->Write();
  FileToWrite->cd();  //Reset to general directory! 
}                  

void getIndividualDirObjects(TDirectory *dir){
  outputFile->cd();

  std::stringstream ssNkeys; ssNkeys << dir->GetNkeys(); std::string sNkeys = ssNkeys.str();
  //if(!(dir->GetNkeys() != NrEvts) ){
  if( dir->GetNkeys() > 0 ){
    
    TList *list_dir = dir->GetListOfKeys();
    TIter next(list_dir); 
    TObject* object_dir = 0;
    
    //Store all the histograms prior to the application of cuts in the same histogram
    TDirectory* OrigDir = outputFile->GetDirectory("OriginalDistributions");
    if(!OrigDir)
      OrigDir = outputFile->mkdir("OriginalDistributions");

    //Convert name of directory (const char*) to a string:
    std::string dirName(dir->GetName(), 0, 100);

    TH1F *h_LLSum = 0;
    TH1F *h_FitSum = 0, *h_LLSumAll = 0;
    TH1F *h_ChiSq = 0;
   //Identify which directory is being studied --> Containing TF1 or TH1F & type of normalisation!
    std::string NormType = "", NormTypeDir = "", FitNr = "", NrUsedPoints = "";
    bool lookAtFits = false;
    if(dirName.find("Acc") != std::string::npos)     {NormType = "acceptance"; NormTypeDir = "_Acc";}
    else if(dirName.find("XS") != std::string::npos) {NormType = "XS";         NormTypeDir = "_XS"; }
    else                                             {NormType = "no";         NormTypeDir = "";    }

    if(dirName.find("FirstPolynomial") != std::string::npos){                   FitNr = "FirstPol";           NrUsedPoints = "all";    }
    else if(dirName.find("SecondPolynomial") != std::string::npos){             FitNr = "SecondPol";          NrUsedPoints = "reduced";}

    if (dirName.find("LL") != std::string::npos)
      h_LLSum = new TH1F(("LL"+NormTypeDir+"_Summed").c_str(), ("Distribution of "+dirName+" after summing over all "+sNkeys.c_str()+" events").c_str(), xBin, xLow, xHigh);
    else{
      lookAtFits = true;
      h_FitSum = new TH1F((FitNr+""+NormTypeDir+"_Summed").c_str(),  "title",xFitBin,xLow,xHigh);
      h_ChiSq = new TH1F( ("ChiSq"+NormTypeDir+"_"+FitNr).c_str(),   "title",200,0,    0.0005);
      h_LLSumAll = new TH1F( ("Hist_"+FitNr+""+NormTypeDir+"_All").c_str(),"title",xBin,xLow,xHigh);
    }
   
    //Initialize counters and histograms:
    int allEvts = 0;
    double VarBinValue[xFitBin] = {0.};
 
    std::string fitName = "";
    while ((object_dir = next())) {         //Keep it general to TObject because both TF1 and TH1F are stored in this ROOT file!
      if(lookAtFits == false){  //Looking at TH1F objects!
        TH1F *h_iLL = (TH1F*) dir->Get(object_dir->GetName());
        h_LLSum->Add(h_iLL);
      }
      else if(lookAtFits == true){ //Looking at TF1 objects!
        TF1 *fitFunc = (TF1*) dir->Get(object_dir->GetName());
        h_FitSum->Add(fitFunc);
        allEvts++;
        h_ChiSq->Fill(fitFunc->GetChisquare());

        //Access the event number of the current fit (first convert name of object (const char*) to a string):
        std::string objectName(object_dir->GetName(), 0, 100);
        std::string EvtNrFit = objectName.substr(objectName.find("Evt"));
        //Get the original LL distribution for the corresponding event number
        TH1F *h_OrigLL = (TH1F*) inFile->GetDirectory(("OriginalLL"+NormTypeDir).c_str())->Get(("LnLik"+NormTypeDir+"_"+EvtNrFit).c_str());
        h_LLSumAll->Add(h_OrigLL);                      //Not possible to access the TH1F object at this point ??

        for(int iBin = 0; iBin < xFitBin; iBin++){
          double VarValue = xLow+(xHigh-xLow)/(2*xFitBin)+((xHigh-xLow)/xFitBin)*iBin;
          for(int iPar = 0; iPar < fitFunc->GetNpar(); iPar++){
            VarBinValue[iBin] += fitFunc->GetParameter(iPar)*pow(VarValue,iPar);
          }
        }
      }
    }

    if(lookAtFits == true){
      //Convert the event counters into strings!
      std::stringstream ssAllEvts; ssAllEvts << allEvts; std::string sAllEvts = ssAllEvts.str();

      //Set the counter information in the titles:
      h_FitSum->SetTitle(("Distribution of "+FitNr+" after summing over "+NrUsedPoints+" points ("+sAllEvts+" evts)").c_str());
      h_ChiSq->SetTitle(("ChiSquared distribution of "+FitNr+" ("+NormType+" norm -- "+sAllEvts+" evts)").c_str());
      h_LLSumAll->SetTitle("LL distribution");

      //TDirectory *dirFit = outputFile->GetDirectory((FitNr+"FitDistributions").c_str());
      //if (!dirFit)
      //  dirFit = outputFile->mkdir((FitNr+"FitDistributions").c_str());
      //dirFit->cd();

      for(int iBin = 0; iBin < xFitBin; iBin++){      
        if(abs((h_FitSum->GetBinContent(iBin+1)) - VarBinValue[iBin]) > 1){             //Why are they all different? (but printout gives same values ..) -->Something went wrong here ...
          h_FitSum->SetBinContent(iBin+1, VarBinValue[iBin]);
        }
      }

      outputFile->cd();

      //Now write away the histograms
      OrigDir->cd(); 
      h_FitSum->Write();
      outputFile->cd();
      PaintOverflow(h_ChiSq,outputFile,"ChiSqDistributions");

      //Also save the fit-distributions in the SplittedCanvasses directory!	
      if(storeSplittedCanvas == true){
        TCanvas* FitFuncCanv = new TCanvas("FitFunc","FitFunc"); FitFuncCanv->cd(); h_FitSum->Draw(); FitFuncCanv->Print((SplittedDir+"/"+FitNr+""+NormTypeDir+".pdf").c_str()); delete FitFuncCanv;
      }
    }
    else{ 
      OrigDir->cd(); h_LLSum->Write(); outputFile->cd();
      if(storeSplittedCanvas == true){
        TCanvas* canv_LLSum = new TCanvas("LLSum","LLSum"); canv_LLSum->cd(); h_LLSum->Draw(); canv_LLSum->Print((SplittedDir+"/TotalLnLik"+NormTypeDir+".pdf").c_str()); delete canv_LLSum;
      }
    }
  }
  else{
    std::cout << "        ### Skipped the empty directory " << dir->GetName() << " ### " << std::endl;
  }
}

void PerformFitOptimization(){

  std::cout << " Looking at TFile " << inFile->GetName() << " (which has " << inFile->GetNkeys() << " directories) " << std::endl;
  TList *list_file = inFile->GetListOfKeys();
  TIter next(list_file);
  TObject* object_file = 0;
  TDirectory* dir_file = 0;
  while( (object_file = next() ) ){
    std::cout << "  ** Name of directory : " << object_file->GetName() << std::endl;

    //Convert name of directory (const char*) to a string:
    std::string objectName(object_file->GetName(), 0, 100);
    if (objectName != "SplitCanvasses" && objectName != "FitResults"){                  //Do not add the histograms in the SplitCanvasses and FitResults directories for the moment!
      dir_file = inFile->GetDirectory(object_file->GetName());
      getIndividualDirObjects(dir_file);
    }
  }

  inFile->Close();
  outputFile->Close();
}
