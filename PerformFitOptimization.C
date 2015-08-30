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

TFile *inFile = new TFile("Events_RecoTest/RecoFirstRun_50000Evts_DblGausTF/FitDistributions_RecoFirstRun_50000Evts_DblGausTF_25000Evts.root","READ"); 
TFile *outputFile = new TFile("Events_RecoTest/RecoFirstRun_50000Evts_DblGausTF/FitOptimizations_RecoFirstRun_50000Evts_DblGausTF_25000Evts.root","RECREATE"); 

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

    //Specify which chi-sq cut values should be considered:
    float ChiSqCutsFstPol[4] = {0.0002, 0.001, 0.0005, 0.005}; 
    float ChiSqCutsScdPol[4] = {0.0002, 0.001, 0.0005, 0.005}; 
    const int NrChiSqCutsFstPol = sizeof(ChiSqCutsFstPol)/sizeof(float);
    const int NrChiSqCutsScdPol = sizeof(ChiSqCutsScdPol)/sizeof(float);
    const int NrChiSqCuts = NrChiSqCutsFstPol;
    //if(NrChiSqCutsFstPol > NrChiSqCutsScdPol) NrChiSqCuts = (const int) NrChiSqCutsScdPol;  //Take the lowest number as the number of cuts!

    TH1F *h_LLSum = 0;
    TH1F *h_FitSum = 0, *h_LLSumAll = 0, *h_FitSumPosSlope = 0;
    TH1F *h_ChiSq = 0, *h_FitSumSmallChiSq[NrChiSqCuts] = {0}, *h_LLSumSmallChiSq[NrChiSqCuts] = {0}, *h_SlopeFit = 0, *h_FitSumSmallChiSqPosSlope[NrChiSqCuts] = {0};
    //Identify which directory is being studied --> Containing TF1 or TH1F & type of normalisation!
    std::string NormType = "", NormTypeDir = "", FitNr = "", NrUsedPoints = "";
    bool lookAtFits = false;
    if(dirName.find("Acc") != std::string::npos)     {NormType = "Acc"; NormTypeDir = "_Acc";}
    else if(dirName.find("XS") != std::string::npos) {NormType = "XS";  NormTypeDir = "_XS"; }
    else                                             {NormType = "";    NormTypeDir = "";    }

    if(dirName.find("FirstPolynomial") != std::string::npos){                   FitNr = "FirstPol";           NrUsedPoints = "all";    }
    else if(dirName.find("SecondPolynomial") != std::string::npos){             FitNr = "SecondPol";          NrUsedPoints = "reduced";}
    else if(dirName.find("FirstFit_Acc_RejectedEvents") != std::string::npos){  FitNr = "FirstPol_RejEvts";   NrUsedPoints = "rejected";}
    else if(dirName.find("FirstFit_Acc_FitDevCutLoose") != std::string::npos){  FitNr = "FirstPol_FitDevCutLoose"; NrUsedPoints = "fitdevCutLoose";}
    else if(dirName.find("FirstFit_Acc_FitDevCutMedium") != std::string::npos){ FitNr = "FirstPol_FitDevCutMedium"; NrUsedPoints = "fitdevCutMedium";}
    else if(dirName.find("FirstFit_Acc_FitDevCutTight") != std::string::npos){  FitNr = "FirstPol_FitDevCutTight"; NrUsedPoints = "fitdevCutTight";}
    else if(dirName.find("FirstFit_Acc_MTop172") != std::string::npos){         FitNr = "FirstPol_MTop172"; NrUsedPoints = "mTop172";}
    else if(dirName.find("FirstFit_Acc_MTop174") != std::string::npos){         FitNr = "FirstPol_MTop174"; NrUsedPoints = "mTop174";}
    else if(dirName.find("FirstFit_Acc_PosScdDerAndSlope") != std::string::npos){ FitNr = "FirstPol_PosScdDerAndSlope"; NrUsedPoints = "all";}
    else if(dirName.find("FirstFit_Acc_PosScdDerAndBothSlope") != std::string::npos){ FitNr = "FirstPol_PosScdDerAndBothSlope"; NrUsedPoints = "all";}

    if (dirName.find("LL") != std::string::npos)
      h_LLSum = new TH1F(("LL"+NormType+"_Summed").c_str(), ("Distribution of "+dirName+" after summing over all "+sNkeys.c_str()+" events").c_str(), xBin, xLow, xHigh);
    else{
      lookAtFits = true;
      h_FitSum = new TH1F((FitNr+""+NormType+"_Summed").c_str(),  "title",xFitBin,xLow,xHigh);
      h_SlopeFit = new TH1F(("Slope_"+FitNr+""+NormType).c_str(), "title",200,-5.0, 5.0  );
      h_ChiSq = new TH1F( ("ChiSq"+NormType+"_"+FitNr).c_str(),   "title",200,0,    0.0005);
      h_LLSumAll = new TH1F( ("Hist_"+FitNr+""+NormType+"_All").c_str(),"title",xBin,xLow,xHigh);
      h_FitSumPosSlope = new TH1F((FitNr+""+NormType+"_PosSlopeSummed").c_str(),"title",xFitBin,xLow,xHigh);
    }
   
    //Initialize counters and histograms:
    int allEvts = 0;
    int evtsSmallChiSq[NrChiSqCuts] = {0}, evtsSmallChiSqPosSlope[NrChiSqCuts] = {0};
    double VarBinValue[xFitBin] = {0.}, VarBinValuePosSlope[xFitBin] = {0.};
    double VarBinValueSmallChiSq[NrChiSqCuts][xFitBin] = {{0.}};
    //double VarBinValueSmallChiSqPosSlope[NrChiSqCuts][xFitBin] = {{0.}};
    for(int ii = 0; ii < NrChiSqCuts; ii++){
      std::stringstream ssii; ssii << ii; std::string sii = ssii.str();
      h_FitSumSmallChiSq[ii]         = new TH1F((FitNr+""+NormType+"_SmallChiSq"+sii).c_str(),        "title",xFitBin,xLow,xHigh);
      h_LLSumSmallChiSq[ii]          = new TH1F(("Hist_"+FitNr+""+NormType+"_SmallChiSq"+sii).c_str(),"title",xBin,xLow,xHigh);
      //h_FitSumSmallChiSqPosSlope[ii] = new TH1F((FitNr+""+NormType+"_SmallChiSqPosSlope"+sii).c_str(),"title",xFitBin,xLow,xHigh);
    }
 
    std::string fitName = "";
    while ((object_dir = next())) {         //Keep it general to TObject because both TF1 and TH1F are stored in this ROOT file!
      if(lookAtFits == false){  //Looking at TH1F objects!
        TH1F *h_iLL = (TH1F*) dir->Get(object_dir->GetName());
        h_LLSum->Add(h_iLL);
      }
      else if(lookAtFits == true){ //Looking at TF1 objects!
        TF1 *fitFunc = (TF1*) dir->Get(object_dir->GetName());
        h_FitSum->Add(fitFunc);
        if(fitFunc->GetParameter(2) > 0) h_FitSumPosSlope->Add(fitFunc);
        h_SlopeFit->Fill(fitFunc->GetParameter(2));
        allEvts++;
        h_ChiSq->Fill(fitFunc->GetChisquare());

        //Access the event number of the current fit (first convert name of object (const char*) to a string):
        std::string objectName(object_dir->GetName(), 0, 100);
        std::string EvtNrFit = objectName.substr(objectName.find("Evt"));
        //Get the original LL distribution for the corresponding event number
        TH1F *h_OrigLL = (TH1F*) inFile->GetDirectory(("OriginalLL"+NormTypeDir).c_str())->Get(("LnLik"+NormType+"_"+EvtNrFit).c_str());
        h_LLSumAll->Add(h_OrigLL);                      //Not possible to access the TH1F object at this point ??

        for(int iBin = 0; iBin < xFitBin; iBin++){
          double VarValue = xLow+(xHigh-xLow)/(2*xFitBin)+((xHigh-xLow)/xFitBin)*iBin;
          for(int iPar = 0; iPar < fitFunc->GetNpar(); iPar++){
            VarBinValue[iBin] += fitFunc->GetParameter(iPar)*pow(VarValue,iPar);
            if(fitFunc->GetParameter(2) > 0) VarBinValuePosSlope[iBin] += fitFunc->GetParameter(iPar)*pow(VarValue,iPar);
          }
        }

        //Apply the chi-sq cuts for each of the values requested:
        for(int iChiSq = 0; iChiSq < NrChiSqCuts; iChiSq++){
          if( (FitNr == "FirstPol" && fitFunc->GetChisquare() < ChiSqCutsFstPol[iChiSq]) || (FitNr == "SecondPol" && fitFunc->GetChisquare() < ChiSqCutsScdPol[iChiSq]) ){
            h_FitSumSmallChiSq[iChiSq]->Add(fitFunc); //h_LLSumSmallChiSq[iChiSq]->Add(h_OrigLL); 
            evtsSmallChiSq[iChiSq]++;
            for(int iBin = 0; iBin < xFitBin; iBin++){
              double VarValue = xLow+(xHigh-xLow)/(2*xFitBin)+((xHigh-xLow)/xFitBin)*iBin;
              for(int iPar = 0; iPar < fitFunc->GetNpar(); iPar++)
                VarBinValueSmallChiSq[iChiSq][iBin] += fitFunc->GetParameter(iPar)*pow(VarValue,iPar);
            }

            /*if(fitFunc->GetParameter(2) > 0){                               //Does this requirement make any sense when a pol4 fit is applied ???
              h_FitSumSmallChiSqPosSlope[iChiSq]->Add(fitFunc);
              evtsSmallChiSqPosSlope[iChiSq]++;
              for(int iBin = 0; iBin < xFitBin; iBin++){
                double VarValue = xLow+(xHigh-xLow)/(2*xFitBin)+((xHigh-xLow)/xFitBin)*iBin;
                for(int iPar = 0; iPar < fitFunc->GetNpar(); iPar++)
                  VarBinValueSmallChiSqPosSlope[iChiSq][iBin] += fitFunc->GetParameter(iPar)*pow(VarValue,iPar);
              }
            }*/
          }
        }
      }
    }

    if(lookAtFits == true){
      //Convert the event counters into strings!
      std::stringstream ssAllEvts; ssAllEvts << allEvts; std::string sAllEvts = ssAllEvts.str();

      //Set the counter information in the titles:
      h_FitSum->SetTitle(("Distribution of "+FitNr+" after summing over "+NrUsedPoints+" points ("+sAllEvts+" evts)").c_str());
      h_SlopeFit->SetTitle(("Distribution of slope of "+FitNr+" ("+sAllEvts+" evts)").c_str());
      h_ChiSq->SetTitle(("ChiSquared distribution of "+FitNr+" (norm = "+NormType+" -- "+sAllEvts+" evts)").c_str());
      h_LLSumAll->SetTitle("LL distribution");

      TDirectory *dir = outputFile->GetDirectory((FitNr+"FitDistributions").c_str());
      if (!dir)
        dir = outputFile->mkdir((FitNr+"FitDistributions").c_str());
      dir->cd();

      for(int iBin = 0; iBin < xFitBin; iBin++){      
        if(abs((h_FitSum->GetBinContent(iBin+1)) - VarBinValue[iBin]) > 1){             //Why are they all different? (but printout gives same values ..) -->Something went wrong here ...
          //std::cout << iBin+1 << ") Difference between bincontent and calculated value (FitSum) is : " << h_FitSum->GetBinContent(iBin+1) << " <--> " << VarBinValue[iBin] << std::endl;
          h_FitSum->SetBinContent(iBin+1, VarBinValue[iBin]);
        }
        if(abs((h_FitSumPosSlope->GetBinContent(iBin+1)) - VarBinValuePosSlope[iBin]) > 1)
          h_FitSumPosSlope->SetBinContent(iBin+1, VarBinValuePosSlope[iBin]);
      }

      for(int ii = 0; ii < NrChiSqCuts; ii++){
        for(int iBin = 0; iBin < xFitBin; iBin++){
          if(abs((h_FitSumSmallChiSq[ii]->GetBinContent(iBin+1)) - VarBinValueSmallChiSq[ii][iBin]) > 1){  
            //cout << iBin+1 << " - " << ii << ") BinContent != VarBinValue (SmallChiSq) : " << h_FitSumSmallChiSq[ii]->GetBinContent(iBin+1) << " vs " << VarBinValueSmallChiSq[ii][iBin] << endl;
            h_FitSumSmallChiSq[ii]->SetBinContent(iBin+1, VarBinValueSmallChiSq[ii][iBin]);
          }

          /*if(abs((h_FitSumSmallChiSqPosSlope[ii]->GetBinContent(iBin+1)) - VarBinValueSmallChiSqPosSlope[ii][iBin]) > 1){  
            //cout << iBin+1 << " - " << ii << ") BinContent != VarBinValue (SmallChiSqPosSlope) : " << h_FitSumSmallChiSqPosSlope[ii]->GetBinContent(iBin+1) << " vs " << VarBinValueSmallChiSqPosSlope[ii][iBin] << endl;
            h_FitSumSmallChiSqPosSlope[ii]->SetBinContent(iBin+1, VarBinValueSmallChiSqPosSlope[ii][iBin]);
          }*/
        }

        std::string sChiSqCut = "";
        if( FitNr == "FirstPol")      { std::stringstream ssChiSqCut; ssChiSqCut << ChiSqCutsFstPol[ii]; sChiSqCut = ssChiSqCut.str();}
        else if (FitNr == "SecondPol"){ std::stringstream ssChiSqCut; ssChiSqCut << ChiSqCutsScdPol[ii]; sChiSqCut = ssChiSqCut.str();}
        std::stringstream ssEvtsSmallChiSq;         ssEvtsSmallChiSq << evtsSmallChiSq[ii];                 std::string sEvtsSmallChiSq = ssEvtsSmallChiSq.str();
        //std::stringstream ssEvtsSmallChiSqPosSlope; ssEvtsSmallChiSqPosSlope << evtsSmallChiSqPosSlope[ii]; std::string sEvtsSmallChiSqPosSlope = ssEvtsSmallChiSqPosSlope.str();

        h_FitSumSmallChiSq[ii]->SetTitle(("Distribution of "+FitNr+" after summing over "+NrUsedPoints+" points (#chi^{2} <"+sChiSqCut+" -- "+sEvtsSmallChiSq+"/"+sAllEvts+" evts)").c_str());
        h_LLSumSmallChiSq[ii]->SetTitle(("Distribution of originalLL when applying cuts on "+FitNr+" (#chi^{2} <"+sChiSqCut+" -- "+sEvtsSmallChiSq+"/"+sAllEvts+" evts)").c_str());
        //h_FitSumSmallChiSqPosSlope[ii]->SetTitle(("Distribution of "+FitNr+" after summing over "+NrUsedPoints+" points (#chi^{2} < "+sChiSqCut+", slope > 0 -- "+sEvtsSmallChiSqPosSlope+"/"+sAllEvts+" evts)").c_str());
        h_FitSumSmallChiSq[ii]->Write();
        //h_LLSumSmallChiSq[ii]->Write();
        //h_FitSumSmallChiSqPosSlope[ii]->Write();
      }
      outputFile->cd();

      //Now write away the histograms
      OrigDir->cd(); 
      h_FitSum->Write(); //h_LLSumAll->Write(); 
      h_FitSumPosSlope->Write();
      outputFile->cd();
      PaintOverflow(h_SlopeFit,outputFile,"SlopeDistributions");
      PaintOverflow(h_ChiSq,outputFile,"ChiSqDistributions");

      //Also save the fit-distributions in the SplittedCanvasses directory!	
      if(storeSplittedCanvas == true){
        TCanvas* FitFuncCanv = new TCanvas("FitFunc","FitFunc"); FitFuncCanv->cd(); h_FitSum->Draw(); FitFuncCanv->Print((SplittedDir+"/"+FitNr+""+NormType+".pdf").c_str()); delete FitFuncCanv;
      }
    }
    else{ 
      OrigDir->cd(); h_LLSum->Write(); outputFile->cd();
      if(storeSplittedCanvas == true){
        TCanvas* canv_LLSum = new TCanvas("LLSum","LLSum"); canv_LLSum->cd(); h_LLSum->Draw(); canv_LLSum->Print((SplittedDir+"/TotalLnLik"+NormType+".pdf").c_str()); delete canv_LLSum;
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
  //while( (object_file = next() ) ){
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
