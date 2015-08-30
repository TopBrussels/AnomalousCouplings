#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

Double_t Minimum(TH1F* histo)
{
  Double_t minimum = histo->GetBinContent(1);
  for(int ii = 0; ii < histo->GetNbinsX(); ii++){
    if( histo->GetBinContent(ii+1) < minimum && histo->GetBinContent(ii+1) > 0.0)
      minimum = histo->GetBinContent(ii+1);
  }
  return minimum;
}

void fitExcludeEmptyBins() {

  //Get the histogram from the ROOT file   --> Name gets changed by pythonScript!!
  TFile* InputFile = new TFile(" Events_ContainingMTopScansForRecLevelEvts/TopMass_Gen_ManyEvts/FitDeviation_Gen_MTop.root ","r");
  TFile* OutputFile = new TFile(" Events_ContainingMTopScansForRecLevelEvts/TopMass_Gen_ManyEvts/LimitedFitResult_Gen_MTop_ScdDerInnerCutApplied.root ","RECREATE");
  
  bool GenLevel = true; 
  
  TH1F* LL    = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPosScdDerInner"); 
  TH1F* LLXS  = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLXSPosScdDerInner"); 
  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLAccPosScdDerInner"); 

  //-- Store the zoomed information in a new histogram which will be used for fitting! --//
  int nBins = 5;         //Initialize nrBins 
  float fitStart = 170.5;  //Initialize fitStart 
  float fitEnd = 175.5;    //Initialize fitEnd 
  TH1F* LLFit    = new TH1F("LLFit",   "LL histogram for fit",   nBins,fitStart,fitEnd);
  TH1F* LLXSFit  = new TH1F("LLXSFit", "LLXS histogram for fit", nBins,fitStart,fitEnd);
  TH1F* LLAccFit = new TH1F("LLAccFit","LLAcc histogram for fit",nBins,fitStart,fitEnd);

  float firstBin = (fitEnd-fitStart)/(nBins*2)+fitStart;
  for(int ii = 0; ii < nBins; ii++){
    LLFit->SetBinContent(ii+1,LL->GetBinContent(LL->FindBin(firstBin+ii))); 
    LLXSFit->SetBinContent(ii+1,LLXS->GetBinContent(LLXS->FindBin(firstBin+ii))); 
    if(GenLevel == false) LLAccFit->SetBinContent(ii+1,LLAcc->GetBinContent(LLAcc->FindBin(firstBin+ii))); 
  }
  LLFit->SetMarkerStyle(LL->GetMarkerStyle());       LLFit->SetMarkerColor(LL->GetMarkerColor());     
  LLXSFit->SetMarkerStyle(LLXS->GetMarkerStyle());   LLXSFit->SetMarkerColor(LLXS->GetMarkerColor()); 
  if(GenLevel == false){ LLAccFit->SetMarkerStyle(LLAcc->GetMarkerStyle()); LLAccFit->SetMarkerColor(LLAcc->GetMarkerColor());} 

  //-- Fit the three different histograms! --//
  TF1 *func = new TF1("func","pol2",fitStart,fitEnd);     func->SetName("fit_LL");     LLFit->Fit(func);
  TF1 *funcXS = new TF1("funcXS","pol2",fitStart,fitEnd); funcXS->SetName("fit_LLXS"); LLXSFit->Fit(funcXS);
  if(GenLevel == false){TF1 *funcAcc = new TF1("funcAcc","pol2",fitStart,fitEnd); funcAcc->SetName("fit_LLAcc"); LLAccFit->Fit(funcAcc);}

  //-- Make graph of fit function! --//
  Double_t x[100], y[100], yXS[100], yAcc[100];
  for(int i = 0; i <LL->GetNbinsX(); i++){
    x[i] = LL->GetBinCenter(i+1);
    y[i]    = func->GetParameter(0)   +func->GetParameter(1)*x[i]+   func->GetParameter(2)*x[i]*x[i];
    yXS[i]  = funcXS->GetParameter(0) +funcXS->GetParameter(1)*x[i]+ funcXS->GetParameter(2)*x[i]*x[i];
    if(GenLevel == false) yAcc[i] = funcAcc->GetParameter(0)+funcAcc->GetParameter(1)*x[i]+funcAcc->GetParameter(2)*x[i]*x[i];
  }
  gr    = new TGraph(LL->GetNbinsX(),x,y);       gr->SetLineColor(2);    gr->SetLineStyle(2);
  grXS  = new TGraph(LLXS->GetNbinsX(),x,yXS);   grXS->SetLineColor(2);  grXS->SetLineStyle(2);
  if(GenLevel == false){grAcc = new TGraph(LLAcc->GetNbinsX(),x,yAcc); grAcc->SetLineColor(2); grAcc->SetLineStyle(2);}

  //-- Draw all the histograms in a single canvas! --//
  TCanvas* LLCanvas = new TCanvas("LLCanvas","Distribution of -ln(likelihood) (no normalisation)"); LLCanvas->cd();
  LL->SetMinimum(Minimum(LL)*0.999); LL->SetMarkerColor(1);
  LL->Draw("P"); gr->Draw("C"); func->Draw("same"); LLCanvas->Write();
  func->SetTitle(LL->GetTitle()); func->Write();
  
  TCanvas* LLXSCanvas = new TCanvas("LLXSCanvas","Distribution of -ln(likelihood) (XS normalisation)"); LLXSCanvas->cd();
  LLXS->SetMinimum(Minimum(LLXS)*0.99); LLXS->SetMarkerColor(1);
  LLXS->Draw("P"); grXS->Draw("C"); funcXS->Draw("same"); LLXSCanvas->Write();
  funcXS->SetTitle(LLXS->GetTitle()); funcXS->Write();
  
  if(GenLevel == false){
    TCanvas* LLAccCanvas = new TCanvas("LLAccCanvas","Distribution of -ln(likelihood) (Acc normalisation)"); LLAccCanvas->cd();
    LLAcc->SetMinimum(Minimum(LLAcc)*0.99); LLAcc->SetMarkerColor(1);
    LLAcc->Draw("P"); grAcc->Draw("C"); funcAcc->Draw("same"); LLAccCanvas->Write();
    funcAcc->SetTitle(LLAcc->GetTitle()); funcAcc->Write();
  }
  OutputFile->Close();
}
