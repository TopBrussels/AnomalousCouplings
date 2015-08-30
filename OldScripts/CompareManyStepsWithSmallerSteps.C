#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

void CompareManyStepsWithSmallerSteps() {
  TFile* InputFileManySteps = new TFile(" Events/RVR_RecoCorrect_SingleGausTF_10000Evts_ManySteps/Likelihood_CorrectReco_RVR.root ","r");
  TFile* InputFileMGAccSmallerSteps = new TFile(" Events/RVR_RecoCorrect_SingleGausTF_10000Evts_SmallerSteps/Likelihood_CorrectReco_RVR.root ","r");
  TFile* InputFileMAAccSmallerSteps = new TFile(" Events/RVR_RecoCorrect_SingleGausTF_10000Evts_SmallerSteps/Likelihood_CorrectReco_RVR_MadAnalysisAcceptance.root ","r");
  TFile* OutputFile = new TFile("ComparisonManyStepsWithSmallerSteps_CorrectReco_RVR.root","RECREATE");

  TH1F* LLManySteps = (TH1F*) InputFileManySteps->Get("LL"); 
  TH1F* LLAccManySteps = (TH1F*) InputFileManySteps->Get("LL_Acc"); 
  TH1F* LLSmallerSteps = (TH1F*) InputFileMGAccSmallerSteps->Get("LL"); 
  TH1F* LLMGAccSmallerSteps = (TH1F*) InputFileMGAccSmallerSteps->Get("LL_Acc"); 
  TH1F* LLMAAccSmallerSteps = (TH1F*) InputFileMAAccSmallerSteps->Get("LL_Acc"); 

  double xMany[21] = {-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1}; 
  double xSmaller[9] = {-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3};

  double yMany[100], ySmaller[100];
  double yManyAcc[100], ySmallerMGAcc[100], ySmallerMAAcc[100];
  for(int ii = 0; ii < 21; ii++){
    yMany[ii] = LLManySteps->GetBinContent(LLManySteps->FindBin(xMany[ii]));
    yManyAcc[ii] = LLAccManySteps->GetBinContent(LLAccManySteps->FindBin(xMany[ii]));
  }
  for(int jj = 0; jj < 9; jj++){
    ySmaller[jj] = LLSmallerSteps->GetBinContent(LLSmallerSteps->FindBin(xSmaller[jj]));
    ySmallerMGAcc[jj] = LLMGAccSmallerSteps->GetBinContent(LLMGAccSmallerSteps->FindBin(xSmaller[jj]));
    ySmallerMAAcc[jj] = LLMAAccSmallerSteps->GetBinContent(LLMAAccSmallerSteps->FindBin(xSmaller[jj]));
  }

  ManyGraph = new TGraph(21,xMany, yMany);          ManyGraph->SetMarkerColor(4);    ManyGraph->SetLineColor(4);
  SmallerGraph = new TGraph(9, xSmaller, ySmaller); SmallerGraph->SetMarkerColor(1); SmallerGraph->SetLineColor(1); SmallerGraph->SetMarkerStyle(21);
  ManyAccGraph = new TGraph(21,xMany, yManyAcc);          ManyAccGraph->SetMarkerColor(4);    ManyAccGraph->SetLineColor(4);
  SmallerMGAccGraph = new TGraph(9, xSmaller, ySmallerMGAcc); SmallerMGAccGraph->SetMarkerColor(1); SmallerMGAccGraph->SetLineColor(1); SmallerMGAccGraph->SetMarkerStyle(21);
  SmallerMAAccGraph = new TGraph(9, xSmaller, ySmallerMAAcc); SmallerMAAccGraph->SetMarkerColor(2); SmallerMAAccGraph->SetLineColor(2); SmallerMAAccGraph->SetMarkerStyle(22);

  TCanvas *c3 = new TCanvas("ManyVsSmallerComp","Comparison of -ln(L) distribution (no normalisation -- Correct Reco)");
  c3->cd();
  ManyGraph->SetTitle(c3->GetTitle());
  ManyGraph->Draw("AC*");
  SmallerGraph->Draw("CP");
  c3->Write();

  TCanvas *c1 = new TCanvas("ManyVsSmallerAccComp","Comparison of -ln(L) distribution (Acc normalisation -- Correct Reco)");
  c1->cd();
  ManyAccGraph->SetTitle(c1->GetTitle());
  ManyAccGraph->Draw("AC*");
  SmallerMGAccGraph->Draw("CP");
  SmallerMAAccGraph->Draw("CP");
  legend = new TLegend(0.3,0.8,0.45,0.9);
  legend->AddEntry(ManyAccGraph,"-ln(L) for fine range (MG Acceptance)","p");
  legend->AddEntry(SmallerMGAccGraph,"-ln(L) for wide range (MG Acceptance)","p");
  legend->AddEntry(SmallerMAAccGraph,"-ln(L) for wide range (MA Acceptance)","p");
  legend->Draw();
  c1->Write();
}
