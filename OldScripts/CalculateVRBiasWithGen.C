#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include <iostream>

void createHistograms(TH1F* h_gen, TH1F* h_recoC, TH1F* h_recoW, std::string ChiSqNr, std::string ChiSqTitle, TFile* file){

  h_gen->SetName(("SecondPolAcc_Gen_"+ChiSqTitle).c_str());
  h_recoC->SetName(("SecondPolAcc_RecoCorrect_"+ChiSqTitle).c_str());
  h_recoW->SetName(("SecondPolAcc_RecoWrong_"+ChiSqTitle).c_str());
  
  int nBins = h_gen->GetNbinsX();
  TH1F* h_RecoCGen = new TH1F(("RecoCorrVSGen_"+ChiSqTitle).c_str(),("Correct-reco and gen result bias (min-normalized ln(L) -- #chi^{2} < "+ChiSqNr+")").c_str(),nBins,h_gen->GetXaxis()->GetXmin(),h_gen->GetXaxis()->GetXmax());
  TH1F* h_RecoWGen = new TH1F(("RecoWrongVSGen_"+ChiSqTitle).c_str(),("Wrong-reco and gen bias (min-normalized ln(L) -- #chi^{2} < "+ChiSqNr+")").c_str(),nBins, h_gen->GetXaxis()->GetXmin(), h_gen->GetXaxis()->GetXmax());

  //Normalize the likelihoods by substracting the minimum!
  double min_gen = h_gen->GetMinimum();
  double min_recoC = h_recoC->GetMinimum();
  double min_recoW = h_recoW->GetMinimum();

  for(int iBin = 0; iBin < nBins; iBin++){

    h_gen->SetBinContent(iBin+1, h_gen->GetBinContent(iBin+1)-min_gen);
    h_recoC->SetBinContent(iBin+1, h_recoC->GetBinContent(iBin+1)-min_recoC);
    h_recoW->SetBinContent(iBin+1, h_recoW->GetBinContent(iBin+1)-min_recoW);

    h_RecoCGen->SetBinContent(iBin+1, h_gen->GetBinContent(iBin+1) - h_recoC->GetBinContent(iBin+1));
    h_RecoWGen->SetBinContent(iBin+1, h_gen->GetBinContent(iBin+1) - h_recoW->GetBinContent(iBin+1));
  }

  //Save the original LL distributions:
  h_gen->Write();
  h_recoC->Write();
  h_recoW->Write();

  TCanvas* canv_RecoCGen = new TCanvas("canv_RecoCGen","canv_RecoCGen");
  canv_RecoCGen->cd();
  h_RecoCGen->GetXaxis()->SetTitle("RVR coefficient");
  h_RecoCGen->Draw();
  canv_RecoCGen->SaveAs(("RecoCorrectGenComparison_"+ChiSqTitle+"_RVR.pdf").c_str());
  h_RecoCGen->Write();
  delete canv_RecoCGen;

  TCanvas* canv_RecoWGen = new TCanvas("canv_RecoWGen","RecoWGen");
  canv_RecoWGen->cd();
  h_RecoWGen->GetXaxis()->SetTitle("RVR coefficient");
  h_RecoWGen->Draw();
  canv_RecoWGen->SaveAs(("RecoWrongGenComparison_"+ChiSqTitle+"_RVR.pdf").c_str());
  h_RecoWGen->Write();
  delete canv_RecoWGen;

  delete h_RecoCGen;
  delete h_RecoWGen;
}

void CalculateVRBiasWithGen(){

  TFile* outputFile = new TFile("BiasTest.root","RECREATE");
  
//  TFile* file_MG = new TFile(" .root","READ");
  TFile* file_Gen = new TFile("Events_RVRTests_UsingNewScript/RVR_Gen_SingleGausTFHalfWidth_20000Evts/FitOptimizations_Gen_RVR_10000Evts.root","READ");
  TFile* file_RecoCorrect = new TFile("Events_RVRTests_UsingNewScript/RVR_RecoCorrect_SingleGausTF_10000Evts_WideRange/FitOptimizations_CorrectReco_RVR_9995Evts.root","READ");
  TFile* file_RecoWrong = new TFile("Events_RVRTests_UsingNewScript/RVR_RecoWrong_SingleGausTF_10000Evts_WideRange/FitOptimizations_WrongReco_RVR_9376Evts.root","READ");

  outputFile->cd();

  //Get the original distributions
  TH1F* h_origGen = (TH1F*) file_Gen->Get("OriginalDistributions/SecondPolAcc_Summed");
  TH1F* h_origRecoCorrect = (TH1F*) file_RecoCorrect->Get("OriginalDistributions/SecondPolAcc_Summed");
  TH1F* h_origRecoWrong = (TH1F*) file_RecoWrong->Get("OriginalDistributions/SecondPolAcc_Summed");
/*  h_origGen->Write();
  h_origRecoCorrect->Write();
  h_origRecoWrong->Write();
*/
  createHistograms(h_origGen, h_origRecoCorrect, h_origRecoWrong, "xx", "NoChiSqCut", outputFile);
//  delete h_origGen; delete h_origRecoCorrect; delete h_origRecoWrong;

  //Get the ChiSq0 distributions (chiSq < 0.0002)
  TH1F* h_chiSq0Gen = (TH1F*) file_Gen->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq0");
  TH1F* h_chiSq0RecoCorrect = (TH1F*) file_RecoCorrect->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq0");
  TH1F* h_chiSq0RecoWrong = (TH1F*) file_RecoWrong->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq0");
/*  h_chiSq0Gen->Write();
  h_chiSq0RecoCorrect->Write();
  h_chiSq0RecoWrong->Write();
*/
  createHistograms(h_chiSq0Gen, h_chiSq0RecoCorrect, h_chiSq0RecoWrong, "0.0002", "ChiSqCut0002", outputFile);
//  delete h_chiSq0Gen; delete h_chiSq0RecoCorrect; delete h_chiSq0RecoWrong;

  //Get the ChiSq1 distributions (chiSq < 0.001
  TH1F* h_chiSq1Gen = (TH1F*) file_Gen->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq1");
  TH1F* h_chiSq1RecoCorrect = (TH1F*) file_RecoCorrect->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq1");
  TH1F* h_chiSq1RecoWrong = (TH1F*) file_RecoWrong->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq1");
/*  h_chiSq1Gen->Write();
  h_chiSq1RecoCorrect->Write();
  h_chiSq1RecoWrong->Write();
*/
  createHistograms(h_chiSq1Gen, h_chiSq1RecoCorrect, h_chiSq1RecoWrong, "0.001", "ChiSqCut001", outputFile);
//  delete h_chiSq1Gen; delete h_chiSq1RecoCorrect; delete h_chiSq1RecoWrong;

  //Get the ChiSq2 distributions (chiSq < 0.0005)
  TH1F* h_chiSq2Gen = (TH1F*) file_Gen->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq2");
  TH1F* h_chiSq2RecoCorrect = (TH1F*) file_RecoCorrect->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq2");
  TH1F* h_chiSq2RecoWrong = (TH1F*) file_RecoWrong->Get("SecondPolFitDistributions/SecondPolAcc_SmallChiSq2");
/*  h_chiSq2Gen->Write();
  h_chiSq2RecoCorrect->Write();
  h_chiSq2RecoWrong->Write();
*/
  createHistograms(h_chiSq2Gen, h_chiSq2RecoCorrect, h_chiSq2RecoWrong, "0.0005", "ChiSqCut0005", outputFile);
  delete h_chiSq2Gen; delete h_chiSq2RecoCorrect; delete h_chiSq2RecoWrong;

}
