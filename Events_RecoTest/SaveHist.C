{ 
  TFile* file = new TFile("Reco_CorrectEvts_DblGausTF_LeptDelta_NonBinned_AllEvts_ISR1_28Oct/FitDistributions_Reco_CorrectEvts_DblGausTF_LeptDelta_NonBinned_AllEvts_ISR1_28Oct_CheckedEvts_OuterBinsExclForFit_117978Evts.root","READ"); 
  TCanvas* canv = new TCanvas("canv","canv"); 
  canv->cd(); 

  TH1F* hist_Sum = (TH1F*)file->Get("FitSums/SummedHist_Acc"); 
  TH1F* SMLik = (TH1F*)file->Get("FitResults/SMLikelihoodValue"); 
  SMLik->Draw();
  canv->SaveAs("SMLikelihoodValue_CorrectEvts.pdf"); 

  hist_Sum->Draw(); 
  canv->SaveAs("SummedHist_CorrectEvts_OuterBinsExclForFit_NoLikCut.pdf"); 

}
