{ 
  TFile* file = new TFile("Data_SemiMu_RgR_AllEvents_LatestEvtSel_Nov24/FitDistributions_Data_SemiMu_RgR_AllEvents_LatestEvtSel_Nov24_16943Evts.root","READ"); 
  TCanvas* canv = new TCanvas("canv","canv"); 
  canv->cd(); 

  TH1F* hist_Sum = (TH1F*)file->Get("FitSums/SummedHist_Acc"); 
  TH1F* SMLik = (TH1F*)file->Get("FitResults/SMLikelihoodValue"); 
  SMLik->Draw();
  canv->SaveAs("SMLikelihoodValue_DataEvts.pdf"); 

  hist_Sum->Draw(); 
  canv->SaveAs("SummedHist_DataEvts_NoLikCut.pdf"); 

}
