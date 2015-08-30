#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include <sstream>
#include <iostream>

void saveInCanvas(TH1F* hist, TCanvas* canv, int counter){

  //std::cout << " Looking at histogram : " << hist->GetName() << std::endl;
  //TCanvas* canvasHist = new TCanvas("Canv","Canv");
  //canvasHist->cd();
  //hist->Draw();
  //canvasHist->SaveAs(hist->GetName());
  //delete canvasHist;

  double min_Hist = hist->GetMinimum();
  //std::cout << " Found minimum is : " << min_Hist << std::endl;
  int nBins = hist->GetNbinsX();
  for(int iBin = 0; iBin < nBins; iBin++){
    hist->SetBinContent(iBin+1, hist->GetBinContent(iBin+1) - min_Hist);
  }

  canv->cd();
  if(counter == 0){ hist->SetMinimum(-0.5); hist->Draw();}
  else             hist->Draw("same");
  canv->Update();
}

void DrawNormalizedLL(){
  TFile* inFile = new TFile("Events_RgRScan/RgR_RECO_SingleGausTF_10000Evts_FullRange/FitOptimizations_RgR_RECO_SingleGausTF_10000Evts_FullRange_8602Evts.root","READ"); 

  //Canvas where all the different distributions will be shown together
  TCanvas* canvas_FirstFit  = new TCanvas("canvas_FirstFit", "canvas_FirstFit");
  TCanvas* canvas_SecondFit = new TCanvas("canvas_SecondFit","canvas_SecondFit");
 
  //First get the original distribution and save it on the canvas:
  TH1F* h_OriginalFstPol = (TH1F*) inFile->Get("OriginalDistributions/FirstPolAcc_Summed");
  TH1F* h_OriginalScdPol = (TH1F*) inFile->Get("OriginalDistributions/SecondPolAcc_Summed");

  //Set the maximum of the combined canvas to be not higher than 50!
  if(h_OriginalFstPol->GetMaximum() > 20)
    h_OriginalFstPol->SetMaximum(20);
  if(h_OriginalScdPol->GetMaximum() > 20)
    h_OriginalScdPol->SetMaximum(20);

  //Set the axis labels and title correct for the first distribution
  h_OriginalFstPol->GetYaxis()->SetTitle("-ln(L) - (-ln(L)_{min})");
  h_OriginalFstPol->GetXaxis()->SetTitle("RgRcoefficient"); 
  h_OriginalFstPol->SetTitle("-Likelihood shapes for different #chi^{2} cuts (Acc norm -- RgR_RECO_SingleGausTF_10000Evts_FullRange -- 8602 evts)"); 
  h_OriginalFstPol->SetLineColor(6);
  saveInCanvas(h_OriginalFstPol, canvas_FirstFit, 0);
  h_OriginalScdPol->GetYaxis()->SetTitle("-ln(L) - (-ln(L)_{min})");
  h_OriginalScdPol->GetXaxis()->SetTitle("RgRcoefficient"); 
  h_OriginalScdPol->SetTitle("-Likelihood shapes for different #chi^{2} cuts (Acc norm -- RgR_RECO_SingleGausTF_10000Evts_FullRange -- 8602 evts)"); 
  h_OriginalScdPol->SetLineColor(6);
  saveInCanvas(h_OriginalScdPol, canvas_SecondFit, 0);

  //Add a legend
  TLegend* legend_Fst = new TLegend(0.1,0.7,0.48,0.9);
  legend_Fst->AddEntry(h_OriginalScdPol,"Original distribution prior to cuts","l");
  TLegend* legend_Scd = new TLegend(0.1,0.7,0.48,0.9);
  legend_Scd->AddEntry(h_OriginalScdPol,"Original distribution prior to cuts","l");

  //Now get all the ones where the chi-sq constraint has been applied!
  float ChiSqCutsFstPol[4] = {0.0002, 0.001, 0.0005, 0.005}; 
  float ChiSqCutsScdPol[4] = {0.0002, 0.001, 0.0005, 0.005}; 
  const int NrChiSqCuts = sizeof(ChiSqCutsFstPol)/sizeof(float);

  for(int iCut = 0; iCut < NrChiSqCuts; iCut++){
    stringstream ssCut; ssCut << iCut; string sCut = ssCut.str();
    stringstream ssCutValueFstPol; ssCutValueFstPol << ChiSqCutsFstPol[iCut]; std::string sCutValueFstPol = ssCutValueFstPol.str();
    stringstream ssCutValueScdPol; ssCutValueScdPol << ChiSqCutsScdPol[iCut]; std::string sCutValueScdPol = ssCutValueScdPol.str();

    TH1F* h_FstPolSmallChiSq = (TH1F*) inFile->Get(("FirstPolFitDistributions/FirstPolAcc_SmallChiSq"+sCut).c_str());
    TH1F* h_ScdPolSmallChiSq = (TH1F*) inFile->Get(("SecondPolFitDistributions/SecondPolAcc_SmallChiSq"+sCut).c_str());
    h_FstPolSmallChiSq->SetLineColor(iCut+1);
    h_ScdPolSmallChiSq->SetLineColor(iCut+1);
    legend_Fst->AddEntry(h_FstPolSmallChiSq,("Distribution when #chi^{2} < "+sCutValueFstPol).c_str(),"l");
    legend_Scd->AddEntry(h_ScdPolSmallChiSq,("Distribution when #chi^{2} < "+sCutValueScdPol).c_str(),"l");
    saveInCanvas(h_FstPolSmallChiSq, canvas_FirstFit,  iCut+1);
    saveInCanvas(h_ScdPolSmallChiSq, canvas_SecondFit, iCut+1);
  } 

  canvas_FirstFit->cd();
  legend_Fst->Draw();
  canvas_FirstFit->SaveAs("Events_RgRScan/RgR_RECO_SingleGausTF_10000Evts_FullRange/LLFirstFitComparison_ChiSqCutsRgR_RECO_SingleGausTF_10000Evts_FullRange_8602Evts.pdf"); 
  
  canvas_SecondFit->cd();
  legend_Scd->Draw();
  canvas_SecondFit->SaveAs("Events_RgRScan/RgR_RECO_SingleGausTF_10000Evts_FullRange/LLSecondFitComparison_ChiSqCutsRgR_RECO_SingleGausTF_10000Evts_FullRange_8602Evts.pdf"); 
}
