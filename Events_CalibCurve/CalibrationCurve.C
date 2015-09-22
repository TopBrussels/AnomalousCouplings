//--------------------------------------------------------------------------------------------------------------------//
//  Script which will compare the measured minimum position with the variable used for creating the MadGraph sample   //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` CalibrationCurve.C -o CalibrationCurve    //
//               ./CalibrationCurve                                                                                   //
//--------------------------------------------------------------------------------------------------------------------//

#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[]){

  // ------------------------------------------------------- //
  // --- Use the input variables to specify some choices --- //
  // ---     * Are cuts applied (y/n)                    --- //
  // ---     * Should cos theta reweighting be applied?  --- //
  // ------------------------------------------------------- //
  bool applyAcc = false, applyCosTheta = false, savePlots = false;
  std::string saveExt = "";
  std::string YesOptions[5] = {"y","yes","Y","Yes","YES"};
  if(argc >= 2){
    if(string(argv[1]) != "n" && string(argv[1]) != "N"){ savePlots = true; saveExt = string(argv[1]);}

    if(argc >= 3 && (string(argv[2]) == "y" || string(argv[2]) == "Y") ) applyAcc = true;

    if(argc >= 4 && (string(argv[3]) == "y" || string(argv[3]) == "Y") && applyAcc) applyCosTheta = true; 
    else if(argc >= 4 && !applyAcc){
      std::cout << "\n ERROR: Cannot apply cos theta reweighting when acceptance is not applied !! " << std::endl;
      std::cout << "    --> Will not apply anything ! \n " << std::endl;
      applyAcc = false; applyCosTheta = false;
    }
  }
  std::cout << " *** Summary of applied corrections :   *** " << std::endl;
  std::cout << "     - Acceptance normalisation : " << applyAcc << std::endl;
  std::cout << "     - Cos theta reweighting    : " << applyCosTheta << " \n" << std::endl;

  std::string CosTheta = "";
  std::string CosThTitle = "no cos #theta applied";
  if(applyCosTheta){ CosTheta = "_CosTheta"; CosThTitle = "cos #theta reweighting applied";}

  std::string Cuts = "NoCuts";
  std::string CutsTitle = "no cuts applied";
  if(applyAcc){ Cuts = "CutsAlsoOnMET"; CutsTitle = "event selection applied";}

  // ----------------------------------------------- //
  // --- Now loop over the different directories --- //
  // ----------------------------------------------- //
  double xAxisRange[11] = {-0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3};
  std::string ConsPars[11] = {"Neg03","Neg02","Neg015","Neg01","Neg005","SM","Pos005","Pos01","Pos015","Pos02","Pos03"};
  TGraph *gr_StraightLine = new TGraph(sizeof(xAxisRange)/sizeof(xAxisRange[0]), xAxisRange, xAxisRange);
  gr_StraightLine->SetLineStyle(2); gr_StraightLine->SetLineColor(13);

  //Store the minimum values and their uncertainties in a vector
  double Min_HistSum[11], Err_HistSum[11];
  double Min_FstFit[11], Err_FstFit[11];
  double Min_ScdFit[11], Err_ScdFit[11];
  for(int iDir = 0; iDir < sizeof(ConsPars)/sizeof(ConsPars[0]); iDir++){
    std::string whichDir = "CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample"+ConsPars[iDir]+"_20000Evts_"+Cuts;
    std::string whichFile = "FitDistributions_CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample"+ConsPars[iDir]+"_20000Evts_"+Cuts+""+CosTheta+"_20000Evts.root";

    std::cout << " - Will be looking at coefficient : " << ConsPars[iDir];
   
    //Get the histograms from the root file:
    TFile* file = new TFile((whichDir+'/'+whichFile).c_str(),"READ");
    
    //Now get the relevant histograms from this file:
    TH1F* h_SummedWeights = (TH1F*) file->Get("FitSums/SummedHist_Acc");
    TF1* pol_FstFit = (TF1*) file->Get("FitSums/SummedFit_FirstFit_Acc");
    TF1* pol_ScdFit = (TF1*) file->Get("FitSums/SummedFit_SecondFit_Acc");

    //The TF1 fitted over the TH1F can be accessed with GetFunction:
    TF1* pol_Sum = h_SummedWeights->GetFunction("polFit_Acc_SummedHist");
    std::cout << " with minimum = " << pol_Sum->GetMinimumX() << std::endl;

    Min_HistSum[iDir] = pol_Sum->GetMinimumX(); 
    Min_FstFit[iDir] = pol_FstFit->GetMinimumX();
    Min_ScdFit[iDir] = pol_ScdFit->GetMinimumX();
    Err_HistSum[iDir] = pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, pol_Sum->GetMinimumX(), 0.4) - pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, -0.4, pol_Sum->GetMinimumX());
    Err_FstFit[iDir] = pol_FstFit->GetX(pol_FstFit->GetMinimum()+0.5, pol_FstFit->GetMinimumX(), 0.4) - pol_FstFit->GetX(pol_FstFit->GetMinimum()+0.5, -0.4, pol_FstFit->GetMinimumX());
    Err_ScdFit[iDir] = pol_FstFit->GetX(pol_FstFit->GetMinimum()+0.5, pol_FstFit->GetMinimumX(), 0.4) - pol_FstFit->GetX(pol_FstFit->GetMinimum()+0.5, -0.4, pol_FstFit->GetMinimumX());
  }

  //Now store everything into a TGraphError
  TGraphErrors *gr_HistSum = new TGraphErrors(sizeof(xAxisRange)/sizeof(xAxisRange[0]), xAxisRange, Min_HistSum, 0, Err_HistSum);
  gr_HistSum->SetTitle(("Calibration curve obtained from summed histogram ("+CutsTitle+" -- "+CosThTitle+")").c_str());
  gr_HistSum->SetMarkerStyle(1); gr_HistSum->SetMarkerColor(1); gr_HistSum->SetLineColor(1);

  TGraphErrors *gr_FstFit = new TGraphErrors(sizeof(xAxisRange)/sizeof(xAxisRange[0]), xAxisRange, Min_FstFit, 0, Err_FstFit);
  gr_FstFit->SetTitle(("Calibration curve obtained from first polynomial fit ("+CutsTitle+" -- "+CosThTitle+")").c_str());
  gr_FstFit->SetMarkerStyle(1); gr_FstFit->SetMarkerColor(2); gr_FstFit->SetLineColor(2);

  TGraphErrors *gr_ScdFit = new TGraphErrors(sizeof(xAxisRange)/sizeof(xAxisRange[0]), xAxisRange, Min_ScdFit, 0, Err_ScdFit);
  gr_ScdFit->SetTitle(("Calibration curve obtained from second polynomial fit ("+CutsTitle+" -- "+CosThTitle+")").c_str());
  gr_ScdFit->SetMarkerStyle(1); gr_ScdFit->SetMarkerColor(4); gr_ScdFit->SetLineColor(4);

  //Fit the calibration curve with a straight line (slope should be as close as possible to 1)
  TF1* polFit_Line = new TF1("polFitLine","pol1",-0.4,0.4); 

  //Store the individual calibration curves in a histogram, and also the three together for comparison!
  TFile *outputFile = new TFile(("Results_CalibrationCurve_"+Cuts+""+CosTheta+".root").c_str(), "RECREATE");
  outputFile->cd();
  TCanvas *canv = new TCanvas("canv","canv");
  TLatex t;
  t.SetTextSize(0.025);    

  //Now draw the individual TGraphs, each fitted with a straight line!
  canv->cd();
  gr_HistSum->Draw("AP");
  polFit_Line->SetLineColor(1);  
  gr_HistSum->Fit(polFit_Line,"Q","",-0.25,0.25);
  stringstream slope_HistSum, slopeErr_HistSum, inter_HistSum, intErr_HistSum, chiSq_HistSum, ndf_HistSum;
  slope_HistSum << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(1); slopeErr_HistSum << std::fixed << std::setprecision(4) << polFit_Line->GetParError(1);
  inter_HistSum << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(0); intErr_HistSum << std::fixed << std::setprecision(4) << polFit_Line->GetParError(0);
  chiSq_HistSum << std::fixed << std::setprecision(3) << polFit_Line->GetChisquare();  ndf_HistSum << polFit_Line->GetNDF();
  t.DrawLatex(0.,-0.25,"Curve fitted with straight line:");
  t.DrawLatex(0.,-0.3,(" - slope = "+slope_HistSum.str()+" #pm "+slopeErr_HistSum.str()).c_str());
  t.DrawLatex(0.,-0.35,(" - intercept = "+inter_HistSum.str()+" #pm "+intErr_HistSum.str()).c_str());
  t.DrawLatex(0.,-0.4,(" - #chi^{2}/ndf = "+chiSq_HistSum.str()+"/"+ndf_HistSum.str()).c_str());
  if(savePlots) canv->SaveAs(("CalibrationCurve_HistSum_"+Cuts+""+CosTheta+"."+saveExt).c_str());
  canv->Write(("CalibCurve_HistSum_"+Cuts+""+CosTheta).c_str());

  canv->cd();
  gr_FstFit->Draw("AP");
  polFit_Line->SetLineColor(2);
  gr_FstFit->Fit(polFit_Line,"Q","",-0.25,0.25); 
  stringstream slope_FstFit, slopeErr_FstFit, inter_FstFit, intErr_FstFit, chiSq_FstFit, ndf_FstFit;
  slope_FstFit << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(1); slopeErr_FstFit << std::fixed << std::setprecision(4) << polFit_Line->GetParError(1);
  inter_FstFit << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(0); intErr_FstFit << std::fixed << std::setprecision(4) << polFit_Line->GetParError(0);
  chiSq_FstFit << std::fixed << std::setprecision(3) << polFit_Line->GetChisquare();  ndf_FstFit << polFit_Line->GetNDF();
  t.DrawLatex(0.,-0.25,"Curve fitted with straight line:");
  t.DrawLatex(0.,-0.3,(" - slope = "+slope_FstFit.str()+" #pm "+slopeErr_FstFit.str()).c_str());
  t.DrawLatex(0.,-0.35,(" - intercept = "+inter_FstFit.str()+" #pm "+intErr_FstFit.str()).c_str());
  t.DrawLatex(0.,-0.4,(" - #chi^{2}/ndf = "+chiSq_FstFit.str()+"/"+ndf_FstFit.str()).c_str());
  if(savePlots) canv->SaveAs(("CalibrationCurve_FirstPolFit_"+Cuts+""+CosTheta+"."+saveExt).c_str());
  canv->Write(("CalibCurve_FirstPolFit_"+Cuts+""+CosTheta).c_str());

  canv->cd();
  gr_ScdFit->Draw("AP");
  polFit_Line->SetLineColor(4);
  gr_ScdFit->Fit(polFit_Line,"Q","",-0.25,0.25);
  stringstream slope_ScdFit, slopeErr_ScdFit, inter_ScdFit, intErr_ScdFit, chiSq_ScdFit, ndf_ScdFit;
  slope_ScdFit << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(1); slopeErr_ScdFit << std::fixed << std::setprecision(4) << polFit_Line->GetParError(1);
  inter_ScdFit << std::fixed << std::setprecision(4) << polFit_Line->GetParameter(0); intErr_ScdFit << std::fixed << std::setprecision(4) << polFit_Line->GetParError(0);
  chiSq_ScdFit << std::fixed << std::setprecision(3) << polFit_Line->GetChisquare();  ndf_ScdFit << polFit_Line->GetNDF();
  t.DrawLatex(0.,-0.25,"Curve fitted with straight line:");
  t.DrawLatex(0.,-0.3,(" - slope = "+slope_ScdFit.str()+" #pm "+slopeErr_ScdFit.str()).c_str());
  t.DrawLatex(0.,-0.35,(" - intercept = "+inter_ScdFit.str()+" #pm "+intErr_ScdFit.str()).c_str());
  t.DrawLatex(0.,-0.4,(" - #chi^{2}/ndf = "+chiSq_ScdFit.str()+"/"+ndf_ScdFit.str()).c_str());
  if(savePlots) canv->SaveAs(("CalibrationCurve_SecondPolFit_"+Cuts+""+CosTheta+"."+saveExt).c_str());
  canv->Write(("CalibCurve_SecondPolFit_"+Cuts+""+CosTheta).c_str());

  canv->cd();
  gr_HistSum->Draw("AP");
  gr_FstFit->Draw("P");
  gr_ScdFit->Draw("P");
  gr_StraightLine->Draw("C");
  TLegend *leg = new TLegend(0.1, 0.7, 0.5, 0.9);
  leg->AddEntry(gr_HistSum,("Summed hist (slope = "+slope_HistSum.str()+" #pm "+slopeErr_HistSum.str()+")").c_str(),"lep");
  leg->AddEntry(gr_FstFit,("First pol fit (slope = "+slope_FstFit.str()+" #pm "+slopeErr_FstFit.str()+")").c_str(),"lep");
  leg->AddEntry(gr_ScdFit,("Second pol fit (slope = "+slope_ScdFit.str()+" #pm "+slopeErr_ScdFit.str()+")").c_str(),"lep");
  leg->AddEntry(gr_StraightLine,"Expected curve","l");
  leg->Draw();
  if(savePlots) canv->SaveAs(("ComparisonDifferentFits_"+Cuts+""+CosTheta+"."+saveExt).c_str());
  canv->Write(("DifferentFits_"+Cuts+""+CosTheta).c_str());

  outputFile->Close();
  return 0;
}
