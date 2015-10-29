//-----------------------------------------------------------------------------------------------------------------------//
//  Script which will compare the different likelihood cut values and help to decide on the ideal cut value!             //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` LikelihoodCutEffect.C -o LikelihoodCutEffect //
//-----------------------------------------------------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"

using namespace std;

int main(int argc, char *argv[]){

  //-------------------------------------------------//
  // Arguments will be given from the python script! //
  //   1) Which directory                            //
  //   2) Which ROOT file                            //
  //   3) Latex file to store information            //
  //   4) Counter of considered files                //
  //   5) Applied cut value                          //
  //-------------------------------------------------//
  std::string whichDir = "", LatexTitle = "";
  std::cout << "Number of arguments given : " << argc << std::endl;
  vector<std::string> RootFiles, CutValues;
  RootFiles.clear(); CutValues.clear();
  if( argc >= 1){
    whichDir = string(argv[1]);
    LatexTitle = string(argv[2]);
    for(int iArg = 3; iArg < argc; iArg++){
      RootFiles.push_back((whichDir+string(argv[iArg])).c_str());
      CutValues.push_back(string(argv[iArg]).substr(string(argv[iArg]).find("LikelihoodCut")+13, string(argv[iArg]).find(".root")-string(argv[iArg]).find("LikelihoodCut")-13));
      std::cout << " Looking at cut-value : " << CutValues[CutValues.size()-1] << std::endl;
    }
  }
  
  ofstream myLatex;
  myLatex.open(LatexTitle.c_str());
  myLatex << "\\documentclass[a4paper]{article} \n\\usepackage{graphicx} \n ";
  myLatex << "\\usepackage{geometry} \n \n";
  //myLatex << "\\usepackage[top=.5in, bottom=1.25in, left=.5in, right=.5in,landscape]{geometry} \n \n";
  myLatex << "\n \\begin{document} \n";
  myLatex << "  \\begin{abstract} \n";
  myLatex << "    File where all the different cut-values on the Likelihood distribution will be compared \\\\ \n";
  myLatex << "    Cut has to be applied in order to remove the observed tail in the likelihood distribution at gR = 0. \n";
  myLatex << "    This behaviour is not present for the Generator-level events and is related with the scanned kinematics region. \n";
  myLatex << "  \\end{abstract} \n \n";

  //Initialize all the histograms and variables
  double CutValue[10] = {0}, minPos[10] = {0}, minPos_E[10] = {0};

  //Store the -ln(L) distribution at gR = 0.0
  myLatex << "  \\begin{figure}[h!t] \n   \\centering \n";
  myLatex << "   \\includegraphics[width = 0.99 \\textwidth]{" << whichDir << "Plots/SMLikelihoodValue.pdf} \n";
  myLatex << "   \\caption{Distribution of the -ln($\\mathcal{L}$) value at gR = 0.0, which shows a clear tail.} \n";
  myLatex << "  \\end{figure} \n \n";

  //Store the overview of the minimum-position shift
  myLatex << "  \\newpage \n";
  myLatex << "  \\begin{figure}[h!t] \n   \\centering \n";
  myLatex << "   \\includegraphics[width = 0.99 \\textwidth]{" << whichDir << "Plots/MinPosition.pdf} \n";
  myLatex << "   \\caption{Effect of the considered -ln($\\mathcal{L}$) cut-value on the retrieved minimum value} \n";
  myLatex << "  \\end{figure} \n \n";

  //Now loop over the different Root files and extract all the information!
  for(int iFile = 0; iFile < RootFiles.size(); iFile++){
    TFile* file = new TFile(RootFiles[iFile].c_str(),"READ");

    //Store the likelihood distribution plot (only need this once since it remains the same for all cuts)
    if(iFile == 0){
      TH1F* h_SMLikelihoodValue = (TH1F*) file->Get("FitResults/SMLikelihoodValue");
      TCanvas* canv_smLik = new TCanvas("smLik","smLik");
      canv_smLik->cd();
      h_SMLikelihoodValue->Draw();
      h_SMLikelihoodValue->SetTitle("");
      canv_smLik->SaveAs((whichDir+"Plots/SMLikelihoodValue.pdf").c_str());
    }

    TH1F* h_SummedWeights = (TH1F*) file->Get("FitSums/SummedHist_Acc");
    TF1* pol_Sum = h_SummedWeights->GetFunction("polFit_Acc_SummedHist");

    //Store this in a map and then order wrt the first element such that the individual histograms are in the correct order!
    minPos[iFile] = pol_Sum->GetMinimumX();
    minPos_E[iFile] = pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, pol_Sum->GetMinimumX(), 0.2) - pol_Sum->GetX(pol_Sum->GetMinimum()+0.5, -0.2, pol_Sum->GetMinimumX());
    CutValue[iFile] = atoi(CutValues[iFile].c_str());

    //Also draw the individual histograms for each considered cut-value!
    // --> Need to lower the y-value such that it becomes a chi-sq
    for(int iBin = 0; iBin <= h_SummedWeights->GetNbinsX(); iBin++){
      h_SummedWeights->SetBinContent(iBin, h_SummedWeights->GetBinContent(iBin)-pol_Sum->GetMinimum());
    }
    pol_Sum->SetParameter(0, (pol_Sum->GetParameter(1)*pol_Sum->GetParameter(1))/(4*pol_Sum->GetParameter(2)));
    TCanvas* canv_indivMin = new TCanvas(("canv_indivMin"+CutValues[iFile]).c_str(),("canv_indivMin"+CutValues[iFile]).c_str());
    canv_indivMin->cd();
    h_SummedWeights->Draw();
    h_SummedWeights->SetTitle(("MadWeight $\\chi^{2}$ distribution for cut-value "+CutValues[iFile]).c_str());
    h_SummedWeights->GetXaxis()->SetTitle("gR coefficient");
    h_SummedWeights->GetYaxis()->SetTitle("$\\chi^{2}$ value");
    pol_Sum->Draw("same");
    canv_indivMin->SaveAs((whichDir+"Plots/IndividualChiSqDistr_Cut"+CutValues[iFile]+".pdf").c_str());

    if(iFile == 0){
      myLatex << "  \\newpage \n";
      myLatex << "  \\paragraph{Individual distributions for the different cut-values} \\\\ \n";
      myLatex << "  \\begin{figure}[h!b] \n   \\centering \n";
    }
    myLatex << "   \\includegraphics[width = 0.49 \\textwidth]{" << whichDir << "Plots/IndividualChiSqDistr_Cut" << CutValues[iFile] << ".pdf} \n";
  }
  myLatex << "   \\caption{$\\chi^{2}$ distributions for the different cut-values applied} \n";
  myLatex << "  \\end{figure} \n";

  //Store all the histograms:
  TGraphErrors* gr_MinPosition = new TGraphErrors( RootFiles.size(), CutValue, minPos, 0, minPos_E);
  TCanvas* canv_minPos = new TCanvas("canv_minPos","canv_minPos");
  canv_minPos->cd();
  gr_MinPosition->Draw("APE");
  gr_MinPosition->GetXaxis()->SetTitle("Considered cut-value on MadWeight -ln(Likelihood) at gR = 0.0");
  gr_MinPosition->GetYaxis()->SetTitle("Minimum value retrieved");
  gr_MinPosition->SetTitle("");
  canv_minPos->SaveAs((whichDir+"Plots/MinPosition.pdf").c_str());

  myLatex << "\n\\end{document} \n";
  myLatex.close();

  return 0;
}
