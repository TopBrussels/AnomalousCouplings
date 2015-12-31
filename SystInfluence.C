//-----------------------------------------------------------------//
//  Use this script to visualize the influence of the systematics  // 
//   --> Compare the shape of the SMLik distribution!              //
//-----------------------------------------------------------------//
// Run with: g++ -m64 -g -I `root-config --incdir` `root-config --libs` SystInfluence.C -o SystInfluence

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "TStyle.h"

//Convert integer to string!
template <typename T>
std::string NumberToString ( T Number ){
    std::stringstream ss; ss << Number;
    return ss.str();
}

int main(int argc, char *argv[]){

  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);

  //Necessary input parameters:
  // 1) Which systematic
  // 2) Which cut-value
  // 3) Which MC samples

  std::string whichSyst = "";       //Possible input variables are : JES, JER, bTag & misTag, Matching & Scaling
  std::string systDirUp = "", systDirDown = "";
  std::string fileUp = "", fileDown = "";
  if(argc >= 2){
    whichSyst = std::string(argv[1]);
    if( (whichSyst.find("Matching") >= 0 && whichSyst.find("Matching") < whichSyst.size()) || 
        (whichSyst.find("Scaling") >= 0 && whichSyst.find("Scaling") < whichSyst.size()) ){
      systDirUp = "MatchingAndScaling"; systDirDown = "MatchingAndScaling";
      fileUp = whichSyst+"Up"; fileDown = whichSyst+"Down";
    }
    else if( whichSyst.find("Tag") >= 0 && whichSyst.find("Tag") < whichSyst.size() ){
      systDirUp = "BTagSF"; systDirDown = "BTagSF";
      fileUp = whichSyst+"Plus"; fileDown = whichSyst+"Minus";
    }
    else if( whichSyst.find("Lumi") >= 0 && whichSyst.find("Lumi") < whichSyst.size() ){
      systDirUp = "PileUp"; systDirDown = "PileUp";
      fileUp = whichSyst+"Down"; fileDown = whichSyst+"Up";
    }
    else{
      systDirUp = whichSyst+"Plus"; systDirDown = whichSyst+"Minus";
      fileUp = systDirUp; fileDown = systDirDown;
    }
  }
  std::cout << " * Will consider the systematics : " << whichSyst << std::endl;

  std::string LikCut = "63.869";
  if(argc >= 3){
    if(std::string(argv[2]) != "opt")
      LikCut = NumberToString(atof(argv[2]));
  }
  std::cout << " * Will consider the likelihood cut-value of : " << LikCut << std::endl;

  std::string whichMC = "";         //Possible input variables are: MC, systMC, AllTT, Signal
  if(argc >= 4){
    whichMC = std::string(argv[3]);
  }
  std::cout << " * Will consider the simulation samples of : " << whichMC << std::endl;

  //Now get the correct histograms from the different ROOT files!
  //1) Nominal comparison:
  TFile* file_Nom = new TFile(("Events_Nom/OutFile_Nom_LikelihoodCut"+LikCut+"_"+whichMC+".root").c_str());
  TH1D* h_Nom_SMLik = (TH1D*)file_Nom->Get("SMLikelihoodValue_Sum");
  h_Nom_SMLik->SetLineColor(1);

  //2) Both systematic contributions:
  TFile* file_Down = new TFile(("Events_"+systDirDown+"/OutFile_"+fileDown+"_LikelihoodCut"+LikCut+"_"+whichMC+".root").c_str());
  TH1D* h_Down_SMLik = (TH1D*)file_Down->Get("SMLikelihoodValue_Sum");
  h_Down_SMLik->SetLineColor(4);
  //h_Down_SMLik->SetLineStyle(2);

  TFile* file_Up = new TFile(("Events_"+systDirUp+"/OutFile_"+fileUp+ "_LikelihoodCut"+LikCut+"_"+whichMC+".root").c_str());
  TH1D* h_Up_SMLik = (TH1D*)file_Up->Get("SMLikelihoodValue_Sum");
  h_Up_SMLik->SetLineColor(6);
  //h_Up_SMLik->SetLineStyle(8);

  TCanvas* canv_Comp = new TCanvas("SystComparison","SystComparison");
  canv_Comp->cd();
  h_Up_SMLik->SetTitle("");
  h_Up_SMLik->GetXaxis()->SetTitle("-ln(L_{MEM}) value at g_{R} = 0"); h_Up_SMLik->GetXaxis()->SetTitleOffset(1.2);
  h_Up_SMLik->GetYaxis()->SetTitle("Number of entries"); h_Up_SMLik->GetYaxis()->SetTitleOffset(1.5);
  h_Up_SMLik->SetMaximum(h_Nom_SMLik->GetMaximum()*1.05);
  h_Up_SMLik->Draw();
  h_Nom_SMLik->Draw("same");
  h_Down_SMLik->Draw("same");
  
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(h_Nom_SMLik,"Nominal result","l");
  leg->AddEntry(h_Up_SMLik,(" "+whichSyst+" up result").c_str(),"l");
  leg->AddEntry(h_Down_SMLik,(" "+whichSyst+" down result").c_str(),"l");
  leg->Draw();

  canv_Comp->SaveAs(("SystComparison_"+whichSyst+".pdf").c_str());
}

