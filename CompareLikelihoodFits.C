//-----------------------------------------------------------------------------------------------------------------------------//
//  Scripts which compares the different root files and directories and directly saves the histogram!                          //
//  Run command: g++ -m64 -g -I `root-config --incdir` `root-config --libs` CompareLikelihoodFits.C -o CompareLikelihoodFits   //
//               ./CompareLikelihoodFits                                                                                       //
//-----------------------------------------------------------------------------------------------------------------------------//


#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TAxis.h"
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[]){

  if(argc >= 2){

    //Try to find what the difference is between the files:
    std::string whichWeightFile[3] = {"_CheckedEvts","_ExtraCuts","_DeletedByExtraCuts"};
    std::string weightTitle[3] = {" no extra cuts"," extra event selection cuts"," rejected by extra event sel cuts"};
    std::string dirOptions[3] = {"Correct","Wrong","Unmatched"};
    std::string dsOptions[5] = {"RecoTest","WJets","ZJets","SingleTop","Data"};
    vector<std::string> whichDir, whichDS;
    std::string CanvTitle = "";

    //Open as many ROOT files as there are input parameters!
    TFile* file_ROOT[argc-1];
    std::string dirDiff[argc-1], weightInfo[argc-1];
    for(int iArg = 0; iArg < argc-1; iArg++){
      std::cout << " ROOT file which will be opened is : " << argv[iArg+1] << endl;
      file_ROOT[iArg] = new TFile(argv[iArg+1], "READ");
      whichDir.push_back(string(argv[iArg+1]).substr(string(argv[iArg+1]).find("/")+1,(string(argv[iArg+1]).find("/Fit")-string(argv[iArg+1]).find("/")-1)));
      whichDS.push_back(string(argv[iArg+1]).substr(0,string(argv[iArg+1]).find("/")));
      dirDiff[iArg] = ""; weightInfo[iArg] = "";
    
      //Get the correct information for the legend title and the stored canvas title:
      for(int iWeight = 0; iWeight < 3; iWeight++){
        if( string(argv[iArg+1]).find(whichWeightFile[iWeight].c_str()) <=  string(argv[iArg+1]).size() )
          weightInfo[iArg] = whichWeightFile[iWeight];
      }
    }

    //Now check whether the entries saved in the whichDir are identical, otherwise their difference should be identified as well for the legend!
    bool sameDir = true, sameDS = true, sameWeight = true; 
    for(int iDir = 0; iDir < whichDir.size()-1; iDir++){
      if( whichDir[iDir] != whichDir[iDir+1]) sameDir = false;
      if( whichDS[iDir] != whichDS[iDir+1]) sameDS = false;
      if( weightInfo[iDir] != weightInfo[iDir+1] && weightInfo[iDir] != "" && weightInfo[iDir+1] != "") sameWeight = false;
    }

    if(sameDS) CanvTitle += whichDS[0]+"/";
    if(sameDir) CanvTitle += whichDir[0]+"/";
    CanvTitle += "Comparison";

    //If dataset names are not identical, store the respective dataset names
    if(!sameDS){
      CanvTitle += "_DataSets_";
      for(int iDS = 0; iDS < whichDS.size(); iDS++){
        for(int iOpt = 0; iOpt < sizeof(dsOptions)/sizeof(dsOptions[0]); iOpt++){
          if( whichDS[iDS].find( dsOptions[iOpt].c_str()) < whichDS[iDS].size() ){
            dirDiff[iDS] += dsOptions[iOpt]+" dataset";
            if(iDS != whichDS.size()-1 ) CanvTitle += dsOptions[iOpt]+"VS";
            else                         CanvTitle += dsOptions[iOpt];
          }
        }
      }
    }

    //If directory names are not identical, check whether it corresponds to Correct/Wrong/Unmatched
    if(!sameDir){
      CanvTitle+= "_Types_";
      for(int iDir = 0; iDir < whichDir.size(); iDir++){
        for(int iOpt = 0; iOpt < sizeof(dirOptions)/sizeof(dirOptions[0]); iOpt++){
          if( whichDir[iDir].find( dirOptions[iOpt].c_str()) < whichDir[iDir].size() ){
            if(dirDiff[iDir] != "") dirDiff[iDir] += ", "; 
            dirDiff[iDir] += dirOptions[iOpt]+" events";
            if( iDir != whichDir.size()-1 ) CanvTitle += dirOptions[iOpt]+"VS";
            else                            CanvTitle += dirOptions[iOpt];
          }
        }
      } 
    }

    CanvTitle += "_EventTypes";
    if(!sameWeight){
      for(int iWeight = 0; iWeight < sizeof(weightInfo)/sizeof(weightInfo[0]); iWeight++){
        if(iWeight != sizeof(weightInfo)/sizeof(weightInfo[0])-1) CanvTitle += weightInfo[iWeight]+"VS";
        else                                                      CanvTitle += weightInfo[iWeight];
      }
    }
    else
      CanvTitle += weightInfo[0];

    //Canvas where the compared histograms will be stored:
    TCanvas *canv = new TCanvas("test","test");
    canv->cd();

    //Legend to keep track of what is being compared:
    TLegend *leg = new TLegend(0.3,0.7,0.9,0.9);
    std::string whichFit[2] = {"all points","best points"};

    //Color and line style
    int lineStyle[2] = {2,1};

    //Now get the TF1's from this ROOT file:
    TF1 *fit[(argc-1)*2];
    for( int iFile = 0; iFile < (argc-1)*2; iFile++){
      if(iFile % 2 == 0) fit[iFile] = (TF1*) file_ROOT[(int)(iFile/2)]->Get("FitSums/SummedFit_FirstFit_Acc");
      else               fit[iFile] = (TF1*) file_ROOT[(int)(iFile/2)]->Get("FitSums/SummedFit_SecondFit_Acc");

      //Shift the fits down such that their minimum equals 0!
      fit[iFile]->SetParameter(0, (fit[iFile]->GetParameter(0)-fit[iFile]->GetMinimum()));
  
      //Transform the likelihood into a chi-squared
      fit[iFile]->SetParameter(0, 2*fit[iFile]->GetParameter(0));
      fit[iFile]->SetParameter(1, 2*fit[iFile]->GetParameter(1));
      fit[iFile]->SetParameter(2, 2*fit[iFile]->GetParameter(2));
    
      fit[iFile]->SetLineColor((int)(iFile/2)+1);
      fit[iFile]->SetLineStyle(lineStyle[iFile%2]);
      if(iFile ==0){        
        fit[iFile]->Draw(); //First histogram should be drawn without "same" option
        fit[iFile]->GetXaxis()->SetTitle("gR coefficient");
        fit[iFile]->GetYaxis()->SetTitle("#chi^{2} value from MadWeight ( 2*[L - L_{min}] )");
        fit[iFile]->SetTitle("Comparison of #chi^{2} distribution for different MadWeight configurations");
      }
      else{
        fit[iFile]->Draw("same");
      }

      //Set the title of the legend entry!
      std::string extraInfoTitle = "";
      for(int iWeight = 0; iWeight < sizeof(whichWeightFile)/sizeof(whichWeightFile[0]); iWeight++){
        if( weightInfo[(int)(iFile/2)] == whichWeightFile[iWeight] )
          extraInfoTitle += weightTitle[iWeight];
      }

      if( string(argv[(int)((iFile/2)+1)]).find("_CosTheta_") <= string(argv[(int)((iFile/2)+1)]).find(".root") ){
        extraInfoTitle += ", cos theta applied";
      }
      if( extraInfoTitle != "" && dirDiff[(int)(iFile/2)] != "") extraInfoTitle += ", ";

      //Now convert the fit information such that it can be stored in the legend!
      double fitSigma = fit[iFile]->GetX(fit[iFile]->GetMinimum()+0.5, fit[iFile]->GetMinimumX(), 0.2) - fit[iFile]->GetX(fit[iFile]->GetMinimum()+0.5, -0.2, fit[iFile]->GetMinimumX());
      std::ostringstream ssFitMin; ssFitMin << std::fixed << setprecision(5) << fit[iFile]->GetMinimumX(); std::string sFitMin = ssFitMin.str();
      std::ostringstream ssFitSig; ssFitSig << std::fixed << setprecision(5) << fitSigma; std::string sFitSigma = ssFitSig.str();
      std::string gRMin = sFitMin+" #pm "+sFitSigma;
      leg->AddEntry(fit[iFile],("Fit on "+whichFit[iFile%2]+" ("+extraInfoTitle+""+dirDiff[(int)(iFile/2)]+") -- gR = "+gRMin).c_str(),"l");

    }

    fit[0]->SetMaximum(500);
    fit[0]->SetMinimum(-5);
    leg->Draw();
    if(argc == 2) CanvTitle += "_FirstVSSecondFit";
    canv->SaveAs((CanvTitle+".pdf").c_str());
  }
  else
    std::cout << " Need to specify at least one directory ! " << endl;

  return 0;
}
