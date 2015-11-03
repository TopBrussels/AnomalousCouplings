#include "../interface/ExtraEvtSelCuts.h"

ExtraEvtSelCuts::ExtraEvtSelCuts(float TopMassHadr, float sTopMassHadr, float WMassHadr, float sWMassHadr, bool oneBTag, int chiSqCut, int massWindow){

  //Initialize the counters for the chi-sq and mW/mTop cuts!
  float ChiSqCutValues[5] = {50, 30, 20, 10, 5};
  float MassWindowSigmas[5] = {5, 3, 2, 1.5, 1}; 
  oneBTag_ = oneBTag;

  nrCuts_ = 5;
  chosenChiSq_ = chiSqCut; chosenMassWindow_ = massWindow;
  WrongEvts_Mu = 0; WrongEvts_El = 0; CorrEvts_Mu = 0; CorrEvts_El = 0;

  //Now transfer this to variables existing throughout the entire class:
  for(int iCut = 0; iCut < nrCuts_; iCut++){
    ChiSqCutVal_[iCut] = ChiSqCutValues[iCut];
    MassWindow_Sigmas_[iCut] = MassWindowSigmas[iCut];

    MTop_Down_[iCut] = TopMassHadr - MassWindow_Sigmas_[iCut]*sTopMassHadr; MTop_Up_[iCut] = TopMassHadr + MassWindow_Sigmas_[iCut]*sTopMassHadr; 
    MW_Down_[iCut] = WMassHadr - MassWindow_Sigmas_[iCut]*sWMassHadr; MW_Up_[iCut] = WMassHadr + MassWindow_Sigmas_[iCut]*sWMassHadr; 

    for(int iCWU = 0; iCWU < 3; iCWU++){
      NrEvts_ChiSq_[iCWU][iCut] = 0;  NrEvts_MT_[iCWU][iCut] = 0; NrEvts_MW_[iCWU][iCut] = 0; NrEvts_MComb_[iCWU][iCut] = 0;
      if(iCut == 0){ NrEvts_AllChosenCuts_[iCWU] = 0; OriginalEvts_[iCWU] = 0;}
    }
  }

}

ExtraEvtSelCuts::~ExtraEvtSelCuts(){

}

void ExtraEvtSelCuts::Initialize(std::string bTagTitle, std::string dataSetName){

  histTitle_ = bTagTitle+"_"+dataSetName;
  //Get the Mlb-Mqqb chiSq for correct and wrong events to decide on a cut-value
  histo1D["LowestChiSq_CorrectEvents_"+histTitle_]   = new TH1F(("LowestChiSq_CorrectEvents_"+histTitle_).c_str(),   ("#chi^{2} distribution for chosen jet-combination (correct events -- "+histTitle_+")").c_str(),   50,0,80);
  histo1D["LowestChiSq_WrongEvents_"+histTitle_]     = new TH1F(("LowestChiSq_WrongEvents_"+histTitle_).c_str(),     ("#chi^{2} distribution for chosen jet-combination (wrong events -- "+histTitle_+")").c_str(),     50,0,80);
  histo1D["LowestChiSq_UnmatchedEvents_"+histTitle_] = new TH1F(("LowestChiSq_UnmatchedEvents_"+histTitle_).c_str(), ("#chi^{2} distribution for chosen jet-combination (unmatched events -- "+histTitle_+")").c_str(), 50,0,80);

  //Reset the counters!
  for(int iCut = 0; iCut < nrCuts_; iCut++){
    for(int iCWU = 0; iCWU < 3; iCWU++){
      NrEvts_ChiSq_[iCWU][iCut] = 0;  NrEvts_MT_[iCWU][iCut] = 0; NrEvts_MW_[iCWU][iCut] = 0; NrEvts_MComb_[iCWU][iCut] = 0;
      if(iCut == 0){ NrEvts_AllChosenCuts_[iCWU] = 0; OriginalEvts_[iCWU] = 0;}
    }
  }
}

bool ExtraEvtSelCuts::KeepEvent(vector<int> jetComb, TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int> selJetComb, float chiSq, int CWUIndex, int decayCh){

  bool keepEvt = false;  

  //Should only be doing all this when one b-tag should be considered!    
  if(oneBTag_){

    //Booleans in order to combine the different cuts!
    bool ChiSqSurv[5] = {false}, MTSurv[5] = {false}, MWSurv[5] = {false};
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      if(chiSq < ChiSqCutVal_[iCut] ) ChiSqSurv[iCut] = true;
      if((Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MTop_Down_[iCut] && (Jets[selJetComb[1]]+Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MTop_Up_[iCut]) MTSurv[iCut] = true;
      if((Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() > MW_Down_[iCut] && (Jets[selJetComb[2]]+Jets[selJetComb[3]]).M() < MW_Up_[iCut] ) MWSurv[iCut] = true;
    }
    if(ChiSqSurv[chosenChiSq_] && MTSurv[chosenMassWindow_] && MWSurv[chosenMassWindow_]) keepEvt = true;

    //Now get the numbers for the correct/wrong/unmatched categories!
    if(CWUIndex == 0)      histo1D["LowestChiSq_CorrectEvents_"+histTitle_]->Fill(chiSq);  
    else if(CWUIndex == 1) histo1D["LowestChiSq_WrongEvents_"+histTitle_]->Fill(chiSq);    
    else if(CWUIndex == 2) histo1D["LowestChiSq_UnmatchedEvents_"+histTitle_]->Fill(chiSq);
 
    OriginalEvts_[CWUIndex]++;

    for(int iCut = 0; iCut < nrCuts_; iCut++){
      if( ChiSqSurv[iCut] )              NrEvts_ChiSq_[CWUIndex][iCut]++;
      if( MTSurv[iCut] )                 NrEvts_MT_[CWUIndex][iCut]++;
      if( MWSurv[iCut] )                 NrEvts_MW_[CWUIndex][iCut]++;
      if( MTSurv[iCut] && MWSurv[iCut] ) NrEvts_MComb_[CWUIndex][iCut]++;
    }
    if(keepEvt) NrEvts_AllChosenCuts_[CWUIndex]++; 
  }

  return keepEvt;
}

void ExtraEvtSelCuts::StoreCutInfluence(TFile* outfile){

  if(oneBTag_){

//    std::cout << " Number of correct muon events (cuts) : " << CorrEvts_Mu << endl;
//    std::cout << " Number of correct electron events (cuts) : " << CorrEvts_El << endl;
//    std::cout << " Number of wrong muon events (cuts): " << WrongEvts_Mu << endl;
//    std::cout << " Number of wrong electron events (cuts): " << WrongEvts_El << endl;

    //Count the total number of events for each of the different categories:
    int TotalEvts = 0, TotalEvts_ChiSq[5] = {0}, TotalEvts_MT[5] = {0}, TotalEvts_MW[5] = {0}, TotalEvts_MComb[5] = {0}, TotalEvts_AllChosenCuts = 0;
    for(int iCWU = 0; iCWU < 3; iCWU++){
      TotalEvts += OriginalEvts_[iCWU];
      TotalEvts_AllChosenCuts += NrEvts_AllChosenCuts_[iCWU];
      for(int iCut = 0; iCut< nrCuts_; iCut++){
        TotalEvts_ChiSq[iCut] += NrEvts_ChiSq_[iCWU][iCut];
        TotalEvts_MT[iCut] += NrEvts_MT_[iCWU][iCut];
        TotalEvts_MW[iCut] += NrEvts_MW_[iCWU][iCut];
        TotalEvts_MComb[iCut] += NrEvts_MComb_[iCWU][iCut];
      }
    }

    ofstream CutInfl;
    CutInfl.open("EventSelectionResults/AnalyzerOutput/ExtraEventSelectionCuts.tex");

    //Initialize the tex document
    CutInfl << "\\documentclass{article} \n\\usepackage[margin=0.5in]{geometry} \n\\begin{document} \n" << endl;

    CutInfl << " \\begin{abstract} \n " << endl;
    CutInfl << "   The tables in this document represent the influence of the additional event selection cuts that were applied in order to reduce the number of selected events for CPU reasons. \\\\ " << endl;
    CutInfl << "   The considered cuts are rather basic and are merely developed to reduce the number of so-called wrong events \\\\ " << endl;
    CutInfl << "   \\begin{itemize} \n     \\item Cut on Mlb-Mqqb $\\chi^{2}$ distribution \n     \\item Cut on top and W-mass window \n   \\end{itemize} \n " << endl;
    CutInfl << "   \\textbf{Created on :} \\today " << endl;
    CutInfl << " \\end{abstract} \n " << endl;
 
    //Statistics for all events 
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Original number of events, before any additional event selection cuts are applied} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c} " << endl;
    CutInfl << "     Correct evts    & Wrong evts     & Unmatched evts      &  Total evts & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    CutInfl << "     " << OriginalEvts_[0] << "   &   " << OriginalEvts_[1] << "  & " << OriginalEvts_[2] << " & " << TotalEvts << "  & " << OriginalEvts_[0]*100.0/TotalEvts << " \n " << endl;
    CutInfl << "   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after chi-sq cut
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the Mlb-Mqqb $\\chi^{2}$.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c} " << endl;
    CutInfl << "     $\\chi^{2}$ cut-value    & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << ChiSqCutVal_[iCut] << "  &   ";
      for(int iCWU = 0; iCWU < 3; iCWU++) CutInfl << NrEvts_ChiSq_[iCWU][iCut] << "  &  "; 
      CutInfl << 100.0-(TotalEvts_ChiSq[iCut]*100.0/TotalEvts) << " & " << NrEvts_ChiSq_[0][iCut]*100.0/TotalEvts_ChiSq[iCut];
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after mT mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the $M_{top}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
    CutInfl << "     Nr $\\sigma$'s & Mass-window ($m_{top}$)    & Correct evts    & Wrong evts     & Unmatched evts   & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & " << MTop_Down_[iCut] << " - " << MTop_Up_[iCut] << "  &   ";
      for(int iCWU = 0; iCWU < 3; iCWU++) CutInfl << NrEvts_MT_[iCWU][iCut] << "  &  ";
      CutInfl << 100.0-(TotalEvts_MT[iCut]*100.0/TotalEvts) << " &  " << NrEvts_MT_[0][iCut]*100.0/TotalEvts_MT[iCut];
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after mW mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the $M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
    CutInfl << "     Nr $\\sigma$'s & Mass-window ($m_{top}$)    & Correct evts    & Wrong evts     & Unmatched evts   & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & " << MW_Down_[iCut] << " - " << MW_Up_[iCut] << "  &   ";
      for(int iCWU = 0; iCWU < 3; iCWU++) CutInfl << NrEvts_MW_[iCWU][iCut] << "  &  ";
      CutInfl << 100.0-(TotalEvts_MW[iCut]*100.0/TotalEvts) << " &  " << NrEvts_MW_[0][iCut]*100.0/TotalEvts_MW[iCut];
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after combined mT-mW mass-window
    CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on the combined $M_{top}-M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c} " << endl;
    CutInfl << "     Mass-window $\\sigma$'s & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)   & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
    for(int iCut = 0; iCut < nrCuts_; iCut++){
      CutInfl << "     " << MassWindow_Sigmas_[iCut] << " & ";
      for(int iCWU = 0; iCWU < 3; iCWU++) CutInfl << NrEvts_MComb_[iCWU][iCut] << "  &  ";
      CutInfl << 100.0-(TotalEvts_MComb[iCut]*100.0/TotalEvts) << " &  " << NrEvts_MComb_[0][iCut]*100.0/TotalEvts_MComb[iCut];
      if(iCut < nrCuts_-1) CutInfl << " \\\\ " << endl;
    }
    CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;

    //Numbers after all cuts combined
    if(chosenChiSq_ != -1 && chosenMassWindow_ != -1){
      CutInfl << " \\begin{table}[h!t] \n  \\caption{Remaining number of events, after applying cuts on both the Mlb-Mqqb $\\chi^{2}$ and the combined $M_{top}-M_{W}$ mass-window.} \n  \\centering \n   \\begin{tabular}{c|c|c|c|c|c|c} " << endl;
      CutInfl << "     Mass-window $\\sigma$'s & $\\chi^{2}$ cut-value    & Correct evts    & Wrong evts     & Unmatched evts  & Evt reduction ($\\%$)    & s/b ($\\%$)     \\\\ \n     \\hline" << endl;
      CutInfl << "     " << MassWindow_Sigmas_[chosenMassWindow_] << " & " << ChiSqCutVal_[chosenChiSq_] << " & ";
      for(int iCWU = 0; iCWU < 3; iCWU++) CutInfl << NrEvts_AllChosenCuts_[iCWU] << " & ";
      CutInfl << 100.0-(TotalEvts_AllChosenCuts*100.0/TotalEvts)  << " & " << NrEvts_AllChosenCuts_[0]*100.0/TotalEvts_AllChosenCuts << endl;
      CutInfl << "\n   \\end{tabular} \n \\end{table} \n " << endl;
    }
    else{
      cout << " No cut-values given for chi-sq and mass-window " << endl;
    }
  
    CutInfl << "\\end{document} " << endl;
    CutInfl.close();

    //Also store the plots!!
    outfile->cd();
    if(histo1D.size() > 0){
      TDirectory* th1dir = outfile->GetDirectory("1D_histograms_ExtraEvtSelCuts");   //Check whether directory already exists ..
      if(!th1dir) th1dir = outfile->mkdir("1D_histograms_ExtraEvtSelCuts");          // .. and otherwise create it!
      th1dir->cd();
      for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
      }
    }
   
    if(histo2D.size() > 0){ 
      TDirectory* th2dir = outfile->GetDirectory("2D_histograms_ExtraEvtSelCuts");
      if(!th2dir) th2dir = outfile->mkdir("2D_histograms_ExtraEvtSelCuts");
      th2dir->cd();
      for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
        TH2F *temp = it->second;
        temp->Write();
      }
    }
    outfile->cd(); 
  }
}
