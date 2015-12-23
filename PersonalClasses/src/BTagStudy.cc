#include "../interface/BTagStudy.h"

BTagStudy::BTagStudy(int outputVerbose, vector<Dataset*> datasets, bool oneWP, int whichCombi, float Masslb, float sMasslb, float Massqqb, float sMassqqb){

  verbose = outputVerbose;
  singleWP_ = oneWP;
  chosenBTag_ = whichCombi;
  use5Jets_ = true;
  nrDatasets_ = datasets.size();

  //Start to initialize all variables!!
  for(int ii = 0; ii < 2; ii++){
    for(int jj = 0; jj < 9; jj++){
      NotReconstructedEvent[ii][jj]=0;
      for(int kk = 0; kk < 3; kk++){
        CorrectlyMatched[ii][kk][jj]=0;
        atLeastOneWrongMatch[ii][kk][jj] = 0;
      }
    }
  }

  //Try to include this in a more cleaner way ....
  float bjetWP[9]      = { 0.244, 0.244, 0.679, 0.679, 0.679, 0.898, 0.898, 0.898, 0.898};
  float lightjetWP[9]  = {     1, 0.244,     1, 0.679, 0.244,     1, 0.898, 0.679, 0.244};
  string optionName[9] = {"  2 L b-tags             ","  2 L b-tags, light L-veto","  2 M b-tags             ","  2 M b-tags, light M-veto","  2 M b-tags, light L-veto","  2 T b-tags             ","  2 T b-tags, light T-veto","  2 T b-tags, light M-veto","  2 T b-tags, light L-veto"};
  string bName[9] = {"(2 L b-tags)","(2L b-tags, light L veto)","(2 M b-tags)","(2 M b-tags, light M-veto)","(2M b-tags, light L-veto)","(2 T b-tags)","(2 T b-tags, light T-veto)","(2 T b-tags, light M-veto)","(2 T b-tags, light L-veto)"};
  string bTitle[9] = {"LooseTags","LooseTagsLVeto","MediumTags","MediumTagsMVeto","MediumTagsLVeto","TightTags","TightTagsTVeto","TightTagsMVeto","TightTagsLVeto"};

  if(singleWP_) nrBTags_ = 1;
  else          nrBTags_ = 9;
  for(int itBTag = 0; itBTag < nrBTags_; itBTag++){

    //In case only one b-tag option should be considered, only the first element should be filled!
    if(singleWP_){BJetWP[0] = bjetWP[chosenBTag_]; LightJetWP[0] = lightjetWP[chosenBTag_]; OptionName[0] = optionName[chosenBTag_]; BName[0] = bName[chosenBTag_]; BTitle[0] = bTitle[chosenBTag_];}
    else         {BJetWP[itBTag] = bjetWP[itBTag]; LightJetWP[itBTag] = lightjetWP[itBTag]; OptionName[itBTag] = optionName[itBTag]; BName[itBTag] = bName[itBTag]; BTitle[itBTag] = bTitle[itBTag];}

    //---  Create the MSPlots (general for all datasets)  ---//
    //MSPlot["Mlb_MassDistribution_"+BTitle[itBTag]] =new MultiSamplePlot(datasets,("Mlb_MassDistribution_"+BTitle[itBTag]).c_str(),150,50,180,("Mass of lepton and leptonic b-jet "+BName[itBTag]).c_str());
    //MSPlot["Mlb_ChiSqDistribution_"+BTitle[itBTag]] = new MultiSamplePlot(datasets,("Mlb_ChiSqDistribution_"+BTitle[itBTag]).c_str(),150,0,10,("#chi^{2}_{mlb} for b-jet choice "+BName[itBTag]).c_str());
  }

  //---  ChiSq Mlb and Mqqb information ---//
  Mlb  = Masslb; S_Mlb  = sMasslb;
  Mqqb = Massqqb; S_Mqqb = sMassqqb;
}

BTagStudy::~BTagStudy(){
    evtSelOutput[0].close();
    evtSelOutput[1].close();
}

void BTagStudy::InitializeDataSet(std::string datasetName){

  //In case a histogram specific for each dataset should be created it should be initialized here!!
  dataSetName_ = datasetName;

  for(int itBTag = 0; itBTag < nrBTags_; itBTag++){

    //---  Create the histograms  ---//
    std::string TitleInfo = BTitle[itBTag]+"_"+dataSetName_;
    histo1D["Mlb_CorrectCombiChiSq_"+TitleInfo] = new TH1F(("Mlb_CorrectCombiChiSq_"+TitleInfo).c_str(),("#chi^{2}_{mlb} for correct b-jet choice -- "+TitleInfo).c_str(),150,0,10);
    histo1D["Mlb_WrongCombiChiSq_"+TitleInfo]   = new TH1F(("Mlb_WrongCombiChiSq_"+TitleInfo).c_str(),("#chi^{2}_{mlb} for wrong b-jet choice -- "+TitleInfo).c_str(),150,0,10);

    histo1D["LowestChiSq_4Jets_"+TitleInfo] = new TH1F(("LowestChiSq_4Jets_"+TitleInfo).c_str(), ("#chi^{2} distribution for chosen jet-combination -- "+TitleInfo).c_str(), 50,0,25);
    histo1D["WrongChiSq_4Jets_"+TitleInfo] = new TH1F(("WrongChiSq_4Jets_"+TitleInfo).c_str(), ("#chi^{2} distribution for the non-chosen jet combination -- "+TitleInfo).c_str(), 50,0,25);
    histo1D["LowestChiSq_"+TitleInfo] = new TH1F(("LowestChiSq_"+TitleInfo).c_str(), ("#chi^{2} distribution for chosen jet-combination -- "+TitleInfo).c_str(), 50,0,25);
    histo1D["WrongChiSq_"+TitleInfo] = new TH1F(("WrongChiSq_"+TitleInfo).c_str(), ("#chi^{2} distribution for the non-chosen jet combinations -- "+TitleInfo).c_str(), 50,0,25);
  }
}   

void BTagStudy::ResetEventArrays(){

  for(int ii = 0; ii < 9; ii++){
    bTagJetNr[ii].clear();
    NonbTagJetNr[ii].clear();
    LightJetNr[ii].clear();
    jetIndices[ii].clear();
  }
}

void BTagStudy::CalculateJets(vector<TLorentzVector> Jets, vector<float> CSVbTagValues, vector<int> jetCombi, TLorentzVector lepton, Dataset* dataset, float weight){

  ResetEventArrays();

  for (int bTagOption = 0; bTagOption < nrBTags_; bTagOption++){

    for(unsigned int iJet = 0; iJet<Jets.size();iJet++){
      
      //Get the b-tagged and light jets:
      if(CSVbTagValues[iJet] >= BJetWP[bTagOption]){
        bTagJetNr[bTagOption].push_back(iJet);
        if(bTagJetNr[bTagOption].size() > 2 && CSVbTagValues[iJet] < LightJetWP[bTagOption]) LightJetNr[bTagOption].push_back(iJet);
      }
      else if(CSVbTagValues[iJet] < LightJetWP[bTagOption])
        LightJetNr[bTagOption].push_back(iJet);
    }
  
    //Copy the Nonbtagged vector into the light one in case the two b-tags are the same!	
    //if(BJetWP[bTagOption] == LightJetWP[bTagOption])
    //  LightJetNr[bTagOption] = NonbTagJetNr[bTagOption];     

    /*//Count how often an event has three or two light jets (for comparison of 4- and 5-jet case)    --> Maybe move this 4/5 jet comparison to a separate function!
    if(LightJetNr[bTagOption].size() >1 ){
      EventWithTwoLightJets[bTagOption]++;
      if(bTagJetNr[bTagOption].size() >1)
	EventWithTwoLightJetsAndBTagged[bTagOption]++;
    }
    if(LightJetNr[bTagOption].size() >2 ){
      EventWithThreeLightJets[bTagOption]++;
      if(bTagJetNr[bTagOption].size() >1)
	EventWithThreeLightJetsAndBTagged[bTagOption]++;
    } 
    */

    if(bTagJetNr[bTagOption].size() >= 2 && LightJetNr[bTagOption].size() >=2 ){

      //--- General case where Mqqb-Mlb is calculated simultaneously ---//
      //--   --> Lowest chi-sq is correct configuration!              --//
      int chosenCombi = getLowestMlbMqqbChiSquared( bTagOption, Jets, lepton);
      if(chosenCombi%2 == 0){ bLeptIndex[bTagOption] = (bTagJetNr[bTagOption])[0]; bHadrIndex[bTagOption] = (bTagJetNr[bTagOption])[1]; }
      else                  { bLeptIndex[bTagOption] = (bTagJetNr[bTagOption])[1]; bHadrIndex[bTagOption] = (bTagJetNr[bTagOption])[0]; }

      //Now decide on the light jets (depends whether more than 4 jets should be considered)
      light1Index[bTagOption] = (LightJetNr[bTagOption])[0]; light2Index[bTagOption] = (LightJetNr[bTagOption])[1];
      if( use5Jets_ == true && LightJetNr[bTagOption].size() > 2 && chosenCombi > 1){
        if(chosenCombi < 4){ light1Index[bTagOption] = (LightJetNr[bTagOption])[0]; light2Index[bTagOption] = (LightJetNr[bTagOption])[2];}
        else               { light1Index[bTagOption] = (LightJetNr[bTagOption])[1]; light2Index[bTagOption] = (LightJetNr[bTagOption])[2];}
      }

      //--- Now use the two chosen light jets and go to CorrectJetCombi calculation! ---//
      if(dataset->Name().find("TTbarJets_SemiLept") == 0) CompareJetCombi(jetCombi, bTagOption, use5Jets_, light1Index[bTagOption], light2Index[bTagOption]); 
    }
    else{
      if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;
      bHadrIndex[bTagOption] = 999;
      bLeptIndex[bTagOption] = 999;
    }

    jetIndices[bTagOption].push_back(bLeptIndex[bTagOption]);
    jetIndices[bTagOption].push_back(bHadrIndex[bTagOption]);
    jetIndices[bTagOption].push_back(light1Index[bTagOption]);
    jetIndices[bTagOption].push_back(light2Index[bTagOption]);
      
    //--- Additional output for debugging --//
    if(verbose > 3)
      cout<<"(BTagStudy class) -- Size of bTaggedJets: "<<bTagJetNr[bTagOption].size()<<", of NonbTaggedJets: "<<NonbTagJetNr[bTagOption].size()<<" & of lightJets: "<<LightJetNr[bTagOption].size()<<endl;
  }//Loop over all btag options!
}

int BTagStudy::getLowestMlbMqqbChiSquared(int bTagNr, vector<TLorentzVector> Jets, TLorentzVector lept){

  //Call this vector for every b-tag combination
  // --> The distinction between the different b-tag options will be made in the main class where the combiNr is stored!
  // ==> Only the lowest combiNr is stored, so in case distributions for other combi's are wanted they should be created here!!
  std::vector< std::pair<int,double> > ChiSq;
  ChiSq.clear();

  ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[0]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[1]]+Jets[(LightJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[1]]).M())/S_Mqqb,2)));
  ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[1]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[1]]).M())/S_Mqqb,2)));

  //Sort the elements of the pair on ascending order (using the second value)!
  //--> This way the first value still remains the original index and is not changed!
  std::sort(ChiSq.begin(), ChiSq.end(), sort_pred());
  histo1D["LowestChiSq_4Jets_"+BTitle[bTagNr]+"_"+dataSetName_]->Fill(ChiSq[0].second);
  histo1D["WrongChiSq_4Jets_"+BTitle[bTagNr]+"_"+dataSetName_]->Fill(ChiSq[1].second);

  //Continue adding elements in the case more than 4 jets should be considered (and they are actually present!)
  if(use5Jets_ && LightJetNr[bTagNr].size() > 2){
    ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[0]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[1]]+Jets[(LightJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[2]]).M())/S_Mqqb,2)));
    ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[1]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[2]]).M())/S_Mqqb,2)));
    ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[0]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[1]]+Jets[(LightJetNr[bTagNr])[1]]+Jets[(LightJetNr[bTagNr])[2]]).M())/S_Mqqb,2)));
    ChiSq.push_back(std::make_pair(ChiSq.size(), pow((Mlb-(lept+Jets[(bTagJetNr[bTagNr])[1]]).M())/S_Mlb,2)+pow((Mqqb-(Jets[(bTagJetNr[bTagNr])[0]]+Jets[(LightJetNr[bTagNr])[1]]+Jets[(LightJetNr[bTagNr])[2]]).M())/S_Mqqb,2)));
  }
     
  //Sort the elements of the pair on ascending order (using the second value)!
  //--> This way the first value still remains the original index and is not changed!
  std::sort(ChiSq.begin(), ChiSq.end(), sort_pred());

  //Add some information in histograms
  histo1D["LowestChiSq_"+BTitle[bTagNr]+"_"+dataSetName_]->Fill(ChiSq[0].second);
  for(unsigned int ii = 1; ii < ChiSq.size(); ii++)
    histo1D["WrongChiSq_"+BTitle[bTagNr]+"_"+dataSetName_]->Fill(ChiSq[ii].second);

  //Store the lowest chi-sq value such that it can be passed on to main class and used for additional event selection
  LowestChiSq[bTagNr] = ChiSq[0].second;

  return ChiSq[0].first;
}

void BTagStudy::CompareJetCombi(vector<int> jetCombi, int OptionNr, bool fifthJet, int lightOneIndex, int lightTwoIndex){

  int nrOptions = 2;
  if(fifthJet == false){nrOptions = 1;}

  for(int iOpt = 0; iOpt < nrOptions; iOpt++){
    int lightOne = lightOneIndex, lightTwo = lightTwoIndex;
    if(iOpt == 0){ lightOne = (LightJetNr[OptionNr])[0]; lightTwo = (LightJetNr[OptionNr])[1];}  //In case only 4 jets are considered the two light jets have index 0 & 1!
    
    //jet combi order is : (0 = BLeptonic, 1 = BHadronic, 2 = Quark1 & 3 =  Quark2)
    if(jetCombi[1] == 9999 || jetCombi[0] == 9999 || jetCombi[2] == 9999 || jetCombi[3] == 9999){
      NotReconstructedEvent[iOpt][OptionNr]++;
    }
    else{
      if( jetCombi[0] == bLeptIndex[OptionNr]  && jetCombi[1] == bHadrIndex[OptionNr]   &&
         (jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
        (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) ){
          CorrectlyMatched[iOpt][0][OptionNr]++;
      }
      if( jetCombi[0] != bLeptIndex[OptionNr]  || jetCombi[1] != bHadrIndex[OptionNr]   ||
        (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
        (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) ){
	  atLeastOneWrongMatch[iOpt][0][OptionNr]++;
      }
      if(jetCombi[0] == bLeptIndex[OptionNr] && jetCombi[1] == bHadrIndex[OptionNr])
        CorrectlyMatched[iOpt][1][OptionNr]++;

      if(jetCombi[0] != bLeptIndex[OptionNr] || jetCombi[1] != bHadrIndex[OptionNr] )
        atLeastOneWrongMatch[iOpt][1][OptionNr]++;

      if((jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
        (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) )
	  CorrectlyMatched[iOpt][2][OptionNr]++;
	
      if( (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
          (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) )
	    atLeastOneWrongMatch[iOpt][2][OptionNr]++;	
    }
  }//If 5-jet case is considered, this loop is done twice!!
}

void BTagStudy::CreateHistograms(TFile* outfile, bool savePDF, std::string pathPNG, int dataSetNr){
  //--- Use this function to create ChiSq histograms ---//
  outfile->cd();
  if(verbose > 3) std::cout << " Inside CreateHistograms function of BTagStudy class ! \n  Histograms will be filled in file : " << outfile->GetName() << " ********************************" << std::endl;

  if(dataSetNr == nrDatasets_-1 && MSPlot.size() > 0){
    TDirectory* msdir = outfile->GetDirectory("MSPlots_BTagStudy");
    if(!msdir) msdir = outfile->mkdir("MSPlots_BTagStudy");
    msdir->cd(); 
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
      MultiSamplePlot *temp = it->second;
      string name = it->first;
      temp->Draw(name, 0, false, false, false, 1);     //string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSSignal 
      temp->Write(outfile, name, savePDF, (pathPNG+"/MSPlots_BTagStudy/").c_str(), "pdf", "MSPlots_BTagStudy");
    }
  }

  //TCanvas* canv_Low = new TCanvas("canv_LowestChiSq","canv_LowestChiSq");
  //canv_Low->cd();
  //histo1D["LowestChiSq_4Jets_TightTags_TTbarJets_SemiLept_Nominal"]->SetLineColor(1);
  //histo1D["LowestChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->SetLineColor(3);
  //histo1D["LowestChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->Draw();
  //histo1D["LowestChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->GetXaxis()->SetTitle("#chi^{2} value");
  //histo1D["LowestChiSq_4Jets_TightTags_TTbarJets_SemiLept_Nominal"]->Draw("same");
  //canv_Low->SaveAs("LowestChiSqComparison.pdf");

  //TCanvas* canv_Wrong = new TCanvas("canv_WrongChiSq","canv_WrongChiSq");
  //canv_Wrong->cd();
  //histo1D["WrongChiSq_4Jets_TightTags_TTbarJets_SemiLept_Nominal"]->SetLineColor(1);
  //histo1D["WrongChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->SetLineColor(3);
  //histo1D["WrongChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->DrawNormalized();
  //histo1D["WrongChiSq_TightTags_TTbarJets_SemiLept_Nominal"]->GetXaxis()->SetTitle("#chi^{2} value");
  //histo1D["WrongChiSq_4Jets_TightTags_TTbarJets_SemiLept_Nominal"]->DrawNormalized("same");
  //canv_Wrong->SaveAs("WrongChiSqComparison.pdf");

  //Histo1D's
  if(histo1D.size() > 0){
    TDirectory* th1dir = outfile->GetDirectory("1D_histograms_BTagStudy");
    if(!th1dir) th1dir = outfile->mkdir("1D_histograms_BTagStudy");
    th1dir->cd();
    for(std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++){
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
      temp->Write();
    }
    histo1D.clear();  //Clear the map such that it can be filled with only the histograms of the considered dataset!!
  }

  //Histo2D's
  if(histo2D.size() > 0){
    TDirectory* th2dir = outfile->GetDirectory("2D_histograms_BTagStudy");
    if(!th2dir) th2dir = outfile->mkdir("2D_histograms_BTagStudy");
    th2dir->cd();
    for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++){    
      TH2F *temp = it->second;
      temp->Write();
    }
    histo2D.clear(); 
  }
  outfile->cd(); 
}

void BTagStudy::ReturnBTagTable(){ 

  std::string TexFileTitle = "AllBTags";
  if(nrBTags_ == 1) TexFileTitle = BTitle[0];	
  evtSelOutput[0].open(("EventSelectionResults/AnalyzerOutput/evtSelChoice_4JetCase_"+TexFileTitle+".tex").c_str());
  evtSelOutput[1].open(("EventSelectionResults/AnalyzerOutput/evtSelChoice_5JetCase_"+TexFileTitle+".tex").c_str());

  string Title[3]= {"  \\textbf{Option} & all 4 correct & $\\geq$ 1 wrong       & correct ($\\%$)       & $\\frac{s}{b}$ & non-matched \\\\", 
		    "  \\textbf{Option} & 2 b's correct & $\\geq$ 1 b wrong     & b's correct ($\\%$)   & $\\frac{s}{b}$ & non-matched \\\\", 
		    "  \\textbf{Option} & 2 light good  & $\\geq$ 1 light wrong & light correct ($\\%$) & $\\frac{s}{b}$ & non-matched \\\\"};

  string Caption[3] = {" \\caption{Overview of correct and wrong reconstructed events for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied", 
		       " \\caption{Overview of correct and wrong reconstructed b-jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied", 
		       " \\caption{Overview of correct and wrong reconstructed light jets for the different b-tags (no $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ applied"};

  for(int itJetCase = 0; itJetCase < 2; itJetCase++){
    for(int itWhichJets = 0; itWhichJets < 3; itWhichJets++){      //Specifying whether all jets, only b-jets or only light jets are considered!

      //--- Set the Title, Caption and variables correct ---//
      string CaptionNrJets;
      if(itJetCase == 0){      CaptionNrJets = Caption[itWhichJets]+")}"; }
      else if(itJetCase == 1){ CaptionNrJets = Caption[itWhichJets]+", 5 jets considered)}"; }
      evtSelOutput[itJetCase] << "\\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c|c} \n" << Title[itWhichJets] << " \\hline " << endl;
         
      for(int itBTag = 0; itBTag < nrBTags_; itBTag++){

	//float SOverSqrtB =  (float)(CorrectlyMatched[itJetCase][itWhichJets][itBTag])/(float)(sqrt(atLeastOneWrongMatch[itJetCase][itWhichJets][itBTag]));
        float SOverB =      (float)(CorrectlyMatched[itJetCase][itWhichJets][itBTag])/(float)(atLeastOneWrongMatch[itJetCase][itWhichJets][itBTag]);
        float CorrectPerc = (float)(CorrectlyMatched[itJetCase][itWhichJets][itBTag]*100.0)/(float)(CorrectlyMatched[itJetCase][itWhichJets][itBTag]+atLeastOneWrongMatch[itJetCase][itWhichJets][itBTag]);
	
        evtSelOutput[itJetCase] << OptionName[itBTag]                 << 
        " & " << CorrectlyMatched[itJetCase][itWhichJets][itBTag]     << 
        " & " << atLeastOneWrongMatch[itJetCase][itWhichJets][itBTag] << 
        " & " << CorrectPerc                                          << 
        " & " << SOverB                                               << 
        " & " << NotReconstructedEvent[itJetCase][itBTag]             <<
        "\\\\ " << endl;

      } //itBTag
      evtSelOutput[itJetCase] << " \\end{tabular} \n" << CaptionNrJets << "\n" << "\\end{table} \n " << endl;
    }//itWhichJets
  }//itJetCase
}

/*void BTagStudy::ReturnThirdJetTable(){
  //--- Create a separate output file for 3rd jet efficiencies! ---//
	//ThirdJetPercentage[jj] = (float)(Correct3rdJet[jj]*100.0)/(float)(CorrectOnes5Jets[jj]);
        //" & " << Correct3rdJet[itWhichJets]          << 
        //" & " << ThirdJetPercentage[itWhichJets]     << 
}*/
