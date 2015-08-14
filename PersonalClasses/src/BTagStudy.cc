#include "../interface/BTagStudy.h"

BTagStudy::BTagStudy(int outputVerbose, vector<Dataset*> datasets){
  BTagStudy::InitializeBegin(datasets);
  verbose = outputVerbose;
  nrDatasets_ = datasets.size();  
    //evtSelOutput[0].open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/EventSelectionResults/AnalyzerOutput/eventSelectionChoiceTables4JetCase.tex");
    //evtSelOutput[1].open("/user/aolbrech/GitTopTree_Feb2014/TopBrussels/AnomalousCouplings/EventSelectionResults/AnalyzerOutput/eventSelectionChoiceTables5JetCase.tex");
}

BTagStudy::~BTagStudy(){
    evtSelOutput[0].close();
    evtSelOutput[1].close();
}

void BTagStudy::InitializeBegin(vector<Dataset*> datasets){

  for(int ii = 0; ii < 6; ii++){
    for(int jj = 0; jj < 2; jj++){
      NotReconstructedEvent[ii][jj]=0;
      for(int kk = 0; kk < 2; kk++){
        CorrectlyMatched[kk][ii][jj]=0;
        atLeastOneWrongMatch[kk][ii][jj] = 0;
      }
    }
  }

  //Try to include this in a more cleaner way ....
  float bjetWP[6]      = {   0.244,                      0.679,                      0.679,                       0.898,                      0.898,                       0.898                    };
  float lightjetWP[6]  = {   0.244,                      0.679,                                  0.244,           0.898,                                  0.679,                       0.244        };
  string optionName[6] = {"  2 L b-tags             ","  2 M b-tags             ","  2 M b-tags, light L-veto","  2 T b-tags             ","  2 T b-tags, light M-veto","  2 T b-tags, light L-veto"};
  string bName[6] = {"(2 L b-tags)","(2 M b-tags)","(2 M b-tags, light L-veto)","(2 T b-tags)","(2 T b-tags, light M-veto)","(2 T b-tags, light L-veto)"};
  string bTitle[6] = {"LooseTags","MediumTags","MediumTagsLVeto","TightTags","TightTagsMVeto","TightTagsLVeto"};
  for(int ii = 0; ii < 6; ii++){
    BJetWP[ii] = bjetWP[ii];
    LightJetWP[ii] = lightjetWP[ii];
    OptionName[ii] = optionName[ii];
    BName[ii] = bName[ii];
    BTitle[ii] = bTitle[ii];
  }

  //---  ChiSq Mlb and Mqqb information ---//
  Mlb = 103.286;
  SigmaMlb = 26.7764;
  Mqqb = 178.722;
  SigmaMqqb = 18.1385;

  //---  Create the histograms  ---//
  for(int itBTag = 0; itBTag < 6; itBTag++){
    MSPlot["Mlb_MassDistribution_"+BTitle[itBTag]] = new MultiSamplePlot(datasets,("Mlb_MassDistribution_"+BTitle[itBTag]).c_str(),150,50,180, ("Mass of lepton and leptonic b-jet "+BName[itBTag]).c_str());
    MSPlot["Mlb_ChiSqDistribution_"+BTitle[itBTag]] = new MultiSamplePlot(datasets,("Mlb_ChiSqDistribution_"+BTitle[itBTag]).c_str(),150, 0, 10, ("#chi^{2}_{mlb} for b-jet choice "+BName[itBTag]).c_str());
    histo1D["Mlb_CorrectCombiChiSq_"+BTitle[itBTag]] = new TH1F(("Mlb_CorrectCombiChiSq_"+BTitle[itBTag]).c_str(),("#chi^{2}_{mlb} for correct b-jet choice "+BName[itBTag]).c_str(),150,0,10);
    histo1D["Mlb_WrongCombiChiSq_"+BTitle[itBTag]]   = new TH1F(("Mlb_WrongCombiChiSq_"+BTitle[itBTag]).c_str(),("#chi^{2}_{mlb} for wrong b-jet choice "+BName[itBTag]).c_str(),150,0,10);
  }
}

void BTagStudy::ResetEventArrays(){

  for(int ii = 0; ii < 6; ii++){
    bTaggedJetNr[ii].clear();
    NonbTaggedJetNr[ii].clear();
    LightJetNr[ii].clear();
    ChiSquaredMlb[ii].clear();
    ChiSquaredMqqb[ii].clear();
  }
}

void BTagStudy::CalculateJets(vector<TLorentzVector> Jets, vector<float> CSVbTagValues, vector<int> jetCombi, TLorentzVector lepton, Dataset* dataset, float weight){

  ResetEventArrays();

  for (unsigned int bTagOption = 0; bTagOption < 6; bTagOption++){
    for(unsigned int iJet = 0; iJet<Jets.size();iJet++){
      if(CSVbTagValues[iJet] >= BJetWP[bTagOption])
        bTaggedJetNr[bTagOption].push_back(iJet);            
      else
        NonbTaggedJetNr[bTagOption].push_back(iJet);            

      //Calculate the light jets when an additional working point is required for these:
      if(BJetWP[bTagOption] != LightJetWP[bTagOption] && CSVbTagValues[iJet] < LightJetWP[bTagOption])
	LightJetNr[bTagOption].push_back(iJet);
    }

    //Copy the Nonbtagged vector into the light one in case the two b-tags are the same!	
    if(BJetWP[bTagOption] == LightJetWP[bTagOption])
      LightJetNr[bTagOption] = NonbTaggedJetNr[bTagOption];     

    /*//Count how often an event has three or two light jets (for comparison of 4- and 5-jet case)    --> Maybe move this 4/5 jet comparison to a separate function!
    if(LightJetNr[bTagOption].size() >1 ){
      EventWithTwoLightJets[bTagOption]++;
      if(bTaggedJetNr[bTagOption].size() >1)
	EventWithTwoLightJetsAndBTagged[bTagOption]++;
    }
    if(LightJetNr[bTagOption].size() >2 ){
      EventWithThreeLightJets[bTagOption]++;
      if(bTaggedJetNr[bTagOption].size() >1)
	EventWithThreeLightJetsAndBTagged[bTagOption]++;
    } 
    */

    if(bTaggedJetNr[bTagOption].size() >= 2 && LightJetNr[bTagOption].size() >=2 ){

      //--- Use Mlb chi-sq to select the b-jets ---//
      CalculateMlbChiSq(bTagOption, lepton, Jets, jetCombi, dataset, weight);

      //--- For 4-jet case jet combi numbers can be calculated directly ---//
      light1Index4Jets[bTagOption] = (LightJetNr[bTagOption])[0];
      light2Index4Jets[bTagOption] = (LightJetNr[bTagOption])[1];
      CompareJetCombi(jetCombi, bTagOption, 0, light1Index4Jets[bTagOption], light2Index4Jets[bTagOption]);      //4-jet case

      //--- In case of 5 jets, select the two correct light ones ---//
      light1Index5Jets[bTagOption] = (LightJetNr[bTagOption])[0];
      light2Index5Jets[bTagOption] = (LightJetNr[bTagOption])[1];
      if(LightJetNr[bTagOption].size() > 2){ 
        vector<int> lightMqqbIndex = CalculateMqqbChiSq(bTagOption, Jets);
        light1Index5Jets[bTagOption] = lightMqqbIndex[0];
        light2Index5Jets[bTagOption] = lightMqqbIndex[1];
      }

      //--- Now use the two chosen light jets and go to CorrectJetCombi calculation! ---//
      CompareJetCombi(jetCombi, bTagOption, 1, light1Index5Jets[bTagOption], light2Index5Jets[bTagOption]); //4+5-jet case
    }
    else{
      if(verbose > 3) std::cout << " Event doesn't have two b-tagged jets and/or two light jets ! " << std::endl;
      bHadrIndex[bTagOption] = 999;
      bLeptIndex[bTagOption] = 999;
      LowestChiSqMlb[bTagOption] = 999;
      LowestChiSqMqqb[bTagOption] = 999;
      for(int ii = 0; ii < 2; ii++) 
        ChiSquaredMlb[bTagOption].push_back(999);  //Also something similar for ChiSquaredMqqb?
    }
      
    //--- Additional output for debugging --//
    if(verbose > 3)
      cout<<"(BTagStudy class) -- Size of bTaggedJets: "<<bTaggedJetNr[bTagOption].size()<<", of NonbTaggedJets: "<<NonbTaggedJetNr[bTagOption].size()<<" & of lightJets: "<<LightJetNr[bTagOption].size()<<endl;
  }//Loop over all btag options!
}

void BTagStudy::CalculateMlbChiSq(int bTagNr, TLorentzVector lepton, vector<TLorentzVector> Jets, vector<int> correctJetCombi, Dataset* dataSet, float msPlotScale){

  //---  This will distinguish the leptonic and hadronic b-jets  ---//
  //---     ==> Identical for 4- and 5-jet case                  ---//
  float MlbValues[2]  = {(lepton+Jets[(bTaggedJetNr[bTagNr])[0]]).M(), (lepton+Jets[(bTaggedJetNr[bTagNr])[1]]).M()};
 
  LowestChiSqMlb[bTagNr] = 0;
  bHadrIndex[bTagNr] = (bTaggedJetNr[bTagNr])[1];
  for(int ii = 0; ii<2; ii++){
    ChiSquaredMlb[bTagNr].push_back((float)(((Mlb-(float)(MlbValues[ii]))/SigmaMlb)*((Mlb-(float)(MlbValues[ii]))/SigmaMlb))); 
    if((ChiSquaredMlb[bTagNr])[ii] < (ChiSquaredMlb[bTagNr])[LowestChiSqMlb[bTagNr]]){
      LowestChiSqMlb[bTagNr] = ii;
      bHadrIndex[bTagNr] = (bTaggedJetNr[bTagNr])[0];
    }
  }
  bLeptIndex[bTagNr] = (bTaggedJetNr[bTagNr])[LowestChiSqMlb[bTagNr]];
  MSPlot["Mlb_MassDistribution_"+BTitle[bTagNr]]->Fill( (lepton+Jets[bLeptIndex[bTagNr]]).M(), dataSet, true, msPlotScale);
  MSPlot["Mlb_ChiSqDistribution_"+BTitle[bTagNr]]->Fill( (ChiSquaredMlb[bTagNr])[LowestChiSqMlb[bTagNr]], dataSet, true, msPlotScale );

  //--- Save information in histograms ---//
  //jet combi order is : 0 = BLeptonic & 1 = BHadronic
  if(correctJetCombi[0] != 9999 && correctJetCombi[1] != 9999){
    if(bLeptIndex[bTagNr] == correctJetCombi[0] && bHadrIndex[bTagNr] == correctJetCombi[1]) histo1D["Mlb_CorrectCombiChiSq_"+BTitle[bTagNr]]->Fill( (ChiSquaredMlb[bTagNr])[LowestChiSqMlb[bTagNr]]);
    if(bLeptIndex[bTagNr] != correctJetCombi[0] || bHadrIndex[bTagNr] != correctJetCombi[1]) histo1D["Mlb_WrongCombiChiSq_"+BTitle[bTagNr]]->Fill( (ChiSquaredMlb[bTagNr])[LowestChiSqMlb[bTagNr]]);
  }
}

vector<int> BTagStudy::CalculateMqqbChiSq(int bTagNr, vector<TLorentzVector> Jets){
  //---  This option chooses 2 out of three light jets  ---//
  //---     ==> Only considered in 5-jet case           ---//
  vector<int> lightIndices;
  lightIndices.clear();

  float MqqbValues[3] = {(Jets[bHadrIndex[bTagNr]] + Jets[(LightJetNr[bTagNr])[0]] + Jets[(LightJetNr[bTagNr])[1]]).M(), 
                         (Jets[bHadrIndex[bTagNr]] + Jets[(LightJetNr[bTagNr])[0]] + Jets[(LightJetNr[bTagNr])[2]]).M(),
                         (Jets[bHadrIndex[bTagNr]] + Jets[(LightJetNr[bTagNr])[1]] + Jets[(LightJetNr[bTagNr])[2]]).M()};
  
  LowestChiSqMqqb[bTagNr] = 0;
  for(int ii = 0; ii<3; ii++){
    ChiSquaredMqqb[bTagNr].push_back((float)(((Mqqb-(float)(MqqbValues[ii]))/SigmaMqqb)*((Mqqb-(float)(MqqbValues[ii]))/SigmaMqqb))); 
    if((ChiSquaredMqqb[bTagNr])[ii] < (ChiSquaredMqqb[bTagNr])[LowestChiSqMqqb[bTagNr]]){
      LowestChiSqMqqb[bTagNr] = ii;
    }
  }
  
  //Set the light jet indices:
  if(     LowestChiSqMqqb[bTagNr] == 0){ lightIndices.push_back( (LightJetNr[bTagNr])[0]); lightIndices.push_back((LightJetNr[bTagNr])[1]);}
  else if(LowestChiSqMqqb[bTagNr] == 1){ lightIndices.push_back( (LightJetNr[bTagNr])[0]); lightIndices.push_back((LightJetNr[bTagNr])[2]);}
  else if(LowestChiSqMqqb[bTagNr] == 2){ lightIndices.push_back( (LightJetNr[bTagNr])[1]); lightIndices.push_back((LightJetNr[bTagNr])[2]);} 

  return lightIndices;
}

void BTagStudy::CreateHistograms(TFile* outfile, bool savePDF, std::string pathPNG, int datasetNr){
  //--- Use this function to create ChiSq histograms ---//
  outfile->cd();
  if(verbose > 3) std::cout << " Inside CreateHistograms function of BTagStudy class ! \n  Histograms will be filled in file : " << outfile->GetName() << " ********************************" << std::endl;

  //Only write out the MSPlots if the datasetNr == nrDatasets_-1!
  std::cout << " Does MSPlots directory exist ?? --> " << outfile->GetDirectory("MSPlots") << std::endl;
  if(datasetNr == nrDatasets_-1 && MSPlot.size() > 0){
    TDirectory* msdir = outfile->GetDirectory("MSPlots_BTagStudy");
    if(!msdir) msdir = outfile->mkdir("MSPlots_BTagStudy");
    std::cout << " Number of MultiSamplePlots is : " <<MSPlot.size() << std::endl; 
    msdir->cd(); 
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){    
      MultiSamplePlot *temp = it->second;
      string name = it->first;
      temp->Draw(name, 0, false, false, false, 1);     //string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSSignal 
      temp->Write(outfile, name, savePDF, (pathPNG+"/MSPlots_BTagStudy/").c_str(), "pdf", "MSPlots_BTagStudy");
    }
    std::cout << " Does MSPlots directory exists after saving the bTag MSPlots ?? --> " << outfile->GetDirectory("MSPlots") << std::endl;
  }

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
  }
  outfile->cd(); 
}

void BTagStudy::CompareJetCombi(vector<int> jetCombi, int OptionNr, int NrJets, int lightOne, int lightTwo){

  //jet combi order is : (0 = BLeptonic, 1 = BHadronic, 2 = Quark1 & 3 =  Quark2)
  if(jetCombi[1] == 9999 || jetCombi[0] == 9999 || jetCombi[2] == 9999 || jetCombi[3] == 9999){
    NotReconstructedEvent[NrJets][OptionNr]++;
  }
  else{
    if( jetCombi[0] == bLeptIndex[OptionNr]  && jetCombi[1] == bHadrIndex[OptionNr]   &&
       (jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
       (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) ){
        CorrectlyMatched[NrJets][0][OptionNr]++;
    }
    if( jetCombi[0] != bLeptIndex[OptionNr]  || jetCombi[1] != bHadrIndex[OptionNr]   ||
       (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
       (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) ){
	atLeastOneWrongMatch[NrJets][0][OptionNr]++;
    }
    if(jetCombi[0] == bLeptIndex[OptionNr] && jetCombi[1] == bHadrIndex[OptionNr])
      CorrectlyMatched[NrJets][1][OptionNr]++;

    if(jetCombi[0] != bLeptIndex[OptionNr] || jetCombi[1] != bHadrIndex[OptionNr] )
      atLeastOneWrongMatch[NrJets][1][OptionNr]++;

    if((jetCombi[2] == lightOne || jetCombi[2] == lightTwo) &&
       (jetCombi[3] == lightOne || jetCombi[3] == lightTwo) )
	CorrectlyMatched[NrJets][2][OptionNr]++;
	
    if( (jetCombi[2] != lightOne && jetCombi[2] != lightTwo) ||
        (jetCombi[3] != lightOne && jetCombi[3] != lightTwo) )
	atLeastOneWrongMatch[NrJets][2][OptionNr]++;	
  }
}

void BTagStudy::ReturnBTagTable(std::string dataSetName){ 

  evtSelOutput[0].open(("EventSelectionResults/AnalyzerOutput/evtSelChoice_4JetCase_"+dataSetName+".tex").c_str());
  evtSelOutput[1].open(("EventSelectionResults/AnalyzerOutput/evtSelChoice_5JetCase_"+dataSetName+".tex").c_str());

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
         
      for(int itBTag = 0; itBTag < 6; itBTag++){

	float SOverSqrtB =  (float)(CorrectlyMatched[itJetCase][itWhichJets][itBTag])/(float)(sqrt(atLeastOneWrongMatch[itJetCase][itWhichJets][itBTag]));
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

void BTagStudy::ReturnThirdJetTable(){
  //--- Create a separate output file for 3rd jet efficiencies! ---//
	//ThirdJetPercentage[jj] = (float)(Correct3rdJet[jj]*100.0)/(float)(CorrectOnes5Jets[jj]);
        //" & " << Correct3rdJet[itWhichJets]          << 
        //" & " << ThirdJetPercentage[itWhichJets]     << 
}
