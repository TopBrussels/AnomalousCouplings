{

  std::string HistoNames[27] = {"YPlusGausTest","YPlusGausTestXS", "YPlusGausTestAcc", "YPlus", "YPlusXS", "YPlusAcc", "YPlusPlus", "YPlusPlusXS", "YPlusPlusAcc", "YRelPlus", "YRelPlusXS", "YRelPlusAcc", "YMin", "YMinXS", "YMinAcc", "YMinMin", "YMinMinXS", "YMinMinAcc", "YRelMin", "YRelMinXS", "YRelMinAcc", "SecondDerivative01", "SecondDerivativeXS01", "SecondDerivativeAcc01", "SecondDerivative02", "SecondDerivativeXS02", "SecondDerivativeAcc02"};

  TFile *RecoCorrectFile = new TFile("Events/RVR_RecoCorrect_SingleGausTF_10000Evts_SecondRun/FitDeviation.root");
  TFile *RecoWrongFile = new TFile("Events/RVR_RecoWrong_SingleGausTF_10000Evts/FitDeviation.root");
  TFile *outputFile = new TFile("CorrectAndWrongReco.root","RECREATE");

  for(int ii= 0; ii < 27; ii++){
    TH1F* CorrectHisto = (TH1F*) RecoCorrectFile.Get(HistoNames[ii].c_str());
    TH1F* WrongHisto = (TH1F*) RecoWrongFile.Get(HistoNames[ii].c_str());

    //TCanvas* canvas = new TCanvas(istoNames[ii].c_str()+"Combined","Combined view for correct and wrong combinations of variable "+HistoNames[ii]);
    TCanvas* canvas = new TCanvas(("Combined_"+HistoNames[ii]).c_str(),"Combined view for correct and wrong combinations");
    canvas.cd();
    CorrectHisto.SetLineColor(kGreen);
    WrongHisto.SetLineColor(kRed);
    CorrectHisto.Draw();
    WrongHisto.Draw("same");
    canvas.Write();
  }
  
  outputFile.Close();
}
