{
  //----------------------------------//
  //   Only draw specific histograms  //
  //----------------------------------//
  bool drawIndivHistos = false;
  bool drawColorHistos = true;
  bool drawDoubleGaus = false;

  //-------------------------------------------------------------------------------------//
  //  Values which need to be changed when a different directory should be accessed !!   //
  //-------------------------------------------------------------------------------------//
  //ROOT file which should be considered:
  string Directory = "FromTree";
  //Parameters:
  int NrBins = 11;  //Also update on loop over bins!
  string HistoBin[NrBins] = {"1","2","3","4","5","6","7","8","9","10","11"};
  const int NrParams = 6;                                                   //Looking at doubleGaussian with 6 parameters!
  const int NrDblGaus = 15;                                                 //Number of different E-bins have been considered for creating separate double Gaussian distributions!
  string GenE[NrDblGaus] = {"10","15","20","30","40","55","70","85","100","115","130","145","160","180","200"};   
  string GenInvE[NrDblGaus] = {"0.1","0.0667", "0.05", "0.033", "0.025", "0.01818", "0.014", "0.0118", "0.01", "0.00869", "0.00769", "0.00689", "0.00625", "0.00556", "0.005"};
  string Param[NrParams] = {"a1","a2","a3","a4","a5","a6"};                 //Names of the fit parameters  --> Check whether this is also correct for the ROOT class!
  bool GetMkdirOutput = false;                                              //Output for directories which need to be created wanted ??
  const int NrFitHistos = 12;                                               //Number of histograms which need to be considered!
  if(NrFitHistos == 1) int ChosenHisto = 1;                                 //Choose 1 specific histogram!
  //-------------------------------------------------------------------------------------//
  
  TFile* histoFile = new TFile( ("../CreatedTFFromDistributions_"+Directory+".root").c_str() ,"READ");
  std::cout << endl << " //--------------------------------------------------------------------------------------------------------" << endl;
  std::cout <<         " //-----------  Used ROOT file : " << histoFile->GetName() << endl;
  std::cout <<         " //--------------------------------------------------------------------------------------------------------" << endl;

  //Store the EtaBin Title and Name for the histograms!
  nEtaBins = 4;
  std::string EtaBin[5], EtaTitle[5];
  float EtaValues[6];
  std::string EtaValuesString[6];
  EtaBin[0] = ""; EtaTitle[0] = " -- for all |#eta| values";
  if(nEtaBins == 4){         
        EtaValues[1] = 0.; EtaValues[2] = 0.375; EtaValues[3] = 0.750; EtaValues[4] = 1.450; EtaValues[5] = 2.5;
        EtaValuesString[1] = "0"; EtaValuesString[2] = "0.375"; EtaValuesString[3] = "0.75"; EtaValuesString[4] = "1.45"; EtaValuesString[5] = "2.5";
        for(int ii = 1; ii <= nEtaBins; ii++){
            EtaBin[ii] = "_Eta_"+EtaValuesString[ii]+"_"+EtaValuesString[ii+1];
            EtaTitle[ii] = " -- "+EtaValuesString[ii]+" < |#eta| #leq "+EtaValuesString[ii+1];
        }
  }

  //***********************//
  int usedEta = 0;
  //***********************//
  Directory = Directory+EtaBin[usedEta];

  std::string Titles[12] ={"Light_DiffEVsGenE"+EtaBin[usedEta],      // Number 0
			   "Light_DiffThetaVsGenE"+EtaBin[usedEta],   // Number 4
			   "Light_DiffPhiVsGenE"+EtaBin[usedEta],     // Number 8
			   "BJet_DiffEVsGenE"+EtaBin[usedEta],       // Number 1
			   "BJet_DiffThetaVsGenE"+EtaBin[usedEta],    // Number 5
			   "BJet_DiffPhiVsGenE"+EtaBin[usedEta],      // Number 9
			   "El_DiffEVsGenE"+EtaBin[usedEta],         // Number 2
			   "El_DiffThetaVsGenE"+EtaBin[usedEta],      // Number 6
			   "El_DiffPhiVsGenE"+EtaBin[usedEta],        // Number 10
			   "Mu_DiffInvEVsGenInvE"+EtaBin[usedEta],   // Number 3
			   "Mu_DiffThetaVsGenInvE"+EtaBin[usedEta],   // Number 7
			   "Mu_DiffPhiVsGenInvE"+EtaBin[usedEta]};    // Number 11
  string HistoAxisName[5] ={"Difference (gen-reco) of E","Difference (gen-reco) of inverse E ", "Difference (gen-reco) of transverse p","Difference (gen-reco) of #phi","Difference (gen-reco) of #theta"};
  string ChiXAxisName[3] = {"Energy of generator parton (GeV)","Inverse energy of generator parton (1/GeV)","Transverse momentum of generator parton (GeV)"};
  string HistoTitle[4] =   {"light jets","b-jets","muon","electron"};

  if(GetMkdirOutput == true){
    //--------------------------------------------//
    //  First time: Create needed directories !!  //    
    //--------------------------------------------//
    std::cout << endl << " *** Directories which should be created ! *** " << endl << endl;
    std::cout << " mkdir " << Directory << " ";
    for(int ii = 0; ii < NrFitHistos; ii++) std::cout << Directory << "/" << Titles[ii] << " ";
    std::cout << " " << endl;
    //--------------------------------------------//
    std::cout << endl << " --> If histograms need to be created, change value of boolean 'GetMkdirOutput' to false !! " << endl;
  }
  else{
    
    //////////////////////////////////
    // Save the fitting histograms! //
    //////////////////////////////////
    TCanvas* twoDCanvas = new TCanvas("fitCanvas","2-dimensional distributions used as fitting input");
    twoDCanvas->Divide(3,4);
    TH1D* twoDHisto;
   
    for(int iHisto = 0; iHisto < NrFitHistos; iHisto++){    //Loop over each of the fit histograms
      if(NrFitHistos == 1) iHisto = ChosenHisto;
     
      string Title = Titles[iHisto];
      string histoName, HistoName, histoTitle, histoAxisName, chiXAxisName;
      TH1D* ProjHisto;
      
      //Set the correct X-axis name for the chi2 distribution
      if(Title.find("GenE") <= Title.size())          chiXAxisName = ChiXAxisName[0];   
      else if(Title.find("GenInvE") <= Title.size())  chiXAxisName = ChiXAxisName[1];
      else if(Title.find("GenPt") <= Title.size())    chiXAxisName = ChiXAxisName[2];
	  
      //Set the title for the correct particle
      if(Title.find("Light_") == 0)      histoTitle = HistoTitle[0];         
      else if(Title.find("BJet_") == 0)  histoTitle = HistoTitle[1];
      else if(Title.find("Mu_") == 0)    histoTitle = HistoTitle[2];
      else if(Title.find("El_") == 0)    histoTitle = HistoTitle[3];
      
      //Set the correct axis variable
      if(Title.find("DiffE") <= Title.size())         histoAxisName = HistoAxisName[0];
      else if(Title.find("DiffInvE") <= Title.size())  histoAxisName = HistoAxisName[1];
      else if(Title.find("DiffPt") <= Title.size())     histoAxisName = HistoAxisName[2];
      else if(Title.find("DiffPhi") <= Title.size())   histoAxisName = HistoAxisName[3];
      else if(Title.find("DiffTheta") <= Title.size()) histoAxisName = HistoAxisName[4];      
      std::cout << " --- Found strings : (chiXAxis, histoTitle & histoAxisName) = " << chiXAxisName << ", " << histoTitle << " & " << histoAxisName << std::endl;
 
      ///////////////////////////////
      //  Draw fitting histograms  //
      ///////////////////////////////
      twoDHisto = (TH1D*) histoFile->Get( ("2D_histograms_graphs/"+Title).c_str() );
    std::cout << " Looking at histogram : " << ("2D_histograms_graphs/"+Title).c_str() <<std::endl;
      twoDHisto->GetXaxis()->SetTitle( chiXAxisName.c_str() );
      twoDHisto->GetYaxis()->SetTitle( histoAxisName.c_str() );
      twoDHisto->SetTitle( ("Color plot for "+histoTitle+EtaTitle[usedEta]).c_str() );
      twoDHisto->GetXaxis()->SetTitleOffset(0.85);
      twoDHisto->GetYaxis()->SetTitleOffset(0.75);
      twoDCanvas->cd(iHisto+1);
      twoDHisto->Draw("colz");

    if(drawIndivHistos == true && drawColorHistos == true){      
      TCanvas* indiv2DCanvas =  new TCanvas("indiv2DCanvas","title");
      indiv2DCanvas->cd();
      twoDHisto->Draw("colz");
      //indiv2DCanvas->SaveAs( (Directory+"/"+Title+".png").c_str() );
      indiv2DCanvas->SaveAs( (Directory+"/"+Title+".pdf").c_str() );
      delete indiv2DCanvas;
    }
       
      ///////////////////////////////////////
      //  Save the ProjectionY histograms  //
      ///////////////////////////////////////
      TCanvas* projCanvas = new TCanvas("projCanvas","Fit results of projectionY for each bin");
      projCanvas->Divide(3,4);
    
      if(Title == "Mu_DiffInvEVsGenInvE") NrBins = 10;
      else NrBins = 11;
      for(int iBin = 0; iBin <= NrBins; iBin++){
	//Set name for .root file and for saving the histograms!
	if(iBin < NrBins){
	  histoName = ("sliceYbin"+HistoBin[iBin]).c_str();          
	  HistoName = ("ProjectionY_Bin"+HistoBin[iBin]).c_str();
	}
	else{
	  histoName = "chi2";
	  HistoName = "Chi2Distribution";
	}
        projHisto = (TH1D*) histoFile->Get( (Title+"/"+Title+"_"+histoName).c_str() );
	if(projHisto == 0) continue;

	//Set the Marker to a small dot
	projHisto->SetMarkerStyle(8); //Only use a small dot!
	projHisto->SetMarkerSize(0.3);
	projHisto->GetXaxis()->SetTitleOffset(0.85);
	projHisto->GetYaxis()->SetTitleOffset(0.80);

	//Draw the histograms in the different pads of the created canvas and set Title & axis label!
	projCanvas->cd(iBin+1);
	if(iBin < NrBins)
	  projHisto->Draw("e");
	else{
	  projCanvas_12->SetLogy();          //Last canvas (chi2) should be in logY 
	  projHisto->Draw();                 //No error bars for the chi2 histogram!
	}	
	
	//Save each histogram separately as well!
      if(drawIndivHistos == true){
	TCanvas* projIndivCanvas = new TCanvas("projIndivCanvas","title");;
	projIndivCanvas->cd();
	if(iBin < NrBins){
	  projHisto->SetTitle( ("Fit result for "+histoTitle+" (bin "+HistoBin[iBin]+")").c_str() );
	  projHisto->GetXaxis()->SetTitle( histoAxisName.c_str() );	  
	}
	else{
	  projHisto->SetTitle(  (histoAxisName+": Fit result").c_str());     
	  projHisto->GetYaxis()->SetTitle("#chi^{2} for double Gaussian fit");
	  projIndivCanvas->SetLogy();
	  projHisto->GetXaxis()->SetTitle(chiXAxisName.c_str() );	       
	}
	projHisto->Draw();	
	//projIndivCanvas->SaveAs( (Directory+"/"+Title+"/"+HistoName+".png").c_str() );    //Individual histograms!
	projIndivCanvas->SaveAs( (Directory+"/"+Title+"/"+HistoName+".pdf").c_str() );
	delete projIndivCanvas;
      }
	
      }//End of loop over NrBins

      //projCanvas->SaveAs( (Directory+"/"+Title+"/Overview_FitDistributions/.png").c_str() ); //Histograms together in one canvas!
      projCanvas->SaveAs( (Directory+"/"+Title+"/Overview_FitDistributions.pdf").c_str() );      
      delete projCanvas;
      delete projHisto;      

      ////////////////////////////////////////////
      //  Save the fit parameter distributions  //
      ////////////////////////////////////////////
      TCanvas* paramCanvas = new TCanvas("paramCanvas","Fit results for the different parameters");
      paramCanvas->Divide(2,3);

      TH1D* paramHisto;
      for(int iPar = 0; iPar < NrParams; iPar++){
	paramHisto = (TH1D*) histoFile->Get( (Title+"/"+Title+"_"+Param[iPar]+"_PointsAndFit").c_str() );
	paramHisto->SetTitle( ("Fitted value of parameter "+Param[iPar]+EtaTitle[usedEta]).c_str());
	paramHisto->GetXaxis()->SetTitle( ChiXAxisName.c_str() );
	paramHisto->GetXaxis()->SetTitleOffset(0.85);
	paramHisto->SetMarkerStyle(8); //Only use a small dot!
	paramHisto->SetMarkerSize(0.3);
	
	//Draw the histograms in the different pads of the created canvas!
	paramCanvas->cd(iPar+1);
	paramHisto->Draw("e");

	//Also draw the histograms individually
        if(drawIndivHistos == true){
	TCanvas* paramIndivCanvas = new TCanvas("paramIndivCanvas","title");
	paramIndivCanvas->cd();
	paramHisto->Draw("e");
	string ParamHistoName = ("FitParameter_"+Param[iPar]).c_str();
	//paramIndivCanvas->SaveAs( (Directory+"/"+Title+"/"+ParamHistoName+".png").c_str() ); //Parameter histograms individual!
	paramIndivCanvas->SaveAs( (Directory+"/"+Title+"/"+ParamHistoName+".pdf").c_str() );
	delete paramIndivCanvas;
        }
	
      }//End of loop over NrParams

      //paramCanvas->SaveAs( (Directory+"/"+Title+"/Overview_FitParameters.png").c_str() ); //Parameter histograms together in one canvs!
      paramCanvas->SaveAs( (Directory+"/"+Title+"/Overview_FitParameters.pdf").c_str() );
      delete paramCanvas;
      delete paramHisto;

      ///////////////////////////////////////////////////////////
      //  Save the DblGaus distributions using the fit params  //
      ///////////////////////////////////////////////////////////
      if(drawDoubleGaus == true){
      TCanvas* dblGausCanvas = new TCanvas("DblGausCanvas","Double Gaussian distribution using the fitted parameters");
      dblGausCanvas->Divide(4,4);

      TH1D* dblGausHisto;
      for(int iGaus = 0; iGaus < NrDblGaus; iGaus++){

        if(Title.find("GenInvE") <= Title.size()){
            dblGausHisto = (TH1D*) histoFile->Get( (Title+"/"+Title+"_DblGausPlot_GenE"+GenInvE[iGaus]).c_str());
            dblGausHisto->SetTitle( ("Double Gaussian distribution ("+histoTitle+") for InvGenE "+GenInvE[iGaus]).c_str());
        }
        else{
            dblGausHisto = (TH1D*) histoFile->Get( (Title+"/"+Title+"_DblGausPlot_GenE"+GenE[iGaus]).c_str());
            dblGausHisto->SetTitle( ("Double Gaussian distribution ("+histoTitle+") for GenE "+GenE[iGaus]).c_str());
        }
        dblGausHisto->GetXaxis()->SetTitle( histoAxisName.c_str() );
        dblGausHisto->GetXaxis()->SetTitleOffset(0.85);
        dblGausHisto->SetMarkerStyle(8);
        dblGausHisto->SetMarkerSize(0.3);

        //Draw the histograms in the different pads of the created canvas!
        dblGausCanvas->cd(iGaus+1);
        dblGausHisto->Draw("e");
      }
      dblGausCanvas->SaveAs( (Directory+"/"+Title+"/Overview_DblGausDistribution.pdf").c_str() );
      delete dblGausCanvas;
      delete dblGausHisto;
      }

    }//End of loop over NrFitHistos
    
    //twoDCanvas->SaveAs( (Directory+"/"+"ColorPlots.png").c_str() );
    twoDCanvas->SaveAs( (Directory+"/"+"ColorPlots.pdf").c_str() ); 
    delete twoDCanvas;
    delete twoDHisto;    
  }
}
