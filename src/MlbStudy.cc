#include "../interface/MlbStudy.h"

//MlbStudy::MlbStudy(){
//}

//MlbStudy::~MlbStudy(){

//}

void MlbStudy::initializePerEvent(){
 
  for(int ii = 0; ii < 6; ii++){
     ChiSquared[ii] = 0.;
     MlbValues[ii] = 0.;
     MqqbValues[ii] = 0.;
  }
}

void MlbStudy::initializeBegin(){

}

void MlbStudy::calculateChiSquared(vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets, TLorentzVector* lepton, vector<TRootJet*> Jets, float MassMlb, float SigmaMlb, float MassMqqb, float SigmaMqqb){

   if(bTaggedJets.size() > 1 && lightJets.size() > 1){                                            //Event needs to have two b-tagged jets and two light jets!

       float MlbValuesLocal[6]  = {(*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M(),
			      (*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M(),
			      (*lepton+*Jets[bTaggedJets[0]]).M(), (*lepton+*Jets[bTaggedJets[1]]).M()};

       float MqqbValuesLocal[6] = {(*Jets[bTaggedJets[1]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(), (*Jets[bTaggedJets[0]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(),
	         	       10000., 10000.,10000., 10000.};
       if(lightJets.size() > 2){  //5-jet case!
         float MqqbValuesLocal[6] = {(*Jets[bTaggedJets[1]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(), (*Jets[bTaggedJets[0]] + *Jets[lightJets[0]] + *Jets[lightJets[1]]).M(),
	         	             (*Jets[bTaggedJets[1]] + *Jets[lightJets[0]] + *Jets[lightJets[2]]).M(), (*Jets[bTaggedJets[0]] + *Jets[lightJets[0]] + *Jets[lightJets[2]]).M(),
				     (*Jets[bTaggedJets[1]] + *Jets[lightJets[1]] + *Jets[lightJets[2]]).M(), (*Jets[bTaggedJets[0]] + *Jets[lightJets[1]] + *Jets[lightJets[2]]).M()};
       }

       //Copy events to general MlbValues and MqqbValues array!
       for(int ii = 0; ii < 6; ii++){
	  MlbValues[ii] = MlbValuesLocal[ii];
	  MqqbValues[ii]= MqqbValuesLocal[ii];
       }

       LowestChiSq = 0;

	std::cout << " MqqbValues in class : " << MqqbValues[0] << ", " << MqqbValues[1] << ", " << MqqbValues[2] << ", " << MqqbValues[3] << ", " << MqqbValues[4] << ", " << MqqbValues[5] << std::endl;

       for(int ii = 0; ii<6; ii++){
           ChiSquared[ii] = ((MassMlb-MlbValues[ii])/SigmaMlb)*((MassMlb-MlbValues[ii])/SigmaMlb) + ((MassMqqb - MqqbValues[ii])/SigmaMqqb)*((MassMqqb - MqqbValues[ii])/SigmaMqqb);
	   if(ii == 0) LowestChiSq = ChiSquared[0];
	   if( ii > 0 && ChiSquared[ii] < LowestChiSq)
	      LowestChiSq = ii; 
       }

       //Match the correct indices to the quarks!
       chosenBLept=0;           //This means that the hightest Pt jet is considered as the leptonic b-jet! 
       chosenBHadr=1;
       int chosenQuarks[2] = {0,1};

       //if( LowestChiSq == 0 ||2 || 4


   }
   else{
	chosenBLept = 999;
	chosenBHadr = 999;
	ChiSquared[0] = 999;
	ChiSquared[1] = 999;
   }
}

void MlbStudy::calculateEfficiency(int option, vector<int> CorrectValues, vector<int> bTaggedJets, vector<int> lightJets){


   if(  CorrectValues[0] != 9999 && CorrectValues[1] != 9999 && CorrectValues[2] != 9999 && CorrectValues[3] != 9999 &&   //Event has been matched!
	bTaggedJets.size() > 1 && lightJets.size() > 1){                                                                  //Event has the correct amount of b-tagged and light jets!

	NumberMatchedEvents[option]++;

        //How often is the correct option chosen?
        if( (CorrectValues[0] == bTaggedJets[0] || CorrectValues[0] == bTaggedJets[1]) &&
            (CorrectValues[1] == bTaggedJets[1] || CorrectValues[1] == bTaggedJets[0])   ){  //Correct option is available!
		CorrectOptionAvailable[option]++;

		if( CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr] ){
			CorrectOptionChosen[option]++;
			h_CorrectOptionChiSq.push_back(ChiSquared[chosenBLept]);
		}
		else{
			WrongOptionChosen[option]++;
			h_WrongOptionChiSq.push_back(ChiSquared[chosenBHadr]);
		}
	}

       //How often is the full event found?
       if( CorrectValues[0] == bTaggedJets[chosenBLept] && CorrectValues[1] == bTaggedJets[chosenBHadr] &&
           (CorrectValues[2] == lightJets[0]   || CorrectValues[2] == lightJets[1])   &&
           (CorrectValues[3] == lightJets[0]   || CorrectValues[3] == lightJets[1])     ){
                CorrectEventMlbMqqb[option]++;
       }
       else{
                WrongEventMlbMqqb[option]++;
       }
   }
   else{
	NumberNotMatchedEvents[option]++;
   }
}

void MlbStudy::saveNumbers(std::string NameOfOption[6], int WhichJets, int NrOptionsConsidered, ofstream &output){

   std::string Title[2] = {" \\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & all 4 correct      & $\\geq$ 1 wrong  & $\\frac{s}{\\sqrt{b}}$ & $\\frac{s}{b}$ & non-matched \\\\", 
			   " \\textbf{Option} (with $\\chi^{2}$ $m_{lb}$) & Correct b's chosen & Wrong b-s chosen & $\\frac{s}{\\sqrt{b}}$ & $\\frac{s}{b}$ & Events with correct option available \\\\"};

   std::string Caption[2] = {" \\caption{Overview of correct and wrong reconstructed events for the different b-tags when a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method is applied} ", 
			      " \\caption{Overview of the number of times the correct b-jet combination is chosen when using a $\\chi^{2}$ $m_{lb}$ - $m_{qqb}$ method} "};

   output << " \\begin{table}[!h] \n \\begin{tabular}{c|c|c|c|c|c} " << endl;
   output << Title[WhichJets] << " \\hline " << endl;
   for(int ii = 0; ii < NrOptionsConsidered; ii++){

       int CorrectOnes[2] = {CorrectEventMlbMqqb[ii],   CorrectOptionChosen[ii] };
       int WrongOnes[2]   = {WrongEventMlbMqqb[ii],     WrongOptionChosen[ii]};
       int LastOnes[2]    = {NumberNotMatchedEvents[ii],CorrectOptionAvailable[ii]}; 

       output << NameOfOption[ii]                                   << " & " <<
       CorrectOnes[WhichJets]                                          << " & " <<
       WrongOnes[WhichJets]                                            << " & " <<
       float(CorrectOnes[WhichJets])/float(sqrt(WrongOnes[WhichJets])) << " & " <<
       float(CorrectOnes[WhichJets])/float(WrongOnes[WhichJets])       << " & " <<
       LastOnes[WhichJets]                                       << " \\\\ " << endl;
     }
    output << " \\end{tabular} " << endl;
    output << Caption[WhichJets] << endl;
    output << " \\end{table} \n " << endl;

}

