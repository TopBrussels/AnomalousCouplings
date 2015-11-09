# ------------------------------------------------------------------------------------ #
#  Script that will summarise the influence of the cut on the likelihood distribution  #
#   --> Will plot the distribution for both Correct, Unmatched and Wrong events        #
#   --> Will compare the cuts 60, 62, 65, 70 and no cut                                #
#   --> Will plot the minimum of the three together for each cut                       #
#   --> Will plot the likelihood value at gR = 0.0 for each type of events             #
#                                                                                      #
#  All output will automatically be stored in a PDF file!                              #
# ------------------------------------------------------------------------------------ #

import os
import sys

# Correct syntax is:
if len(sys.argv) <= 1:
    print "Will run all options"
    print "Additional syntax can be : python LikCutInfluence_RecoEvts.py (opt - int/-1)CutValue (opt - y/n/-1)ExclOuterBins (opt - Correct/Wrong/Unmatched/-1)whichReco"

LikCuts = [60, 62, 65, 70, 100]
ExclOuterBinsFit = ['y', 'n']
ExclOuterBinsFitStr = ["_OuterBinsExclForFit", ""]
whichReco = ["CorrectEvts", "WrongEvts", "UnmatchedEvts"]
if len(sys.argv) > 1 and str(sys.argv[1]) != '-1':
    LikCuts = [int(sys.argv[1])]
    print "Will only look at cut-value = ", sys.argv[1]
if len(sys.argv) > 2 and str(sys.argv[2]) != '-1':
    ExclOuterBinsFit = [str(sys.argv[2])]
    print "Will exclude outer bins from fit: ", sys.argv[2]
if len(sys.argv) > 3 and str(sys.argv[3]) != '-1':
    whichReco = [str(sys.argv[3])]
    print "Will only look at ", sys.argv[3], " events"

# **  Run over the different ROOT files and store the necessary histos  ** #
for iCut in range(len(LikCuts)):
    for iFit in range(len(ExclOuterBinsFit)):
        for iReco in range(len(whichReco)):

            # Go to the correct directory
            for recoDir in os.listdir('.'):
                if whichReco[iReco] in recoDir and not recoDir.endswith(".pdf"):

                    # Now get the correct ROOT file
                    for File in os.listdir(recoDir):
                        if File.startswith("FitDistributions") and File.endswith(".root"):
                            if ( (str(LikCuts[iCut]) != "100" and 'LikelihoodCut'+str(LikCuts[iCut]) in File) or (str(LikCuts[iCut]) == "100" and not 'LikelihoodCut' in File)) and ((ExclOuterBinsFit[iFit] == 'y' and 'OuterBinsExclForFit' in File) or (ExclOuterBinsFit[iFit] == 'n' and not 'OuterBinsExclForFit' in File)):
                                ROOTFile = recoDir + '/' + File
                                # print "Looking at ROOT file : ", ROOTFile, "for iCut = ", LikCuts[iCut], ", iFit = ", ExclOuterBinsFit[iFit], ' and whichReco = ', whichReco[iReco]

                                # Now open the ROOT file and store the histograms
                                SaveHist = open('SaveHist.C', 'w')
                                SaveHist.write('{ \n  TFile* file = new TFile("'+ROOTFile+'","READ"); \n')
                                SaveHist.write('  TCanvas* canv = new TCanvas("canv","canv"); \n  canv->cd(); \n\n')
                                SaveHist.write('  TH1F* hist_Sum = (TH1F*)file->Get("FitSums/SummedHist_Acc"); \n')
                                if str(LikCuts[iCut]) == '100' or len(LikCuts) == 1:
                                    SaveHist.write('  hist_Sum->Draw(); \n  canv->SaveAs("SummedHist_'+str(whichReco[iReco])+str(ExclOuterBinsFitStr[iFit])+'_LikCut'+str(LikCuts[iCut])+'.pdf"); \n\n')
                                    SaveHist.write('  TH1F* SMLik = (TH1F*)file->Get("FitResults/SMLikelihoodValue"); \n')
                                    SaveHist.write('  SMLik->Draw();\n  canv->SaveAs("SMLikelihoodValue_'+str(whichReco[iReco])+'.pdf"); \n\n')
                                else:
                                    SaveHist.write('  hist_Sum->Draw(); \n  canv->SaveAs("SummedHist_'+str(whichReco[iReco])+str(ExclOuterBinsFitStr[iFit])+'_LikCut'+str(LikCuts[iCut])+'.pdf"); \n\n')
                                SaveHist.write('}\n')
                                SaveHist.close()

                                os.system('root -l -b -q SaveHist.C+')

# LatexOutput = open('LikCutInfluence_RecoEvts.tex', 'w')
