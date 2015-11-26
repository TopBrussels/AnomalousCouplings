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

LikCuts = [60, 62, 64, 65, 66, 67, 68, 69, 70, 100]
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
rootFileList = ""
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
                                if whichReco[iReco] == 'CorrectEvts' and ExclOuterBinsFitStr[iFit] == "_OuterBinsExclForFit":
                                    rootFileList += str(ROOTFile) + ' '

                                # Now open the ROOT file and store the histograms
                                SaveHist = open('SaveHist.C', 'w')
                                SaveHist.write('{ \n  TFile* file = new TFile("'+ROOTFile+'","READ"); \n')
                                SaveHist.write('  TCanvas* canv = new TCanvas("canv","canv"); \n  canv->cd(); \n\n')
                                SaveHist.write('  TH1F* hist_Sum = (TH1F*)file->Get("FitSums/SummedHist_Acc"); \n')
                                if str(LikCuts[iCut]) == '100' or len(LikCuts) == 1:
                                    SaveHist.write('  TH1F* SMLik = (TH1F*)file->Get("FitResults/SMLikelihoodValue"); \n')
                                    SaveHist.write('  SMLik->Draw();\n  canv->SaveAs("SMLikelihoodValue_'+str(whichReco[iReco])+'.pdf"); \n\n')

                                if str(LikCuts[iCut]) == '100':
                                    SaveHist.write('  hist_Sum->Draw(); \n  canv->SaveAs("SummedHist_'+str(whichReco[iReco])+str(ExclOuterBinsFitStr[iFit])+'_NoLikCut.pdf"); \n\n')
                                else:
                                    SaveHist.write('  hist_Sum->Draw(); \n  canv->SaveAs("SummedHist_'+str(whichReco[iReco])+str(ExclOuterBinsFitStr[iFit])+'_LikCut'+str(LikCuts[iCut])+'.pdf"); \n\n')
                                SaveHist.write('}\n')
                                SaveHist.close()

                                os.system('root -l -b -q SaveHist.C+')

print "Full list of Root files is : ", rootFileList

consCuts = ' '
for cuts in range(len(LikCuts)):
    consCuts += str(LikCuts[cuts]) + ' '

os.system('g++ -m64 -g -I `root-config --incdir` `root-config --libs` LikCut_MinimumEvolution.C -o LikCut_MinimumEvolution')
os.system('./LikCut_MinimumEvolution ' + str(len(LikCuts)) + ' ' + rootFileList + ' ' + consCuts)

# Output ROOT file
LatexOutput = open('LikCutInfluence_RecoEvts.tex', 'w')
LatexOutput.write('\\documentclass[a4paper]{article} \n')
LatexOutput.write('\usepackage{graphicx} \n\usepackage{geometry}\n')
LatexOutput.write('\\begin{document} \n ')
LatexOutput.write(' \\begin{abstract} \n  Summary of the influence of the cut on the MadWeight likelihood distribution for reco-level events.\n')
LatexOutput.write('  Will compare the difference in obtained minimum for both correct, wrong and unmatched events. \n')
LatexOutput.write(' \\end{abstract}\n')

LatexOutput.write(' \\section{Original result before any cut} \n')
LatexOutput.write('   Original result for the reconstructed events indicated that the shape of the MadWeight likelihood (-ln($\\mathcal{L}$)) did not correspond with the expectation of a minimum close to $g_R$ = 0.0. \n')
LatexOutput.write('   The distributions were closest related to a descending straight line with minimum at the border of the considered range, which does not match at all with the expected distribution. \n')
LatexOutput.write('   These distributions can be found in the following figure (Fig \\ref{fig::OrigResult} and \\ref{fig::OrigResult_AllBins}) for both the correct, wrong and unmatched reco-events. \n')
LatexOutput.write('   The only difference between the two figures is the exclusion of the outer bin or not.')
LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_CorrectEvts_OuterBinsExclForFit_NoLikCut.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_WrongEvts_OuterBinsExclForFit_NoLikCut.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_UnmatchedEvts_OuterBinsExclForFit_NoLikCut.pdf} \n')
LatexOutput.write('    \\caption{MadWeight likelihood distributions (-ln($\\mathcal{L}$) for correct, wrong and unmatched reco-level events when the outer bins are excluded from the fit. This because the last bin lies rather far away from the expectations and is mainly added to avoid having an important value on the edge of the considered range.} \n')
LatexOutput.write('    \\label{fig::OrigResult} \n')
LatexOutput.write('   \\end{figure} \n \n')

LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_CorrectEvts_NoLikCut.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_WrongEvts_NoLikCut.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_UnmatchedEvts_NoLikCut.pdf} \n')
LatexOutput.write('    \\caption{MadWeight likelihood distributions (-ln($\\mathcal{L}$) for correct, wrong and unmatched reco-level events when also the outer bins are used in the fit.} \n')
LatexOutput.write('    \\label{fig::OrigResult_AllBins} \n')
LatexOutput.write('   \\end{figure} \n \n')

LatexOutput.write(' \\section{Issue: Shape dominated by badly reconstructed events}')
LatexOutput.write('   Studying the individual events clearly showed that there were a rather large number of events for which the MadWeight likelihood value was higher than average. \n')
LatexOutput.write('   Since this corresponds to calculated weights which are very small, it is rather likely that these are actually events which had trouble to be correctly determined. \n')
LatexOutput.write('   However due to their high Likelihood value they contribute excessively to the overall shape when summing all the individual likelihoods. \n')
LatexOutput.write('   Therefore the influence of excluding events with such an extremely high likelihood value has been studied, especially since these does not seem to occur for generator-level events explaining the straightforward minimum determination for the latter type of events. \n')
LatexOutput.write('   The MadWeight Likelihood distribution evaluated in the point $g_R$ = 0.0 for both correct, wrong and unmatched events can be found in Figure \\ref{fig::SMLik}. \n')
LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SMLikelihoodValue_CorrectEvts.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SMLikelihoodValue_WrongEvts.pdf} \n')
LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SMLikelihoodValue_UnmatchedEvts.pdf} \n')
LatexOutput.write('    \\caption{Distribution of the MadWeight likelihood for $g_R$ = 0.0, which gives an idea of the smallness of the individually calculated weight for the considered event since the MadWeight likelihood is in general rather flat. Events with a very high likelihood value actually had a very small weight in the calculations and are therefore rather likely to be miscalculated.} \n')
LatexOutput.write('    \\label{fig::SMLik} \n')
LatexOutput.write('   \\end{figure} \n \n')

LatexOutput.write('  In the following subsections the shape of the MadWeight likelihood distribution will be given for the different considered cut-values (60, 62, 65 and 70). \n')
LatexOutput.write('  As was the case for the originally obtained distributions the first one will contain the result retrieved when the outer bins are excluded from the polynomial fit. \n')

for likCut in range(len(LikCuts)):
    if str(LikCuts[likCut]) != '100':
        LatexOutput.write('  \\newpage \n')
        LatexOutput.write('  \\subsection{Result after cut = '+str(LikCuts[likCut])+'} \n')
        LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_CorrectEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_WrongEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_UnmatchedEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\caption{MadWeight likelihood distributions (-ln($\\mathcal{L}$)) for correct, wrong and unmatched reco-level events after a likelihood cut of ' + str(LikCuts[likCut]) + ' (outer bins excluded from fit).} \n')
        LatexOutput.write('    \\label{fig::Cut' + str(LikCuts[likCut]) + '} \n')
        LatexOutput.write('   \\end{figure} \n \n')

        LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_CorrectEvts_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_WrongEvts_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\includegraphics[width = 0.32 \\textwidth]{SummedHist_UnmatchedEvts_LikCut' + str(LikCuts[likCut]) + '.pdf} \n')
        LatexOutput.write('    \\caption{MadWeight likelihood distributions (-ln($\\mathcal{L}$)) for correct, wrong and unmatched reco-level events after a likelihood cut of ' + str(LikCuts[likCut]) + ' (outer bins used in  fit).} \n')
        LatexOutput.write('    \\label{fig::Cut' + str(LikCuts[likCut]) + '_AllBins} \n')
        LatexOutput.write('   \\end{figure} \n \n')

LatexOutput.write(' \\section{Summary} \n')
LatexOutput.write('   \\begin{figure}[h!t] \n    \\centering \n')
LatexOutput.write('    \\includegraphics[width = 0.75 \\textwidth]{MinimumPosition.pdf} \n')
LatexOutput.write('    \\caption{Position of the minimum for the different cut-values (correct reco events with outer bins excluded from fit).} \n')
LatexOutput.write('    \\label{fig::MinPos} \n')
LatexOutput.write('   \\end{figure} \n \n')

LatexOutput.write('\\end{document} \n')
LatexOutput.close()

os.system('pdflatex -interaction=batchmode '+ LatexOutput.name)
os.system('rm '+ LatexOutput.name[:-4] + ".aux")
os.system('rm '+ LatexOutput.name[:-4] + ".log")

for LikCut in range(len(LikCuts)):
    if str(LikCuts[LikCut]) != '100':
        os.system('rm SummedHist_CorrectEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[LikCut]) + '.pdf')
        os.system('rm SummedHist_WrongEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[LikCut]) + '.pdf')
        os.system('rm SummedHist_UnmatchedEvts_OuterBinsExclForFit_LikCut' + str(LikCuts[LikCut]) + '.pdf')

        os.system('rm SummedHist_CorrectEvts_LikCut' + str(LikCuts[LikCut]) + '.pdf')
        os.system('rm SummedHist_WrongEvts_LikCut' + str(LikCuts[LikCut]) + '.pdf')
        os.system('rm SummedHist_UnmatchedEvts_LikCut' + str(LikCuts[LikCut]) + '.pdf')