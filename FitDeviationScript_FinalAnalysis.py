# ########################################################
#                                                       ##
#  Python macro which takes care of fitting the -ln(L)  ##
#  Uses the following ROOT macros:                      ##
#   - fitDeviationMacro.C                               ##
#   - PerformFitOptimization.C                          ##
#                                                       ##
#  ** First one performs 2 consecutive fits             ##
#       * First one using all points                    ##
#       * Second one use 66% of best points             ##
#  ** Second one allows for quick evt sel tests         ##
#       * Uses the TF1's created in first macro         ##
#       * Allows to test different chi-sq cuts          ##
#                                                       ##
# #########################################################

# ! python
import os
import sys
import re
import shutil
from array import array

# Get all the input from the command line:
if len(sys.argv) <= 1:
    print "Need to give the directory of interest, the type of sample and the number of events in command line !"
    print " Correct syntax : python FitDeviation.py Events/blabla/ MGorRECO(MG/RECO)"
    sys.exit()
elif len(sys.argv) == 2:
    print "Need to specify the type of sample (MG or RECO)"
    print " Correct syntax : python FitDeviation.py Events/blabla/ MGorRECO(MG/RECO)"
    sys.exit()

whichDir = sys.argv[1]
MGorRECO = sys.argv[2]


# -----------------------------------------#
#  Set all the optional parameters correct #
# -----------------------------------------#
# Full syntax is : python FitDeviation.py Events/blabla/ MGorRECO #evts(opt) Range(opt) weightsFile(opt) likCutValue(opt) excludeOuterBinsFit(y/n - opt)  applyAcc(y/n - opt) TexWanted(y/n - opt)
YesOptions = ['y', 'yes', 'Y', 'YES', 'Yes']

nEvts = "-1"
if len(sys.argv) > 3:
    if sys.argv[3] != "-1":
        nEvts = int(sys.argv[3])
print "Number of events should be : ", nEvts

VarWindowGiven = False
if len(sys.argv) > 4:
    VarWindowGiven = True
    if sys.argv[4] == "Normal":
        VarWindow = "1"
    elif sys.argv[4] == "Wide":
        VarWindow = "2"
    elif sys.argv[4] == "Middle":
        VarWindow = "3"

WeightsFileGiven = False
if len(sys.argv) > 5:
    WeightsFileGiven = True
    WeightsFileName = sys.argv[5]
    WeightsFile = open(os.path.join(WeightsFileName), 'r')

# Give the cut-value on the -ln(Likelihood) as an input
if len(sys.argv) > 6:
    LikCut = sys.argv[6]
else:
    LikCut = 100

# Specify whether the cuts should exclude the two outer bins!
excludeOuterBinsFit = False
if len(sys.argv) > 7:
    if sys.argv[7] in YesOptions:
        excludeOuterBinsFit = True

applyAccNorm = "n"
if MGorRECO == "RECO" or MGorRECO == "DATA":
    applyAccNorm = "y"
    print "\n To notice: For RECO and DATA events acceptance is always applied!!    **********"
elif MGorRECO == "MG":
    if len(sys.argv) > 8:
        applyAccNorm = sys.argv[8]

# Set the 'CreateTexFile' correctly:
CreateTexFile = False
if len(sys.argv) > 9:
    if sys.argv[9] in YesOptions:
        CreateTexFile = True

# Use this boolean to calculate the acceptance XS by scaling the MadGraph XS with the SM-acceptance ratio
UseScaledXS = False
if len(sys.argv) > 10:
    if sys.argv[10] in YesOptions:
        UseScaledXS = True

applyCosTheta = "n"
#if len(sys.argv) > 4:
#    applyCosTheta = sys.argv[4]

# --------------------------------------------------------------#
#  Now continue by selecting the desired configs form the range #
# --------------------------------------------------------------#
ValuesToDelete = []
if not VarWindowGiven:
    print '\n ** Choose the correct RgR-window corresponding to the studied file ** \n** Different options are: '
    print ' 1) Normal : [-0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20] '
    print ' 2) Wide   : [-0.4, -0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4] '
    print ' 3) Middle : [-0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3]'
    VarWindow = raw_input(' --> Choose the correct number : ')
Var = array('d', [-0.4, -0.3, -0.2, -0.15, -0.1, -0.05, -0.025, 0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4])

xMinValue = [4, 6, 5]
KinVar = "Re(g_{R})"
NrPointsToRemove = [2, 4, 2]

if VarWindow == "1":
    ValuesToDelete = [-0.4, -0.3, -0.025, 0.025, 0.3, 0.4]
    NrBins = 8

elif VarWindow == "2":
    ValuesToDelete = [-0.025, 0.025]
    NrBins = 16

elif VarWindow == "3":
    ValuesToDelete = [-0.4, -0.025, 0.025, 0.4]
    NrBins = 12

OnlyXSPossible = True
if MGorRECO == "MG":
    MGXS = array('d',      [3.95248, 5.5612,  8.24066,  10.09161, 12.39876, 15.16567, 0.0, 18.54042, 0.0, 22.52782, 27.3093, 32.9283, 39.4799, 55.9507,  77.7365])
    MGXSCut = array('d',   [0.93159, 1.27966, 1.825208, 2.194079, 2.6393,   3.17698,  0.0, 3.80921,  0.0, 4.5645,   5.45665, 6.47791, 7.66805, 10.63243, 14.46786])  # Also MET cuts!
elif MGorRECO == "RECO" or MGorRECO == "DATA":
    MGXS = array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #MGXSCut = array('d',[0.0, 0.657583, 0.948253, 1.13713, 1.36514, 1.63988, 0.0, 1.96892, 0.0, 2.35987, 2.8203,  3.3578,  3.97995, 5.50855, 0.0])  # Up to 1 additional jet
    MGXSCut = array('d',[0.0, 0.65669,  0.947244, 1.13624, 1.36448, 1.63952, 0.0, 1.96892, 0.0, 2.36027, 2.82111, 3.35903, 3.98157, 5.51083, 0.0])  # Up to 2 additional jets
    OnlyXSPossible = False
else:
    print "\n Should specify whether MG or RECO is desired! "
    sys.exit()

# Now delete the values stored in array 'ValuesToDelete'
for iVar in range(len(ValuesToDelete)):
    MGXS.pop(int(Var.index(ValuesToDelete[iVar])))
    MGXSCut.pop(Var.index(ValuesToDelete[iVar]))
    if UseScaledXS: MGXSScale.pop(Var.index(ValuesToDelete[iVar]))
    Var.pop(Var.index(ValuesToDelete[iVar]))

# Set all the parameters accordingly
NrConfigs = len(Var)
xMin = xMinValue[int(VarWindow) - 1]
NumberOfPointsToRemove = NrPointsToRemove[int(VarWindow) - 1]
FitMin, FitMax = Var[0], Var[NrConfigs - 1]
if excludeOuterBinsFit:
    FitMin, FitMax = Var[1], Var[NrConfigs - 2]

# ---------------------------------------#
#  Open the correct weight file          #
#   --> Only if WeightsFileGiven = false #
# ---------------------------------------#
if not WeightsFileGiven:
    list_dir = os.listdir(whichDir)
    WeightsFileArray, weightsFileCounter = [], 0
    for file in list_dir:
        if (applyCosTheta == "n" and not "CosTheta" in file or applyCosTheta == "y" and "_CosTheta" in file) and file.startswith("weights") and file.endswith(".out"):
            weightsFileCounter += 1
            WeightsFileArray.append(file)

    if int(weightsFileCounter) == 0:
        print "\n ERROR: No weights file found in this directory !"
        sys.exit()
    elif int(weightsFileCounter) == 1:
        WeightsFileName = str(whichDir) + '' + str(WeightsFileArray[0])
    elif int(weightsFileCounter) > 1:
        print "\n** Still need to choose the weights file of interest! ** \n** Different options are:"
        for ii in range(len(WeightsFileArray)):
            print " ", ii, " ) ", WeightsFileArray[ii]
        fileNr = raw_input(' --> Choose the number of the desired file! : ')
        WeightsFileName = str(whichDir) + '' + str(WeightsFileArray[int(fileNr)])

# If _SFAdded is in the name of the used weights file, this should also be applied in the analyzer!
applySF = False
if '_SFAdded' in str(WeightsFileName) and MGorRECO == "RECO":
    applySF = True

# Open the selected files!
WeightsFile, LikelihoodFile = open(WeightsFileName, 'r'), open(WeightsFileName, 'r')

# Set the title using the directory name and the considered weights file
title = str(whichDir[whichDir.find("/") + 1:-1])  # Only need the part after the "/"!!
title += WeightsFileName[WeightsFileName.find("weights") + 7: -4]

# Count the number of events actually present in the WeightsFile as the maximum numbers which can actually be used:
# --> Check whether the file has enough events, otherwise use the maximum number
maxNrEvts = os.popen('grep " 1 1 " ' + str(WeightsFileName) + ' | wc -l').read()
if nEvts == "-1" or int(maxNrEvts) < int(nEvts):
    nEvts = int(maxNrEvts)

print "\n----- Summary ------ "
print " Will be using file : ", WeightsFileName
print " Will consider ", nEvts, " events ! "
if str(LikCut) != '100':
    print " Will apply cut on -ln(Likelihood) of ", str(LikCut)
else:
    print ' Will apply no cut on -ln(Likelihood) '
print " Will fit between ", FitMin, " and ", FitMax, " (full range = [", Var[0], ",", Var[NrConfigs-1], "] ) "

# --------------------------------------------------------------------------------------#
#   Set the MGXSCut values based on whether acceptance normalisation should be applied  #
# --------------------------------------------------------------------------------------#
if applyAccNorm == "n" or applyAccNorm == "no" or applyAccNorm == "No":
    MGXSCut = MGXS
    print " No acceptance normalisation will be applied! "
    if (len(whichDir) >= whichDir.find("RECO") > 0) or (len(whichDir) >= whichDir.find("Reco") > 0):
        print " \n ************************* ERROR **************************** "
        print " -----> Applying no Acceptance normalisation for reco-level events .... \n"
else:
    print " Will apply acceptance normalisation! "
    if (len(whichDir) >= whichDir.find("MGSample") > 0 > whichDir.find("Cut")) or (
                    len(whichDir) >= whichDir.find("GEN") > 0):
        print " \n ************************* ERROR **************************** "
        print " -----> Applying Acceptance normalisation to generator-level events .... \n"
if UseScaledXS and MGorRECO =="MG":
    MGXSCut = MGXSScale

# ---------------------------#
#   Setting output title     #
# ---------------------------#
# *** Indicating that the cos theta* reweighting has been applied
if applyCosTheta == "y" or applyCosTheta == "yes" or applyCosTheta == "Y":
    if nEvts != int(maxNrEvts):
        print "\n *** ERROR: Since cos theta weights have to be normalized, this can only be applied if ALL events are considered ***"
        cont = raw_input(' *** Decide whether the script should be continued, but then using all events (y).\n *** If (n) is chosen script will be terminated. ')
        if cont == "y":
            nEvts = int(maxNrEvts)
            print " --> Will now be using ", nEvts
        else:
            sys.exit()

# *** Indicating whether cut on Likelihood is applied!
if str(LikCut) != '100':
    title = title + '_LikelihoodCut' + str(LikCut)

# *** Indicating whether fit range has been reduced!
if excludeOuterBinsFit:
    title += '_OuterBinsExclForFit'

# *** Indicating that XS-values are changed!
if len(whichDir) >= whichDir.find("ChangingXS") > 0:
    print "\n ******************** IMPORTANT : XS VALUE WILL BE SCALED !!!!  ************************** \n"
    MGXSScaled, MGXSCutScaled = [], []
    ScaleValue, TitleChange = 0.05, "_XSScaledWithPos005"
    title += TitleChange
    for iVar in range(len(MGXS)):
        MGXSScaled.append(MGXS[iVar] * (1 + Var[iVar] * ScaleValue))
        MGXSCutScaled.append(MGXSCut[iVar] * (1 + Var[iVar] * ScaleValue))
    MGXS, MGXSCut = MGXSScaled, MGXSCutScaled

# *** Indicating that Pt-cuts have been applied!
if len(WeightsFileName) >= WeightsFileName.find("NoLowPt") >= 0:
    title = title + "_" + "NoLowPtEvts" + WeightsFileName[WeightsFileName.find("_Cut"):-4]

# *** Indicating that scaled acceptance XS have been used
if UseScaledXS and MGorRECO =="MG":
    title = title + "_" + "XSScaledWithSMAcc"

# -------------------------------------------------#
# --  Pass on all variables to the ROOT macro !  --#
# -------------------------------------------------#
RootAnalyzer, NewRootAnalyzer = open('doublePolFitMacro.C', 'r'), open('output.C', 'w')

NormType = ["no", "XS", "acceptance"]
NormTypeName = ["", "_XS", "_Acc"]
if not OnlyXSPossible:
    NormType.pop(1), NormTypeName.pop(1)
NormTypeC = 'std::string NormType[' + str(len(NormType)) + '] = {'
NormTypeNameC = 'std::string NormTypeName[' + str(len(NormTypeName)) + '] = {'
for ii in range(len(NormType)):
    NormTypeC += '"' + NormType[ii] + '"'
    NormTypeNameC += '"' + NormTypeName[ii] + '"'
    if ii < len(NormType) - 1:
        NormTypeC += ','
        NormTypeNameC += ','
    else:
        NormTypeC += '}; \n'
        NormTypeNameC += '}; \n'

VarLine = 'double Var[] = {'
MGXSLine, MGXSCutLine = 'double MGXS[] = {', 'double MGXSCut[] = {'
for ii in range(len(Var)):
    VarLine, MGXSLine, MGXSCutLine = VarLine + str(Var[ii]), MGXSLine + str(MGXS[ii]), MGXSCutLine + str(MGXSCut[ii])
    if ii < len(Var) - 1:
        VarLine, MGXSLine, MGXSCutLine = VarLine + ',', MGXSLine + ',', MGXSCutLine + ','
    else:
        VarLine, MGXSLine, MGXSCutLine = VarLine + '};\n', MGXSLine + '};\n', MGXSCutLine + '};\n'

# -->Create the directory SplittedCanvasses if needed (otherwise delete the created pdf files ...)!
if CreateTexFile:
    if not (os.path.exists(os.path.join(whichDir + 'SplittedCanvasses/'))):
        os.makedirs(os.path.join(whichDir + 'SplittedCanvasses/'))
    else:
        print "Deleting the existing SplitCanvas files ! "
        os.system('rm ' + whichDir + 'SplittedCanvasses/SplitCanvasLL*.pdf')

for RootLine in RootAnalyzer:
    RootWord = RootLine.split()

    if re.search(r"int nEvts", RootLine):
        NewRootAnalyzer.write('const int nEvts = ' + str(nEvts) + '; \n')
    elif re.search(r"int NrConfigs", RootLine):
        NewRootAnalyzer.write('const int NrConfigs = ' + str(NrConfigs) + '; \n')
    elif re.search(r"double Var", RootLine):
        NewRootAnalyzer.write(VarLine)
    elif re.search(r"double MGXSCut", RootLine):
        NewRootAnalyzer.write(MGXSCutLine)
    elif re.search(r"double MGXS", RootLine):
        NewRootAnalyzer.write(MGXSLine)
    elif re.search(r"int NrBins", RootLine):
        NewRootAnalyzer.write('int NrBins = ' + str(NrBins) + '; \n')
    elif re.search(r"double FitMin", RootLine):
        NewRootAnalyzer.write('double FitMin = ' + str(FitMin) + ', FitMax = ' + str(FitMax) + '; \n')
    elif re.search(r"double LikCut", RootLine):
        NewRootAnalyzer.write('double LikCut = ' + str(LikCut) + '; \n')
    elif re.search(r"string KinVar", RootLine):
        NewRootAnalyzer.write('std::string KinVar = "' + str(KinVar) + '"; \n')
    elif re.search(r"int NrToDel", RootLine):
        NewRootAnalyzer.write('const unsigned int NrToDel = ' + str(NumberOfPointsToRemove) + '; \n')
    elif re.search(r"std::string NormTypeName", RootLine):
        NewRootAnalyzer.write(NormTypeNameC)
    elif re.search(r"std::string NormType", RootLine):
        NewRootAnalyzer.write(NormTypeC)
    elif re.search(r"std::ifstream ifs", RootLine):
        NewRootAnalyzer.write('  std::ifstream ifs ("' + str(WeightsFileName) + '", std::ifstream::in); \n')
    elif re.search(r"std::string title", RootLine):
        NewRootAnalyzer.write('std::string title = "' + str(title) + '"; \n')
    elif re.search(r"string SplittedDir", RootLine):
        NewRootAnalyzer.write('std::string SplittedDir = "' + str(whichDir) + 'SplittedCanvasses"; \n')
    elif re.search(r"iss >> evt", RootLine):
        if applyCosTheta == "y" or applyCosTheta == "Y":
            NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> weightUnc >> CosThetaCorr ){ \n')
            print " Cos theta* reweighting will be applied! \n"
        elif applySF == True and MGorRECO == "RECO":
            NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> weightUnc >> MCScaleFactor >> Luminosity >> NormFactor >> SampleName){ \n')
            print " MC Scale Factor will be applied! \n"
        elif MGorRECO == "DATA":
            NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> weightUnc >> SampleName){ \n')
        else:
            NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> weightUnc){ \n')
    elif re.search(r"bool storeSplittedCanvas", RootLine):
        if CreateTexFile:
            NewRootAnalyzer.write('bool storeSplittedCanvas = true; \n')
        else:
            NewRootAnalyzer.write('bool storeSplittedCanvas = false; \n')
    elif re.search(r"new TFile", RootLine) and re.search(r"file_FitDist", RootLine):
            NewRootAnalyzer.write('TFile* file_FitDist = new TFile("' + str(whichDir) + 'FitDistributions_' + str(title) + '_' + str(nEvts) + 'Evts.root","RECREATE"); \n')
    else:
        NewRootAnalyzer.write(RootLine)
NewRootAnalyzer.close(), RootAnalyzer.close()

# Run the root macro!
os.rename('output.C', 'doublePolFitMacro.C'), os.system("root -l -b -q doublePolFitMacro.C+")

# -- Now store the stacked canvasses in a .txt file --#
if CreateTexFile:

    # Change the working directory of the script!
    os.chdir(os.path.join(whichDir + 'SplittedCanvasses/'))
    # CanvasOutputFile_NoNorm = open(os.path.join('FitDeviationSplitCanvas_' + str(title) + '_' + str(nEvts) + 'Evts_NoNorm.tex'), 'w')
    # CanvasOutputFile_XSNorm = open(os.path.join('FitDeviationSplitCanvas_' + str(title) + '_' + str(nEvts) + 'Evts_XSNorm.tex'), 'w')
    # CanvasOutputFile_AccNorm = open(os.path.join('FitDeviationSplitCanvas_' + str(title) + '_' + str(nEvts) + 'Evts_AccNorm.tex'), 'w')

    print " Current working directory is : ", os.getcwd()
    # CanvasOutputFile = [CanvasOutputFile_NoNorm, CanvasOutputFile_XSNorm, CanvasOutputFile_AccNorm]
    # if not OnlyXSPossible:
    #    CanvasOutputFile.pop(1)

    # Store the information in the correct directory:
    Canvaslist_dir = os.listdir('.')  # os.path.join(whichDir+'SplittedCanvasses/'))

    OtherNorms = []
    for iNormType in range(len(NormType)):
        for ii in range(len(NormType)):
            if ii != iNormType:
                OtherNorms.append(ii)

        CanvasOutputFile = open(os.path.join('FitDeviationSplitCanvas_' + str(title) + '_' + str(nEvts) + 'Evts'+NormTypeName[iNormType]+'.tex'), 'w')

        # Check whether these output files already exist, otherwise delete them !
        if os.path.isfile(os.path.join("../" + CanvasOutputFile.name[:-4] + ".pdf")):
            os.system('rm ../' + CanvasOutputFile.name[:-4] + '.pdf')

        CanvasOutputFile.write('\\documentclass[a4paper,landscape]{article} \n')
        CanvasOutputFile.write('\\usepackage{graphicx} \n ')
        CanvasOutputFile.write('\\usepackage[top=.5in, bottom=1.25in, left=.5in, right=.5in,landscape]{geometry} \n \n')
        CanvasOutputFile.write('\\begin{document} \n')

        CanvasOutputFile.write('\\section{Distributions of -ln(L) when ' + NormType[iNormType] + ' normalisation is applied} \n')
        CanvasOutputFile.write('\\begin{abstract} \n')
        CanvasOutputFile.write('  Looking at directory : $'+ whichDir +'$')
        CanvasOutputFile.write('\\end{abstract} \n')
        CanvasOutputFile.write('\n \\centering \n')

        # Include the overall likelihood and secondPol distribution (together one 1 page!):
        # if ("Summed" + NormTypeName[iNormType] + ".pdf") in Canvaslist_dir and ("SummedFit_FirstFit" + NormTypeName[iNormType] + ".pdf") in Canvaslist_dir:
        CanvasOutputFile.write('\\includegraphics[width = 0.32 \\textwidth]{SummedHist' + NormTypeName[iNormType] + '.pdf} \n')
        CanvasOutputFile.write('\\includegraphics[width = 0.32 \\textwidth]{SummedFit_FirstFit' + NormTypeName[iNormType] + '.pdf} \n')
        CanvasOutputFile.write('\\includegraphics[width = 0.32 \\textwidth]{SummedFit_SecondFit' + NormTypeName[iNormType] + '.pdf} \n')

        for File in Canvaslist_dir:
            # Include the stacked canvasses:
            if File.endswith(".pdf") and ("LL" + NormTypeName[iNormType] + "_Nr" in File):
                CanvasOutputFile.write('\\includegraphics[width = 0.9 \\textwidth]{' + File + '} \n')
                # else:
                #    CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.9 \\textwidth]{' + File + '} \n')
        # Clear the contents of the OtherNorms array
        OtherNorms[:] = []

        # Write the \end document sentence:
        CanvasOutputFile.write('\\end{document} \n')
        CanvasOutputFile.close()

        # Now create the PDF document using pdflatex!
        os.system('pdflatex -interaction=batchmode ' + CanvasOutputFile.name)
        print "Want to move file with name : ", (CanvasOutputFile.name[:-4] + ".pdf \n")
        shutil.move(CanvasOutputFile.name[:-4] + ".pdf", "..")

print " --> All information is stored in : ", whichDir
