# --------------------------------------------------------------------------------------------------------------- #
#  This script will take care of the creation of all necessary files and plots for the Calibration-curve studies  #
#  At first it will be used to create the different ROOT files where the likelihoods are stored                   #
#    --> Here all the different cuts can be specified that need to be considered                                  #
#  Secondly the calibration-curves itself will be created for the cut-values                                      #
#                                                                                                                 #
#  Remark: Script should be able to work for individual updates as well ...                                       #
# --------------------------------------------------------------------------------------------------------------- #

import os
import sys
import shutil

# Correct syntax is : python RunCalibrationForCuts.py
if len(sys.argv) <= 1:
    print "Need to give the directory of interest !"
    print "Correct syntax is : python RunCalibrationForCuts.py whichDir doCalc(y/n) (opt - int/-1)CutValue (opt - y/n/-1)ExclOuterBins (opt - NoCuts/CutsAlsoOnMet/-1)cutType"
    sys.exit()
whichDir = sys.argv[1]
doCalc = False
if str(sys.argv[2]) == 'y':
    doCalc = True

LikelihoodCuts = [52, 53, 54, 56, 58, 100]     # 100 corresponds to no cut applied!
ExclOuterBinsFit = ['y', 'n']
cutType = ['NoCuts', 'CutsAlsoOnMET']

# *** Other arguments are optional and can be used to avoid running over all possible cases!  *** #
# Which likelihood cut:
if len(sys.argv) > 3 and str(sys.argv[3]) != '-1':
    LikelihoodCuts = [int(sys.argv[3])]
if len(sys.argv) > 4 and str(sys.argv[4]) != '-1':
    ExclOuterBinsFit = [str(sys.argv[4])]
if len(sys.argv) > 5 and str(sys.argv[5]) != '-1':
    cutType = [str(sys.argv[5])]

# ******* Now start with the actual calculations  ******* #
#   --> Will call two different files:                    #
#       1) FitDeviationScript_FinalAnalysis.py            #
#       2) CalibrationCurve.C                             #
# ------------------------------------------------------- #

list_dir = os.listdir(sys.argv[1])
if doCalc:
    for CalibDir in list_dir:
        if CalibDir.startswith('CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample'):
            if (len(cutType) == 1 and CalibDir.endswith(cutType[0])) or len(cutType) > 1:
                consPath = os.path.join(whichDir + CalibDir)
                for iCut in range(len(LikelihoodCuts)):
                    for iFit in range(len(ExclOuterBinsFit)):
                        os.system('python FitDeviationScript_FinalAnalysis.py ' + str(consPath) + '/ MG -1 n Wide ' + str(consPath) + '/weights.out y n n ' + str(LikelihoodCuts[iCut]) + ' ' + str(ExclOuterBinsFit[iFit]))

# Now use this same script to create the different calibration curves if requested!
# --> Will need to make sure the script can also be ran only if one curve needs te be recreated
os.system('g++ -m64 -g -I `root-config --incdir` `root-config --libs` CalibrationCurve.C -o CalibrationCurve')

# Make a list containing all the ROOT files that need to be opened for this specific calibration curve!
HistSum, FirstPolFit, SecondPolFit, FitComp = ' \n ', ' \n ', ' \n ', ' \n '
for iCut in range(len(cutType)):
    for iExcl in range(len(ExclOuterBinsFit)):
        for iLik in range(len(LikelihoodCuts)):

            # Set the name of the output ROOT file correctly
            OutputTitle = 'CalibrationCurve'
            if str(LikelihoodCuts[iLik]) != '100':
                OutputTitle += '_LikelihoodCut'+str(LikelihoodCuts[iLik])
            if ExclOuterBinsFit[iExcl] == 'y':
                OutputTitle += '_OuterBinsExclForFit'
            OutputTitle += '_' + cutType[iCut] + '.root'
            print "\nOutput Title is : ", OutputTitle

            nrFiles = 0
            ROOTFiles = ''

            # For this specific configuration store the relevant ROOT files for all variables!
            for calibDir in os.listdir(sys.argv[1]):
                if calibDir.startswith('CalibCurve_SemiMu_RgR_AllDeltaTF_MGSample') and calibDir.endswith(cutType[iCut]):

                    # Now check all the files existing in this variable directory and select the desired ROOT file
                    for rFile in os.listdir(os.path.join(whichDir + calibDir)):
                        if ((ExclOuterBinsFit[iExcl] == 'y' and 'OuterBinsExclForFit' in rFile) or (ExclOuterBinsFit[iExcl] == 'n' and not 'OuterBinsExclForFit' in rFile)) and ((str(LikelihoodCuts[iLik]) != '100' and 'LikelihoodCut'+str(LikelihoodCuts[iLik]) in rFile) or (str(LikelihoodCuts[iLik]) == '100' and not 'LikelihoodCut' in rFile and rFile.endswith('.root'))) and not 'CosTheta' in rFile:
                            ROOTFiles += whichDir+calibDir+'/'+rFile+' '
                            nrFiles += 1
            # print "Stored list of RootFiles is : ", ROOTFiles
            os.system('./CalibrationCurve ' + str(whichDir) + str(OutputTitle) + ' ' + str(nrFiles) + ' ' + ROOTFiles)
            HistSum += '  \\includegraphics[width = 0.5 \\textwidth]{' + whichDir + 'Plots/' + OutputTitle[:-5] + '_HistSum.pdf' + '} \n'
            FirstPolFit += '  \\includegraphics[width = 0.5 \\textwidth]{' + whichDir + 'Plots/' + OutputTitle[:-5] + '_FirstPolFit.pdf' + '} \n'
            SecondPolFit += '  \\includegraphics[width = 0.5 \\textwidth]{' + whichDir + 'Plots/' + OutputTitle[:-5] + '_SecondPolFit.pdf' + '} \n'
            FitComp += '  \\includegraphics[width = 0.5 \\textwidth]{' + whichDir + 'Plots/' + OutputTitle[:-5] + '_FitComp.pdf' + '} \n'

# Now store everything into a latex file and create a pdf!

# Check whether these output files already exist, otherwise delete them !
latexOutputName = 'CalibrationCurveInfo.tex'
if os.path.isfile(os.path.join(whichDir + latexOutputName[:-4] + ".pdf")):
    os.system('rm ' + whichDir + latexOutputName[:-4] + '.pdf')
if os.path.isfile(os.path.join(whichDir + latexOutputName)):
    os.system('rm ' + whichDir + latexOutputName)

print " \n *** Saving as pdf file ***"
latexOutput = open('CalibrationCurveInfo.tex', 'w')
latexOutput.write('\\documentclass[a4paper, landscape]{article} \n\\usepackage{graphicx} \n ')
latexOutput.write('\\usepackage{geometry} \n \n')
latexOutput.write('\\usepackage[top=.5in, bottom=1.25in, left=.5in, right=.5in,landscape]{geometry} \n \n')
latexOutput.write('\n \\begin{document} \n')
latexOutput.write(' \\begin{abstract} \n')
latexOutput.write('   File where all the calibration curves are compared for different likelihood cuts \\\\ \n')
latexOutput.write(' \\end{abstract} \n \n')

latexOutput.write(' \\section{Original summed histogram} ')
latexOutput.write(' ' + HistSum)
latexOutput.write(' \\section{First polynomial fit} ')
latexOutput.write(' ' + FirstPolFit)
latexOutput.write(' \\section{Second polynomial fit} ')
latexOutput.write(' ' + SecondPolFit)
latexOutput.write(' \\section{Fit comparison}')
latexOutput.write(' ' + FitComp)

latexOutput.write(' \\end{document}')
latexOutput.close()

os.system('pdflatex -interaction=batchmode '+ latexOutput.name)
shutil.move(latexOutput.name, whichDir)
shutil.move(latexOutput.name[:-4] + ".pdf", whichDir)
os.system('rm '+ latexOutput.name[:-4] + ".aux")
os.system('rm '+ latexOutput.name[:-4] + ".log")