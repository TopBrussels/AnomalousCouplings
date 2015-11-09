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
    print "Correct syntax is : python RunCalibrationForCuts.py whichDir (opt - int/-1)CutValue (opt - y/n/-1)ExclOuterBins"
    sys.exit()
whichDir = sys.argv[1]

LikelihoodCuts = [60, 62, 65, 70, 100]       # 100 corresponds to no cut applied! (Reco cuts)
ExclOuterBinsFit = ['y', 'n']

# *** Other arguments are optional and can be used to avoid running over all possible cases!  *** #
# Which likelihood cut:
if len(sys.argv) > 2 and str(sys.argv[2]) != '-1':
    LikelihoodCuts = [int(sys.argv[2])]
if len(sys.argv) > 3 and str(sys.argv[3]) != '-1':
    ExclOuterBinsFit = [str(sys.argv[3])]

# ******* Now start with the actual calculations  ******* #
#   --> Will call the file                                #
#          FitDeviationScript_FinalAnalysis.py            #
# ------------------------------------------------------- #

list_dir = os.listdir(sys.argv[1])
for iCut in range(len(LikelihoodCuts)):
    for iFit in range(len(ExclOuterBinsFit)):
        if os.path.exists(os.path.join(whichDir+'weights_CheckedEvts.out')):
            os.system('python FitDeviationScript_FinalAnalysis.py ' + str(whichDir) + ' RECO -1 n Normal ' + str(whichDir) + 'weights_CheckedEvts.out y n n ' + str(LikelihoodCuts[iCut]) + ' ' + str(ExclOuterBinsFit[iFit]))
        else:
            os.system('python FitDeviationScript_FinalAnalysis.py ' + str(whichDir) + ' RECO -1 n Normal ' + str(whichDir) + 'weights.out y n n ' + str(LikelihoodCuts[iCut]) + ' ' + str(ExclOuterBinsFit[iFit]))
