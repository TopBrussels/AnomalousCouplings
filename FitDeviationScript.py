###########################################################
##                                                       ##
##  Python macro which takes care of fitting the -ln(L)  ##
##  Uses the following ROOT macros:                      ##
##   - fitDeviationMacro.C                               ##
##   - PerformFitOptimization.C                          ##
##                                                       ##
##  ** First one performs 2 consecutive fits             ##
##       * First one using all points                    ##
##       * Second one use 66% of best points             ##
##  ** Second one allows for quick evt sel tests         ##
##       * Uses the TF1's created in first macro         ##
##       * Allows to test different chi-sq cuts          ##
##                                                       ##
###########################################################

#! python
import os
import sys
import re
import shutil
from array import array

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and the number of events in the command line !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
  sys.exit()
elif len(sys.argv) == 3:
  print "Need to give the number of considered events !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
  sys.exit()
elif len(sys.argv) == 4:
  print "Need to specify whether the performed fits have to be updated (hence run fitDeviationMacro.C)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
  sys.exit()
elif len(sys.argv) == 5:
  print "Need to mention whether Tex output is wanted (y/n)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
elif len(sys.argv) == 6:
  print "Need to specify which analysis should be executed: fctDeviationMacro or doublePolFitMacro (fctDeviation/doublePolFit)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"
elif len(sys.argv) == 7:
  print "Need to specify which whether the acceptance normalisation should be applied! (y/n)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ Var(MTop/RVR/RgR) #evts ForceFitUpdate(y/n) TexWanted(y/n) WhichAnalysis(fctDeviation.C/doublePolFit.C) applyAccNorm(y/n)"

whichDir = sys.argv[1]
KinVariable = sys.argv[2]
nEvts = sys.argv[3]
whichAnalysis = sys.argv[6]
applyAccNorm = sys.argv[7]

applyCosTheta = "n"   #Initialized to "no"!
WeightsFileGiven, VarWindowGiven = False, False
#Special case when more than the obligatory number of input is given!
#  --> This way it is possible to skip the questions asked by the script
#  --> Correct syntax is : python FitDeviation.py(#0) Events/blabla/(#1) RVR(#2) 10000(#3) y(#4) n(#5) doublePolFit.C(#6) n(#7) n(#8-optional) Events/blabla/weights.out(#9-optional) Range(#10-optional)
#
#-- Cos theta normalisation

if len(sys.argv) >= 9: applyCosTheta = sys.argv[8]
#-- Specific weights file
if len(sys.argv) >= 10:
  WeightsFileGiven = True
  WeightsFileName = sys.argv[9]
  WeightsFile = open(os.path.join(WeightsFileName),'r')
#-- Specific range
if len(sys.argv) >= 11:
  VarWindowGiven = True
  if sys.argv[10] == "Wide" and KinVariable == "RVR": VarWindow = "1"
  elif sys.argv[10] == "Narrow" and KinVariable == "RgR": VarWindow = "1"
  elif sys.argv[10] == "Full" and KinVariable == "RgR": VarWindow = "3"
  elif sys.argv[10] == "CalibCurve" and KinVariable == "RgR": VarWindow = "5"

#Set the 'CreateTexFile' correctly:
if sys.argv[5] == "y" or sys.argv[5] == "yes":  CreateTexFile = True
elif sys.argv[5] == "n" or sys.argv[5] == "no": CreateTexFile = False
else: print "!!!!! Simple yes/no was needed for TexWanted boolean!!!!!! ", sys.exit()

if KinVariable != "MTop" and KinVariable != "RVR" and KinVariable != "RgR":
  print "Need to specify which kinematic variable should be considered (MTop, RVR or RgR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the three : ')

ValuesToDelete = []
if KinVariable == "RVR":
  #Information about the scanned RVR values and the corresponding cross-section
  Var =        array('d',[-1.5,    -1.0,     -0.5,    -0.3,     -0.2,     -0.1,       -0.05,    0.0,       0.05,     0.1,       0.2,       0.3,       0.5,      1.0,     1.5     ])
  MGXS =       array('d',[122.082, 46.4474,  17.9275, 13.3944,  12.06555, 11.25909,   11.02784, 10.90059,  10.88228, 10.97767,  11.49883,  12.49056,  16.1508,  40.8074, 108.249 ])
  MGXSCut =    array('d',[28.507,  10.63436, 3.95435, 2.92922,  2.62439,  2.4352,     2.38608,  2.35285,   2.35117,  2.37359,   2.49101,   2.72632,   3.58445,  9.36921, 25.4672 ])
  MGXSe =      array('d',[1,       1,        0.01231, 0.009958, 0.009346, 0.00836608, 0.0,      0.0082221, 0.0,      0.0084729, 0.0090198, 0.0087468, 0.011308, 1,       1       ])
  Acceptance = array('d',[1,       1,        0.22164, 0.21742,  0.21672,  0.21737,    0.21614,  0.21670,   0.21531,  0.21677,   0.21437,   0.21793,   0.22205,  1,       1       ])

  #Select which window of RVR values was considered!
  if VarWindowGiven == False:
    VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  3) Many   : [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1] \n  4) Very wide : [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] \n 5) Wide & Many : [-0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3] \n 6) Narrow short : [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3] \n --> Choose the correct number : ')

  xMinValue = [4, 4, 10, 3, 12, 3]
  NrPointsToRemove = [3, 3, 7, 2, 7, 1]
  KinVar = "Re(V_{R})"
  FitType = "pol4"

  if VarWindow == "1":
    ValuesToDelete = [-1.5, -1.0, -0.05, 0.05, 1.0, 1.5]
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    ValuesToDelete = [-1.5, -1.0, -0.5, 0.5, 1.0, 1.5]
    xBin, xLow, xHigh = 13, -0.325, 0.325
  elif VarWindow == "3":
    Var = array('d',[-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    MGXS = array('d',[11.25909, 11.19599, 11.14075, 11.10667, 11.06145, 11.02784, 10.99474, 10.9534, 10.93846, 10.91787, 10.90059, 10.89643, 10.86829, 10.87792, 10.87266, 10.88228, 10.87684, 10.89534, 10.91641, 10.9317, 10.97767])
    MGXSCut = array('d',[2.43546, 2.42182, 2.41059, 2.40057, 2.38839, 2.38187, 2.36976, 2.36513, 2.35512, 2.35666, 2.35415, 2.35694, 2.35174, 2.34909, 2.34392, 2.35108, 2.34767, 2.35477, 2.36148, 2.3643, 2.37424])
    xBin, xLow, xHigh = 21, -0.105, 0.105
  elif VarWindow == "4":
    ValuesToDelete = [-0.3, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.3]
    xBin, xLow, xHigh = 7, -1.75, 1.75
  elif VarWindow == "5":
    Var  = array('d',[-0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3])
    MGXS = array('d',[13.3944, 13.037, 12.66011, 12.37463, 12.06555, 11.83271, 11.60956, 11.4194, 11.25909, 11.12321, 11.02784, 10.9524, 10.90059, 10.87549, 10.88228, 10.93437, 10.97767, 11.07142, 11.17366, 11.32792, 11.49883, 11.69063, 11.90668, 12.18904, 12.49056])
    MGXSCut = array('d',[2.92922, 1,      1,        1,        2.62439,  1,        1,        1,       2.4352,   1,        2.38608,  1,       2.35285,  1,        2.35117,  1,        2.37359,  1,        1,        1,        2.49101,  1,        1,        1,        2.72632 ])
    xBin, xLow, xHigh = 25, -0.3025, 0.3025
  elif VarWindow == "6":
    ValuesToDelete = [-1.5, -1.0, -0.5, -0.05, 0.05, 0.5, 1.0, 1.5]
    xBin, xLow, xHigh = 7, -0.35, 0.35
elif KinVariable == "MTop":
  #Information about the scanned MTop values and the corresponding cross-section
  Var =        array('d',[153,      163,      170,       171,       172,      173,        174,        175,       183,       193     ])
  MGXS =       array('d',[8.20916,  9.6299,   10.57123,  10.70485,  10.8257,  10.96469,   11.08428,   11.22448,  12.18068,  13.3046 ])
  MGXSCut =    array('d',[1.35059,  1.85406,  2.21902,   2.27174,   2.32261,  2.38097,    2.43678,    2.49254,   2.93184,   3.50146 ])
  MGXSe =      array('d',[0.006413, 0.007759, 0.0085714, 0.0081430, 0.008789, 0.00816802, 0.00904798, 0.0086538, 0.0093129, 0.010331])
  Acceptance = array('d',[0.16203,  0.19152,  0.21008,   0.21460,   0.21735,  0.21290,    0.21752,    0.22185,   0.23941,   0.26413 ])

  #Select which window of masses was considered!
  if VarWindowGiven == False:
    VarWindow = raw_input('** Choose the correct mass-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [153, 163, 170, 171, 172, 173, 174, 175, 183, 193] \n  2) Normal : [153, 163, 171, 172, 173, 174, 175, 183, 193] \n  3) Narrow : [171, 172, 173, 174, 175] \n --> Choose the correct number : ')

  xMinValue = [5, 4, 2]
  KinVar = "m_{top}"
  NrPointsToRemove = [3, 3, 1]
  FitType = "pol2"

  if VarWindow == "1":
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "2":
    ValuesToDelete = [170]
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "3":
    ValuesToDelete = [153, 163, 170, 183, 193]
    Bin, xLow, xHigh = 5, 170.5, 175.5

elif KinVariable == "RgR":
  Var =     array('d',[-0.5,     -0.4,     -0.3,     -0.2,     -0.15,    -0.1,    -0.05,   -0.025,  0.0,      0.025,    0.05,    0.1,     0.15,    0.2,      0.3,     0.4,     0.5     ])
  MGXS =    array('d',[1.8647,   1.97357,  3.36424,  4.92909,  6.02588,  7.34593, 8.94878, 9.88333, 10.89487, 11.92922, 13.1987, 15.9457, 19.1623, 22.9185,  32.2975, 38.8312, 60.5617 ])
  MGXSCut = array('d',[0.363334, 0.465328, 0.639413, 0.912292, 1.098184, 1.32024, 1.58727, 1.73857, 1.90447,  2.0856,   2.28471, 2.71983, 3.2418,  3.838581, 5.31357, 7.23931, 9.66734 ]) #Also MET cuts!
#  MGXSCut = array('d',[0.462736, 0.807034,  1.14397,   1.37121,   1.6414,    1.96796,            2.35719,              2.80657,    3.34178,   3.96076,    4.67808,  6.42868,   11.61381  ])

  if VarWindowGiven == False:
    VarWindow = raw_input('** Choose the correct RgR-window corresponding to the studied file ** \n** Different options are: \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20] \n  3) Full   : [-0.5, -0.3, -0.2, -0.15, -0.1, -0.05, -0.025, 0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5] \n  4) Wide narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  5) Calibration Curve range : [-0.4, -0.3, -0.2, -0.15, 0.1, 0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4] \n  6) Zoomed : [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.05, 0.1, 0.15, 0.2] \n --> Choose the correct number : ')

  xMinValue = [4, 4, 7, 4, 6, 4]
  KinVar = "Re(g_{R})"
  NrPointsToRemove = [2, 2, 5, 2, 4, 2]
  FitType = "pol2"
 
  if VarWindow == "1":
    ValuesToDelete = [-0.4, -0.15, -0.05, -0.025, 0.025, 0.05, 0.15, 0.4]
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    ValuesToDelete = [-0.5, -0.4, -0.3, -0.025, 0.025, 0.3, 0.4, 0.5]
    xBin, xLow, xHigh = 9, -0.225, 0.225
  elif VarWindow == "3":
    ValuesToDelete = [-0.4, 0.4]
    xBin, xLow, xHigh = 41, -0.5125, 0.5125
  elif VarWindow == "4":
    ValuesToDelete = [-0.5, -0.4, -0.15, -0.025, 0.025, 0.15, 0.4, 0.5]
    xBin, xLow, xHigh = 13, -0.325, 0.325
  elif VarWindow == "5":
    ValuesToDelete = [-0.5, -0.025, 0.025, 0.5]
    xBin, xLow, xHigh = 17, -0.425, 0.425
  elif VarWindow == "6":
    ValuesToDelete = [-0.5, -0.4, -0.3, -0.025, 0.025, 0.3, 0.4, 0.5]
    xBin, xLow, xHigh = 9, -0.225, 0.225

#Now delete the values stored in array 'ValuesToDelete'
for iVar in range(len(ValuesToDelete)):
  MGXS.pop(Var.index(ValuesToDelete[iVar])), MGXSCut.pop(Var.index(ValuesToDelete[iVar])), Var.pop(Var.index(ValuesToDelete[iVar])) 

#---------------------------------------------#
#   Special cases for MGXSCut initialization  #
#---------------------------------------------#
#In the case that MGSample or GEN info is used, the acceptance-normalisation shouldn't be applied
if applyAccNorm == "n"or applyAccNorm == "no" or applyAccNorm == "No":
  MGXSCut = MGXS
  print "  ==> No acceptance normalisation will be applied! \n"
  if (whichDir.find("RECO") <= len(whichDir) and whichDir.find("RECO") > 0) or (whichDir.find("Reco") <= len(whichDir) and whichDir.find("Reco") > 0):
    print " \n ************************* ERROR **************************** \n -----> Applying no Acceptance normalisation for reco-level events .... \n"
else:
  print " Applying acceptance normalisation ! \n"
  if (whichDir.find("MGSample") <= len(whichDir) and whichDir.find("MGSample") > 0 and whichDir.find("Cut") < 0) or (whichDir.find("GEN") <= len(whichDir) and whichDir.find("GEN") > 0):
    print " \n ************************* ERROR **************************** \n -----> Applying Acceptance normalisation to generator-level events .... \n"

#Perform Acc-norm by applying the normalisation for some cases
#if (whichDir.find("AccNormForCuts") <= len(whichDir) and whichDir.find("AccNormForCuts") >= 0) or (FullInfoGiven == True and (applyAccNorm == "y" or applyAccNorm == "Y" or applyAccNorm == "yes")):
#  if WeightsFileName.find("Cut15") <= len(WeightsFileName) and WeightsFileName.find("Cut15") >= 0:
#    MGXSCut = array('d',[13.460863375, 10.01499288, 8.983043286, 8.376425187, 8.118650426, 8.174741519, 8.559729052, 9.324577757, 12.131027388]) 
#    title += "_AccNormApplied"
#  elif WeightsFileName.find("Cut30") <= len(WeightsFileName) and WeightsFileName.find("Cut30") >= 0:
#    MGXSCut = array('d',[4.772121225, 3.471694536, 3.103621427, 2.907659993, 2.789024957, 2.832458413, 2.945885258, 3.251167862, 4.277700888])
#    title += "_AccNormApplied"
#  elif WeightsFileName.find("Cut20") <= len(WeightsFileName) and WeightsFileName.find("Cut20") >= 0:

print " List of considered Var values is : ",Var,"\n"
NrConfigs = len(Var)
xMin = xMinValue[int(VarWindow)-1]
NumberOfPointsToRemove = NrPointsToRemove[int(VarWindow)-1]

#File of interest (only search if WeightsFileGiven is set to false):
if WeightsFileGiven == False:
  list_dir = os.listdir(whichDir)
  WeightsFileArray, weightsFileCounter = [], 0
  for file in list_dir:
    if (applyCosTheta == "n" and file.endswith(".out") or applyCosTheta == "y" and file.endswith("ApplyCosThetaReweighting.out") ) and file.startswith("weights"):
      weightsFileCounter += 1
      WeightsFileArray.append(file)

  if int(weightsFileCounter) == 0: print "No weights file found in this directory !", sys.exit()
  elif int(weightsFileCounter) == 1: WeightsFileName = str(whichDir)+''+str(WeightsFileArray[0])
  elif int(weightsFileCounter) > 1:
    for ii in range(len(WeightsFileArray)):
      print " ",ii," ) ",WeightsFileArray[ii]
    fileNr = raw_input('Choose the number of the file of interest! : ')
    WeightsFileName = str(whichDir)+''+str(WeightsFileArray[int(fileNr)])

#Open the selected files!   
WeightsFile, LikelihoodFile = open(WeightsFileName,'r'), open(WeightsFileName,'r')  

#Count the number of events actually present in the WeightsFile as the maximum numbers which can actually be used:
#Check whether the file has enough events, otherwise use the maximum number
maxNrEvts = os.popen('grep " 1 1 " '+str(WeightsFileName)+' | wc -l').read()
if int(maxNrEvts) < int(nEvts): nEvts = int(maxNrEvts)
print " Will be using file : ",WeightsFileName," with ",nEvts," events ! "

#----
#Special case: Scaling the XS-array with a xx% function:
#----
TitleChange = ""
if whichDir.find("ChangingXS") <= len(whichDir) and whichDir.find("ChangingXS") > 0:
  MGXSScaled, MGXSCutScaled = [], []
  ScaleValue, TitleChange = 0.05, "_XSScaledWithPos005"
  for iVar in range(len(MGXS)):
    MGXSScaled.append(MGXS[iVar]*(1+Var[iVar]*ScaleValue)), MGXSCutScaled.append(MGXSCut[iVar]*(1+Var[iVar]*ScaleValue))
  MGXS, MGXSCut = MGXSScaled, MGXSCutScaled

#----------------------------#
#   Setting output title     #
#   --> Use directory name   #
#   --> Add weightsFileName  #
#----------------------------#
title = str(whichDir[whichDir.find("/")+1:-1])  #Only need the part after the "/"!!
#Special cases!
#*** Indicating that the cos theta* reweighting has been applied
if applyCosTheta == "y" or applyCosTheta == "yes" or applyCosTheta == "Y" or applyCosTheta == "Yes": title += "_CosThetaReweightingApplied"
#*** Indicating that XS-values are changed!
if whichDir.find("ChangingXS") <= len(whichDir) and whichDir.find("ChangingXS") > 0:
  title += TitleChange
#*** Indicating that Pt-cuts have been applied!
if WeightsFileName.find("NoLowPt") <= len(WeightsFileName) and WeightsFileName.find("NoLowPt") >= 0: 
  title = title+"_"+"NoLowPtEvts"+WeightsFileName[WeightsFileName.find("_Cut"):-4]

#Set the 'RunFitMacro' correctly:
if sys.argv[4] == "y" or sys.argv[4] == "yes": 
  RunFitMacro = True
  print " **** Will update the fits ! ***** "
elif sys.argv[4] == "n" or sys.argv[4] == "no":
  RunFitMacro = False
  #However in case the ROOT file of interest doesn't exist, this boolean has to be overwritten since fit has to be done anyway!
  if not (os.path.exists(os.path.join(whichDir+"FitDistributions_"+str(title)+"_"+str(nEvts)+"Evts.root"))):
    RunFitMacro = True
    print "The ROOT file which was needed for the PerformFitOptimization.C macro to run is missing !! "
    print "     ==> Will run the fitDeviationMacro.C with the requested configuration! "

#-------------------------------------------------#
#--  Pass on all variables to the ROOT macro !  --#
#-------------------------------------------------#
if RunFitMacro == True:
  RootAnalyzer, NewRootAnalyzer = open(whichAnalysis,'r'), open('output.C','w')
  
  VarLine = 'double Var[] = {'
  MGXSLine, MGXSCutLine = 'double MGXS[] = {', 'double MGXSCut[] = {'
  for ii in range(len(Var)):
    VarLine, MGXSLine, MGXSCutLine = VarLine+str(Var[ii]), MGXSLine+str(MGXS[ii]), MGXSCutLine+str(MGXSCut[ii])
    if ii < len(Var)-1: VarLine, MGXSLine, MGXSCutLine = VarLine+',',   MGXSLine+',',    MGXSCutLine+','
    else:               VarLine, MGXSLine, MGXSCutLine = VarLine+'};\n', MGXSLine+'};\n', MGXSCutLine+'};\n'

  xMinValueLine = 'int xMinValue[] = {'
  for ii in range(len(xMinValue)):
    if ii < len(xMinValue)-1: xMinValueLine += str(xMinValue[ii])+','
    else:                     xMinValueLine += str(xMinValue[ii])+'}; \n'

  #-->Create the directory SplittedCanvasses if needed (otherwise delete the created pdf files ...)!
  if CreateTexFile == True:
    if not (os.path.exists(os.path.join(whichDir+'SplittedCanvasses/')) ): os.makedirs(os.path.join(whichDir+'SplittedCanvasses/'))
    else:
      print "Deleting the existing SplitCanvas files ! "
      os.system('rm '+whichDir+'SplittedCanvasses/SplitCanvasLL*.pdf')

  for RootLine in RootAnalyzer:
    RootWord = RootLine.split()

    #Changes valid for both files!
    if   re.search( r"int nEvts",          RootLine): NewRootAnalyzer.write('const int nEvts = '+str(nEvts)+'; \n')
    elif re.search( r"int NrConfigs",      RootLine): NewRootAnalyzer.write('const int NrConfigs = '+str(NrConfigs)+'; \n')
    elif re.search( r"double Var",         RootLine): NewRootAnalyzer.write(VarLine)
    elif re.search( r"double MGXSCut",     RootLine): NewRootAnalyzer.write(MGXSCutLine)
    elif re.search( r"double MGXS",        RootLine): NewRootAnalyzer.write(MGXSLine)
    elif re.search( r"int xBin",           RootLine): NewRootAnalyzer.write('int xBin = '+str(xBin)+'; \n')
    elif re.search( r"float xLow",         RootLine): NewRootAnalyzer.write('float xLow = '+str(xLow)+'; \n')
    elif re.search( r"float xHigh",        RootLine): NewRootAnalyzer.write('float xHigh = '+str(xHigh)+'; \n')
    elif re.search( r"int xMinValue",      RootLine): NewRootAnalyzer.write(xMinValueLine)
    elif re.search( r"string KinVar",      RootLine): NewRootAnalyzer.write('std::string KinVar = "'+str(KinVar)+'"; \n')
    elif re.search( r"int VarWindow",      RootLine): NewRootAnalyzer.write('int VarWindow = '+str(VarWindow)+'; \n')
    elif re.search( r"std::ifstream ifs",  RootLine): NewRootAnalyzer.write('  std::ifstream ifs ("'+str(WeightsFileName)+'", std::ifstream::in); \n')
    elif re.search( r"std::string title",  RootLine): NewRootAnalyzer.write('std::string title = "'+str(title)+'"; \n')
    elif re.search( r"string SplittedDir", RootLine): NewRootAnalyzer.write('std::string SplittedDir = "'+str(whichDir)+'SplittedCanvasses"; \n')
    elif re.search( r"iss >> evt",         RootLine):
      if applyCosTheta == "y" or applyCosTheta == "Y":
        NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> CosThetaCorr >> weightUnc){ \n')
        print " Cos theta* reweighting will be applied! \n"
      else:
        NewRootAnalyzer.write('    if( iss >> evt >> config >> tf >> weight >> weightUnc){ \n') 
    elif re.search( r"bool storeSplittedCanvas", RootLine):
      if   CreateTexFile == True: NewRootAnalyzer.write('bool storeSplittedCanvas = true; \n')
      else:                       NewRootAnalyzer.write('bool storeSplittedCanvas = false; \n')
    #Changes specific for any of the two files
    elif whichAnalysis == "doublePolFitMacro.C" and re.search( r"int NrToDel", RootLine): 
      NewRootAnalyzer.write('const unsigned int NrToDel = '+str(NumberOfPointsToRemove)+'; \n')
    elif whichAnalysis == "doublePolFitMacro.C" and re.search( r"new TF1",     RootLine):
      print "\n ---> Fit will go between ",Var[1]," and ",Var[NrConfigs-2],"\n"
      if   re.search( r"AllPoints",     RootLine): NewRootAnalyzer.write('  polFit_AllPoints = new TF1(("polFit"+Type+"_AllPoints_Evt"+EvtNumber).c_str(),"'+str(FitType)+'",Var[1],Var[NrConfigs-2]); \n')
      elif re.search( r"ReducedPoints", RootLine): NewRootAnalyzer.write('  polFit_ReducedPoints = new TF1(("polFit"+Type+"_"+sNrRemaining+"ReducedPoints_Evt"+EvtNumber).c_str(),"'+str(FitType)+'",Var[1],Var[NrConfigs-2]); \n')
    elif whichAnalysis == "doublePolFitMacro.C" and re.search( r"new TFile", RootLine) and re.search( r"file_FitDist", RootLine):
      NewRootAnalyzer.write('TFile* file_FitDist = new TFile("'+str(whichDir)+'FitDistributions_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n')
    elif whichAnalysis == "fctDeviationMacro.C" and re.search( r"new TFile", RootLine): 
      NewRootAnalyzer.write('TFile* Tfile = new TFile("'+str(whichDir)+'FctDeviation_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n')
    else:                                             NewRootAnalyzer.write(RootLine)
  NewRootAnalyzer.close(), RootAnalyzer.close()


  #Run the root macro!  
  os.rename('output.C',whichAnalysis), os.system("root -l -b -q "+whichAnalysis+"+")

#Now perform the fitOptimizations! 
if whichAnalysis == "doublePolFitMacro.C": 
  PerformFitOptAnalyzer, NewPerformFitOptAnalyzer = open('PerformFitOptimization.C','r'), open('fitOptimization.C','w')

  for FitOptLine in PerformFitOptAnalyzer:
    FitOptWord = FitOptLine.split()
    if   re.search( r"int nEvts",         FitOptLine): NewPerformFitOptAnalyzer.write('const int nEvts = '+str(nEvts)+'; \n')
    elif re.search( r"int xBin",          FitOptLine): NewPerformFitOptAnalyzer.write('const int xBin = '+str(xBin)+'; \n')
    elif re.search( r"float xLow",        FitOptLine): NewPerformFitOptAnalyzer.write('float xLow = '+str(xLow)+'; \n')
    elif re.search( r"float xHigh",       FitOptLine): NewPerformFitOptAnalyzer.write('float xHigh = '+str(xHigh)+'; \n')
    elif re.search( r"string SplittedDir", FitOptLine): NewPerformFitOptAnalyzer.write('std::string SplittedDir = "'+str(whichDir)+'SplittedCanvasses"; \n')
    elif re.search( r"new TFile", FitOptLine):
      if   re.search( r"inFile",     FitOptLine): NewPerformFitOptAnalyzer.write('TFile *inFile = new TFile("'+str(whichDir)+'FitDistributions_'+str(title)+'_'+str(nEvts)+'Evts.root","READ"); \n')
      elif re.search( r"outputFile", FitOptLine): NewPerformFitOptAnalyzer.write('TFile *outputFile = new TFile("'+str(whichDir)+'FitOptimizations_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n')
    elif re.search( r"bool storeSplittedCanvas", FitOptLine):
      if   CreateTexFile == True: NewPerformFitOptAnalyzer.write('bool storeSplittedCanvas = true; \n')
      else:                       NewPerformFitOptAnalyzer.write('bool storeSplittedCanvas = false; \n')
    elif re.search( r"float ChiSqCutsFstPol", FitOptLine):
      if (title.find("Gen") <= len(title) and title.find("Gen") >= 0) or (title.find("MGSample") <= len(title) and title.find("MGSample") >= 0): 
        NewPerformFitOptAnalyzer.write('    float ChiSqCutsFstPol[4] = {0.0005, 0.00008, 0.00005, 0.00001}; \n')
      else:
        NewPerformFitOptAnalyzer.write('    float ChiSqCutsFstPol[4] = {0.0002, 0.001, 0.0005, 0.005}; \n')
    elif re.search( r"float ChiSqCutsScdPol", FitOptLine):
      if (title.find("Gen") <= len(title) and title.find("Gen") >= 0) or (title.find("MGSample") <= len(title) and title.find("MGSample") >= 0): 
        NewPerformFitOptAnalyzer.write('    float ChiSqCutsScdPol[4] = {0.0002, 0.00005, 0.00001, 0.000005}; \n')
      else:
        NewPerformFitOptAnalyzer.write('    float ChiSqCutsScdPol[4] = {0.0002, 0.001, 0.0005, 0.005}; \n')
    else:                                        NewPerformFitOptAnalyzer.write(FitOptLine)
  NewPerformFitOptAnalyzer.close(), PerformFitOptAnalyzer.close()
  os.rename('fitOptimization.C','PerformFitOptimization.C'), os.system("root -l -b -q PerformFitOptimization.C+") 

  #-- Now plot the comparison of the different chisq-cut LL's
#  DrawLLAnalyzer, NewDrawAnalyzer = open('DrawNormalizedLL.C','r'), open('drawLL.C','w')
#
#  for DrawLine in DrawLLAnalyzer:
#    if   re.search( r"new TFile",                 DrawLine): NewDrawAnalyzer.write('  TFile* inFile = new TFile("'+str(whichDir)+'FitOptimizations_'+str(title)+'_'+str(nEvts)+'Evts.root","READ"); \n')
#    elif re.search( r"coefficient",               DrawLine): 
#      if   re.search( r"FstPol", DrawLine): NewDrawAnalyzer.write('  h_OriginalFstPol->GetXaxis()->SetTitle("'+str(KinVariable)+'coefficient"); \n')
#      elif re.search( r"ScdPol", DrawLine): NewDrawAnalyzer.write('  h_OriginalScdPol->GetXaxis()->SetTitle("'+str(KinVariable)+'coefficient"); \n')
#    elif re.search( r"LLFirstFit",                DrawLine): NewDrawAnalyzer.write('  canvas_FirstFit->SaveAs("'+str(whichDir)+'LLFirstFitComparison_ChiSqCuts'+str(title)+'_'+str(nEvts)+'Evts.pdf"); \n')
#    elif re.search( r"LLSecondFit",               DrawLine): NewDrawAnalyzer.write('  canvas_SecondFit->SaveAs("'+str(whichDir)+'LLSecondFitComparison_ChiSqCuts'+str(title)+'_'+str(nEvts)+'Evts.pdf");\n')
#    elif re.search( r"Likelihood shape",          DrawLine): 
#      if   re.search( r"FstPol", DrawLine): NewDrawAnalyzer.write('  h_OriginalFstPol->SetTitle("-Likelihood shapes for different #chi^{2} cuts (Acc norm -- '+str(title)+' -- '+str(nEvts)+' evts)"); \n')
#      elif re.search( r"ScdPol", DrawLine): NewDrawAnalyzer.write('  h_OriginalScdPol->SetTitle("-Likelihood shapes for different #chi^{2} cuts (Acc norm -- '+str(title)+' -- '+str(nEvts)+' evts)"); \n')
#    elif re.search( r"float ChiSqCutsFstPol", DrawLine):
#      if (title.find("Gen") <= len(title) and title.find("Gen") >= 0) or (title.find("MGSample") <= len(title) and title.find("MGSample") >= 0): 
#        NewDrawAnalyzer.write('  float ChiSqCutsFstPol[4] = {0.0005, 0.00008, 0.00005, 0.00001}; \n')
#      else:
#        NewDrawAnalyzer.write('  float ChiSqCutsFstPol[4] = {0.0002, 0.001, 0.0005, 0.005}; \n')
#    elif re.search( r"float ChiSqCutsScdPol", DrawLine):
#      if (title.find("Gen") <= len(title) and title.find("Gen") >= 0) or (title.find("MGSample") <= len(title) and title.find("MGSample") >= 0): 
#        NewDrawAnalyzer.write('  float ChiSqCutsScdPol[4] = {0.0002, 0.00005, 0.00001, 0.000005}; \n')
#      else:
#        NewDrawAnalyzer.write('  float ChiSqCutsScdPol[4] = {0.0002, 0.001, 0.0005, 0.005}; \n')
#    else:                                                    NewDrawAnalyzer.write(DrawLine) 
#    
#  NewDrawAnalyzer.close(), DrawLLAnalyzer.close()
#  os.rename('drawLL.C','DrawNormalizedLL.C'), os.system("root -l -b -q DrawNormalizedLL.C+")

#-- Now store the stacked canvasses in a .txt file --#
if CreateTexFile == True and RunFitMacro == True:

  #Change the working directory of the script!
  os.chdir(os.path.join(whichDir+'SplittedCanvasses/'))
  CanvasOutputFile_NoNorm  = open(os.path.join('FitDeviationSplitCanvas_'+str(title)+'_'+str(nEvts)+'Evts_NoNorm.tex'),'w')
  CanvasOutputFile_XSNorm  = open(os.path.join('FitDeviationSplitCanvas_'+str(title)+'_'+str(nEvts)+'Evts_XSNorm.tex'),'w')
  CanvasOutputFile_AccNorm = open(os.path.join('FitDeviationSplitCanvas_'+str(title)+'_'+str(nEvts)+'Evts_AccNorm.tex'),'w')

  print " Current working directory is : ", os.getcwd()
  NormType = ["no","XS","acceptance"]
  NormTypeName = ["","XS","Acc"]
  CanvasOutputFile = [CanvasOutputFile_NoNorm, CanvasOutputFile_XSNorm, CanvasOutputFile_AccNorm]
  
  #Store the information in the correct directory:
  Canvaslist_dir = os.listdir('.') #os.path.join(whichDir+'SplittedCanvasses/'))

  OtherNorms = []
  for iNormType in range(len(NormType)):
    for ii in range(len(NormType)):
      if ii != iNormType:
        OtherNorms.append(ii)

    #Check whether these output files already exist, otherwise delete them !	
    if os.path.isfile(os.path.join("../"+CanvasOutputFile[iNormType].name[:-4]+".pdf") ): os.system('rm ../'+CanvasOutputFile[iNormType].name[:-4]+'.pdf')

    CanvasOutputFile[iNormType].write('\\documentclass[a4paper,landscape]{article} \n')
    CanvasOutputFile[iNormType].write('\\usepackage{graphicx} \n ')
    CanvasOutputFile[iNormType].write('\\usepackage[top=.5in, bottom=1.25in, left=.5in, right=.5in,landscape]{geometry} \n \n')
    CanvasOutputFile[iNormType].write('\\begin{document} \n')

    CanvasOutputFile[iNormType].write('\\section{Distributions of -ln(L) when '+NormType[iNormType]+' normalisation is applied} \n')
    CanvasOutputFile[iNormType].write('\n \\centering \n')

    #Include the overall likelihood and secondPol distribution (together one 1 page!):
    if ("TotalLnLik"+NormTypeName[iNormType]+".pdf") in Canvaslist_dir and ("SecondPol"+NormTypeName[iNormType]+".pdf") in Canvaslist_dir:
      CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.32 \\textwidth]{TotalLnLik'+NormTypeName[iNormType]+'.pdf} \n')
      CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.32 \\textwidth]{FirstPol'+NormTypeName[iNormType]+'.pdf} \n')
      CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.32 \\textwidth]{SecondPol'+NormTypeName[iNormType]+'.pdf} \n')

    for file in Canvaslist_dir:
      #Include the stacked canvasses:
      if file.endswith(".pdf") and ("LL"+NormTypeName[iNormType] in file):
        if NormType[iNormType] == "no":
          if not (NormTypeName[OtherNorms[0]] in file) and not (NormTypeName[OtherNorms[1]] in file): # eg: '.txt'
            CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.9 \\textwidth]{'+file+'} \n')
        else:
          CanvasOutputFile[iNormType].write('\\includegraphics[width = 0.9 \\textwidth]{'+file+'} \n')
    #Clear the contents of the OtherNorms array
    OtherNorms[:] = []
    
    #Write the \end document sentence:
    CanvasOutputFile[iNormType].write('\\end{document} \n')
    CanvasOutputFile[iNormType].close()

    #Now create the PDF document using pdflatex!
    os.system('pdflatex -interaction=batchmode '+CanvasOutputFile[iNormType].name)
    print "Want to move file with name : ",(CanvasOutputFile[iNormType].name[:-4]+".pdf")
    shutil.move(CanvasOutputFile[iNormType].name[:-4]+".pdf","..");
    
print "\n --> All information is stored in : ", whichDir
