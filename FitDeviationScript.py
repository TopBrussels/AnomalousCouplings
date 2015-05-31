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
from array import array

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and the number of events in the command line !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop #evts ForceFitUpdate(y/n) TexWanted(y/n)"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop #evts ForceFitUpdate(y/n) TexWanted(y/n)"
  sys.exit()
elif len(sys.argv) == 3:
  print "Need to give the number of considered events !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop #evts ForceFitUpdate(y/n) TexWanted(y/n)"
  sys.exit()
elif len(sys.argv) == 4:
  print "Need to specify whether the performed fits have to be updated (hence run fitDeviationMacro.C)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop #evts ForceFitUpdate(y/n) TexWanted(y/n)"
  sys.exit()
elif len(sys.argv) == 5:
  print "Need to mention whether Tex output is wanted (y/n)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop #evts ForceFitUpdate(y/n) TexWanted(y/n)"

whichDir = sys.argv[1]
KinVariable = sys.argv[2]
nEvts = sys.argv[3]

#Set the 'CreateTexFile' correctly:
if sys.argv[5] == "y" or sys.argv[5] == "yes":
  CreateTexFile = True
  print " ***** Will also create the TexFile *****"
elif sys.argv[5] == "n" or sys.argv[5] == "no":
  CreateTexFile = False
  print " ***** Will only update the Root macro, not create the TexFile ***** "
else:
  print "!!!!! Simple yes/no was needed for TexWanted boolean!!!!!! "
  sys.exit()

if KinVariable != "MTop" and KinVariable != "RVR" and KinVariable != "RgR":
  print "Need to specify which kinematic variable should be considered (MTop, RVR or RgR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the three : ')
print "Interested in directory :",whichDir," (using ",nEvts," events)"

if KinVariable == "RVR":
  #Information about the scanned RVR values and the corresponding cross-section
  VarValues =  ["Re(V_{R}) = -1.5","Re(V_{R}) = -1.0","Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5","Re(V_{R}) = 1.0","Re(V_{R}) = 1.5"]
  Var =        array('d',[-1.5,    -1.0,     -0.5,    -0.3,     -0.2,     -0.1,       -0.05,    0.0,       0.05,     0.1,       0.2,       0.3,       0.5,      1.0,     1.5     ])
  MGXS =       array('d',[122.082, 46.4474,  17.9275, 13.3944,  12.06555, 11.25909,   11.02784, 10.90059,  10.88228, 10.97767,  11.49883,  12.49056,  16.1508,  40.8074, 108.249 ])
  MGXSCut =    array('d',[28.507,  10.63436, 3.95435, 2.92922,  2.62439,  2.4352,     2.38608,  2.35285,   2.35117,  2.37359,   2.49101,   2.72632,   3.58445,  9.36921, 25.4672 ])
  MGXSe =      array('d',[1,       1,        0.01231, 0.009958, 0.009346, 0.00836608, 0.0,      0.0082221, 0.0,      0.0084729, 0.0090198, 0.0087468, 0.011308, 1,       1       ])
  Acceptance = array('d',[1,       1,        0.22164, 0.21742,  0.21672,  0.21737,    0.21614,  0.21670,   0.21531,  0.21677,   0.21437,   0.21793,   0.22205,  1,       1       ])

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  3) Many   : [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1] \n  4) Very wide : [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] \n 5) Wide & Many : [-0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3] \n --> Choose the correct number : ')

  xMinValue = [4, 4, 10, 3, 12]
  NrPointsToRemove = [3, 3, 7, 2, 7]
  KinVar = "Re(V_{R})"
  FitType = "pol4"

  if VarWindow == "1":
    VarValues.pop(14), Var.pop(14), MGXS.pop(14), MGXSe.pop(14), MGXSCut.pop(14), Acceptance.pop(14)
    VarValues.pop(13), Var.pop(13), MGXS.pop(13), MGXSe.pop(13), MGXSCut.pop(13), Acceptance.pop(13)
    VarValues.pop(8),  Var.pop(8),  MGXS.pop(8),  MGXSe.pop(8),  MGXSCut.pop(8),  Acceptance.pop(8)
    VarValues.pop(6),  Var.pop(6),  MGXS.pop(6),  MGXSe.pop(6),  MGXSCut.pop(6),  Acceptance.pop(6)
    VarValues.pop(1),  Var.pop(1),  MGXS.pop(1),  MGXSe.pop(1),  MGXSCut.pop(1),  Acceptance.pop(1)
    VarValues.pop(0),  Var.pop(0),  MGXS.pop(0),  MGXSe.pop(0),  MGXSCut.pop(0),  Acceptance.pop(0)
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    VarValues.pop(14), Var.pop(14), MGXS.pop(14), MGXSe.pop(14), MGXSCut.pop(14), Acceptance.pop(14)
    VarValues.pop(13), Var.pop(13), MGXS.pop(13), MGXSe.pop(13), MGXSCut.pop(13), Acceptance.pop(13)
    VarValues.pop(12), Var.pop(12), MGXS.pop(12), MGXSe.pop(12), MGXSCut.pop(12), Acceptance.pop(12)
    VarValues.pop(2),  Var.pop(2),  MGXS.pop(2),  MGXSe.pop(2),  MGXSCut.pop(2),  Acceptance.pop(2)
    VarValues.pop(1),  Var.pop(1),  MGXS.pop(1),  MGXSe.pop(1),  MGXSCut.pop(1),  Acceptance.pop(1)
    VarValues.pop(0),  Var.pop(0),  MGXS.pop(0),  MGXSe.pop(0),  MGXSCut.pop(0),  Acceptance.pop(0)
    xBin, xLow, xHigh = 13, -0.325, 0.325
  elif VarWindow == "3":
    VarValues =["Re(V_{R}) = -0.1","Re(V_{R}) = -0.09","Re(VR) = -0.08","Re(V_{R}) = -0.07","Re(VR) = -0.06","Re(V_{R}) = -0.05","Re(V_{R}) = -0.04","Re(V_{R}) = -0.03","Re(V_{R}) = -0.02","Re(V_{R}) = -0.01","Re(VR) = 0.0","Re(V_{R}) = 0.01","Re(VR) = 0.02","Re(V_{R}) = 0.03","Re(V_{R}) = 0.04","Re(V_{R}) = 0.05","Re(V_{R}) = 0.06","Re(V_{R}) = 0.07","Re(V_{R}) = 0.08","Re(V_{R}) = 0.09","Re(V_{R}) = 0.1"]
    Var = array('d',[-0.1,     -0.09,    -0.08,    -0.07,    -0.06,    -0.05,    -0.04,    -0.03,   -0.02,    -0.01,    0.0,      0.01,     0.02,     0.03,     0.04,     0.05,     0.06,     0.07,     0.08,     0.09,    0.1])
    MGXS = array('d',[11.25909, 11.19599, 11.14075, 11.10667, 11.06145, 11.02784, 10.99474, 10.9534, 10.93846, 10.91787, 10.90059, 10.89643, 10.86829, 10.87792, 10.87266, 10.88228, 10.87684, 10.89534, 10.91641, 10.9317, 10.97767])
    MGXSCut = array('d',[2.43546,  2.42182,  2.41059,  2.40057,  2.38839,  2.38187,  2.36976,  2.36513, 2.35512,  2.35666,  2.35415,  2.35694,  2.35174,  2.34909,  2.34392,  2.35108,  2.34767,  2.35477,  2.36148,  2.3643,  2.37424])
    xBin, xLow, xHigh = 21, -0.105, 0.105
  elif VarWindow == "4":
    VarValues.pop(11), Var.pop(11), MGXS.pop(11), MGXSe.pop(11), Acceptance.pop(11)
    VarValues.pop(10), Var.pop(10), MGXS.pop(10), MGXSe.pop(10), Acceptance.pop(10)
    VarValues.pop(9),  Var.pop(9),  MGXS.pop(9),  MGXSe.pop(9),  Acceptance.pop(9)
    VarValues.pop(8),  Var.pop(8),  MGXS.pop(8),  MGXSe.pop(8),  Acceptance.pop(8)
    VarValues.pop(6),  Var.pop(6),  MGXS.pop(6),  MGXSe.pop(6),  Acceptance.pop(6)
    VarValues.pop(5),  Var.pop(5),  MGXS.pop(5),  MGXSe.pop(5),  Acceptance.pop(5)
    VarValues.pop(4),  Var.pop(4),  MGXS.pop(4),  MGXSe.pop(4),  Acceptance.pop(4)
    VarValues.pop(3),  Var.pop(3),  MGXS.pop(3),  MGXSe.pop(3),  Acceptance.pop(3)
    xBin, xLow, xHigh = 7, -1.75, 1.75
  elif VarWindow == "5":
    VarValues = ["Re(V_{R}) = -0.3","Re(V_{R}) = -0.275","Re(V_{R}) = -0.25","Re(V_{R}) = -0.225","Re(V_{R}) = -0.2","Re(V_{R}) = -0.175","Re(V_{R}) = -0.15","Re(V_{R}) = -0.125","Re(V_{R}) = -0.1","Re(V_{R}) = -0.075","Re(V_{R}) = -0.05","Re(V_{R}) = -0.025","Re(V_{R}) = 0.0","Re(V_{R}) = 0.025","Re(V_{R}) = 0.05","Re(V_{R}) = 0.075","Re(V_{R}) = 0.1","Re(V_{R}) = 0.125","Re(V_{R}) = 0.15","Re(V_{R}) = 0.175","Re(V_{R}) = 0.2","Re(V_{R}) = 0.225","Re(V_{R}) = 0.25","Re(V_{R}) = 0.275","Re(V_{R}) = 0.3"]
    Var     = array('d',[-0.3,    -0.275, -0.25,    -0.225,   -0.2,     -0.175,   -0.15,    -0.125,  -0.1,     -0.075,   -0.05,    -0.025,  0.0,      0.025,    0.05,     0.075,    0.1,      0.125,    0.15,     0.175,    0.2,      0.225,    0.25,     0.275,    0.3     ])
    MGXS    = array('d',[13.3944, 13.037, 12.66011, 12.37463, 12.06555, 11.83271, 11.60956, 11.4194, 11.25909, 11.12321, 11.02784, 10.9524, 10.90059, 10.87549, 10.88228, 10.93437, 10.97767, 11.07142, 11.17366, 11.32792, 11.49883, 11.69063, 11.90668, 12.18904, 12.49056])
    MGXSCut = array('d',[2.92922, 1,      1,        1,        2.62439,  1,        1,        1,       2.4352,   1,        2.38608,  1,       2.35285,  1,        2.35117,  1,        2.37359,  1,        1,        1,        2.49101,  1,        1,        1,        2.72632 ])
    xBin, xLow, xHigh = 25, -0.3125, 0.3125

elif KinVariable == "MTop":
  #Information about the scanned MTop values and the corresponding cross-section
  VarValues = ["m_{top} = 153","m_{top} = 163","m_{top} = 170","m_{top} = 171","m_{top} = 172","m_{top} = 173","m_{top} = 174","m_{top} = 175","m_{top} = 183","m_{top} = 193"]
  Var =        array('d',[153,      163,      170,       171,       172,      173,        174,        175,       183,       193     ])
  MGXS =       array('d',[8.20916,  9.6299,   10.57123,  10.70485,  10.8257,  10.96469,   11.08428,   11.22448,  12.18068,  13.3046 ])
  MGXSCut =    array('d',[1.35059,  1.85406,  2.21902,   2.27174,   2.32261,  2.38097,    2.43678,    2.49254,   2.93184,   3.50146 ])
  MGXSe =      array('d',[0.006413, 0.007759, 0.0085714, 0.0081430, 0.008789, 0.00816802, 0.00904798, 0.0086538, 0.0093129, 0.010331])
  Acceptance = array('d',[0.16203,  0.19152,  0.21008,   0.21460,   0.21735,  0.21290,    0.21752,    0.22185,   0.23941,   0.26413 ])

  #Select which window of masses was considered!
  VarWindow = raw_input('** Choose the correct mass-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [153, 163, 170, 171, 172, 173, 174, 175, 183, 193] \n  2) Normal : [153, 163, 171, 172, 173, 174, 175, 183, 193] \n  3) Narrow : [171, 172, 173, 174, 175] \n --> Choose the correct number : ')

  xMinValue = [5, 4, 2]
  KinVar = "m_{top}"
  NrPointsToRemove = [3, 3, 1]
  FitType = "pol2"

  if VarWindow == "1":
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "2":
    VarValues.pop(2), Var.pop(2), MGXS.pop(2), MGXSe.pop(2), MGXSCut.pop(2), Acceptance.pop(2)
    xBin, xLow, xHigh = 41, 152.5, 193.5
  elif VarWindow == "3":
    VarValues.pop(9), Var.pop(9), MGXS.pop(9), MGXSe.pop(9), MGXSCut.pop(9), Acceptance.pop(9)
    VarValues.pop(8), Var.pop(8), MGXS.pop(8), MGXSe.pop(8), MGXSCut.pop(8), Acceptance.pop(8)
    VarValues.pop(2), Var.pop(2), MGXS.pop(2), MGXSe.pop(2), MGXSCut.pop(2), Acceptance.pop(2)
    VarValues.pop(1), Var.pop(1), MGXS.pop(1), MGXSe.pop(1), MGXSCut.pop(1), Acceptance.pop(1)
    VarValues.pop(0), Var.pop(0), MGXS.pop(0), MGXSe.pop(0), MGXSCut.pop(0), Acceptance.pop(0)
    xBin, xLow, xHigh = 5, 170.5, 175.5

elif KinVariable == "RgR":
  VarValues = ["Re(g_{R}) = -0.5","Re(g_{R}) = -0.3","Re(g_{R}) = -0.2","Re(g_{R}) = -0.15","Re(g_{R}) = -0.1","Re(g_{R}) = -0.05","Re(g_{R}) = 0.0","Re(g_{R}) = 0.05","Re(g_{R}) = 0.1","Re(g_{R}) = 0.15","Re(g_{R}) = 0.2","Re(g_{R}) = 0.3","Re(g_{R}) = 0.5"]
  Var =     array('d',[-0.5,     -0.3,     -0.2,    -0.15,   -0.1,    -0.05,   0.0,      0.05,    0.1,     0.15,    0.2,     0.3,     0.5     ])
  MGXS =    array('d',[1.8647,   3.36424,  4.92909, 6.02588, 7.34593, 8.94878, 10.89487, 13.1987, 15.9457, 19.1623, 22.9185, 32.2975, 60.5617 ])
  MGXSCut = array('d',[0.462736, 0.807034, 1.14397, 1,       1.6414,  1,       2.35719,  1,       3.34178, 1,       4.67808, 6.42868, 11.61381])

  VarWindow = raw_input('** Choose the correct RgR-window corresponding to the studied file ** \n** Different optiosn are: \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20] \n --> Choose the correct number : ')

  xMinValue = [4, 4]
  KinVar = "Re(g_{R})"
  NrPointsToRemove = [2, 2]
  FitType = "pol2"
 
  if VarWindow == "1": 
    Var.pop(9), VarValues.pop(9), MGXS.pop(9), MGXSCut.pop(9)
    Var.pop(7), VarValues.pop(7), MGXS.pop(7), MGXSCut.pop(7)
    Var.pop(5), VarValues.pop(5), MGXS.pop(5), MGXSCut.pop(5)
    Var.pop(3), VarValues.pop(3), MGXS.pop(3), MGXSCut.pop(3)
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    Var.pop(12), VarValues.pop(12), MGXS.pop(12), MGXSCut.pop(12)
    Var.pop(11), VarValues.pop(11), MGXS.pop(11), MGXSCut.pop(11)
    Var.pop(1),  VarValues.pop(1),  MGXS.pop(1),  MGXSCut.pop(1)
    Var.pop(0),  VarValues.pop(0),  MGXS.pop(0),  MGXSCut.pop(0)
    xBin, xLow, xHigh = 9, -0.225, 0.225

if KinVariable == "RVR" and VarWindow == "3":
  xPos = [15, 20]
  xNeg = [5, 0]
else:
  xPos = [xMinValue[int(VarWindow)-1]+1,xMinValue[int(VarWindow)-1]+2]
  xNeg = [xMinValue[int(VarWindow)-1]-1,xMinValue[int(VarWindow)-1]-2]

print " List of considered Var values is : ",Var    
NrConfigs = len(Var)
xMin = xMinValue[int(VarWindow)-1]
NumberOfPointsToRemove = NrPointsToRemove[int(VarWindow)-1]
xStep = [Var[xNeg[0]]-Var[xNeg[1]], Var[xMin]-Var[xNeg[0]], Var[xPos[0]]-Var[xMin], Var[xPos[1]]-Var[xPos[0]] ]
print "Step size = ",xStep

#File of interest:
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".out") and file.startswith("weights"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
  WeightsFileName = str(whichDir)+''+str(WeightsFileArray[0])
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
  WeightsFileName = str(whichDir)+''+str(WeightsFileArray[int(fileNr)])
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
print "Will be using file : ",WeightsFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"
elif whichDir.find("MGSample") <= len(whichDir)  and whichDir.find("MGSample") > 0:  
  title = "MGSample"
  MGXSCut = MGXS

if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:
  title += "Gen"
  MGXSCut = MGXS
title = title+"_"+KinVariable

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
if not (os.path.exists(os.path.join(whichDir+'StackedCanvasses/')) ) and CreateTexFile == True and RunFitMacro == True:
  os.makedirs(os.path.join(whichDir+'StackedCanvasses/'))

if RunFitMacro == True:
  RootAnalyzer = open('fitDeviationMacro.C','r')
  NewRootAnalyzer = open('output.C','w')
  
  VarValuesLine = 'std::string VarValues[] = {"'
  VarLine = 'double Var[] = {'
  MGXSLine = 'double MGXS[] = {'
  MGXSCutLine = 'double MGXSCut[] = {'
  for ii in range(len(VarValues)):
    if ii < len(VarValues)-1: 
      VarValuesLine = VarValuesLine+str(VarValues[ii])+'","' 
      VarLine = VarLine+str(Var[ii])+','
      MGXSLine = MGXSLine+str(MGXS[ii])+','
      MGXSCutLine = MGXSCutLine+str(MGXSCut[ii])+','
    else:
      VarValuesLine = VarValuesLine+str(VarValues[ii])+'"}; \n'
      VarLine+= str(Var[ii])+'}; \n'
      MGXSLine += str(MGXS[ii])+'}; \n'
      MGXSCutLine += str(MGXSCut[ii])+'}; \n'

  xMinValueLine = 'int xMinValue[] = {'
  for ii in range(len(xMinValue)):
    if ii < len(xMinValue)-1: xMinValueLine += str(xMinValue[ii])+','
    else:                     xMinValueLine += str(xMinValue[ii])+'}; \n'

  xPosLine = 'int xPos[] = {' 
  xNegLine = 'int xNeg[] = {'
  for ii in range(len(xPos)):
    if ii < len(xPos)-1:
      xPosLine += str(xPos[ii])+','
      xNegLine += str(xNeg[ii])+','
    else:
      xPosLine += str(xPos[ii])+'}; \n'
      xNegLine += str(xNeg[ii])+'}; \n'

  for RootLine in RootAnalyzer:
    RootWord = RootLine.split()
    if   re.search( r"int nEvts",             RootLine): NewRootAnalyzer.write('const int nEvts = '+str(nEvts)+'; \n')
    elif re.search( r"int NrConfigs",         RootLine): NewRootAnalyzer.write('const int NrConfigs = '+str(NrConfigs)+'; \n')
    elif re.search( r"int NrToDel",           RootLine): NewRootAnalyzer.write('const int NrToDel = '+str(NumberOfPointsToRemove)+'; \n') 
    elif re.search( r"std::string VarValues", RootLine): NewRootAnalyzer.write(VarValuesLine)
    elif re.search( r"double Var",            RootLine): NewRootAnalyzer.write(VarLine)
    elif re.search( r"double MGXSCut",        RootLine): NewRootAnalyzer.write(MGXSCutLine)
    elif re.search( r"double MGXS",           RootLine): NewRootAnalyzer.write(MGXSLine)
    elif re.search( r"int xBin",              RootLine): NewRootAnalyzer.write('int xBin = '+str(xBin)+'; \n')
    elif re.search( r"float xLow",            RootLine): NewRootAnalyzer.write('float xLow = '+str(xLow)+'; \n')
    elif re.search( r"float xHigh",           RootLine): NewRootAnalyzer.write('float xHigh = '+str(xHigh)+'; \n')
    elif re.search( r"int xMinValue",         RootLine): NewRootAnalyzer.write(xMinValueLine)
    elif re.search( r"std::string KinVar",    RootLine): NewRootAnalyzer.write('std::string KinVar = "'+str(KinVar)+'"; \n')
    elif re.search( r"int VarWindow",         RootLine): NewRootAnalyzer.write('int VarWindow = '+str(VarWindow)+'; \n')
    elif re.search( r"int xPos",              RootLine): NewRootAnalyzer.write(xPosLine)
    elif re.search( r"int xNeg",              RootLine): NewRootAnalyzer.write(xNegLine)
    elif re.search( r"std::ifstream ifs",     RootLine): NewRootAnalyzer.write('  std::ifstream ifs ("'+str(WeightsFileName)+'", std::ifstream::in); \n')
    elif re.search( r"std::string title",     RootLine): NewRootAnalyzer.write('std::string title = "'+str(title)+'"; \n')
    elif re.search( r"string StackedDir",     RootLine): NewRootAnalyzer.write('std::string StackedDir = "'+str(whichDir)+'StackedCanvasses"; \n')
    elif re.search( r"new TF1",               RootLine): 
      if   re.search( r"AllPoints",     RootLine): NewRootAnalyzer.write('  polFit_AllPoints = new TF1(("polFit"+Type+"_AllPoints_Evt"+EvtNumber).c_str(),"'+str(FitType)+'",Var[0],Var[NrConfigs-1]); \n')
      elif re.search( r"ReducedPoints", RootLine): NewRootAnalyzer.write('  polFit_ReducedPoints = new TF1(("polFit"+Type+"_"+sNrRemaining+"ReducedPoints_Evt"+EvtNumber).c_str(),"'+str(FitType)+'",Var[0],Var[NrConfigs-1]); \n')
    elif re.search( r"new TFile",             RootLine):
      if   re.search( r"Tfile",        RootLine): NewRootAnalyzer.write('TFile* Tfile = new TFile("'+str(whichDir)+'FitDeviation_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n')
      elif re.search( r"file_FitDist", RootLine): NewRootAnalyzer.write('TFile* file_FitDist = new TFile("'+str(whichDir)+'FitDistributions_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n') 
    elif re.search( r"bool storeStackedCanvas", RootLine):
      if   CreateTexFile == True: NewRootAnalyzer.write('bool storeStackedCanvas = true; \n')
      else:                       NewRootAnalyzer.write('bool storeStackedCanvas = false; \n')
    else:                                                NewRootAnalyzer.write(RootLine)
  NewRootAnalyzer.close()
  RootAnalyzer.close()
    
  os.rename('output.C','fitDeviationMacro.C'), os.system("root -l -b -q fitDeviationMacro.C+")
  
PerformFitOptAnalyzer = open('PerformFitOptimization.C','r')
NewPerformFitOptAnalyzer = open('fitOptimization.C','w')

for FitOptLine in PerformFitOptAnalyzer:
  FitOptWord = FitOptLine.split()
  if   re.search( r"int nEvts",         FitOptLine): NewPerformFitOptAnalyzer.write('const int nEvts = '+str(nEvts)+'; \n')
  elif re.search( r"int xBin",          FitOptLine): NewPerformFitOptAnalyzer.write('const int xBin = '+str(xBin)+'; \n')
  elif re.search( r"float xLow",        FitOptLine): NewPerformFitOptAnalyzer.write('float xLow = '+str(xLow)+'; \n')
  elif re.search( r"float xHigh",       FitOptLine): NewPerformFitOptAnalyzer.write('float xHigh = '+str(xHigh)+'; \n')
  elif re.search( r"string StackedDir", FitOptLine): NewPerformFitOptAnalyzer.write('std::string StackedDir = "'+str(whichDir)+'StackedCanvasses"; \n')
  elif re.search( r"new TFile",   FitOptLine):
    if   re.search( r"inFile",     FitOptLine): NewPerformFitOptAnalyzer.write('TFile *inFile = new TFile("'+str(whichDir)+'FitDistributions_'+str(title)+'_'+str(nEvts)+'Evts.root","READ"); \n')
    elif re.search( r"outputFile", FitOptLine): NewPerformFitOptAnalyzer.write('TFile *outputFile = new TFile("'+str(whichDir)+'FitOptimizations_'+str(title)+'_'+str(nEvts)+'Evts.root","RECREATE"); \n')
  elif re.search( r"bool storeStackedCanvas", FitOptLine):
    if   CreateTexFile == True: NewPerformFitOptAnalyzer.write('bool storeStackedCanvas = true; \n')
    else:                       NewPerformFitOptAnalyzer.write('bool storeStackedCanvas = false; \n')
  else:                                        NewPerformFitOptAnalyzer.write(FitOptLine)
NewPerformFitOptAnalyzer.close()
PerformFitOptAnalyzer.close()
os.rename('fitOptimization.C','PerformFitOptimization.C'), os.system("root -l -b -q PerformFitOptimization.C+") 

#-- Now store the stacked canvasses in a .text file --#
if CreateTexFile == True and RunFitMacro == True:
  CanvasOutputFile = open(os.path.join(whichDir+'StackedCanvasses/FitDeviationStackCanvas_'+str(title)+'_'+str(nEvts)+'Evts.tex'),'w')

  CanvasOutputFile.write('\\documentclass[a4paper,landscape]{article} \n')
  CanvasOutputFile.write('\\usepackage{graphicx} \n ')
  CanvasOutputFile.write('\\usepackage[top=.5in, bottom=1.25in, left=.5in, right=.5in,landscape]{geometry} \n \n')
  CanvasOutputFile.write('\\begin{document} \n')

  NormType = ["no","XS","acceptance"]
  NormTypeName = ["","XS","Acc"]
  
  #Store the information in the correct directory:
  Canvaslist_dir = os.listdir(os.path.join(whichDir+'StackedCanvasses/'))

  OtherNorms = []
  for iNormType in range(len(NormType)):
    for ii in range(len(NormType)):
      if ii != iNormType:
        OtherNorms.append(ii)
    CanvasOutputFile.write('\\section{Distributions of -ln(L) when '+NormType[iNormType]+' normalisation is applied} \n')
    CanvasOutputFile.write('\n \\centering \n')

    #Include the overall likelihood and secondPol distribution (together one 1 page!):
    if ("TotalLnLik"+NormTypeName[iNormType]+".pdf") in Canvaslist_dir and ("SecondPol"+NormTypeName[iNormType]+".pdf") in Canvaslist_dir:
      CanvasOutputFile.write('\\includegraphics[width = 0.49 \\textwidth]{TotalLnLik'+NormTypeName[iNormType]+'.pdf} \n')
      CanvasOutputFile.write('\\includegraphics[width = 0.49 \\textwidth]{SecondPol'+NormTypeName[iNormType]+'.pdf} \n')

    for file in Canvaslist_dir:
      #Include the stacked canvasses:
      if file.endswith(".pdf") and ("LL"+NormTypeName[iNormType] in file):
        if NormType[iNormType] == "no":
          if not (NormTypeName[OtherNorms[0]] in file) and not (NormTypeName[OtherNorms[1]] in file): # eg: '.txt'
            CanvasOutputFile.write('\\includegraphics[width = 0.9 \\textwidth]{'+file+'} \n')
        else:
          CanvasOutputFile.write('\\includegraphics[width = 0.9 \\textwidth]{'+file+'} \n')
    #Clear the contents of the OtherNorms array
    OtherNorms[:] = []
  
  CanvasOutputFile.write('\\end{document} \n')

print "\n --> All information is stored in : ", whichDir