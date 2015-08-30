#! python
import os
import sys
from math import log,sqrt,pow
from array import array
import ROOT
from ROOT import TH1F,TH2F,TFile,TCanvas,TLegend,gStyle,gROOT,TGraph
ROOT.gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and the number of events in the command line !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()
elif len(sys.argv) == 3:
  print "Need to give the number of considered events !"
  print " Correct syntax is : python FitDeviation.py Events/blabla/ MTop 100"
  sys.exit()

whichDir = sys.argv[1]
KinVariable = sys.argv[2]
nEvts = sys.argv[3]

if KinVariable != "MTop" and KinVariable != "RVR":
  print "Need to specify which kinematic variable should be considered (MTop or RVR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the two : ')
print "Interested in directory :",whichDir," (using ",nEvts," events)"

if KinVariable == "RVR":
  #Information about the scanned RVR values and the corresponding cross-section
  VarValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
  Var =        array('d',[-0.5,    -0.3,     -0.2,     -0.1,       -0.05,    0.0,       0.05,     0.1,       0.2,       0.3,       0.5            ])
  MGXS =       array('d',[17.9275, 13.3944,  12.06555, 11.25909,   11.02784, 10.90059,  10.88228, 10.97767,  11.49883,  12.49056,  16.1508        ])
  MGXSCut =    array('d',[3.95435, 2.92922,  2.62439,  2.4352,     2.38608,  2.35285,   2.35117,  2.37359,   2.49101,   2.72632,   3.58445        ])
  MGXSe =      array('d',[0.01231, 0.009958, 0.009346, 0.00836608, 0.0,      0.0082221, 0.0,      0.0084729, 0.0090198, 0.0087468, 0.0113081652137])
  Acceptance = array('d',[0.22164, 0.21742,  0.21672,  0.21737,    0.21614,  0.21670,   0.21531,  0.21677,   0.21437,   0.21793,   0.22205        ])

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  3) Many   : [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1] \n --> Choose the correct number : ')

  xMinValue = [4, 4, 10]
  KinVar = "Re(V_{R})"

  if VarWindow == "1":
    VarValues.pop(6), Var.pop(6), MGXS.pop(6), MGXSe.pop(6), MGXSCut.pop(6), Acceptance.pop(6)
    VarValues.pop(4), Var.pop(4), MGXS.pop(4), MGXSe.pop(4), MGXSCut.pop(4), Acceptance.pop(4)
    xBin, xLow, xHigh = 11, -0.55, 0.55
  elif VarWindow == "2":
    VarValues.pop(10), Var.pop(10), MGXS.pop(10), MGXSe.pop(10), MGXSCut.pop(10), Acceptance.pop(10)
    VarValues.pop(0),  Var.pop(0),  MGXS.pop(0),  MGXSe.pop(0),  MGXSCut.pop(0),  Acceptance.pop(0)
    xBin, xLow, xHigh = 13, -0.325, 0.325
  elif VarWindow == "3":
    VarValues =["Re(V_{R}) = -0.1","Re(V_{R}) = -0.09","Re(VR) = -0.08","Re(V_{R}) = -0.07","Re(VR) = -0.06","Re(V_{R}) = -0.05","Re(V_{R}) = -0.04","Re(V_{R}) = -0.03","Re(V_{R}) = -0.02","Re(V_{R}) = -0.01","Re(VR) = 0.0","Re(V_{R}) = 0.01","Re(VR) = 0.02","Re(V_{R}) = 0.03","Re(V_{R}) = 0.04","Re(V_{R}) = 0.05","Re(V_{R}) = 0.06","Re(V_{R}) = 0.07","Re(V_{R}) = 0.08","Re(V_{R}) = 0.09","Re(V_{R}) = 0.1"]
    Var = array('d',[-0.1,     -0.09,    -0.08,    -0.07,    -0.06,    -0.05,    -0.04,    -0.03,   -0.02,    -0.01,    0.0,      0.01,     0.02,     0.03,     0.04,     0.05,     0.06,     0.07,     0.08,     0.09,    0.1])
    MGXS = array('d',[11.25909, 11.19599, 11.14075, 11.10667, 11.06145, 11.02784, 10.99474, 10.9534, 10.93846, 10.91787, 10.90059, 10.89643, 10.86829, 10.87792, 10.87266, 10.88228, 10.87684, 10.89534, 10.91641, 10.9317, 10.97767])
    MGXSCut = array('d',[2.43546,  2.42182,  2.41059,  2.40057,  2.38839,  2.38187,  2.36976,  2.36513, 2.35512,  2.35666,  2.35415,  2.35694,  2.35174,  2.34909,  2.34392,  2.35108,  2.34767,  2.35477,  2.36148,  2.3643,  2.37424])
    xBin, xLow, xHigh = 21, -0.105, 0.105
  
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

if KinVariable == "RVR" and VarWindow == "3":
  xPos = [15, 20]
  xNeg = [5, 0]
else:
  xPos = [xMinValue[int(VarWindow)-1]+1,xMinValue[int(VarWindow)-1]+2]
  xNeg = [xMinValue[int(VarWindow)-1]-1,xMinValue[int(VarWindow)-1]-2]

print " List of considered Var values is : ",Var    
NrConfigs = len(Var)
xMin = xMinValue[int(VarWindow)-1]
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
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
  LikelihoodFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
print "Will be using file : ",WeightsFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:
  title += "Gen"
  MGXSCut = MGXS
title = title+"_"+KinVariable

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"FitDeviation_"+title+"_"+nEvts+"Evts.root"),'recreate')
gStyle.SetOptStat(0)

LnLikAll    = TH1F('LnLikAll',   'Distribution of -ln(L) using all events (no normalisation -- '+title+' evts)', xBin,xLow,xHigh)
LnLikXSAll  = TH1F('LnLikXSAll', 'Distribution of -ln(L) using all events (XS normalisation -- '+title+' evts)', xBin,xLow,xHigh)
LnLikAccAll = TH1F('LnLikAccAll','Distribution of -ln(L) using all events (Acc normalisation -- '+title+' evts)',xBin,xLow,xHigh)

YPlusGausTest    = TH1F('YPlusGausTest',   'Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestXS  = TH1F('YPlusGausTestXS', 'Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestAcc = TH1F('YPlusGausTestAcc','Comparison of fit deviation in '+KinVar+' = '+str(Var[xPos[1]])+' and '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.15,0.15)
YPlusGausTestPosScdDer    = TH1F('YPlusGausTestPosScdDer',   'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (no norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSPosScdDer  = TH1F('YPlusGausTestXSPosScdDer', 'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (XS norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccPosScdDer = TH1F('YPlusGausTestAccPosScdDer','Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (Acc norm -- outer 2nd der > 0 -- '+title+' evts)',250,-0.1,0.1)
YPlusGausTestNegScdDer    = TH1F('YPlusGausTestNegScdDer',   'Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (no norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSNegScdDer  = TH1F('YPlusGausTestXSNegScdDer', 'Gaussing fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (XS norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccNegScdDer = TH1F('YPlusGausTestAccNegScdDer','Gaussian fit deviation using '+KinVar+' = '+str(Var[xPos[1]])+' & '+str(Var[xPos[0]])+' (Acc norm -- outer 2nd der < 0 -- '+title+' evts)',250,-0.1,0.1)

YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YPlusPosScdDer    = TH1F('YPlusPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSPosScdDer  = TH1F('YPlusXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccPosScdDer = TH1F('YPlusAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YPlusNegScdDer    = TH1F('YPlusNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSNegScdDer  = TH1F('YPlusXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccNegScdDer = TH1F('YPlusAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YPlusPlus    = TH1F('YPlusPlus',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusXS  = TH1F('YPlusPlusXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusAcc = TH1F('YPlusPlusAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YPlusPlusPosScdDer    = TH1F('YPlusPlusPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSPosScdDer  = TH1F('YPlusPlusXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccPosScdDer = TH1F('YPlusPlusAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YPlusPlusNegScdDer    = TH1F('YPlusPlusNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSNegScdDer  = TH1F('YPlusPlusXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccNegScdDer = TH1F('YPlusPlusAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xPos[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YMin    = TH1F('YMin',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YMinPosScdDer    = TH1F('YMinPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSPosScdDer  = TH1F('YMinXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccPosScdDer = TH1F('YMinAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YMinNegScdDer    = TH1F('YMinNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSNegScdDer  = TH1F('YMinXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccNegScdDer = TH1F('YMinAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YMinMin    = TH1F('YMinMin',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinXS  = TH1F('YMinMinXS', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinAcc = TH1F('YMinMinAcc','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YMinMinPosScdDer    = TH1F('YMinMinPosScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSPosScdDer  = TH1F('YMinMinXSPosScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccPosScdDer = TH1F('YMinMinAccPosScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YMinMinNegScdDer    = TH1F('YMinMinNegScdDer',   'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSNegScdDer  = TH1F('YMinMinXSNegScdDer', 'Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccNegScdDer = TH1F('YMinMinAccNegScdDer','Deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YRelPlus    = TH1F('YRelPlus',   'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusXS  = TH1F('YRelPlusXS', 'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusAcc = TH1F('YRelPlusAcc','Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.1,0.1)

YRelMin    = TH1F('YRelMin',   'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinXS  = TH1F('YRelMinXS', 'Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinAcc = TH1F('YRelMinAcc','Relative deviation from parabolic fit for '+KinVar+' = '+str(Var[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.1,0.1)

FstDer    = TH1F('FirstDer',   'First derivative of -ln(likelihood) distribution -- '+title+' evts', 4,0,4)
FstDer.GetXaxis().SetBinLabel(1,'(y_{DATA}(x='+str(Var[xNeg[1]])+') - y_{DATA}(x='+str(Var[xNeg[0]])+'))/'+str(xStep[0]))
FstDer.GetXaxis().SetBinLabel(2,'(y_{DATA}(x='+str(Var[xNeg[0]])+') - y_{DATA}(x='+str(Var[xMin])+'))/'+str(xStep[1]))
FstDer.GetXaxis().SetBinLabel(3,'(y_{DATA}(x='+str(Var[xMin])+') - y_{DATA}(x='+str(Var[xPos[0]])+'))/'+str(xStep[2]))
FstDer.GetXaxis().SetBinLabel(4,'(y_{DATA}(x='+str(Var[xPos[0]])+') - y_{DATA}(x='+str(Var[xPos[1]])+'))/'+str(xStep[3]))
FstDerXS  = TH1F('FirstDerivativeXS', 'First derivative of -ln(likelihood) distribution (XS normalisation -- '+title+' evts)', 5,-0.25,0.25)
FstDerAcc = TH1F('FirstDerivativeAcc','First derivative of -ln(likelihood) distribution (Acc normalisation -- '+title+' evts)',5,-0.25,0.25)
ScdDerInner    = TH1F('SecondDerivativeInner',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using inner points -- '+title+' evts)', 250,-5,5)
ScdDerXSInner  = TH1F('SecondDerivativeXSInner', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using inner points -- '+title+' evts)', 250,-5,5)
ScdDerAccInner = TH1F('SecondDerivativeAccInner','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using inner points -- '+title+' evts)',250,-5,5)
ScdDerOuter    = TH1F('SecondDerivativeOuter',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using outer points -- '+title+' evts)', 250,-5,5)
ScdDerXSOuter  = TH1F('SecondDerivativeXSOuter', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using outer points -- '+title+' evts)', 250,-5,5)
ScdDerAccOuter = TH1F('SecondDerivativeAccOuter','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using outer points -- '+title+' evts)',250,-5,5)
ScdDerScatter    = TH2F('ScdDerScatterPlot',   'Second derivative of -ln(L) using inner points versus using outer points (no normalisation -- '+title+' evts)', 250,-5,5,250,-5,5)
ScdDerXSScatter  = TH2F('ScdDerXSScatterPlot', 'Second derivative of -ln(L) using inner points versus using outer points (XS normalisation -- '+title+' evts)', 250,-5,5,250,-5,5)
ScdDerAccScatter = TH2F('ScdDerAccScatterPlot','Second derivative of -ln(L) using inner points versus using outer points (Acc normalisation -- '+title+' evts)',250,-5,5,250,-5,5)

FstDerInnerPlusRelToUnc = TH1F('FirstDerInner_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (no norm -- pos inner point)',250,0,25)
FstDerXSInnerPlusRelToUnc = TH1F('FirstDerXSInner_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (XS norm -- pos inner point)',250,0,25)
FstDerAccInnerPlusRelToUnc = TH1F('FirstDerAccInner_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- pos inner point)',250,0,25)
FstDerInnerMinRelToUnc = TH1F('FirstDerInner_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (no norm -- neg inner point)',250,0,25)
FstDerXSInnerMinRelToUnc = TH1F('FirstDerXSInner_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (XS norm -- neg inner point)',250,0,25)
FstDerAccInnerMinRelToUnc = TH1F('FirstDerAccInner_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- neg inner point)',250,0,25)
FstDerOuterPlusRelToUnc = TH1F('FirstDerOuter_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (no norm -- pos outer point)',250,0,25)
FstDerXSOuterPlusRelToUnc = TH1F('FirstDerXSOuter_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (XS norm -- pos outer point)',250,0,25)
FstDerAccOuterPlusRelToUnc = TH1F('FirstDerAccOuter_Plus_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- pos outer point)',250,0,25)
FstDerOuterMinRelToUnc = TH1F('FirstDerOuter_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (no norm -- neg outer point)',250,0,25)
FstDerXSOuterMinRelToUnc = TH1F('FirstDerXSOuter_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (XS norm -- neg outer point)',250,0,25)
FstDerAccOuterMinRelToUnc = TH1F('FirstDerAccouter_Min_RelativeToUnc','First derivative of -ln(L) distribution wrt unc of '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- neg outer point)',250,0,25)

WeightVsUnc = TH2F('WeightVsUnc','Distribution of weight vs uncertainty (ln applied) for '+KinVar+' = '+str(Var[xMin]),250,20,100,250,0,0.2)
WeightVsUncPlus = TH2F('WeightVsUnc_Plus','Distribution of weight vs uncertainty (ln applied) for '+KinVar+' = '+str(Var[xPos[0]]),250,20,100,250,0,0.2)
WeightVsUncMin = TH2F('WeightVsUnc_Min','Distribution of weight vs uncertainty (ln applied) for '+KinVar+' = '+str(Var[xNeg[0]]),250,20,100,250,0,0.2)

LnLikUncDist = TH1F('LnLikUncDist_RVR0','Distribution of uncertainty on -ln(L) of '+KinVar+' = '+str(Var[xMin])+' (all norm -- approx -- '+title+' events)',250,0,0.2)
LnLikUncDistPlus = TH1F('LnLikUncDist_RVRPlus','Distribution of uncertainty on -ln(L) of '+KinVar+' = '+str(Var[xPos[0]])+' (all norm -- approx -- '+title+' events)',250,0,0.2)
LnLikUncDistMin = TH1F('LnLikUncDist_RVRMin','Distribution of uncertainty on -ln(L) of '+KinVar+' = '+str(Var[xNeg[0]])+' (all norm -- approx -- '+title+' events)',250,0,0.2)
LikUncDist    = TH1F('LikUncDist_RVR0','Distribution of uncertainty on Likelihood of '+KinVar+' = '+str(Var[xMin])+' (no norm -- approx -- '+title+' events)',250,0,0.00000000000000000002)
LikXSUncDist  = TH1F('LikXSUncDist_RVR0','Distribution of uncertainty on Likelihood of '+KinVar+' = '+str(Var[xMin])+' (XS norm -- approx -- '+title+' events)',250,0,0.00000000000000000002)
LikAccUncDist = TH1F('LikAccUncDist_RVR0','Distribution of uncertainty on Likelihood of '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- approx -- '+title+' events)',250,0,0.00000000000000000002)
RelLnLikUncDist    = TH1F('RelLnLikUncDist_RVR0','Distribution of #sigma(-ln(L))/-ln(L) for '+KinVar+' = '+str(Var[xMin])+' (no norm -- '+title+' events)',250,0,0.002)
RelLnLikXSUncDist  = TH1F('RelLnLikXSUncDist_RVR0','Distribution of #sigma(-ln(L))/-ln(L) for '+KinVar+' = '+str(Var[xMin])+' (XS norm -- '+title+' events)',250,0,0.002)
RelLnLikAccUncDist = TH1F('RelLnLikAccUncDist_RVR0','Distribution of #sigma(-ln(L))/-ln(L) for '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- '+title+' events)',250,0,0.002)
RelLikUncDist    = TH1F('RelLikUncDist_RVR0','Distribution of #sigma(L)/L for '+KinVar+' = '+str(Var[xMin])+' (no norm -- '+title+' events)',250,0,0.2)
RelLikXSUncDist  = TH1F('RelLikXSUncDist_RVR0','Distribution of #sigma(L)/L for '+KinVar+' = '+str(Var[xMin])+' (XS norm -- '+title+' events)',250,0,0.2)
RelLikAccUncDist = TH1F('RelLikAccUncDist_RVR0','Distribution of #sigma(L)/L for '+KinVar+' = '+str(Var[xMin])+' (Acc norm -- '+title+' events)',250,0,0.2)
#LnLikXSUncDist = TH1F('LnLikUncDist_RVR0','Distribution of uncertainty on -ln(L) of '+KinVar+' = '+str(Var[xMin])+' (no norm -- '+title+' events)',250,0,50)
#LnLikAccUncDist = TH1F('LnLikUncDist_RVR0','Distribution of uncertainty on -ln(L) of '+KinVar+' = '+str(Var[xMin])+' (no norm -- '+title+' events)',250,0,50)

SigmaVarianceDist = TH1F('SigmaVariance','Spread of the uncertainty on the -ln(L) distribution',250,0,0.0002)
AverageSigmaDist = TH1F('AverageSigma','Average value of the uncertainties for each of the considered configurations',250,0,0.2)
LnLikVariationDist = TH1F('LnLikVariation','Difference between the maximum and the minimum value of the -ln(L) distribution',250,0,1)

ProcContrDist = TH1F('ProcentualContribution','Procentual contribution of each of the events to the total likelihood',250,0,1)

TotalFitDevDist    = TH1F('TotalFitDeviation',   'Sum of difference between likelihood and fit value for each point in range (no normalisation -- '+title+' events)', 250,0,5)
TotalFitDevXSDist  = TH1F('TotalFitXSDeviation', 'Sum of difference between likelihood and fit value for each point in range (XS normalisation -- '+title+' events)', 250,0,5)
TotalFitDevAccDist = TH1F('TotalFitAccDeviation','Sum of difference between likelihood and fit value for each point in range (Acc normalisation -- '+title+' events)',250,0,5)
TotalFctDevDist    = TH1F('TotalFctDeviation',   'Sum of difference between likelihood and function value for each point in range (no normalisation -- '+title+' events)', 250,0,5)
TotalFctDevXSDist  = TH1F('TotalFctXSDeviation', 'Sum of difference between likelihood and function value for each point in range (XS normalisation -- '+title+' events)', 250,0,5)
TotalFctDevAccDist = TH1F('TotalFctAccDeviation','Sum of difference between likelihood and function value for each point in range (Acc normalisation -- '+title+' events)',250,0,5)

LnLikDist = TH1F("LnLik","title",xBin,xLow,xHigh)
LnLikDist.SetMarkerStyle(20), LnLikDist.SetLineColor(1), LnLikDist.SetMarkerColor(1), LnLikDist.SetMarkerSize(1.2)
LnLikXSDist = TH1F("LnLikXS","title",xBin,xLow,xHigh)
LnLikXSDist.SetMarkerStyle(21), LnLikXSDist.SetLineColor(3), LnLikXSDist.SetMarkerColor(3), LnLikXSDist.SetMarkerSize(1.2)
LnLikAccDist = TH1F("LnLikAcc","title",xBin,xLow,xHigh)
LnLikAccDist.SetMarkerStyle(22), LnLikAccDist.SetLineColor(4), LnLikAccDist.SetMarkerColor(4), LnLikAccDist.SetMarkerSize(1.2)

FitComp = Tfile.mkdir("FitComparison")
LnLikDir = Tfile.mkdir("LnLikDist")
LnLikXSDir = Tfile.mkdir("LnLikXSDist")
LnLikAccDir = Tfile.mkdir("LnLikAccDist")
LnLikAccDirVarVsUnc = Tfile.mkdir("LnLikAccDist_VarLargerThanAvgUnc")
LnLikAccDirVarVsDUnc = Tfile.mkdir("LnLikAccDist_VarLargerThanTwiceAvgUnc")
FstDerDir = Tfile.mkdir("FirstDerivativeDist")
LnLikFitCanvas = TCanvas('name','title')
print "Name of LnLikFitCanvas is : ",LnLikFitCanvas.GetName()
LnLikFitCanvasName = 'LnLikFitCanvas_Evt'
LnLikFitCanvasTitle = 'Comparing fit function from ROOT fit and algebraic function for event -- '+title+' evts '

#Create the arrays where the likelihood values will be stored
LnLik, LnLikXS, LnLikAcc, Lik, LikXS, LikAcc, LnLikErr, LnLikXSErr, LnLikAccErr, LikErr, LikXSErr, LikAccErr = [], [], [], [], [], [], [], [], [], [], [], []
for ii in range(NrConfigs):
  LnLik.append(0),    LnLikXS.append(0),    LnLikAcc.append(0),    Lik.append(0),    LikXS.append(0),    LikAcc.append(0)
  LnLikErr.append(0), LnLikXSErr.append(0), LnLikAccErr.append(0), LikErr.append(0), LikXSErr.append(0), LikAccErr.append(0)

#---  Create arrays where the events passing specific cuts will be stored  ---#
EvtsPosScdDerInner, EvtsPosScdDerXSInner, EvtsPosScdDerAccInner = [], [], []
EvtsPosScdDerOuter, EvtsPosScdDerXSOuter, EvtsPosScdDerAccOuter = [], [], []
EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc = [], [], []
EvtsWithSmallFctDev = []

nrEvtsWithVarLargerThanAverageUnc, nrEvtsWithVarLargerThanTwiceAverageUnc = 0,0
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WWord[0]) != "#" and str(WWord[3]) != "0.0" :
    for iEvt in range(int(nEvts)):
      if str(WWord[0]) == str(iEvt+1):                                                               #Look at one single event!
        if str(WWord[1]) == "1":
          cHat, cHatXS, cHatAcc = [0,0],[0,0],[0,0]
          bHat, bHatXS, bHatAcc = [0,0],[0,0],[0,0]
          aHat, aHatXS, aHatAcc = [0,0],[0,0],[0,0]
          LnLikDist.SetName("LnLik_Evt"+str(int(iEvt)+1)),       LnLikDist.SetTitle('LnLik distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLikXSDist.SetName("LnLikXS_Evt"+str(int(iEvt)+1)),   LnLikXSDist.SetTitle('LnLikXS distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLikAccDist.SetName("LnLikAcc_Evt"+str(int(iEvt)+1)), LnLikAccDist.SetTitle('LnLikAcc distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')

	Lik[int(WWord[1])-1],    LnLik[int(WWord[1])-1]    = float(WWord[3]),                          -log(float(WWord[3]))
        LikXS[int(WWord[1])-1],  LnLikXS[int(WWord[1])-1]  = float(WWord[3])/MGXS[int(WWord[1])-1],    -log(float(WWord[3])) + log(MGXS[int(WWord[1])-1])
        LikAcc[int(WWord[1])-1], LnLikAcc[int(WWord[1])-1] = float(WWord[3])/MGXSCut[int(WWord[1])-1], -log(float(WWord[3])) + log(MGXSCut[int(WWord[1])-1])
	LikErr[int(WWord[1])-1], LikXSErr[int(WWord[1])-1], LikAccErr[int(WWord[1])-1] = float(WWord[4]), float(WWord[4])/MGXS[int(WWord[1])-1], float(WWord[4])/MGXSCut[int(WWord[1])-1]
	LnLikErr[int(WWord[1])-1] = float(WWord[4])/float(WWord[3])  #Simplification of uncertainty on ln of weight!

        #---  Fill the LnLik histograms for each event and for all events together  ---#
        LnLikDist.SetBinContent(   LnLikDist.FindBin(Var[int(WWord[1])-1]),    LnLik[int(WWord[1])-1]),    LnLikDist.SetBinError(   LnLikDist.FindBin(Var[int(WWord[1])-1]),    LnLikErr[int(WWord[1])-1]) 
        LnLikXSDist.SetBinContent( LnLikXSDist.FindBin(Var[int(WWord[1])-1]),  LnLikXS[int(WWord[1])-1]),  LnLikXSDist.SetBinError( LnLikXSDist.FindBin(Var[int(WWord[1])-1]),  LnLikErr[int(WWord[1])-1])
        LnLikAccDist.SetBinContent(LnLikAccDist.FindBin(Var[int(WWord[1])-1]), LnLikAcc[int(WWord[1])-1]), LnLikAccDist.SetBinError(LnLikAccDist.FindBin(Var[int(WWord[1])-1]), LnLikErr[int(WWord[1])-1])
        LnLikAll.SetBinContent(LnLikAll.FindBin(Var[int(WWord[1])-1]),       LnLikAll.GetBinContent(LnLikAll.FindBin(Var[int(WWord[1])-1]))       + LnLik[int(WWord[1])-1]    )
        LnLikXSAll.SetBinContent(LnLikXSAll.FindBin(Var[int(WWord[1])-1]),   LnLikXSAll.GetBinContent(LnLikXSAll.FindBin(Var[int(WWord[1])-1]))   + LnLikXS[int(WWord[1])-1]  )
        LnLikAccAll.SetBinContent(LnLikAccAll.FindBin(Var[int(WWord[1])-1]), LnLikAccAll.GetBinContent(LnLikAccAll.FindBin(Var[int(WWord[1])-1])) + LnLikAcc[int(WWord[1])-1] )

        #---  Only perform the fit after all configurations are considered!  ---#
        if str(WWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          LnLikFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLikXSFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLikAccFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          FctDevOuter, RelFctDevOuter,FctDevXSOuter, RelFctDevXSOuter , FctDevAccOuter, RelFctDevAccOuter = [],[],[],[],[],[]
          TotalFctDevOuter, TotalRelFctDevOuter, TotalFctDevXSOuter, TotalRelFctDevXSOuter, TotalFctDevAccOuter, TotalRelFctDevAccOuter = 0, 0, 0, 0, 0, 0

          for ii in range(2):
            cHat[ii] = LnLik[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLik[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            cHatXS[ii] = LnLikXS[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikXS[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            cHatAcc[ii] = LnLikAcc[xNeg[ii]]*Var[xPos[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]])) - LnLikAcc[xMin]*Var[xNeg[ii]]*Var[xPos[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xPos[ii]]*Var[xNeg[ii]]*Var[xMin]/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))

            bHat[ii] = LnLik[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLik[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))
            bHatXS[ii] = LnLikXS[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikXS[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))
            bHatAcc[ii] = LnLikAcc[xMin]*(Var[xPos[ii]]+Var[xNeg[ii]])/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xPos[ii]]*(Var[xNeg[ii]]+Var[xMin])/((Var[xPos[ii]]-Var[xMin])*(Var[xPos[ii]]-Var[xNeg[ii]]))-LnLikAcc[xNeg[ii]]*(Var[xMin]+Var[xPos[ii]])/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xMin]-Var[xNeg[ii]]))

            aHat[ii] = LnLik[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLik[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLik[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            aHatXS[ii] = LnLikXS[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikXS[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikXS[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))
            aHatAcc[ii] = LnLikAcc[xPos[ii]]/((Var[xPos[ii]]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) - LnLikAcc[xMin]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xMin])) + LnLikAcc[xNeg[ii]]/((Var[xMin]-Var[xNeg[ii]])*(Var[xPos[ii]]-Var[xNeg[ii]]))

            #-- Now calculate each of the 9 Var configurations using this parabola --#
            for var in range(NrConfigs):
              LnLikFunction[ii].append(    aHat[ii]*Var[var]*Var[var]+bHat[ii]*Var[var]+cHat[ii])
              LnLikXSFunction[ii].append(  aHatXS[ii]*Var[var]*Var[var]+bHatXS[ii]*Var[var]+cHatXS[ii])
              LnLikAccFunction[ii].append( aHatAcc[ii]*Var[var]*Var[var]+bHatAcc[ii]*Var[var]+cHatAcc[ii])

	      if ii == 1:
                FctDevOuter.append(abs(LnLikAccFunction[1][var]-LnLik[var])),       RelFctDevOuter.append(FctDevOuter[var]/LnLik[var])
                FctDevXSOuter.append(abs(LnLikXSFunction[1][var]-LnLikXS[var])),    RelFctDevXSOuter.append(FctDevXSOuter[var]/LnLikXS[var])
                FctDevAccOuter.append(abs(LnLikAccFunction[1][var]-LnLikAcc[var])), RelFctDevAccOuter.append(FctDevAccOuter[var]/LnLikAcc[var])
                TotalFctDevOuter,    TotalRelFctDevOuter    = TotalFctDevOuter+FctDevOuter[var], TotalRelFctDevOuter+RelFctDevOuter[var]
                TotalFctDevXSOuter,  TotalRelFctDevXSOuter  = TotalFctDevXSOuter+FctDevXSOuter[var], TotalRelFctDevXSOuter+RelFctDevXSOuter[var]
                TotalFctDevAccOuter, TotalRelFctDevAccOuter = TotalFctDevAccOuter+FctDevAccOuter[var], TotalRelFctDevAccOuter+RelFctDevAccOuter[var]

          #---  Save the LnLik distributions (event-per-event) together with the function "fit" ---#
          LnLikFctOuter, LnLikXSFctOuter, LnLikAccFctOuter = TGraph(NrConfigs,Var,LnLikFunction[1]), TGraph(NrConfigs,Var,LnLikXSFunction[1]), TGraph(NrConfigs,Var,LnLikAccFunction[1])
          LnLikFctOuter.SetMarkerColor(2),LnLikFctOuter.SetLineColor(2),LnLikXSFctOuter.SetMarkerColor(2),LnLikXSFctOuter.SetLineColor(2),LnLikAccFctOuter.SetMarkerColor(2),LnLikAccFctOuter.SetLineColor(2)
          LnLikFctOuter.SetTitle('LnLik for event '+str(iEvt+1)+' -- Outer points used (Fct deviation is '+str(TotalFctDevOuter)+' -- '+str(TotalRelFctDevOuter)+')')          
          LnLikAccFctOuter.SetTitle('LnLikXS for event '+str(iEvt+1)+' -- Outer points used (Fct deviation is '+str(TotalFctDevXSOuter)+' -- '+str(TotalRelFctDevXSOuter)+')')
          LnLikAccFctOuter.SetTitle('LnLikAcc for event '+str(iEvt+1)+' -- Outer points used (Fct deviation is '+str(TotalFctDevAccOuter)+' -- '+str(TotalRelFctDevAccOuter)+')')
          LnLikCanv, LnLikXSCanv, LnLikAccCanv = TCanvas('LnLikCanv_'+str(int(iEvt)+1),'LnLik'), TCanvas('LnLikXSCanv_'+str(int(iEvt)+1),'LnLikXS'), TCanvas('LnLikAccCanv_'+str(int(iEvt)+1),'LnLikAcc')
          
          LnLikCanv.cd(),    LnLikFctOuter.Draw('AC*'),    LnLikDist.Draw('samep'),    LnLikDir.cd(),    LnLikCanv.Write()
          LnLikXSCanv.cd(),  LnLikXSFctOuter.Draw('AC*'),  LnLikXSDist.Draw('samep'),  LnLikXSDir.cd(),  LnLikXSCanv.Write()
          LnLikAccCanv.cd(), LnLikAccFctOuter.Draw('AC*'), LnLikAccDist.Draw('samep'), LnLikAccDir.cd(), LnLikAccCanv.Write()

          #LnLikGraphInner = TGraph(9, Var, LnLikFunction[0])
          #LnLikGraphInner.SetTitle('Comparison between ROOT fit and algebraic method -- '+title+' evts')
          #LnLikGraphInner.GetYaxis().SetTitle('#Delta ln(Likelihood)')
          #LnLikGraphInner.GetXaxis().SetTitle(KinVar)
          #LnLikGraphInner.SetMarkerColor(6), LnLikGraphInner.SetLineColor(6), LnLikGraphInner.SetMarkerSize(1.2)
          #LnLikGraphOuter = TGraph(9,Var, LnLikFunction[1])
          #LnLikGraphOuter.SetMarkerColor(7), LnLikGraphOuter.SetLineColor(7), LnLikGraphOuter.SetMarkerSize(1.2)	
          #---  Compare with TF1 fit from ROOT  ---#
          #LnLikDist.Fit("pol2","Q")
          #LnLikFitCanvas.SetName(LnLikFitCanvasName+str(int(iEvt)+1)), LnLikFitCanvas.SetTitle(LnLikFitCanvasTitle+str(int(iEvt)+1))
          #legendLnLik = TLegend(0.75,0.7,0.95,0.95)
          #LnLikFitCanvas.cd()
          #legendLnLik.AddEntry(LnLikGraphInner,'Fit result using algebraic function (x = '+str(Var[xNeg[0]])+'/'+str(Var[xMin])+'/'+str(Var[xPos[0]])+')',"p")
          #legendLnLik.AddEntry(LnLikGraphOuter,'Fit result using algebraic function (x = '+str(Var[xNeg[1]])+'/'+str(Var[xMin])+'/'+str(Var[xPos[1]])+')',"p")
          #legendLnLik.AddEntry(LnLikDist,'Obtained LnLik values with ROOT fit',"p")
          #LnLikGraphInner.Draw("AC*"), LnLikGraphOuter.Draw("C*"), LnLikDist.Draw("samep"), legendLnLik.Draw()
          #FitComp.cd(), LnLikFitCanvas.Write()

          #---  Calculate the first derivative distribution and save them in event-by-event plot  ---#
          FstDer.SetName('FirstDerivative_Evt'+str(int(iEvt)+1))
          FstDer.SetTitle('First derivative of -ln(likelihood) distribution for event'+str(int(iEvt)+1)+' (no normalisation) -- '+title+' evts') 
          for ii in range(FstDer.GetNbinsX()):
            FstDer.SetBinContent(ii+1, (LnLik[ii+2] - LnLik[ii+3])/xStep[ii])
          FstDerDir.cd(), FstDer.Write()

          #---  Calculate the positive and negative deviation  ---#
          yPlusPlus = [LnLik[xPos[1]] - LnLikFunction[0][xPos[1]], LnLikXS[xPos[1]] - LnLikXSFunction[0][xPos[1]], LnLikAcc[xPos[1]] - LnLikAccFunction[0][xPos[1]]]
          yMinMin   = [LnLik[xNeg[1]] - LnLikFunction[0][xNeg[1]], LnLikXS[xNeg[1]] - LnLikXSFunction[0][xNeg[1]], LnLikAcc[xNeg[1]] - LnLikAccFunction[0][xNeg[1]]]
          yPlus = [LnLik[xPos[0]] - LnLikFunction[1][xPos[0]], LnLikXS[xPos[0]] - LnLikXSFunction[1][xPos[0]], LnLikAcc[xPos[0]] - LnLikAccFunction[1][xPos[0]]]
          yMin  = [LnLik[xNeg[0]] - LnLikFunction[1][xNeg[0]], LnLikXS[xNeg[0]] - LnLikXSFunction[1][xNeg[0]], LnLikAcc[xNeg[0]] - LnLikAccFunction[1][xNeg[0]]]
          scdDerInner = [LnLik[xNeg[0]]-2*LnLik[xMin]+LnLik[xPos[0]], LnLikXS[xNeg[0]]-2*LnLikXS[xMin]+LnLikXS[xPos[0]], LnLikAcc[xNeg[0]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[0]]]	
          scdDerOuter   = [LnLik[xNeg[1]]-2*LnLik[xMin]+LnLik[xPos[1]], LnLikXS[xNeg[1]]-2*LnLikXS[xMin]+LnLikXS[xPos[1]], LnLikAcc[xNeg[1]]-2*LnLikAcc[xMin]+LnLikAcc[xPos[1]]]
          fstDerInnerPlusRelToUnc = [abs(Lik[xMin]-Lik[xPos[0]])/sqrt(pow(LikErr[xMin],2)+pow(LikErr[xPos[0]],2)),abs(LikXS[xMin]-LikXS[xPos[0]])/sqrt(pow(LikXSErr[xMin],2)+pow(LikXSErr[xPos[0]],2)), abs(LikAcc[xMin]-LikAcc[xPos[0]])/sqrt(pow(LikAccErr[xMin],2)+pow(LikAccErr[xPos[0]],2))]
          fstDerInnerMinRelToUnc =  [abs(Lik[xNeg[0]]-Lik[xMin])/sqrt(pow(LikErr[xMin],2)+pow(LikErr[xNeg[0]],2)), abs(LikXS[xNeg[0]]-LikXS[xMin])/sqrt(pow(LikXSErr[xMin],2)+pow(LikXSErr[xNeg[0]],2)), abs(LikAcc[xNeg[0]]-LikAcc[xMin])/sqrt(pow(LikAccErr[xMin],2)+pow(LikAccErr[xNeg[0]],2))]
          fstDerOuterPlusRelToUnc = [abs(Lik[xMin]-Lik[xPos[1]])/sqrt(pow(LikErr[xMin],2)+pow(LikErr[xPos[1]],2)), abs(LikXS[xMin]-LikXS[xPos[1]])/sqrt(pow(LikXSErr[xMin],2)+pow(LikXSErr[xPos[1]],2)), abs(LikAcc[xMin]-LikAcc[xPos[1]])/sqrt(pow(LikAccErr[xMin],2)+pow(LikAccErr[xPos[1]],2))]
          fstDerOuterMinRelToUnc =  [abs(Lik[xNeg[1]]-Lik[xMin])/sqrt(pow(LikErr[xMin],2)+pow(LikErr[xNeg[1]],2)), abs(LikXS[xNeg[1]]-LikXS[xMin])/sqrt(pow(LikXSErr[xMin],2)+pow(LikXSErr[xNeg[1]],2)), abs(LikAcc[xNeg[1]]-LikAcc[xMin])/sqrt(pow(LikAccErr[xMin],2)+pow(LikAccErr[xNeg[1]],2))]

          #-- Check the spread of the uncertainties for each event ! --#
          AverageSigma = 0
          for ii in range(NrConfigs): AverageSigma += LnLikErr[ii]
          AverageSigma = AverageSigma/NrConfigs	
          SigmaVariance = 0
          for jj in range(NrConfigs): SigmaVariance += pow((LnLikErr[jj]-AverageSigma),2)
          SigmaVariance = SigmaVariance/(NrConfigs-1)

          #-- Calculate the difference between the minimum and maximum of the -ln(L) --#
          maxLnLik = 0
          minLnLik = LnLikAcc[0]
          for ii in range(NrConfigs):
            if LnLikAcc[ii] > maxLnLik: maxLnLik = LnLikAcc[ii]
            if LnLikAcc[ii] < minLnLik: minLnLik = LnLikAcc[ii]
          LnLikVariation = maxLnLik - minLnLik
          if LnLikVariation < 3*AverageSigma:
            LnLikAccDirVarVsUnc.cd(), LnLikAccDist.Write()
            nrEvtsWithVarLargerThanAverageUnc += 1
          if LnLikVariation > 3*AverageSigma and LnLikVariation < 10*AverageSigma and scdDerOuter[2] > 0.0:
            LnLikAccDirVarVsDUnc.cd(), LnLikAccDist.Write()
            nrEvtsWithVarLargerThanTwiceAverageUnc += 1
          #print WWord[0],") Comparison of LnLikVariation and AverageSigma : ", LnLikVariation," vs ",AverageSigma

          #--- Fill the histograms  ---#
          YPlus.Fill(yPlus[0]),                                     YPlusXS.Fill(yPlus[1]),                                     YPlusAcc.Fill(yPlus[2])
          YPlusPlus.Fill(yPlusPlus[0]),                             YPlusPlusXS.Fill(yPlusPlus[1]),                             YPlusPlusAcc.Fill(yPlusPlus[2])
          YPlusGausTest.Fill(yPlus[0] + yPlusPlus[0]/4),            YPlusGausTestXS.Fill(yPlus[1] + yPlusPlus[1]/4),            YPlusGausTestAcc.Fill(yPlus[2] + yPlusPlus[2]/4)
          #YRelPlus.Fill(yPlus[0]/LnLik[6]),                         YRelPlusXS.Fill(yPlus[1]/LnLikXS[6]),                       YRelPlusAcc.Fill(yPlus[2]/LnLikAcc[6])
          YMin.Fill(yMin[0]),                                       YMinXS.Fill(yMin[1]),                                       YMinAcc.Fill(yMin[2])
          YMinMin.Fill(yMinMin[0]),                                 YMinMinXS.Fill(yMinMin[1]),                                 YMinMinAcc.Fill(yMinMin[2])
          #YRelMin.Fill(yMin[0]/LnLik[2]),                           YRelMinXS.Fill(yMin[1]/LnLikXS[2]),                         YRelMinAcc.Fill(yMin[2]/LnLikAcc[2])
          ScdDerInner.Fill(scdDerInner[0]),                         ScdDerXSInner.Fill(scdDerInner[1]),                         ScdDerAccInner.Fill(scdDerInner[2])
          ScdDerOuter.Fill(scdDerOuter[0]),                         ScdDerXSOuter.Fill(scdDerOuter[1]),                         ScdDerAccOuter.Fill(scdDerOuter[2])
          ScdDerScatter.Fill(scdDerOuter[0], scdDerInner[0]),       ScdDerXSScatter.Fill(scdDerOuter[1], scdDerInner[1]),       ScdDerAccScatter.Fill(scdDerOuter[2], scdDerInner[2])
          FstDerInnerPlusRelToUnc.Fill(fstDerInnerPlusRelToUnc[0]), FstDerXSInnerPlusRelToUnc.Fill(fstDerInnerPlusRelToUnc[1]), FstDerAccInnerPlusRelToUnc.Fill(fstDerInnerPlusRelToUnc[2])
          FstDerInnerMinRelToUnc.Fill(fstDerInnerMinRelToUnc[0]),   FstDerXSInnerMinRelToUnc.Fill(fstDerInnerMinRelToUnc[1]),   FstDerAccInnerMinRelToUnc.Fill(fstDerInnerMinRelToUnc[2])
          FstDerOuterPlusRelToUnc.Fill(fstDerOuterPlusRelToUnc[0]), FstDerXSOuterPlusRelToUnc.Fill(fstDerOuterPlusRelToUnc[1]), FstDerAccOuterPlusRelToUnc.Fill(fstDerOuterPlusRelToUnc[2])
          FstDerOuterMinRelToUnc.Fill(fstDerOuterMinRelToUnc[0]),   FstDerXSOuterMinRelToUnc.Fill(fstDerOuterMinRelToUnc[1]),   FstDerAccOuterMinRelToUnc.Fill(fstDerOuterMinRelToUnc[2])
          LnLikUncDist.Fill(LnLikErr[xMin]),                       LnLikUncDistPlus.Fill(LnLikErr[xPos[0]]),                  LnLikUncDistMin.Fill(LnLikErr[xNeg[0]])
          WeightVsUnc.Fill(LnLik[xMin],LnLikErr[xMin]),             WeightVsUncPlus.Fill(LnLik[xPos[0]],LnLikErr[xPos[0]]),     WeightVsUncMin.Fill(LnLik[xNeg[0]],LnLikErr[xNeg[0]])
          LikUncDist.Fill(LikErr[xMin]),                           LikXSUncDist.Fill(LikXSErr[xMin]),                         LikAccUncDist.Fill(LikAccErr[xMin])
          RelLnLikUncDist.Fill(LnLikErr[xMin]/LnLik[xMin]),        RelLnLikXSUncDist.Fill(LnLikErr[xMin]/LnLikXS[xMin]),      RelLnLikAccUncDist.Fill(LnLikErr[xMin]/LnLikAcc[xMin])
          RelLikUncDist.Fill(LikErr[xMin]/Lik[xMin]),              RelLikXSUncDist.Fill(LikXSErr[xMin]/LikXS[xMin]),          RelLikAccUncDist.Fill(LikAccErr[xMin]/LikAcc[xMin])
          SigmaVarianceDist.Fill(SigmaVariance),                   AverageSigmaDist.Fill(AverageSigma),                       LnLikVariationDist.Fill(LnLikVariation)
          TotalFctDevDist.Fill(TotalFctDevAccOuter)
          if TotalFctDevAccOuter > 5: print "Overflow found for TotalFctDevDist : ",TotalFctDevAccOuter

          #-- Check for presence of high overflow! --#
          if LikErr[xMin] > 0.00000000000000000002: print "Overflow found for LikUncDist : ",LikErr[xMin]
          if LnLikErr[xMin]/LnLik[xMin] > 0.002: print "Overflow found for RelLnLikUncDist :", LnLikErr[xMin]/LnLik[xMin]
          if LnLikErr[xMin] > 0.2: print " Overflow found for LnLikUncDist : ",LnLikErr[xMin]
          if LikErr[xMin]/Lik[xMin] > 0.2: print "Overflow found for RelLikUncDist : ",LikErr[xMin]/Lik[xMin]

          #-- Apply cut on YPlusGausTest --#
          if (yPlus[0] + yPlusPlus[0]/4) <= 0.025 and (yPlus[0] + yPlusPlus[0]/4) >= -0.025: EvtsWithYPlusGausSmall.append(iEvt+1)
          if (yPlus[1] + yPlusPlus[1]/4) <= 0.025 and (yPlus[1] + yPlusPlus[1]/4) >= -0.025: EvtsWithYPlusGausSmallXS.append(iEvt+1)
          if (yPlus[2] + yPlusPlus[2]/4) <= 0.025 and (yPlus[2] + yPlusPlus[2]/4) >= -0.025: EvtsWithYPlusGausSmallAcc.append(iEvt+1)
          #-- Apply cut on TotalFctDeviation --#
          if TotalFctDevAccOuter <= 0.5: EvtsWithSmallFctDev.append(iEvt+1)
          #-- Apply cut on ScdDer (using inner Var points) --#
          if scdDerInner[0] > 0.0: EvtsPosScdDerInner.append(iEvt+1)
          if scdDerInner[1] > 0.0: EvtsPosScdDerXSInner.append(iEvt+1)
          if scdDerInner[2] > 0.0: EvtsPosScdDerAccInner.append(iEvt+1)
          #-- Apply cut on ScdDer (using outer Var points) --#
          if scdDerOuter[0] > 0.0:
            EvtsPosScdDerOuter.append(iEvt+1)
            YPlusGausTestPosScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusPosScdDer.Fill(yPlus[0]), YPlusPlusPosScdDer.Fill(yPlusPlus[0])
            YMinPosScdDer.Fill(yMin[0]),   YMinMinPosScdDer.Fill(yMinMin[0])
          else:
            YPlusGausTestNegScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusNegScdDer.Fill(yPlus[0]), YPlusPlusNegScdDer.Fill(yPlusPlus[0])
            YMinNegScdDer.Fill(yMin[0]),   YMinMinNegScdDer.Fill(yMinMin[0])

          if scdDerOuter[1] > 0.0:
            EvtsPosScdDerXSOuter.append(iEvt+1)
            YPlusGausTestXSPosScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSPosScdDer.Fill(yPlus[1]), YPlusPlusXSPosScdDer.Fill(yPlusPlus[1])
            YMinXSPosScdDer.Fill(yMin[1]),   YMinMinXSPosScdDer.Fill(yMinMin[1])
          else:
            YPlusGausTestXSNegScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSNegScdDer.Fill(yPlus[1]), YPlusPlusXSNegScdDer.Fill(yPlusPlus[1])
            YMinXSNegScdDer.Fill(yMin[1]),   YMinMinXSNegScdDer.Fill(yMinMin[1])

          if scdDerOuter[2] > 0.0:
            EvtsPosScdDerAccOuter.append(iEvt+1)
            YPlusGausTestAccPosScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccPosScdDer.Fill(yPlus[2]), YPlusPlusAccPosScdDer.Fill(yPlusPlus[2])
            YMinAccPosScdDer.Fill(yMin[2]),   YMinMinAccPosScdDer.Fill(yMinMin[2])
          else:
            YPlusGausTestAccNegScdDer.Fill(yPlus[2] + yPlusPlus[2]/4)
            YPlusAccNegScdDer.Fill(yPlus[2]), YPlusPlusAccNegScdDer.Fill(yPlusPlus[2])
            YMinAccNegScdDer.Fill(yMin[2]),   YMinMinAccNegScdDer.Fill(yMinMin[2])

print "Nr Evts with variation of -ln(L) > average uncertainty = ", nrEvtsWithVarLargerThanAverageUnc, nrEvtsWithVarLargerThanTwiceAverageUnc

#--- Save all the histograms containing information about all the events! ---#
Tfile.cd()
LnLikAll.Write(),                LnLikXSAll.Write(),                LnLikAccAll.Write()
YPlusGausTest.Write(),           YPlusGausTestXS.Write(),           YPlusGausTestAcc.Write()  
YPlus.Write(),                   YPlusXS.Write(),                   YPlusAcc.Write()
YPlusPlus.Write(),               YPlusPlusXS.Write(),               YPlusPlusAcc.Write()
YRelPlus.Write(),                YRelPlusXS.Write(),                YRelPlusAcc.Write()
YMin.Write(),                    YMinXS.Write(),                    YMinAcc.Write()
YMinMin.Write(),                 YMinMinXS.Write(),                 YMinMinAcc.Write()
YRelMin.Write(),                 YRelMinXS.Write(),                 YRelMinAcc.Write()
ScdDerInner.Write(),             ScdDerXSInner.Write(),             ScdDerAccInner.Write()
ScdDerOuter.Write(),             ScdDerXSOuter.Write(),             ScdDerAccOuter.Write()
ScdDerScatter.Write(),           ScdDerXSScatter.Write(),           ScdDerAccScatter.Write()
FstDerInnerPlusRelToUnc.Write(), FstDerXSInnerPlusRelToUnc.Write(), FstDerAccInnerPlusRelToUnc.Write()
FstDerInnerMinRelToUnc.Write(),  FstDerXSInnerMinRelToUnc.Write(),  FstDerAccInnerMinRelToUnc.Write()
FstDerOuterPlusRelToUnc.Write(), FstDerXSOuterPlusRelToUnc.Write(), FstDerAccOuterPlusRelToUnc.Write()
FstDerOuterMinRelToUnc.Write(),  FstDerXSOuterMinRelToUnc.Write(),  FstDerAccOuterMinRelToUnc.Write()
LnLikUncDist.Write(),            LnLikUncDistPlus.Write(),          LnLikUncDistMin.Write()
LikUncDist.Write(),              LikXSUncDist.Write(),              LikAccUncDist.Write()
RelLnLikUncDist.Write(),         RelLnLikXSUncDist.Write(),         RelLnLikAccUncDist.Write()
RelLikUncDist.Write(),           RelLikXSUncDist.Write(),           RelLikAccUncDist.Write()
TotalFctDevDist.Write(),         TotalFctDevXSDist.Write(),         TotalFctDevAccDist.Write()
WeightVsUnc.Write()
WeightVsUncPlus.Write()
WeightVsUncMin.Write()
SigmaVarianceDist.Write()
AverageSigmaDist.Write()
LnLikVariationDist.Write()

#---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---#
print "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = ",str(Var[xNeg[0]]),"/",str(Var[xMin]),"/",str(Var[xPos[0]]),") :",len(EvtsPosScdDerInner),", ",len(EvtsPosScdDerXSInner)," & ",len(EvtsPosScdDerAccInner)
print "Nr of events with 2nd derivative > 0 (LnLik, LnLikXS & LnLikAcc -- using x = ",str(Var[xNeg[1]]),"/",str(Var[xMin]),"/",str(Var[xPos[1]]),") :",len(EvtsPosScdDerOuter),", ",len(EvtsPosScdDerXSOuter)," & ",len(EvtsPosScdDerAccOuter)
print "Nr of events with Gaussiaanse vergelijking voor + (LnLik, LnLikXS & LnLikAcc) ", len(EvtsWithYPlusGausSmall),", ",len(EvtsWithYPlusGausSmallXS)," & ",len(EvtsWithYPlusGausSmallAcc)
print "Nr of events with total function deviation < 0.5 : ",len(EvtsWithSmallFctDev)

LLPosScdDerInner, LLNegScdDerInner, LLXSPosScdDerInner, LLXSNegScdDerInner, LLAccPosScdDerInner, LLAccNegScdDerInner = [],[],[],[],[],[]
LLPosScdDerOuter, LLNegScdDerOuter, LLXSPosScdDerOuter, LLXSNegScdDerOuter, LLAccPosScdDerOuter, LLAccNegScdDerOuter = [],[],[],[],[],[]
LLPosScdDerBoth,  LLNegScdDerBoth,  LLXSPosScdDerBoth,  LLXSNegScdDerBoth,  LLAccPosScdDerBoth,  LLAccNegScdDerBoth  = [],[],[],[],[],[]
LLSmallFctDev = []
for ii in range(NrConfigs):
  LLPosScdDerInner.append(0), LLNegScdDerInner.append(0), LLXSPosScdDerInner.append(0), LLXSNegScdDerInner.append(0), LLAccPosScdDerInner.append(0), LLAccNegScdDerInner.append(0)
  LLPosScdDerOuter.append(0), LLNegScdDerOuter.append(0), LLXSPosScdDerOuter.append(0), LLXSNegScdDerOuter.append(0), LLAccPosScdDerOuter.append(0), LLAccNegScdDerOuter.append(0)
  LLPosScdDerBoth.append(0),  LLNegScdDerBoth.append(0),  LLXSPosScdDerBoth.append(0),  LLXSNegScdDerBoth.append(0),  LLAccPosScdDerBoth.append(0),  LLAccNegScdDerBoth.append(0)
  LLSmallFctDev.append(0)

EvtsPosScdDerBoth, EvtsPosScdDerXSBoth, EvtsPosScdDerAccBoth = 0, 0, 0
#Loop over all lines in weights file:
for LikelihoodLine in LikelihoodFile:
  LWord = LikelihoodLine.split()
  #Only interested in files starting with a number
  if str(LWord[0]) != "#" and str(LWord[3]) != "0.0" :
    #--- Separate the events with small function deviation ---#
    if int(LWord[0]) in EvtsWithSmallFctDev: LLSmallFctDev[int(LWord[1])-1] = LLSmallFctDev[int(LWord[1])-1] - log(float(LWord[3])) + log(MGXSCut[int(LWord[1])-1])

    #---  Separate the events with positive and negative second derivative (using both inner and outer Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner and int(LWord[0]) in EvtsPosScdDerOuter:
      LLPosScdDerBoth[int(LWord[1])-1] = LLPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))
      if LWord[1] == "1": EvtsPosScdDerBoth += 1
    if int(LWord[0]) in EvtsPosScdDerXSInner and int(LWord[0]) in EvtsPosScdDerXSOuter:
      LLXSPosScdDerBoth[int(LWord[1])-1] = LLXSPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerXSBoth += 1
    if int(LWord[0]) in EvtsPosScdDerAccInner and int(LWord[0]) in EvtsPosScdDerAccOuter:
      LLAccPosScdDerBoth[int(LWord[1])-1] = LLAccPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXSCut[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerAccBoth += 1
    #---  Separate the events with positive and negative second derivative (using inner Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner:    LLPosScdDerInner[int(LWord[1])-1] = LLPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                      LLNegScdDerInner[int(LWord[1])-1] = LLNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSInner:  LLXSPosScdDerInner[int(LWord[1])-1] = LLXSPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                      LLXSNegScdDerInner[int(LWord[1])-1] = LLXSNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccInner: LLAccPosScdDerInner[int(LWord[1])-1] = LLAccPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXSCut[int(LWord[1])-1])
    else:                                      LLAccNegScdDerInner[int(LWord[1])-1] = LLAccNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXSCut[int(LWord[1])-1])
    #---  Separate the events with positive and negative second derivative (using outer Var points) ---#
    if int(LWord[0]) in EvtsPosScdDerOuter:    LLPosScdDerOuter[int(LWord[1])-1] = LLPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                      LLNegScdDerOuter[int(LWord[1])-1] = LLNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSOuter:  LLXSPosScdDerOuter[int(LWord[1])-1] = LLXSPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                      LLXSNegScdDerOuter[int(LWord[1])-1] = LLXSNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccOuter: LLAccPosScdDerOuter[int(LWord[1])-1] = LLAccPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXSCut[int(LWord[1])-1])
    else:                                      LLAccNegScdDerOuter[int(LWord[1])-1] = LLAccNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXSCut[int(LWord[1])-1])

NrEvtsNegScdDerInner, NrEvtsNegScdDerXSInner, NrEvtsNegScdDerAccInner = int(nEvts)-len(EvtsPosScdDerInner), int(nEvts) - len(EvtsPosScdDerXSInner), int(nEvts) - len(EvtsPosScdDerAccInner)
NrEvtsNegScdDerOuter, NrEvtsNegScdDerXSOuter, NrEvtsNegScdDerAccOuter = int(nEvts)-len(EvtsPosScdDerOuter), int(nEvts) - len(EvtsPosScdDerXSOuter), int(nEvts) - len(EvtsPosScdDerAccOuter)
LLPosScdDerDistBoth     = TH1F('LLPosScdDerBoth',    '-ln(L) when both 2nd derivatives > 0 (no norm -- '+str(EvtsPosScdDerBoth)+'/'+nEvts+' evts -- '+title+')',    xBin,xLow,xHigh)
LLXSPosScdDerDistBoth   = TH1F('LLXSPosScdDerBoth',  '-ln(L) when both 2nd derivatives > 0 (XS norm -- '+str(EvtsPosScdDerXSBoth)+'/'+nEvts+' evts -- '+title+')',  xBin,xLow,xHigh)
LLAccPosScdDerDistBoth  = TH1F('LLAccPosScdDerBoth', '-ln(L) when both 2nd derivatives > 0 (Acc norm -- '+str(EvtsPosScdDerAccBoth)+'/'+nEvts+' evts -- '+title+')',xBin,xLow,xHigh)
LLPosScdDerDistInner    = TH1F('LLPosScdDerInner',   '-ln(L) when inner 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerInner))+'/'+nEvts+' evts -- '+title+')',    xBin,xLow,xHigh)
LLXSPosScdDerDistInner  = TH1F('LLXSPosScdDerInner', '-ln(L) when inner 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSInner))+'/'+nEvts+' evts -- '+title+')',  xBin,xLow,xHigh)
LLAccPosScdDerDistInner = TH1F('LLAccPosScdDerInner','-ln(L) when inner 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccInner))+'/'+nEvts+' evts -- '+title+')',xBin,xLow,xHigh)
LLPosScdDerDistOuter    = TH1F('LLPosScdDerOuter',   '-ln(L) when outer 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerOuter))+'/'+nEvts+' evts -- '+title+')',    xBin,xLow,xHigh)
LLXSPosScdDerDistOuter  = TH1F('LLXSPosScdDerOuter', '-ln(L) when outer 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSOuter))+'/'+nEvts+' evts -- '+title+')',  xBin,xLow,xHigh)
LLAccPosScdDerDistOuter = TH1F('LLAccPosScdDerOuter','-ln(L) when outer 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccOuter))+'/'+nEvts+' evts -- '+title+')',xBin,xLow,xHigh)
LLNegScdDerDistInner    = TH1F('LLNegScdDerInner',   '-ln(L) when inner 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerInner)+'/'+nEvts+' evts -- '+title+')',    xBin,xLow,xHigh)
LLXSNegScdDerDistInner  = TH1F('LLXSNegScdDerInner', '-ln(L) when inner 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSInner)+'/'+nEvts+' evts -- '+title+')',  xBin,xLow,xHigh)
LLAccNegScdDerDistInner = TH1F('LLAccNegScdDerInner','-ln(L) when inner 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccInner)+'/'+nEvts+' evts -- '+title+')',xBin,xLow,xHigh)
LLNegScdDerDistOuter    = TH1F('LLNegScdDerOuter',   '-ln(L) when outer 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerOuter)+'/'+nEvts+' evts -- '+title+')',    xBin,xLow,xHigh)
LLXSNegScdDerDistOuter  = TH1F('LLXSNegScdDerOuter', '-ln(L) when outer 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSOuter)+'/'+nEvts+' evts -- '+title+')',  xBin,xLow,xHigh)
LLAccNegScdDerDistOuter = TH1F('LLAccNegScdDerOuter','-ln(L) when outer 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccOuter)+'/'+nEvts+' evts -- '+title+')',xBin,xLow,xHigh)

LLSmallFctDevDist = TH1F('LLSmallFctDeviation','-ln(L) when sum of all deviations between likelihood and function is smaller than 0.5',xBin,xLow,xHigh)

for ii in range(NrConfigs):
  LLPosScdDerDistInner.SetBinContent(LLPosScdDerDistInner.FindBin(Var[ii]), float(LLPosScdDerInner[ii])), 
  LLXSPosScdDerDistInner.SetBinContent(LLXSPosScdDerDistInner.FindBin(Var[ii]),float(LLXSPosScdDerInner[ii]))
  LLAccPosScdDerDistInner.SetBinContent(LLAccPosScdDerDistInner.FindBin(Var[ii]),float(LLAccPosScdDerInner[ii]))
  LLNegScdDerDistInner.SetBinContent(LLNegScdDerDistInner.FindBin(Var[ii]), float(LLNegScdDerInner[ii]))
  LLXSNegScdDerDistInner.SetBinContent(LLXSNegScdDerDistInner.FindBin(Var[ii]),float(LLXSNegScdDerInner[ii]))
  LLAccNegScdDerDistInner.SetBinContent(LLAccNegScdDerDistInner.FindBin(Var[ii]),float(LLAccNegScdDerInner[ii]))

  LLPosScdDerDistOuter.SetBinContent(LLPosScdDerDistOuter.FindBin(Var[ii]), float(LLPosScdDerOuter[ii]))
  LLXSPosScdDerDistOuter.SetBinContent(LLXSPosScdDerDistOuter.FindBin(Var[ii]),float(LLXSPosScdDerOuter[ii]))
  LLAccPosScdDerDistOuter.SetBinContent(LLAccPosScdDerDistOuter.FindBin(Var[ii]),float(LLAccPosScdDerOuter[ii]))
  LLNegScdDerDistOuter.SetBinContent(LLNegScdDerDistOuter.FindBin(Var[ii]), float(LLNegScdDerOuter[ii]))
  LLXSNegScdDerDistOuter.SetBinContent(LLXSNegScdDerDistOuter.FindBin(Var[ii]),float(LLXSNegScdDerOuter[ii]))
  LLAccNegScdDerDistOuter.SetBinContent(LLAccNegScdDerDistOuter.FindBin(Var[ii]),float(LLAccNegScdDerOuter[ii]))

  LLPosScdDerDistBoth.SetBinContent(LLPosScdDerDistBoth.FindBin(Var[ii]), float(LLPosScdDerBoth[ii]))
  LLXSPosScdDerDistBoth.SetBinContent(LLXSPosScdDerDistBoth.FindBin(Var[ii]),float(LLXSPosScdDerBoth[ii]))
  LLAccPosScdDerDistBoth.SetBinContent(LLAccPosScdDerDistBoth.FindBin(Var[ii]),float(LLAccPosScdDerBoth[ii]))

  LLSmallFctDevDist.SetBinContent(LLSmallFctDevDist.FindBin(Var[ii]),float(LLSmallFctDev[ii]))

AppliedCutsDir = Tfile.mkdir("LikelihoodAfterCuts")
SmallFctDevDir = AppliedCutsDir.mkdir("SmallFunctionDeviation")
SmallFctDevDir.cd()
LLSmallFctDevDist.Write()

SignScdDerDir = AppliedCutsDir.mkdir("SignSecondDerivative")
SignScdDerDir.cd()
LLPosScdDerDistInner.Write(), LLXSPosScdDerDistInner.Write(), LLAccPosScdDerDistInner.Write()
LLNegScdDerDistInner.Write(), LLXSNegScdDerDistInner.Write(), LLAccNegScdDerDistInner.Write()
LLPosScdDerDistOuter.Write(), LLXSPosScdDerDistOuter.Write(), LLAccPosScdDerDistOuter.Write()
LLNegScdDerDistOuter.Write(), LLXSNegScdDerDistOuter.Write(),LLAccNegScdDerDistOuter.Write() 
LLPosScdDerDistBoth.Write(),  LLXSPosScdDerDistBoth.Write(), LLAccPosScdDerDistBoth.Write() 

#-- Save the variables separate of scd Der sign  --#
#--   --> Can give hint for good cut!            --#
#- Save the YPlusGausTest variables -#
YPlusGausTestPosScdDer.SetLineColor(3), YPlusGausTestXSPosScdDer.SetLineColor(3), YPlusGausTestAccPosScdDer.SetLineColor(3)
YPlusGausTestNegScdDer.SetLineColor(2), YPlusGausTestXSNegScdDer.SetLineColor(2), YPlusGausTestAccNegScdDer.SetLineColor(2)
YPlusGausTestCanvas = TCanvas('YPlusGausTestCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusGausTestCanvas.cd(), YPlusGausTestPosScdDer.Draw(), YPlusGausTestNegScdDer.Draw("same"), YPlusGausTestCanvas.Write()
YPlusGausTestXSCanvas = TCanvas('YPlusGausTestXSCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusGausTestXSCanvas.cd(), YPlusGausTestXSNegScdDer.Draw(), YPlusGausTestXSPosScdDer.Draw("same"), YPlusGausTestXSCanvas.Write()
YPlusGausTestAccCanvas = TCanvas('YPlusGausTestAccCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusGausTestAccCanvas.cd(), YPlusGausTestAccNegScdDer.Draw(), YPlusGausTestAccPosScdDer.Draw("same"), YPlusGausTestAccCanvas.Write()
#- Save the YPlus variables -#
YPlusPosScdDer.SetLineColor(3), YPlusXSPosScdDer.SetLineColor(3), YPlusAccPosScdDer.SetLineColor(3)
YPlusNegScdDer.SetLineColor(2), YPlusXSNegScdDer.SetLineColor(2), YPlusAccNegScdDer.SetLineColor(2)
YPlusCanvas = TCanvas('YPlusCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusCanvas.cd(), YPlusPosScdDer.Draw(), YPlusNegScdDer.Draw("same"), YPlusCanvas.Write()
YPlusXSCanvas = TCanvas('YPlusXSCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusXSCanvas.cd(), YPlusXSNegScdDer.Draw(), YPlusXSPosScdDer.Draw("same"), YPlusXSCanvas.Write()
YPlusAccCanvas = TCanvas('YPlusAccCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusAccCanvas.cd(), YPlusAccNegScdDer.Draw(), YPlusAccPosScdDer.Draw("same"), YPlusAccCanvas.Write()
#- Save the YPlusPlus variables -#
YPlusPlusPosScdDer.SetLineColor(3), YPlusPlusXSPosScdDer.SetLineColor(3), YPlusPlusAccPosScdDer.SetLineColor(3)
YPlusPlusNegScdDer.SetLineColor(2), YPlusPlusXSNegScdDer.SetLineColor(2), YPlusPlusAccNegScdDer.SetLineColor(2)
YPlusPlusCanvas = TCanvas('YPlusPlusCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YPlusPlusCanvas.cd(), YPlusPlusPosScdDer.Draw(), YPlusPlusNegScdDer.Draw("same"), YPlusPlusCanvas.Write()
YPlusPlusXSCanvas = TCanvas('YPlusPlusXSCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YPlusPlusXSCanvas.cd(), YPlusPlusXSNegScdDer.Draw(), YPlusPlusXSPosScdDer.Draw("same"), YPlusPlusXSCanvas.Write()
YPlusPlusAccCanvas = TCanvas('YPlusPlusAccCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YPlusPlusAccCanvas.cd(), YPlusPlusAccNegScdDer.Draw(), YPlusPlusAccPosScdDer.Draw("same"), YPlusPlusAccCanvas.Write()
#- Save the YMin variables -#
YMinPosScdDer.SetLineColor(3), YMinXSPosScdDer.SetLineColor(3), YMinAccPosScdDer.SetLineColor(3)
YMinNegScdDer.SetLineColor(2), YMinXSNegScdDer.SetLineColor(2), YMinAccNegScdDer.SetLineColor(2)
YMinCanvas = TCanvas('YMinCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YMinCanvas.cd(), YMinPosScdDer.Draw(), YMinNegScdDer.Draw("same"), YMinCanvas.Write()
YMinXSCanvas = TCanvas('YMinXSCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YMinXSCanvas.cd(), YMinXSNegScdDer.Draw(), YMinXSPosScdDer.Draw("same"), YMinXSCanvas.Write()
YMinAccCanvas = TCanvas('YMinAccCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YMinAccCanvas.cd(), YMinAccNegScdDer.Draw(), YMinAccPosScdDer.Draw("same"), YMinAccCanvas.Write()
#- Save the YMinMin variables -#
YMinMinPosScdDer.SetLineColor(3), YMinMinXSPosScdDer.SetLineColor(3), YMinMinAccPosScdDer.SetLineColor(3)
YMinMinNegScdDer.SetLineColor(2), YMinMinXSNegScdDer.SetLineColor(2), YMinMinAccNegScdDer.SetLineColor(2)
YMinMinCanvas = TCanvas('YMinMinCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (no normalisation -- '+title+')')
YMinMinCanvas.cd(), YMinMinPosScdDer.Draw(), YMinMinNegScdDer.Draw("same"), YMinMinCanvas.Write()
YMinMinXSCanvas = TCanvas('YMinMinXSCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (XS normalisation -- '+title+')')
YMinMinXSCanvas.cd(), YMinMinXSNegScdDer.Draw(), YMinMinXSPosScdDer.Draw("same"), YMinMinXSCanvas.Write()
YMinMinAccCanvas = TCanvas('YMinMinAccCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (Acc normalisation -- '+title+')')
YMinMinAccCanvas.cd(), YMinMinAccNegScdDer.Draw(), YMinMinAccPosScdDer.Draw("same"), YMinMinAccCanvas.Write()

