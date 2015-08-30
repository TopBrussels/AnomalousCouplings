#! python
import os
import sys
from math import log
#import ROOT
from array import array
from ROOT import TH1F,TH2F,TFile,TCanvas,TLegend,gStyle,TDirectory,gROOT,TGraph
gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Information about the scanned RVR values and the corresponding cross-section
RVRValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1", "Re(V_{R}) = 0.0","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
RVR =        array('d',[-0.5,              -0.3,              -0.2,              -0.1,               0.0,              0.1,              0.2,              0.3,               0.5              ])
MGXS =       array('d',[17.9275,           13.3944,           12.06555,          11.25909,           10.90059,         10.97767,         11.49883,         12.49056,          16.1508          ])
MGXSe =      array('d',[0.0123100357311,   0.00995808028337,  0.0093464076837,   0.00836607833038,   0.00822214433527, 0.00847293509122, 0.00901976602967, 0.00874682264197,  0.0113081652137  ])
Acceptance = array('d',[0.22164,           0.21742,           0.21672,           0.21737,            0.21670,          0.21677,          0.21437,          0.21793,           0.22205          ])
xPos = [5,6]   #Want minimum around RVR = 0.0, which corresponds to position 4!
xNeg = [3,2]
xMin = 4
xStep = [RVR[xNeg[0]]-RVR[xNeg[1]], RVR[xMin]-RVR[xNeg[0]], RVR[xPos[0]]-RVR[xMin], RVR[xPos[1]]-RVR[xPos[0]] ]
print "Step size = ",xStep

#File of interest:
list_dir = []
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".out") and file.startswith("weights"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
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
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:  title += "Gen"

if (whichDir.find("Narrow") <= len(whichDir) and whichDir.find("Narrow") > 0) or (whichDir.find("SmallerStep") <= len(whichDir) and whichDir.find("SmallerStep") > 0): title += "_NarrowRVR"
else:									     title += "_WideRVR"

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"FitDeviation_"+title+".root"),'recreate')
gStyle.SetOptStat(0)
YPlusGausTest    = TH1F('YPlusGausTest',   'Comparison of fit deviation in Re(VR) = '+str(RVR[xPos[1]])+' and Re(VR) = '+str(RVR[xPos[0]])+' (no normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestXS  = TH1F('YPlusGausTestXS', 'Comparison of fit deviation in Re(VR) = '+str(RVR[xPos[1]])+' and Re(VR) = '+str(RVR[xPos[0]])+' (XS normalisation -- '+title+' evts)', 250,-0.15,0.15)
YPlusGausTestAcc = TH1F('YPlusGausTestAcc','Comparison of fit deviation in Re(VR) = '+str(RVR[xPos[1]])+' and Re(VR) = '+str(RVR[xPos[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.15,0.15)
YPlusGausTestPosScdDer    = TH1F('YPlusGausTestPosScdDer',   'Gaussian fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (no norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSPosScdDer  = TH1F('YPlusGausTestXSPosScdDer', 'Gaussian fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (XS norm -- outer 2nd der > 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccPosScdDer = TH1F('YPlusGausTestAccPosScdDer','Gaussian fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (Acc norm -- outer 2nd der > 0 -- '+title+' evts)',250,-0.1,0.1)
YPlusGausTestNegScdDer    = TH1F('YPlusGausTestNegScdDer',   'Gaussian fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (no norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestXSNegScdDer  = TH1F('YPlusGausTestXSNegScdDer', 'Gaussing fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (XS norm -- outer 2nd der < 0 -- '+title+' evts)', 250,-0.1,0.1)
YPlusGausTestAccNegScdDer = TH1F('YPlusGausTestAccNegScdDer','Gaussian fit deviation using Re(VR) = '+str(RVR[xPos[1]])+' & '+str(RVR[xPos[0]])+' (Acc norm -- outer 2nd der < 0 -- '+title+' evts)',250,-0.1,0.1)

YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YPlusPosScdDer    = TH1F('YPlusPosScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSPosScdDer  = TH1F('YPlusXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccPosScdDer = TH1F('YPlusAccPosScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YPlusNegScdDer    = TH1F('YPlusNegScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusXSNegScdDer  = TH1F('YPlusXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YPlusAccNegScdDer = TH1F('YPlusAccNegScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YPlusPlus    = TH1F('YPlusPlus',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusXS  = TH1F('YPlusPlusXS', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YPlusPlusAcc = TH1F('YPlusPlusAcc','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YPlusPlusPosScdDer    = TH1F('YPlusPlusPosScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSPosScdDer  = TH1F('YPlusPlusXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccPosScdDer = TH1F('YPlusPlusAccPosScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YPlusPlusNegScdDer    = TH1F('YPlusPlusNegScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusXSNegScdDer  = TH1F('YPlusPlusXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YPlusPlusAccNegScdDer = TH1F('YPlusPlusAccNegScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xPos[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YMin    = TH1F('YMin',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.25,0.25)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.25,0.25)
YMinPosScdDer    = TH1F('YMinPosScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSPosScdDer  = TH1F('YMinXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccPosScdDer = TH1F('YMinAccPosScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.2,0.2)
YMinNegScdDer    = TH1F('YMinNegScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinXSNegScdDer  = TH1F('YMinXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.2,0.2)
YMinAccNegScdDer = TH1F('YMinAccNegScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.2,0.2)

YMinMin    = TH1F('YMinMin',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (no normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinXS  = TH1F('YMinMinXS', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (XS normalisation -- '+title+' evts)' ,150,-0.5,0.5)
YMinMinAcc = TH1F('YMinMinAcc','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (Acc normalisation -- '+title+' evts)',150,-0.5,0.5)
YMinMinPosScdDer    = TH1F('YMinMinPosScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (no normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSPosScdDer  = TH1F('YMinMinXSPosScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (XS normalisation -- outer 2nd der > 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccPosScdDer = TH1F('YMinMinAccPosScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (Acc normalisation -- outer 2nd der > 0 -- '+title+' evts)',150,-0.3,0.3)
YMinMinNegScdDer    = TH1F('YMinMinNegScdDer',   'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (no normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinXSNegScdDer  = TH1F('YMinMinXSNegScdDer', 'Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (XS normalisation -- outer 2nd der < 0 -- '+title+' evts)' ,150,-0.3,0.3)
YMinMinAccNegScdDer = TH1F('YMinMinAccNegScdDer','Deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[1]])+' (Acc normalisation -- outer 2nd der < 0 -- '+title+' evts)',150,-0.3,0.3)

YRelPlus    = TH1F('YRelPlus',   'Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusXS  = TH1F('YRelPlusXS', 'Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,250,-0.1,0.1)
YRelPlusAcc = TH1F('YRelPlusAcc','Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',250,-0.1,0.1)

YRelMin    = TH1F('YRelMin',   'Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (no normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinXS  = TH1F('YRelMinXS', 'Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (XS normalisation -- '+title+' evts)' ,150,-0.1,0.1)
YRelMinAcc = TH1F('YRelMinAcc','Relative deviation from parabolic fit for Re(VR) = '+str(RVR[xNeg[0]])+' (Acc normalisation -- '+title+' evts)',150,-0.1,0.1)

FstDer    = TH1F('FirstDer',   'First derivative of -ln(likelihood) distribution -- '+title+' evts', 4,0,4)
FstDer.GetXaxis().SetBinLabel(1,'(y_{DATA}(x='+str(RVR[xNeg[1]])+') - y_{DATA}(x='+str(RVR[xNeg[0]])+'))/'+str(xStep[0]))
FstDer.GetXaxis().SetBinLabel(2,'(y_{DATA}(x='+str(RVR[xNeg[0]])+') - y_{DATA}(x='+str(RVR[xMin])+'))/'+str(xStep[1]))
FstDer.GetXaxis().SetBinLabel(3,'(y_{DATA}(x='+str(RVR[xMin])+') - y_{DATA}(x='+str(RVR[xPos[0]])+'))/'+str(xStep[2]))
FstDer.GetXaxis().SetBinLabel(4,'(y_{DATA}(x='+str(RVR[xPos[0]])+') - y_{DATA}(x='+str(RVR[xPos[1]])+'))/'+str(xStep[3]))
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

LnLogDist = TH1F("LnLog","title",7,-0.35,0.35)
LnLogDist.SetMarkerStyle(20), LnLogDist.SetLineColor(1), LnLogDist.SetMarkerColor(1), LnLogDist.SetMarkerSize(1.2)
LnLogXSDist = TH1F("LnLogXS","title",7,-0.35,0.35)
LnLogXSDist.SetMarkerStyle(21), LnLogXSDist.SetLineColor(3), LnLogXSDist.SetMarkerColor(3), LnLogXSDist.SetMarkerSize(1.2)
LnLogAccDist = TH1F("LnLogAcc","title",7,-0.35,0.35)
LnLogAccDist.SetMarkerStyle(22), LnLogAccDist.SetLineColor(4), LnLogAccDist.SetMarkerColor(4), LnLogAccDist.SetMarkerSize(1.2)

FitComp = Tfile.mkdir("FitComparison")
LnLogDir = Tfile.mkdir("LnLogDistr")
LnLogXSDir = Tfile.mkdir("LnLogXSDistr")
LnLogAccDir = Tfile.mkdir("LnLogAccDistr")
FstDerDir = Tfile.mkdir("FirstDerivativeDistr")
LnLogFitCanvas = TCanvas('name','title')
LnLogFitCanvasName = 'LnLogFitCanvas_Evt'
LnLogFitCanvasTitle = 'Comparing fit function from ROOT fit and algebraic function for event -- '+title+' evts '

#Create the arrays where the likelihood values will be stored
NrConfigs = 9
LnLog, LnLogXS, LnLogAcc = [], [], []
for ii in range(NrConfigs):
  LnLog.append(0), LnLogXS.append(0), LnLogAcc.append(0)

#---  Create arrays where the events passing specific cuts will be stored  ---#
EvtsWithPosScdDerInner, EvtsWithPosScdDerXSInner, EvtsWithPosScdDerAccInner = [], [], []
EvtsWithPosScdDerOuter, EvtsWithPosScdDerXSOuter, EvtsWithPosScdDerAccOuter = [], [], []
EvtsPosScdDerInner, EvtsPosScdDerXSInner, EvtsPosScdDerAccInner = [], [], []
EvtsPosScdDerOuter, EvtsPosScdDerXSOuter, EvtsPosScdDerAccOuter = [], [], []
EvtsWithYPlusGausSmall, EvtsWithYPlusGausSmallXS, EvtsWithYPlusGausSmallAcc = [], [], []


nEvts = 1000
print " **** Will be using ",nEvts," events!!"
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    for iEvt in range(nEvts):
      if str(WeightWord[0]) == str(iEvt+1):                                                               #Look at one single event!
        if str(WeightWord[1]) == "1":
          cHat, cHatXS, cHatAcc = [0,0],[0,0],[0,0]
          bHat, bHatXS, bHatAcc = [0,0],[0,0],[0,0]
          aHat, aHatXS, aHatAcc = [0,0],[0,0],[0,0]
          LnLogDistFit = TH1F('LnLogFit_Evt'+str(int(iEvt)+1),'LnLog distribution for fit (event '+str(int(iEvt)+1)+' -- '+title+' evts)',11,-0.55,0.55)
          LnLogDistFit.SetMarkerStyle(20), LnLogDistFit.SetLineColor(1), LnLogDistFit.SetMarkerColor(1), LnLogDistFit.SetMarkerSize(1.2)
          LnLogDist.SetName("LnLog_Evt"+str(int(iEvt)+1)),       LnLogDist.SetTitle('LnLog distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLogXSDist.SetName("LnLogXS_Evt"+str(int(iEvt)+1)),   LnLogXSDist.SetTitle('LnLogXS distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
          LnLogAccDist.SetName("LnLogAcc_Evt"+str(int(iEvt)+1)), LnLogAccDist.SetTitle('LnLogAcc distribution for event '+str(int(iEvt)+1)+' -- '+title+' evts')
        LnLog[int(WeightWord[1])-1] = -log(float(WeightWord[3]))
        LnLogXS[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
        LnLogAcc[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
        #---  Fill the LnLog histograms for each event  ---#
        LnLogDistFit.SetBinContent(LnLogDistFit.FindBin(RVR[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        if str(WeightWord[1]) != "1" or str(WeightWord[1]) != "9":                                             #Don't fill the bins coresponding to -0.5 & 0.5 (= too wide range)
          LnLogDist.SetBinContent(LnLogDist.FindBin(RVR[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
          LnLogXSDist.SetBinContent(LnLogXSDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogXS[int(WeightWord[1])-1])
          LnLogAccDist.SetBinContent(LnLogAccDist.FindBin(RVR[int(WeightWord[1])-1]), LnLogAcc[int(WeightWord[1])-1])
        #---  Only perform the fit after all 9 configurations are considered!  ---#
        if str(WeightWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          LnLogFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogXSFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogAccFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          for ii in range(2):
            cHat[ii] = LnLog[xNeg[ii]]*RVR[xPos[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]])) - LnLog[xMin]*RVR[xNeg[ii]]*RVR[xPos[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLog[xPos[ii]]*RVR[xNeg[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))
            cHatXS[ii] = LnLogXS[xNeg[ii]]*RVR[xPos[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]])) - LnLogXS[xMin]*RVR[xNeg[ii]]*RVR[xPos[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLogXS[xPos[ii]]*RVR[xNeg[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))
            cHatAcc[ii] = LnLogAcc[xNeg[ii]]*RVR[xPos[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]])) - LnLogAcc[xMin]*RVR[xNeg[ii]]*RVR[xPos[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLogAcc[xPos[ii]]*RVR[xNeg[ii]]*RVR[xMin]/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))

            bHat[ii] = LnLog[xMin]*(RVR[xPos[ii]]+RVR[xNeg[ii]])/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLog[xPos[ii]]*(RVR[xNeg[ii]]+RVR[xMin])/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))-LnLog[xNeg[ii]]*(RVR[xMin]+RVR[xPos[ii]])/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]]))
            bHatXS[ii] = LnLogXS[xMin]*(RVR[xPos[ii]]+RVR[xNeg[ii]])/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLogXS[xPos[ii]]*(RVR[xNeg[ii]]+RVR[xMin])/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))-LnLogXS[xNeg[ii]]*(RVR[xMin]+RVR[xPos[ii]])/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]]))
            bHatAcc[ii] = LnLogAcc[xMin]*(RVR[xPos[ii]]+RVR[xNeg[ii]])/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLogAcc[xPos[ii]]*(RVR[xNeg[ii]]+RVR[xMin])/((RVR[xPos[ii]]-RVR[xMin])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))-LnLogAcc[xNeg[ii]]*(RVR[xMin]+RVR[xPos[ii]])/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xMin]-RVR[xNeg[ii]]))

            aHat[ii] = LnLog[xPos[ii]]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLog[xMin]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLog[xNeg[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))
            aHatXS[ii] = LnLogXS[xPos[ii]]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLogXS[xMin]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLogXS[xNeg[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))
            aHatAcc[ii] = LnLogAcc[xPos[ii]]/((RVR[xPos[ii]]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) - LnLogAcc[xMin]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xMin])) + LnLogAcc[xNeg[ii]]/((RVR[xMin]-RVR[xNeg[ii]])*(RVR[xPos[ii]]-RVR[xNeg[ii]]))

            #-- Now calculate each of the 9 RVR configurations using this parabola --#
            for rvr in range(9):
              LnLogFunction[ii].append(aHat[ii]*RVR[rvr]*RVR[rvr]+bHat[ii]*RVR[rvr]+cHat[ii])
              LnLogXSFunction[ii].append( aHatXS[ii]*RVR[rvr]*RVR[rvr]+bHatXS[ii]*RVR[rvr]+cHatXS[ii])
              LnLogAccFunction[ii].append( aHatAcc[ii]*RVR[rvr]*RVR[rvr]+bHatAcc[ii]*RVR[rvr]+cHatAcc[ii])
          LnLogGraphInner = TGraph(9, RVR, LnLogFunction[0])
          LnLogGraphInner.SetTitle('Comparison between ROOT fit and algebraic method -- '+title+' evts')
          LnLogGraphInner.GetYaxis().SetTitle('#Delta ln(Likelihood)')
          LnLogGraphInner.GetXaxis().SetTitle('Re(VR) value')
          LnLogGraphInner.SetMarkerColor(6), LnLogGraphInner.SetLineColor(6), LnLogGraphInner.SetMarkerSize(1.2)
          LnLogGraphOuter = TGraph(9,RVR, LnLogFunction[1])
          LnLogGraphOuter.SetMarkerColor(7), LnLogGraphOuter.SetLineColor(7), LnLogGraphOuter.SetMarkerSize(1.2)	
          #---  Compare with TF1 fit from ROOT  ---#
          LnLogDistFit.Fit("pol2","Q","",-0.35, 0.35)
          LnLogFitCanvas.SetName(LnLogFitCanvasName+str(int(iEvt)+1)), LnLogFitCanvas.SetTitle(LnLogFitCanvasTitle+str(int(iEvt)+1))
          legendLnLog = TLegend(0.75,0.7,0.95,0.95)
          LnLogFitCanvas.cd()
          legendLnLog.AddEntry(LnLogGraphInner,'Fit result using algebraic function (x = '+str(RVR[xNeg[0]])+'/'+str(RVR[xMin])+'/'+str(RVR[xPos[0]])+')',"p")
          legendLnLog.AddEntry(LnLogGraphOuter,'Fit result using algebraic function (x = '+str(RVR[xNeg[1]])+'/'+str(RVR[xMin])+'/'+str(RVR[xPos[1]])+')',"p")
          legendLnLog.AddEntry(LnLogDistFit,'Obtained LnLog values with ROOT fit',"p")
          LnLogGraphInner.Draw("AC*"), LnLogGraphOuter.Draw("C*"), LnLogDistFit.Draw("samep"), legendLnLog.Draw()
          FitComp.cd(), LnLogFitCanvas.Write()

          #---  Calculate the first derivative distribution and save them in event-by-event plot  ---#
          FstDer.SetName('FirstDerivative_Evt'+str(int(iEvt)+1))
          FstDer.SetTitle('First derivative of -ln(likelihood) distribution for event'+str(int(iEvt)+1)+' (no normalisation) -- '+title+' evts') 
          for ii in range(FstDer.GetNbinsX()):
            FstDer.SetBinContent(ii+1, (LnLog[ii+2] - LnLog[ii+3])/0.1)
          FstDerDir.cd(), FstDer.Write()
          #---  Calculate the positive and negative deviation  ---#
          yPlusPlus = [LnLog[xPos[1]] - LnLogFunction[0][xPos[1]], LnLogXS[xPos[1]] - LnLogXSFunction[0][xPos[1]], LnLogAcc[xPos[1]] - LnLogAccFunction[0][xPos[1]]]
          yMinMin   = [LnLog[xNeg[1]] - LnLogFunction[0][xNeg[1]], LnLogXS[xNeg[1]] - LnLogXSFunction[0][xNeg[1]], LnLogAcc[xNeg[1]] - LnLogAccFunction[0][xNeg[1]]]
          yPlus = [LnLog[xPos[0]] - LnLogFunction[1][xPos[0]], LnLogXS[xPos[0]] - LnLogXSFunction[1][xPos[0]], LnLogAcc[xPos[0]] - LnLogAccFunction[1][xPos[0]]]
          yMin  = [LnLog[xNeg[0]] - LnLogFunction[1][xNeg[0]], LnLogXS[xNeg[0]] - LnLogXSFunction[1][xNeg[0]], LnLogAcc[xNeg[0]] - LnLogAccFunction[1][xNeg[0]]]
          scdDerInner = [LnLog[xNeg[0]]-2*LnLog[xMin]+LnLog[xPos[0]], LnLogXS[xNeg[0]]-2*LnLogXS[xMin]+LnLogXS[xPos[0]], LnLogAcc[xNeg[0]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[0]]]	
          scdDerOuter   = [LnLog[xNeg[1]]-2*LnLog[xMin]+LnLog[xPos[1]], LnLogXS[xNeg[1]]-2*LnLogXS[xMin]+LnLogXS[xPos[1]], LnLogAcc[xNeg[1]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[1]]]	
          #--- Fill the histograms  ---#
          YPlus.Fill(yPlus[0]),                               YPlusXS.Fill(yPlus[1]),                               YPlusAcc.Fill(yPlus[2])
          YPlusPlus.Fill(yPlusPlus[0]),                       YPlusPlusXS.Fill(yPlusPlus[1]),                       YPlusPlusAcc.Fill(yPlusPlus[2])
          YPlusGausTest.Fill(yPlus[0] + yPlusPlus[0]/4),      YPlusGausTestXS.Fill(yPlus[1] + yPlusPlus[1]/4),      YPlusGausTestAcc.Fill(yPlus[2] + yPlusPlus[2]/4)
          #YRelPlus.Fill(yPlus[0]/LnLog[6]),                   YRelPlusXS.Fill(yPlus[1]/LnLogXS[6]),                 YRelPlusAcc.Fill(yPlus[2]/LnLogAcc[6])
          YMin.Fill(yMin[0]),                                 YMinXS.Fill(yMin[1]),                                 YMinAcc.Fill(yMin[2])
          YMinMin.Fill(yMinMin[0]),                           YMinMinXS.Fill(yMinMin[1]),                           YMinMinAcc.Fill(yMinMin[2])
          #YRelMin.Fill(yMin[0]/LnLog[2]),                     YRelMinXS.Fill(yMin[1]/LnLogXS[2]),                   YRelMinAcc.Fill(yMin[2]/LnLogAcc[2])
          ScdDerInner.Fill(scdDerInner[0]),                 ScdDerXSInner.Fill(scdDerInner[1]),                 ScdDerAccInner.Fill(scdDerInner[2])
          ScdDerOuter.Fill(scdDerOuter[0]),                     ScdDerXSOuter.Fill(scdDerOuter[1]),                     ScdDerAccOuter.Fill(scdDerOuter[2])
          ScdDerScatter.Fill(scdDerOuter[0], scdDerInner[0]), ScdDerXSScatter.Fill(scdDerOuter[1], scdDerInner[1]), ScdDerAccScatter.Fill(scdDerOuter[2], scdDerInner[2])
          #---  Save the LnLog distributions (event-per-event) and the Y-deviations in histograms  ---#
          LnLogDir.cd(),    LnLogDist.Write()
          LnLogXSDir.cd(),  LnLogXSDist.Write()
          LnLogAccDir.cd(), LnLogAccDist.Write()
          #-- Apply cut on YPlusGausTest --#
          if (yPlus[0] + yPlusPlus[0]/4) <= 0.025 and (yPlus[0] + yPlusPlus[0]/4) >= -0.025: EvtsWithYPlusGausSmall.append(iEvt+1)
          if (yPlus[1] + yPlusPlus[1]/4) <= 0.025 and (yPlus[1] + yPlusPlus[1]/4) >= -0.025: EvtsWithYPlusGausSmallXS.append(iEvt+1)
          if (yPlus[2] + yPlusPlus[2]/4) <= 0.025 and (yPlus[2] + yPlusPlus[2]/4) >= -0.025: EvtsWithYPlusGausSmallAcc.append(iEvt+1)
          #-- Apply cut on ScdDer (using inner RVR points) --#
          if scdDerInner[0] > 0.0: EvtsPosScdDerInner.append(iEvt+1)
          if scdDerInner[1] > 0.0: EvtsPosScdDerXSInner.append(iEvt+1)
          if scdDerInner[2] > 0.0: EvtsPosScdDerAccInner.append(iEvt+1)
          #-- Apply cut on ScdDer (using outer RVR points) --#
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

#--- Save all the histograms containing information about all the events! ---#
Tfile.cd()
YPlusGausTest.Write(), YPlusGausTestXS.Write(), YPlusGausTestAcc.Write()
YPlus.Write(),         YPlusXS.Write(),         YPlusAcc.Write()
YPlusPlus.Write(),     YPlusPlusXS.Write(),     YPlusPlusAcc.Write()
YRelPlus.Write(),      YRelPlusXS.Write(),      YRelPlusAcc.Write()
YMin.Write(),          YMinXS.Write(),          YMinAcc.Write()
YMinMin.Write(),       YMinMinXS.Write(),       YMinMinAcc.Write()
YRelMin.Write(),       YRelMinXS.Write(),       YRelMinAcc.Write()
ScdDerInner.Write(),   ScdDerXSInner.Write(),   ScdDerAccInner.Write()
ScdDerOuter.Write(),   ScdDerXSOuter.Write(),   ScdDerAccOuter.Write()
ScdDerScatter.Write(), ScdDerXSScatter.Write(), ScdDerAccScatter.Write()

#---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---#
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = ",str(RVR[xNeg[0]]),"/",str(RVR[xMin]),"/",str(RVR[xPos[0]]),") :",len(EvtsPosScdDerInner),", ",len(EvtsPosScdDerXSInner)," & ",len(EvtsPosScdDerAccInner)
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = ",str(RVR[xNeg[1]]),"/",str(RVR[xMin]),"/",str(RVR[xPos[1]]),") :",len(EvtsPosScdDerOuter),", ",len(EvtsPosScdDerXSOuter)," & ",len(EvtsPosScdDerAccOuter)
print "Nr of events with Gaussiaanse vergelijking voor + (LnLog, LnLogXS & LnLogAcc) ", len(EvtsWithYPlusGausSmall),", ",len(EvtsWithYPlusGausSmallXS)," & ",len(EvtsWithYPlusGausSmallAcc)
LLPosScdDerInner,   LLNegScdDerInner,   LLXSPosScdDerInner,   LLXSNegScdDerInner,   LLAccPosScdDerInner,   LLAccNegScdDerInner = [],[],[],[],[],[]
LLPosScdDerOuter,   LLNegScdDerOuter,   LLXSPosScdDerOuter,   LLXSNegScdDerOuter,   LLAccPosScdDerOuter,   LLAccNegScdDerOuter = [],[],[],[],[],[]
LLPosScdDerBoth, LLNegScdDerBoth, LLXSPosScdDerBoth, LLXSNegScdDerBoth, LLAccPosScdDerBoth, LLAccNegScdDerBoth = [],[],[],[],[],[]
for ii in range(len(RVR)):
  LLPosScdDerInner.append(0),   LLNegScdDerInner.append(0),   LLXSPosScdDerInner.append(0),   LLXSNegScdDerInner.append(0),   LLAccPosScdDerInner.append(0),   LLAccNegScdDerInner.append(0)
  LLPosScdDerOuter.append(0),   LLNegScdDerOuter.append(0),   LLXSPosScdDerOuter.append(0),   LLXSNegScdDerOuter.append(0),   LLAccPosScdDerOuter.append(0),   LLAccNegScdDerOuter.append(0)
  LLPosScdDerBoth.append(0), LLNegScdDerBoth.append(0), LLXSPosScdDerBoth.append(0), LLXSNegScdDerBoth.append(0), LLAccPosScdDerBoth.append(0), LLAccNegScdDerBoth.append(0)

EvtsPosScdDerBoth = 0
EvtsPosScdDerXSBoth = 0
EvtsPosScdDerAccBoth = 0
#Loop over all lines in weights file:
for LikelihoodLine in LikelihoodFile:
  LWord = LikelihoodLine.split()
  #Only interested in files starting with a number
  if str(LWord[0]) != "#" and str(LWord[3]) != "0.0" :
    #---  Separate the events with positive and negative second derivative (using both inner and outer RVR points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner and int(LWord[0]) in EvtsPosScdDerOuter:
      LLPosScdDerBoth[int(LWord[1])-1] = LLPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))
      if LWord[1] == "1": EvtsPosScdDerBoth += 1
    if int(LWord[0]) in EvtsPosScdDerXSInner and int(LWord[0]) in EvtsPosScdDerXSOuter:
      LLXSPosScdDerBoth[int(LWord[1])-1] = LLXSPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerXSBoth += 1
    if int(LWord[0]) in EvtsPosScdDerAccInner and int(LWord[0]) in EvtsPosScdDerAccOuter:
      LLAccPosScdDerBoth[int(LWord[1])-1] = LLAccPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
      if LWord[1] == "1": EvtsPosScdDerAccBoth += 1
      ##print EvtsPosScdDerAccBoth,") Looking at event : ",int(LWord[0])
    #---  Separate the events with positive and negative second derivative (using inner RVR points) ---#
    if int(LWord[0]) in EvtsPosScdDerInner:    LLPosScdDerInner[int(LWord[1])-1] = LLPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDerInner[int(LWord[1])-1] = LLNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSInner:  LLXSPosScdDerInner[int(LWord[1])-1] = LLXSPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDerInner[int(LWord[1])-1] = LLXSNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccInner: LLAccPosScdDerInner[int(LWord[1])-1] = LLAccPosScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDerInner[int(LWord[1])-1] = LLAccNegScdDerInner[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    #---  Separate the events with positive and negative second derivative (using outer RVR points) ---#
    if int(LWord[0]) in EvtsPosScdDerOuter:    LLPosScdDerOuter[int(LWord[1])-1] = LLPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDerOuter[int(LWord[1])-1] = LLNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsPosScdDerXSOuter:  LLXSPosScdDerOuter[int(LWord[1])-1] = LLXSPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDerOuter[int(LWord[1])-1] = LLXSNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsPosScdDerAccOuter: LLAccPosScdDerOuter[int(LWord[1])-1] = LLAccPosScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDerOuter[int(LWord[1])-1] = LLAccNegScdDerOuter[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])

NrEvtsNegScdDerInner, NrEvtsNegScdDerXSInner, NrEvtsNegScdDerAccInner = int(nEvts)-len(EvtsPosScdDerInner), nEvts - len(EvtsPosScdDerXSInner), nEvts - len(EvtsPosScdDerAccInner)
NrEvtsNegScdDerOuter, NrEvtsNegScdDerXSOuter, NrEvtsNegScdDerAccOuter = int(nEvts)-len(EvtsPosScdDerOuter), nEvts - len(EvtsPosScdDerXSOuter), nEvts - len(EvtsPosScdDerAccOuter)
LLPosScdDerDistBoth   =TH1F('LLPosScdDerBoth',   '-ln(L) when both 2nd derivatives > 0 (no norm -- '+str(EvtsPosScdDerBoth)+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLXSPosScdDerDistBoth =TH1F('LLXSPosScdDerBoth', '-ln(L) when both 2nd derivatives > 0 (XS norm -- '+str(EvtsPosScdDerXSBoth)+'/'+str(nEvts)+' evts -- '+title+')',  7,-0.35,0.35)
LLAccPosScdDerDistBoth=TH1F('LLAccPosScdDerBoth','-ln(L) when both 2nd derivatives > 0 (Acc norm -- '+str(EvtsPosScdDerAccBoth)+'/'+str(nEvts)+' evts -- '+title+')',7,-0.35,0.35)
LLPosScdDerDistInner   =TH1F('LLPosScdDerInner',   '-ln(L) when inner 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerInner))+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLXSPosScdDerDistInner =TH1F('LLXSPosScdDerInner', '-ln(L) when inner 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSInner))+'/'+str(nEvts)+' evts -- '+title+')',  7,-0.35,0.35)
LLAccPosScdDerDistInner=TH1F('LLAccPosScdDerInner','-ln(L) when inner 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccInner))+'/'+str(nEvts)+' evts -- '+title+')',7,-0.35,0.35)
LLPosScdDerDistOuter   =TH1F('LLPosScdDerOuter',   '-ln(L) when outer 2nd derivative > 0 (no norm -- '+str(len(EvtsPosScdDerOuter))+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLXSPosScdDerDistOuter =TH1F('LLXSPosScdDerOuter', '-ln(L) when outer 2nd derivative > 0 (XS norm -- '+str(len(EvtsPosScdDerXSOuter))+'/'+str(nEvts)+' evts -- '+title+')',  7,-0.35,0.35)
LLAccPosScdDerDistOuter=TH1F('LLAccPosScdDerOuter','-ln(L) when outer 2nd derivative > 0 (Acc norm -- '+str(len(EvtsPosScdDerAccOuter))+'/'+str(nEvts)+' evts -- '+title+')',7,-0.35,0.35)
LLNegScdDerDistInner   =TH1F('LLNegScdDerInner',   '-ln(L) when inner 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerInner)+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLXSNegScdDerDistInner =TH1F('LLXSNegScdDerInner', '-ln(L) when inner 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSInner)+'/'+str(nEvts)+' evts -- '+title+')',  7,-0.35,0.35)
LLAccNegScdDerDistInner=TH1F('LLAccNegScdDerInner','-ln(L) when inner 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccInner)+'/'+str(nEvts)+' evts -- '+title+')',7,-0.35,0.35)
LLNegScdDerDistOuter   =TH1F('LLNegScdDerOuter',   '-ln(L) when outer 2nd derivative < 0 (no norm -- '+str(NrEvtsNegScdDerOuter)+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLXSNegScdDerDistOuter =TH1F('LLXSNegScdDerOuter', '-ln(L) when outer 2nd derivative < 0 (XS norm -- '+str(NrEvtsNegScdDerXSOuter)+'/'+str(nEvts)+' evts -- '+title+')',  7,-0.35,0.35)
LLAccNegScdDerDistOuter=TH1F('LLAccNegScdDerOuter','-ln(L) when outer 2nd derivative < 0 (Acc norm -- '+str(NrEvtsNegScdDerAccOuter)+'/'+str(nEvts)+' evts -- '+title+')',7,-0.35,0.35)

for ii in range(len(RVR)):
  LLPosScdDerDistInner.SetBinContent(LLPosScdDerDistInner.FindBin(RVR[ii]), float(LLPosScdDerInner[ii]))
  LLXSPosScdDerDistInner.SetBinContent(LLXSPosScdDerDistInner.FindBin(RVR[ii]),float(LLXSPosScdDerInner[ii]))
  LLAccPosScdDerDistInner.SetBinContent(LLAccPosScdDerDistInner.FindBin(RVR[ii]),float(LLAccPosScdDerInner[ii]))
  LLNegScdDerDistInner.SetBinContent(LLNegScdDerDistInner.FindBin(RVR[ii]), float(LLNegScdDerInner[ii]))
  LLXSNegScdDerDistInner.SetBinContent(LLXSNegScdDerDistInner.FindBin(RVR[ii]),float(LLXSNegScdDerInner[ii]))
  LLAccNegScdDerDistInner.SetBinContent(LLAccNegScdDerDistInner.FindBin(RVR[ii]),float(LLAccNegScdDerInner[ii]))

  LLPosScdDerDistOuter.SetBinContent(LLPosScdDerDistOuter.FindBin(RVR[ii]), float(LLPosScdDerOuter[ii]))
  LLXSPosScdDerDistOuter.SetBinContent(LLXSPosScdDerDistOuter.FindBin(RVR[ii]),float(LLXSPosScdDerOuter[ii]))
  LLAccPosScdDerDistOuter.SetBinContent(LLAccPosScdDerDistOuter.FindBin(RVR[ii]),float(LLAccPosScdDerOuter[ii]))
  LLNegScdDerDistOuter.SetBinContent(LLNegScdDerDistOuter.FindBin(RVR[ii]), float(LLNegScdDerOuter[ii]))
  LLXSNegScdDerDistOuter.SetBinContent(LLXSNegScdDerDistOuter.FindBin(RVR[ii]),float(LLXSNegScdDerOuter[ii]))
  LLAccNegScdDerDistOuter.SetBinContent(LLAccNegScdDerDistOuter.FindBin(RVR[ii]),float(LLAccNegScdDerOuter[ii]))

  LLPosScdDerDistBoth.SetBinContent(LLPosScdDerDistBoth.FindBin(RVR[ii]), float(LLPosScdDerBoth[ii]))
  LLXSPosScdDerDistBoth.SetBinContent(LLXSPosScdDerDistBoth.FindBin(RVR[ii]),float(LLXSPosScdDerBoth[ii]))
  LLAccPosScdDerDistBoth.SetBinContent(LLAccPosScdDerDistBoth.FindBin(RVR[ii]),float(LLAccPosScdDerBoth[ii]))

AppliedCutsDir = Tfile.mkdir("LikelihoodAfterCuts")
SignScdDerDir = AppliedCutsDir.mkdir("SignSecondDerivative")
SignScdDerDir.cd()
LLPosScdDerDistInner.Write(),   LLXSPosScdDerDistInner.Write(),   LLAccPosScdDerDistInner.Write()
LLNegScdDerDistInner.Write(),   LLXSNegScdDerDistInner.Write(),   LLAccNegScdDerDistInner.Write()
LLPosScdDerDistOuter.Write(),   LLXSPosScdDerDistOuter.Write(),   LLAccPosScdDerDistOuter.Write()
LLNegScdDerDistOuter.Write(),   LLXSNegScdDerDistOuter.Write(),   LLAccNegScdDerDistOuter.Write()
LLPosScdDerDistBoth.Write(), LLXSPosScdDerDistBoth.Write(), LLAccPosScdDerDistBoth.Write()
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

