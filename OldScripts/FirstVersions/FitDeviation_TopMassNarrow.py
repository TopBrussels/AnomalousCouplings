#! python
import os
import sys
from math import log
#import ROOT
from array import array
from ROOT import TH1F,TH2F,TFile,TCanvas,TLegend,gStyle,gROOT,TGraph
gROOT.SetBatch(True)                       #Don't print the histograms (necessary for fit)

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Information about the scanned MTop values and the corresponding cross-section
MTopValues = ["m_{top} = 153",  "m_{top} = 163",  "m_{top} = 170",  "m_{top} = 171",  "m_{top} = 172",  "m_{top} = 173",  "m_{top} = 174",  "m_{top} = 175","m_{top} = 183",  "m_{top} = 193"]
MTop =       array('d',[153,              163,              170,              171,              172,              173,              174,              175,            183,              193            ])
MGXS =       array('d',[8.20916,          9.6299,           10.57123,         10.70485,         10.8257,          10.96469,         11.08428,         11.22448,       12.18068,         13.3046        ])
MGXSe =      array('d',[0.00641266950107, 0.00775899974932, 0.00857143063963, 0.00814303595657, 0.00878899028501, 0.00816801717126, 0.00904797742371, 0.008653800078, 0.00931290317946, 0.0103310001752])
Acceptance = array('d',[0.16203,          0.19152,          0.21008,          0.21460,          0.21735,          0.21290,          0.21752,          0.22185,        0.23941,          0.26413        ])

#Select which window of masses was considered!
MTopWindow = raw_input('** Choose the correct mass-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [153, 163, 170, 171, 172, 173, 174, 175, 183, 193] \n  2) Normal : [153, 163, 171, 172, 173, 174, 175, 183, 193] \n  3) Narrow : [171, 172, 173, 174, 175] \n --> Choose the correct number : ')

xMinValue = [5, 4, 2]
xPos = [xMinValue[int(MTopWindow)-1]+1,xMinValue[int(MTopWindow)-1]+2]   #Want minimum around mtop = 173, which corresponds to position 2!
xNeg = [xMinValue[int(MTopWindow)-1]-1,xMinValue[int(MTopWindow)-1]-2]
xMin = xMinValue[int(MTopWindow)-1]

if MTopWindow == "2":
  MTopValues.pop(0), MTop.pop(0), MGXS.pop(0), MGXSe.pop(0), Acceptance.pop(0)
elif MTopWindow == "3":
  MTopValues.pop(9), MTop.pop(9), MGXS.pop(9), MGXSe.pop(9), Acceptance.pop(9)
  MTopValues.pop(8), MTop.pop(8), MGXS.pop(8), MGXSe.pop(8), Acceptance.pop(8)
  MTopValues.pop(2), MTop.pop(2), MGXS.pop(2), MGXSe.pop(2), Acceptance.pop(2)
  MTopValues.pop(1), MTop.pop(1), MGXS.pop(1), MGXSe.pop(1), Acceptance.pop(1)
  MTopValues.pop(0), MTop.pop(0), MGXS.pop(0), MGXSe.pop(0), Acceptance.pop(0)
NrConfigs = len(MTop)
print "List of studied top-mass values : ",MTopValues

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:  title += "Gen" 

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

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+'FitDeviation_MTop_'+title+'.root'),'recreate')
gStyle.SetOptStat(0)
YPlusGausTest    = TH1F('YPlusGausTest',   'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (no normalisation)', 250,-0.5,0.5)
YPlusGausTestXS  = TH1F('YPlusGausTestXS', 'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (XS normalisation)', 250,-0.5,0.5)
YPlusGausTestAcc = TH1F('YPlusGausTestAcc','Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (Acc normalisation)',250,-0.5,0.5)
YPlusGausTestPosScdDer    = TH1F('YPlusGausTestPosScdDer',   'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (no normalisation -- second der > 0)', 250,-0.2,0.2)
YPlusGausTestXSPosScdDer  = TH1F('YPlusGausTestXSPosScdDer', 'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (XS normalisation -- second der > 0)', 250,-0.2,0.2)
YPlusGausTestAccPosScdDer = TH1F('YPlusGausTestAccPosScdDer','Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (Acc normalisation -- second der > 0)',250,-0.2,0.2)
YPlusGausTestNegScdDer    = TH1F('YPlusGausTestNegScdDer',   'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (no normalisation -- second der < 0)', 250,-0.2,0.2)
YPlusGausTestXSNegScdDer  = TH1F('YPlusGausTestXSNegScdDer', 'Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (XS normalisation -- second der < 0)', 250,-0.2,0.2)
YPlusGausTestAccNegScdDer = TH1F('YPlusGausTestAccNegScdDer','Comparison of fit deviation in m_{top} = '+str(MTop[xPos[2]])+' and m_{top} = 174 (Acc normalisation -- second der < 0)',250,-0.2,0.2)

YPlus    = TH1F('YPlus',   'Deviation from parabolic fit for m_{top} = 174 (no normalisation)' ,150,-0.5,0.5)
YPlusXS  = TH1F('YPlusXS', 'Deviation from parabolic fit for m_{top} = 174 (XS normalisation)' ,150,-0.5,0.5)
YPlusAcc = TH1F('YPlusAcc','Deviation from parabolic fit for m_{top} = 174 (Acc normalisation)',150,-0.5,0.5)
YPlusPosScdDer    = TH1F('YPlusPosScdDer',   'Deviation from parabolic fit for m_{top} = 174 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusXSPosScdDer  = TH1F('YPlusXSPosScdDer', 'Deviation from parabolic fit for m_{top} = 174 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusAccPosScdDer = TH1F('YPlusAccPosScdDer','Deviation from parabolic fit for m_{top} = 174 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YPlusNegScdDer    = TH1F('YPlusNegScdDer',   'Deviation from parabolic fit for m_{top} = 174 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusXSNegScdDer  = TH1F('YPlusXSNegScdDer', 'Deviation from parabolic fit for m_{top} = 174 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusAccNegScdDer = TH1F('YPlusAccNegScdDer','Deviation from parabolic fit for m_{top} = 174 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YPlusPlus    = TH1F('YPlusPlus',   'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (no normalisation)' ,150,-0.5,0.5)
YPlusPlusXS  = TH1F('YPlusPlusXS', 'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (XS normalisation)' ,150,-0.5,0.5)
YPlusPlusAcc = TH1F('YPlusPlusAcc','Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (Acc normalisation)',150,-0.5,0.5)
YPlusPlusPosScdDer    = TH1F('YPlusPlusPosScdDer',   'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusPlusXSPosScdDer  = TH1F('YPlusPlusXSPosScdDer', 'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YPlusPlusAccPosScdDer = TH1F('YPlusPlusAccPosScdDer','Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YPlusPlusNegScdDer    = TH1F('YPlusPlusNegScdDer',   'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusPlusXSNegScdDer  = TH1F('YPlusPlusXSNegScdDer', 'Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YPlusPlusAccNegScdDer = TH1F('YPlusPlusAccNegScdDer','Deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YMin    = TH1F('YMin',   'Deviation from parabolic fit for m_{top} = 172 (no normalisation)' ,150,-0.5,0.5)
YMinXS  = TH1F('YMinXS', 'Deviation from parabolic fit for m_{top} = 172 (XS normalisation)' ,150,-0.5,0.5)
YMinAcc = TH1F('YMinAcc','Deviation from parabolic fit for m_{top} = 172 (Acc normalisation)',150,-0.5,0.5)
YMinPosScdDer    = TH1F('YMinPosScdDer',   'Deviation from parabolic fit for m_{top} = 172 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinXSPosScdDer  = TH1F('YMinXSPosScdDer', 'Deviation from parabolic fit for m_{top} = 172 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinAccPosScdDer = TH1F('YMinAccPosScdDer','Deviation from parabolic fit for m_{top} = 172 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YMinNegScdDer    = TH1F('YMinNegScdDer',   'Deviation from parabolic fit for m_{top} = 172 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinXSNegScdDer  = TH1F('YMinXSNegScdDer', 'Deviation from parabolic fit for m_{top} = 172 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinAccNegScdDer = TH1F('YMinAccNegScdDer','Deviation from parabolic fit for m_{top} = 172 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YMinMin    = TH1F('YMinMin',   'Deviation from parabolic fit for m_{top} = 171 (no normalisation)' ,150,-0.5,0.5)
YMinMinXS  = TH1F('YMinMinXS', 'Deviation from parabolic fit for m_{top} = 171 (XS normalisation)' ,150,-0.5,0.5)
YMinMinAcc = TH1F('YMinMinAcc','Deviation from parabolic fit for m_{top} = 171 (Acc normalisation)',150,-0.5,0.5)
YMinMinPosScdDer    = TH1F('YMinMinPosScdDer',   'Deviation from parabolic fit for m_{top} = 171 (no normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinMinXSPosScdDer  = TH1F('YMinMinXSPosScdDer', 'Deviation from parabolic fit for m_{top} = 171 (XS normalisation -- second der > 0)' ,150,-0.2,0.2)
YMinMinAccPosScdDer = TH1F('YMinMinAccPosScdDer','Deviation from parabolic fit for m_{top} = 171 (Acc normalisation -- second der > 0)',150,-0.2,0.2)
YMinMinNegScdDer    = TH1F('YMinMinNegScdDer',   'Deviation from parabolic fit for m_{top} = 171 (no normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinMinXSNegScdDer  = TH1F('YMinMinXSNegScdDer', 'Deviation from parabolic fit for m_{top} = 171 (XS normalisation -- second der < 0)' ,150,-0.2,0.2)
YMinMinAccNegScdDer = TH1F('YMinMinAccNegScdDer','Deviation from parabolic fit for m_{top} = 171 (Acc normalisation -- second der < 0)',150,-0.2,0.2)

YRelPlus    = TH1F('YRelPlus',   'Relative deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (no normalisation)' ,250,-0.1,0.1)
YRelPlusXS  = TH1F('YRelPlusXS', 'Relative deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (XS normalisation)' ,250,-0.1,0.1)
YRelPlusAcc = TH1F('YRelPlusAcc','Relative deviation from parabolic fit for m_{top} = '+str(MTop[xPos[2]])+' (Acc normalisation)',250,-0.1,0.1)

YRelMin    = TH1F('YRelMin',   'Relative deviation from parabolic fit for m_{top} = 171 (no normalisation)' ,150,-0.1,0.1)
YRelMinXS  = TH1F('YRelMinXS', 'Relative deviation from parabolic fit for m_{top} = 171 (XS normalisation)' ,150,-0.1,0.1)
YRelMinAcc = TH1F('YRelMinAcc','Relative deviation from parabolic fit for m_{top} = 171 (Acc normalisation)',150,-0.1,0.1)

FstDer    = TH1F('FirstDer',   'First derivative of -ln(likelihood) distribution', 4,0,4)
FstDer.GetXaxis().SetBinLabel(1,"y_{DATA}(x=171) - y_{DATA}(x=172)")
FstDer.GetXaxis().SetBinLabel(1,"y_{DATA}(x=172) - y_{DATA}(x=173)")
FstDer.GetXaxis().SetBinLabel(1,"y_{DATA}(x=173) - y_{DATA}(x=174)")
FstDer.GetXaxis().SetBinLabel(1,"y_{DATA}(x=174) - y_{DATA}(x='+str(MTop[xPos[2]])+')")
FstDerXS  = TH1F('FirstDerivativeXS', 'First derivative of -ln(likelihood) distribution (XS normalisation)', 5,-0.25,0.25)
FstDerAcc = TH1F('FirstDerivativeAcc','First derivative of -ln(likelihood) distribution (Acc normalisation)',5,-0.25,0.25)
ScdDerNarrow    = TH1F('SecondDerivativeNarrow',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using x = 172/173/174)', 250,-20,20)
ScdDerXSNarrow  = TH1F('SecondDerivativeXSNarrow', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using x = 172/173/174)', 250,-20,20)
ScdDerAccNarrow = TH1F('SecondDerivativeAccNarrow','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using x = 172/173/174)',250,-20,20)
ScdDerWide    = TH1F('SecondDerivativeWide',   'Second derivative of -ln(likelihood) distribution (no normalisation -- using x = 171/173/'+str(MTop[xPos[2]])+')', 250,-20,20)
ScdDerXSWide  = TH1F('SecondDerivativeXSWide', 'Second derivative of -ln(likelihood) distribution (XS normalisation -- using x = 171/173/'+str(MTop[xPos[2]])+')', 250,-20,20)
ScdDerAccWide = TH1F('SecondDerivativeAccWide','Second derivative of -ln(likelihood) distribution (Acc normalisation -- using x = 171/173/'+str(MTop[xPos[2]])+')',250,-20,20)
ScdDerScatter    = TH2F('ScdDerScatterPlot',   'Second derivative of -ln(L) using Re(VR) = 172/173/174 versus using Re(VR) = 171/173/'+str(MTop[xPos[2]])+' (no normalisation)', 250,-15,15,250,-15,15)
ScdDerXSScatter  = TH2F('ScdDerXSScatterPlot', 'Second derivative of -ln(L) using Re(VR) = 172/173/174 versus using Re(VR) = 171/173/'+str(MTop[xPos[2]])+' (XS normalisation)', 250,-15,15,250,-15,15)
ScdDerAccScatter = TH2F('ScdDerAccScatterPlot','Second derivative of -ln(L) using Re(VR) = 172/173/174 versus using Re(VR) = 171/173/'+str(MTop[xPos[2]])+' (Acc normalisation)',250,-15,15,250,-15,15)

LnLogDist = TH1F("LnLog","title",5,170.5,175.5)
LnLogDist.SetMarkerStyle(20), LnLogDist.SetLineColor(1), LnLogDist.SetMarkerColor(1), LnLogDist.SetMarkerSize(1.2)
LnLogXSDist = TH1F("LnLogXS","title",5,170.5,175.5)
LnLogXSDist.SetMarkerStyle(21), LnLogXSDist.SetLineColor(3), LnLogXSDist.SetMarkerColor(3), LnLogXSDist.SetMarkerSize(1.2)
LnLogAccDist = TH1F("LnLogAcc","title",5,170.5,175.5)
LnLogAccDist.SetMarkerStyle(22), LnLogAccDist.SetLineColor(4), LnLogAccDist.SetMarkerColor(4), LnLogAccDist.SetMarkerSize(1.2)

FitComp = Tfile.mkdir("FitComparison")
LnLogDir = Tfile.mkdir("LnLogDistr")
LnLogXSDir = Tfile.mkdir("LnLogXSDistr")
LnLogAccDir = Tfile.mkdir("LnLogAccDistr")
FstDerDir = Tfile.mkdir("FirstDerivativeDistr")
LnLogFitCanvas = TCanvas('name','title')
LnLogFitCanvasName = "LnLogFitCanvas_Evt"
LnLogFitCanvasTitle = "Comparing fit function from ROOT fit and algebraic function for event "

#Create the arrays where the likelihood values will be stored
LnLog, LnLogXS, LnLogAcc = [], [], []
for ii in range(NrConfigs):
  LnLog.append(0), LnLogXS.append(0), LnLogAcc.append(0)

#---  Create arrays where the events passing specific cuts will be stored  ---#
EvtsWithPosScdDerNarrow, EvtsWithPosScdDerXSNarrow, EvtsWithPosScdDerAccNarrow = [], [], []
EvtsWithPosScdDerWide, EvtsWithPosScdDerXSWide, EvtsWithPosScdDerAccWide = [], [], []
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
          LnLogDistFit = TH1F("LnLogFit_Evt"+str(int(iEvt)+1),"LnLog distribution for fit (event "+str(int(iEvt)+1)+")",5,170.5,175.5)
          LnLogDistFit.SetMarkerStyle(20), LnLogDistFit.SetLineColor(1), LnLogDistFit.SetMarkerColor(1), LnLogDistFit.SetMarkerSize(1.2)
          LnLogDist.SetName("LnLog_Evt"+str(int(iEvt)+1)),       LnLogDist.SetTitle("LnLog distribution for event "+str(int(iEvt)+1))
          LnLogXSDist.SetName("LnLogXS_Evt"+str(int(iEvt)+1)),   LnLogXSDist.SetTitle("LnLogXS distribution for event "+str(int(iEvt)+1))
          LnLogAccDist.SetName("LnLogAcc_Evt"+str(int(iEvt)+1)), LnLogAccDist.SetTitle("LnLogAcc distribution for event "+str(int(iEvt)+1))
	#print "weightword[1] = ",WeightWord[1]," & weightWord[3] = ",WeightWord[3]
        LnLog[int(WeightWord[1])-1] = -log(float(WeightWord[3]))
        LnLogXS[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
        LnLogAcc[int(WeightWord[1])-1] = -log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
        #---  Fill the LnLog histograms for each event  ---#
        LnLogDistFit.SetBinContent(LnLogDistFit.FindBin(MTop[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        LnLogDist.SetBinContent(LnLogDist.FindBin(MTop[int(WeightWord[1])-1]), LnLog[int(WeightWord[1])-1])
        LnLogXSDist.SetBinContent(LnLogXSDist.FindBin(MTop[int(WeightWord[1])-1]), LnLogXS[int(WeightWord[1])-1])
        LnLogAccDist.SetBinContent(LnLogAccDist.FindBin(MTop[int(WeightWord[1])-1]), LnLogAcc[int(WeightWord[1])-1])
        #---  Only perform the fit after all configurations are considered!  ---#
        if str(WeightWord[1]) == str(NrConfigs):
          #---  Calculate the fit parameters (polynomial)  ---#
          LnLogFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogXSFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          LnLogAccFunction = [array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d'),array('d')]
          for ii in range(2):
            cHat[ii] = LnLog[xNeg[ii]]*MTop[xPos[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]])) - LnLog[xMin]*MTop[xNeg[ii]]*MTop[xPos[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLog[xPos[ii]]*MTop[xNeg[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))
            cHatXS[ii] = LnLogXS[xNeg[ii]]*MTop[xPos[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]])) - LnLogXS[xMin]*MTop[xNeg[ii]]*MTop[xPos[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLogXS[xPos[ii]]*MTop[xNeg[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))
            cHatAcc[ii] = LnLogAcc[xNeg[ii]]*MTop[xPos[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]])) - LnLogAcc[xMin]*MTop[xNeg[ii]]*MTop[xPos[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLogAcc[xPos[ii]]*MTop[xNeg[ii]]*MTop[xMin]/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))

            bHat[ii] = LnLog[xMin]*(MTop[xPos[ii]]+MTop[xNeg[ii]])/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLog[xPos[ii]]*(MTop[xNeg[ii]]+MTop[xMin])/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))-LnLog[xNeg[ii]]*(MTop[xMin]+MTop[xPos[ii]])/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]]))
            bHatXS[ii] = LnLogXS[xMin]*(MTop[xPos[ii]]+MTop[xNeg[ii]])/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLogXS[xPos[ii]]*(MTop[xNeg[ii]]+MTop[xMin])/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))-LnLogXS[xNeg[ii]]*(MTop[xMin]+MTop[xPos[ii]])/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]]))
            bHatAcc[ii] = LnLogAcc[xMin]*(MTop[xPos[ii]]+MTop[xNeg[ii]])/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLogAcc[xPos[ii]]*(MTop[xNeg[ii]]+MTop[xMin])/((MTop[xPos[ii]]-MTop[xMin])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))-LnLogAcc[xNeg[ii]]*(MTop[xMin]+MTop[xPos[ii]])/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xMin]-MTop[xNeg[ii]]))

            aHat[ii] = LnLog[xPos[ii]]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLog[xMin]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLog[xNeg[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))
            aHatXS[ii] = LnLogXS[xPos[ii]]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLogXS[xMin]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLogXS[xNeg[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))
            aHatAcc[ii] = LnLogAcc[xPos[ii]]/((MTop[xPos[ii]]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) - LnLogAcc[xMin]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xMin])) + LnLogAcc[xNeg[ii]]/((MTop[xMin]-MTop[xNeg[ii]])*(MTop[xPos[ii]]-MTop[xNeg[ii]]))

            #-- Now calculate each of the MTop configurations using this parabola --#
            for mtop in range(NrConfigs):
              LnLogFunction[ii].append(aHat[ii]*MTop[mtop]*MTop[mtop]+bHat[ii]*MTop[mtop]+cHat[ii])
              LnLogXSFunction[ii].append( aHatXS[ii]*MTop[mtop]*MTop[mtop]+bHatXS[ii]*MTop[mtop]+cHatXS[ii])
              LnLogAccFunction[ii].append( aHatAcc[ii]*MTop[mtop]*MTop[mtop]+bHatAcc[ii]*MTop[mtop]+cHatAcc[ii])
          LnLogGraphNarrow = TGraph(5, MTop, LnLogFunction[0])
          LnLogGraphNarrow.SetTitle('Comparison between ROOT fit and algebraic method')
          LnLogGraphNarrow.GetYaxis().SetTitle('#Delta ln(Likelihood)')
          LnLogGraphNarrow.GetXaxis().SetTitle('Top mass (GeV)')
          LnLogGraphNarrow.SetMarkerColor(6), LnLogGraphNarrow.SetLineColor(6), LnLogGraphNarrow.SetMarkerSize(1.2)
          LnLogGraphWide = TGraph(5,MTop, LnLogFunction[1])
          LnLogGraphWide.SetMarkerColor(7), LnLogGraphWide.SetLineColor(7), LnLogGraphWide.SetMarkerSize(1.2)	
          #---  Compare with TF1 fit from ROOT  ---#
          LnLogDistFit.Fit("pol2","Q")
          LnLogFitCanvas.SetName(LnLogFitCanvasName+str(int(iEvt)+1)), LnLogFitCanvas.SetTitle(LnLogFitCanvasTitle+str(int(iEvt)+1))
          legendLnLog = TLegend(0.75,0.7,0.95,0.95)
          LnLogFitCanvas.cd()
          legendLnLog.AddEntry(LnLogGraphNarrow,'Fit result using algebraic function (x = 172/173/174)',"p")
          legendLnLog.AddEntry(LnLogGraphWide,'Fit result using algebraic function (x = 171/173/'+str(MTop[xPos[2]])+')',"p")
          legendLnLog.AddEntry(LnLogDistFit,'Obtained LnLog values with ROOT fit',"p")
          LnLogGraphNarrow.Draw("AC*"), LnLogGraphWide.Draw("C*"), LnLogDistFit.Draw("samep"), legendLnLog.Draw()
          FitComp.cd(), LnLogFitCanvas.Write()

          #---  Calculate the first derivative distribution and save them in event-by-event plot  ---#
          FstDer.SetName("FirstDerivative_Evt"+str(int(iEvt)+1))
          FstDer.SetTitle("First derivative of -ln(likelihood) distribution for event"+str(int(iEvt)+1)+" (no normalisation)") 
          for ii in range(FstDer.GetNbinsX()):
            FstDer.SetBinContent(ii+1, LnLog[ii] - LnLog[ii+1])
          FstDerDir.cd(), FstDer.Write()
          #---  Calculate the positive and negative deviation  ---#
          yPlusPlus = [LnLog[xPos[1]] - LnLogFunction[0][xPos[1]], LnLogXS[xPos[1]] - LnLogXSFunction[0][xPos[1]], LnLogAcc[xPos[1]] - LnLogAccFunction[0][xPos[1]]]
          yMinMin   = [LnLog[xNeg[1]] - LnLogFunction[0][xNeg[1]], LnLogXS[xNeg[1]] - LnLogXSFunction[0][xNeg[1]], LnLogAcc[xNeg[1]] - LnLogAccFunction[0][xNeg[1]]]
          yPlus = [LnLog[xPos[0]] - LnLogFunction[1][xPos[0]], LnLogXS[xPos[0]] - LnLogXSFunction[1][xPos[0]], LnLogAcc[xPos[0]] - LnLogAccFunction[1][xPos[0]]]
          yMin  = [LnLog[xNeg[0]] - LnLogFunction[1][xNeg[0]], LnLogXS[xNeg[0]] - LnLogXSFunction[1][xNeg[0]], LnLogAcc[xNeg[0]] - LnLogAccFunction[1][xNeg[0]]]
          scdDerNarrow = [LnLog[xNeg[0]]-2*LnLog[xMin]+LnLog[xPos[0]], LnLogXS[xNeg[0]]-2*LnLogXS[xMin]+LnLogXS[xPos[0]], LnLogAcc[xNeg[0]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[0]]]	
          scdDerWide   = [LnLog[xNeg[1]]-2*LnLog[xMin]+LnLog[xPos[1]], LnLogXS[xNeg[1]]-2*LnLogXS[xMin]+LnLogXS[xPos[1]], LnLogAcc[xNeg[1]]-2*LnLogAcc[xMin]+LnLogAcc[xPos[1]]]	
          #--- Fill the histograms  ---#
          YPlus.Fill(yPlus[0]),                               YPlusXS.Fill(yPlus[1]),                               YPlusAcc.Fill(yPlus[2])
          YPlusPlus.Fill(yPlusPlus[0]),                       YPlusPlusXS.Fill(yPlusPlus[1]),                       YPlusPlusAcc.Fill(yPlusPlus[2])
          YPlusGausTest.Fill(yPlus[0] + yPlusPlus[0]/4),      YPlusGausTestXS.Fill(yPlus[1] + yPlusPlus[1]/4),      YPlusGausTestAcc.Fill(yPlus[2] + yPlusPlus[2]/4)
          #YRelPlus.Fill(yPlus[0]/LnLog[6]),                   YRelPlusXS.Fill(yPlus[1]/LnLogXS[6]),                 YRelPlusAcc.Fill(yPlus[2]/LnLogAcc[6])
          YMin.Fill(yMin[0]),                                 YMinXS.Fill(yMin[1]),                                 YMinAcc.Fill(yMin[2])
          YMinMin.Fill(yMinMin[0]),                           YMinMinXS.Fill(yMinMin[1]),                           YMinMinAcc.Fill(yMinMin[2])
          #YRelMin.Fill(yMin[0]/LnLog[2]),                     YRelMinXS.Fill(yMin[1]/LnLogXS[2]),                   YRelMinAcc.Fill(yMin[2]/LnLogAcc[2])
          ScdDerNarrow.Fill(scdDerNarrow[0]),                 ScdDerXSNarrow.Fill(scdDerNarrow[1]),                 ScdDerAccNarrow.Fill(scdDerNarrow[2])
          ScdDerWide.Fill(scdDerWide[0]),                     ScdDerXSWide.Fill(scdDerWide[1]),                     ScdDerAccWide.Fill(scdDerWide[2])
          ScdDerScatter.Fill(scdDerWide[0], scdDerNarrow[0]), ScdDerXSScatter.Fill(scdDerWide[1], scdDerNarrow[1]), ScdDerAccScatter.Fill(scdDerWide[2], scdDerNarrow[2])
          #---  Save the LnLog distributions (event-per-event) and the Y-deviations in histograms  ---#
          LnLogDir.cd(),    LnLogDist.Write()
          LnLogXSDir.cd(),  LnLogXSDist.Write()
          LnLogAccDir.cd(), LnLogAccDist.Write()
          #-- Apply cut on YPlusGausTest --#
          if (yPlus[0] + yPlusPlus[0]/4) <= 0.025 and (yPlus[0] + yPlusPlus[0]/4) >= -0.025: EvtsWithYPlusGausSmall.append(iEvt+1)
          if (yPlus[1] + yPlusPlus[1]/4) <= 0.025 and (yPlus[1] + yPlusPlus[1]/4) >= -0.025: EvtsWithYPlusGausSmallXS.append(iEvt+1)
          if (yPlus[2] + yPlusPlus[2]/4) <= 0.025 and (yPlus[2] + yPlusPlus[2]/4) >= -0.025: EvtsWithYPlusGausSmallAcc.append(iEvt+1)
          #-- Apply cut on ScdDer (using x = -0.1/0.0/0.1) --#
          if scdDerNarrow[0] > 0.0: EvtsWithPosScdDerNarrow.append(iEvt+1)
          if scdDerNarrow[1] > 0.0: EvtsWithPosScdDerXSNarrow.append(iEvt+1)
          if scdDerNarrow[2] > 0.0: EvtsWithPosScdDerAccNarrow.append(iEvt+1)
          #-- Apply cut on ScdDer (using x = -0.2/0.0/0.2) --#
          if scdDerWide[0] > 0.0:
            EvtsWithPosScdDerWide.append(iEvt+1)
            YPlusGausTestPosScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusPosScdDer.Fill(yPlus[0]), YPlusPlusPosScdDer.Fill(yPlusPlus[0])
            YMinPosScdDer.Fill(yMin[0]),   YMinMinPosScdDer.Fill(yMinMin[0])
          else:
            YPlusGausTestNegScdDer.Fill(yPlus[0] + yPlusPlus[0]/4)
            YPlusNegScdDer.Fill(yPlus[0]), YPlusPlusNegScdDer.Fill(yPlusPlus[0])
            YMinNegScdDer.Fill(yMin[0]),   YMinMinNegScdDer.Fill(yMinMin[0])
          if scdDerWide[1] > 0.0:
            EvtsWithPosScdDerXSWide.append(iEvt+1)
            YPlusGausTestXSPosScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSPosScdDer.Fill(yPlus[1]), YPlusPlusXSPosScdDer.Fill(yPlusPlus[1])
            YMinXSPosScdDer.Fill(yMin[1]),   YMinMinXSPosScdDer.Fill(yMinMin[1])
          else:
            YPlusGausTestXSNegScdDer.Fill(yPlus[1] + yPlusPlus[1]/4)
            YPlusXSNegScdDer.Fill(yPlus[1]), YPlusPlusXSNegScdDer.Fill(yPlusPlus[1])
            YMinXSNegScdDer.Fill(yMin[1]),   YMinMinXSNegScdDer.Fill(yMinMin[1])
          if scdDerWide[2] > 0.0:
            EvtsWithPosScdDerAccWide.append(iEvt+1)
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
ScdDerNarrow.Write(),      ScdDerXSNarrow.Write(),      ScdDerAccNarrow.Write()
ScdDerWide.Write(),      ScdDerXSWide.Write(),      ScdDerAccWide.Write()
ScdDerScatter.Write(), ScdDerXSScatter.Write(), ScdDerAccScatter.Write()

#---  Draw the likelihood distribution separately for events surviving and passing the cuts!  ---#
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = 172/173/174) ",len(EvtsWithPosScdDerNarrow),", ",len(EvtsWithPosScdDerXSNarrow)," & ",len(EvtsWithPosScdDerAccNarrow)
print "Nr of events with 2nd derivative > 0 (LnLog, LnLogXS & LnLogAcc -- using x = 171/173/'+str(MTop[xPos[2]])+') ",len(EvtsWithPosScdDerWide),", ",len(EvtsWithPosScdDerXSWide)," & ",len(EvtsWithPosScdDerAccWide)
LLPosScdDerNarrow,   LLNegScdDerNarrow,   LLXSPosScdDerNarrow,   LLXSNegScdDerNarrow,   LLAccPosScdDerNarrow,   LLAccNegScdDerNarrow = [],[],[],[],[],[]
LLPosScdDerWide,   LLNegScdDerWide,   LLXSPosScdDerWide,   LLXSNegScdDerWide,   LLAccPosScdDerWide,   LLAccNegScdDerWide = [],[],[],[],[],[]
LLPosScdDerBoth, LLNegScdDerBoth, LLXSPosScdDerBoth, LLXSNegScdDerBoth, LLAccPosScdDerBoth, LLAccNegScdDerBoth = [],[],[],[],[],[]
for ii in range(len(MTop)):
  LLPosScdDerNarrow.append(0),   LLNegScdDerNarrow.append(0),   LLXSPosScdDerNarrow.append(0),   LLXSNegScdDerNarrow.append(0),   LLAccPosScdDerNarrow.append(0),   LLAccNegScdDerNarrow.append(0)
  LLPosScdDerWide.append(0),   LLNegScdDerWide.append(0),   LLXSPosScdDerWide.append(0),   LLXSNegScdDerWide.append(0),   LLAccPosScdDerWide.append(0),   LLAccNegScdDerWide.append(0)
  LLPosScdDerBoth.append(0), LLNegScdDerBoth.append(0), LLXSPosScdDerBoth.append(0), LLXSNegScdDerBoth.append(0), LLAccPosScdDerBoth.append(0), LLAccNegScdDerBoth.append(0)

EvtsWithPosScdDerBoth = 0
EvtsWithPosScdDerXSBoth = 0
EvtsWithPosScdDerAccBoth = 0
#Loop over all lines in weights file:
for LikelihoodLine in LikelihoodFile:
  LWord = LikelihoodLine.split()
  #Only interested in files starting with a number
  if str(LWord[0]) != "#" and str(LWord[3]) != "0.0" :
    #---  Separate the events with positive and negative second derivative (using both x = -0.1/0.0/0.1 and x = -0.2/0.0/0.2) ---#
    if int(LWord[0]) in EvtsWithPosScdDerNarrow and int(LWord[0]) in EvtsWithPosScdDerWide:
      LLPosScdDerBoth[int(LWord[1])-1] = LLPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))
      if LWord[1] == "1": EvtsWithPosScdDerBoth += 1
    if int(LWord[0]) in EvtsWithPosScdDerXSNarrow and int(LWord[0]) in EvtsWithPosScdDerXSWide:
      LLXSPosScdDerBoth[int(LWord[1])-1] = LLXSPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
      if LWord[1] == "1": EvtsWithPosScdDerXSBoth += 1
    if int(LWord[0]) in EvtsWithPosScdDerAccNarrow and int(LWord[0]) in EvtsWithPosScdDerAccWide:
      LLAccPosScdDerBoth[int(LWord[1])-1] = LLAccPosScdDerBoth[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
      if LWord[1] == "1": EvtsWithPosScdDerAccBoth += 1
      ##print EvtsWithPosScdDerAccBoth,") Looking at event : ",int(LWord[0])
    #---  Separate the events with positive and negative second derivative (using x = -0.1/0.0/0.1) ---#
    if int(LWord[0]) in EvtsWithPosScdDerNarrow:    LLPosScdDerNarrow[int(LWord[1])-1] = LLPosScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDerNarrow[int(LWord[1])-1] = LLNegScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsWithPosScdDerXSNarrow:  LLXSPosScdDerNarrow[int(LWord[1])-1] = LLXSPosScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDerNarrow[int(LWord[1])-1] = LLXSNegScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsWithPosScdDerAccNarrow: LLAccPosScdDerNarrow[int(LWord[1])-1] = LLAccPosScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDerNarrow[int(LWord[1])-1] = LLAccNegScdDerNarrow[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    #---  Separate the events with positive and negative second derivative (using x = -0.2/0.0/0.2) ---#
    if int(LWord[0]) in EvtsWithPosScdDerWide:    LLPosScdDerWide[int(LWord[1])-1] = LLPosScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))
    else:                                       LLNegScdDerWide[int(LWord[1])-1] = LLNegScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))
    if int(LWord[0]) in EvtsWithPosScdDerXSWide:  LLXSPosScdDerWide[int(LWord[1])-1] = LLXSPosScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    else:                                       LLXSNegScdDerWide[int(LWord[1])-1] = LLXSNegScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])
    if int(LWord[0]) in EvtsWithPosScdDerAccWide: LLAccPosScdDerWide[int(LWord[1])-1] = LLAccPosScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])
    else:                                       LLAccNegScdDerWide[int(LWord[1])-1] = LLAccNegScdDerWide[int(LWord[1])-1]-log(float(LWord[3]))+log(MGXS[int(LWord[1])-1])+log(Acceptance[int(LWord[1])-1])

LLPosScdDerDistBoth   =TH1F('LLPosScdDerBoth',   '-ln(L) when both 2nd derivatives > 0 (no norm -- '+str(EvtsPosScdDerBoth)+'/'+str(nEvts)+' evts -- '+title+')',    7,-0.35,0.35)
LLPosScdDerDistBoth   =TH1F('LLPosScdDerBoth',   '-ln(L) when both 2nd derivatives > 0 (no norm -- '+str(EvtsWithPosScdDerBoth)+'/'+str(nEvts)+' evts -- '+title+')',    5,170.5,175.5)
LLXSPosScdDerDistBoth =TH1F('LLXSPosScdDerBoth', '-ln(L) for events with 2nd derivative > 0 (XS norm -- using x = 174 & '+str(MTop[xPos[2]])+' -- '+str(EvtsWithPosScdDerXSBoth)+'/'+str(nEvts)+' evts)',  5,170.5,175.5)
LLAccPosScdDerDistBoth=TH1F('LLAccPosScdDerBoth','-ln(L) for events with 2nd derivative > 0 (Acc norm -- using x = 174 & '+str(MTop[xPos[2]])+' -- '+str(EvtsWithPosScdDerAccBoth)+'/'+str(nEvts)+' evts)',5,170.5,175.5)
LLPosScdDerDistNarrow   =TH1F('LLPosScdDerNarrow',   '-ln(L) for events with 2nd derivative > 0 (no norm -- x = 172/173/174 -- '+str(len(EvtsWithPosScdDerNarrow))+'/'+str(nEvts)+' evts)',    5,170.5,175.5)
LLXSPosScdDerDistNarrow =TH1F('LLXSPosScdDerNarrow', '-ln(L) for events with 2nd derivative > 0 (XS norm -- x = 172/173/174 -- '+str(len(EvtsWithPosScdDerXSNarrow))+'/'+str(nEvts)+' evts)',  5,170.5,175.5)
LLAccPosScdDerDistNarrow=TH1F('LLAccPosScdDerNarrow','-ln(L) for events with 2nd derivative > 0 (Acc norm -- x = 172/173/174 -- '+str(len(EvtsWithPosScdDerAccNarrow))+'/'+str(nEvts)+' evts)',5,170.5,175.5)
LLPosScdDerDistWide   =TH1F('LLPosScdDerWide',   '-ln(L) for events with 2nd derivative > 0 (no norm -- x = 171/173/'+str(MTop[xPos[2]])+' -- '+str(len(EvtsWithPosScdDerWide))+'/'+str(nEvts)+' evts)',    5,170.5,175.5)
LLXSPosScdDerDistWide =TH1F('LLXSPosScdDerWide', '-ln(L) for events with 2nd derivative > 0 (XS norm -- x = 171/173/'+str(MTop[xPos[2]])+' -- '+str(len(EvtsWithPosScdDerXSWide))+'/'+str(nEvts)+' evts)',  5,170.5,175.5)
LLAccPosScdDerDistWide=TH1F('LLAccPosScdDerWide','-ln(L) for events with 2nd derivative > 0 (Acc norm -- x = 171/173/'+str(MTop[xPos[2]])+' -- '+str(len(EvtsWithPosScdDerAccWide))+'/'+str(nEvts)+' evts)',5,170.5,175.5)
LLNegScdDerDistNarrow   =TH1F('LLNegScdDerNarrow',   '-ln(L) when 2nd derivative < 0 (no norm - narrow fit - '+str(int(nEvts)-len(EvtsWithPosScdDerNarrow))+'/'+str(nEvts)+' evts)',    5,170.5,175.5)
LLXSNegScdDerDistNarrow =TH1F('LLXSNegScdDerNarrow', '-ln(L) when 2nd derivative < 0 (XS norm - narrow fit - '+str(int(nEvts)-len(EvtsWithPosScdDerXSNarrow))+'/'+str(nEvts)+' evts)',  5,170.5,175.5)
LLAccNegScdDerDistNarrow=TH1F('LLAccNegScdDerNarrow','-ln(L) when 2nd derivative < 0 (Acc norm - narrow fit - '+str(int(nEvts)-len(EvtsWithPosScdDerAccNarrow))+'/'+str(nEvts)+' evts)',5,170.5,175.5)
LLNegScdDerDistWide   =TH1F('LLNegScdDerWide',   '-ln(L) when 2nd derivative < 0 (no norm - wide fit - '+str(int(nEvts)-len(EvtsWithPosScdDerNarrow))+'/'+str(nEvts)+' evts)',    5,170.5,175.5)
LLXSNegScdDerDistWide =TH1F('LLXSNegScdDerWide', '-ln(L) when 2nd derivative < 0 (XS norm - wide fit - '+str(int(nEvts)-len(EvtsWithPosScdDerXSNarrow))+'/'+str(nEvts)+' evts)',  5,170.5,175.5)
LLAccNegScdDerDistWide=TH1F('LLAccNegScdDerWide','-ln(L) when 2nd derivative < 0 (Acc norm - wide fit - '+str(int(nEvts)-len(EvtsWithPosScdDerAccNarrow))+'/'+str(nEvts)+' evts)',5,170.5,175.5)

for ii in range(len(MTop)):
  LLPosScdDerDistNarrow.SetBinContent(LLPosScdDerDistNarrow.FindBin(MTop[ii]), float(LLPosScdDerNarrow[ii]))
  LLXSPosScdDerDistNarrow.SetBinContent(LLXSPosScdDerDistNarrow.FindBin(MTop[ii]),float(LLXSPosScdDerNarrow[ii]))
  LLAccPosScdDerDistNarrow.SetBinContent(LLAccPosScdDerDistNarrow.FindBin(MTop[ii]),float(LLAccPosScdDerNarrow[ii]))
  LLNegScdDerDistNarrow.SetBinContent(LLNegScdDerDistNarrow.FindBin(MTop[ii]), float(LLNegScdDerNarrow[ii]))
  LLXSNegScdDerDistNarrow.SetBinContent(LLXSNegScdDerDistNarrow.FindBin(MTop[ii]),float(LLXSNegScdDerNarrow[ii]))
  LLAccNegScdDerDistNarrow.SetBinContent(LLAccNegScdDerDistNarrow.FindBin(MTop[ii]),float(LLAccNegScdDerNarrow[ii]))

  LLPosScdDerDistWide.SetBinContent(LLPosScdDerDistWide.FindBin(MTop[ii]), float(LLPosScdDerWide[ii]))
  LLXSPosScdDerDistWide.SetBinContent(LLXSPosScdDerDistWide.FindBin(MTop[ii]),float(LLXSPosScdDerWide[ii]))
  LLAccPosScdDerDistWide.SetBinContent(LLAccPosScdDerDistWide.FindBin(MTop[ii]),float(LLAccPosScdDerWide[ii]))
  LLNegScdDerDistWide.SetBinContent(LLNegScdDerDistWide.FindBin(MTop[ii]), float(LLNegScdDerWide[ii]))
  LLXSNegScdDerDistWide.SetBinContent(LLXSNegScdDerDistWide.FindBin(MTop[ii]),float(LLXSNegScdDerWide[ii]))
  LLAccNegScdDerDistWide.SetBinContent(LLAccNegScdDerDistWide.FindBin(MTop[ii]),float(LLAccNegScdDerWide[ii]))

  LLPosScdDerDistBoth.SetBinContent(LLPosScdDerDistBoth.FindBin(MTop[ii]), float(LLPosScdDerBoth[ii]))
  LLXSPosScdDerDistBoth.SetBinContent(LLXSPosScdDerDistBoth.FindBin(MTop[ii]),float(LLXSPosScdDerBoth[ii]))
  LLAccPosScdDerDistBoth.SetBinContent(LLAccPosScdDerDistBoth.FindBin(MTop[ii]),float(LLAccPosScdDerBoth[ii]))

AppliedCutsDir = Tfile.mkdir("LikelihoodAfterCuts")
SignScdDerDir = AppliedCutsDir.mkdir("SignSecondDerivative")
SignScdDerDir.cd()
LLPosScdDerDistNarrow.Write(), LLXSPosScdDerDistNarrow.Write(), LLAccPosScdDerDistNarrow.Write()
LLNegScdDerDistNarrow.Write(), LLXSNegScdDerDistNarrow.Write(), LLAccNegScdDerDistNarrow.Write()
LLPosScdDerDistWide.Write(),   LLXSPosScdDerDistWide.Write(),   LLAccPosScdDerDistWide.Write()
LLNegScdDerDistWide.Write(),   LLXSNegScdDerDistWide.Write(),   LLAccNegScdDerDistWide.Write()
LLPosScdDerDistBoth.Write(),   LLXSPosScdDerDistBoth.Write(),   LLAccPosScdDerDistBoth.Write()
#-- Save the variables separate of scd Der sign  --#
#--   --> Can give hint for good cut!            --#
#- Save the YPlusGausTest variables -#
YPlusGausTestPosScdDer.SetLineColor(3), YPlusGausTestXSPosScdDer.SetLineColor(3), YPlusGausTestAccPosScdDer.SetLineColor(3)
YPlusGausTestNegScdDer.SetLineColor(2), YPlusGausTestXSNegScdDer.SetLineColor(2), YPlusGausTestAccNegScdDer.SetLineColor(2)
YPlusGausTestCanvas = TCanvas('YPlusGausTestCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (no normalisation)')
YPlusGausTestCanvas.cd(), YPlusGausTestPosScdDer.Draw(), YPlusGausTestNegScdDer.Draw("same"), YPlusGausTestCanvas.Write()
YPlusGausTestXSCanvas = TCanvas('YPlusGausTestXSCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (XS normalisation)')
YPlusGausTestXSCanvas.cd(), YPlusGausTestXSNegScdDer.Draw(), YPlusGausTestXSPosScdDer.Draw("same"), YPlusGausTestXSCanvas.Write()
YPlusGausTestAccCanvas = TCanvas('YPlusGausTestAccCanvas','Combined distribution of YPlusGausTest variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusGausTestAccCanvas.cd(), YPlusGausTestAccNegScdDer.Draw(), YPlusGausTestAccPosScdDer.Draw("same"), YPlusGausTestAccCanvas.Write()
#- Save the YPlus variables -#
YPlusPosScdDer.SetLineColor(3), YPlusXSPosScdDer.SetLineColor(3), YPlusAccPosScdDer.SetLineColor(3)
YPlusNegScdDer.SetLineColor(2), YPlusXSNegScdDer.SetLineColor(2), YPlusAccNegScdDer.SetLineColor(2)
YPlusCanvas = TCanvas('YPlusCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (no normalisation)')
YPlusCanvas.cd(), YPlusPosScdDer.Draw(), YPlusNegScdDer.Draw("same"), YPlusCanvas.Write()
YPlusXSCanvas = TCanvas('YPlusXSCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (XS normalisation)')
YPlusXSCanvas.cd(), YPlusXSNegScdDer.Draw(), YPlusXSPosScdDer.Draw("same"), YPlusXSCanvas.Write()
YPlusAccCanvas = TCanvas('YPlusAccCanvas','Combined distribution of YPlus variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusAccCanvas.cd(), YPlusAccNegScdDer.Draw(), YPlusAccPosScdDer.Draw("same"), YPlusAccCanvas.Write()
#- Save the YPlusPlus variables -#
YPlusPlusPosScdDer.SetLineColor(3), YPlusPlusXSPosScdDer.SetLineColor(3), YPlusPlusAccPosScdDer.SetLineColor(3)
YPlusPlusNegScdDer.SetLineColor(2), YPlusPlusXSNegScdDer.SetLineColor(2), YPlusPlusAccNegScdDer.SetLineColor(2)
YPlusPlusCanvas = TCanvas('YPlusPlusCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (no normalisation)')
YPlusPlusCanvas.cd(), YPlusPlusPosScdDer.Draw(), YPlusPlusNegScdDer.Draw("same"), YPlusPlusCanvas.Write()
YPlusPlusXSCanvas = TCanvas('YPlusPlusXSCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (XS normalisation)')
YPlusPlusXSCanvas.cd(), YPlusPlusXSNegScdDer.Draw(), YPlusPlusXSPosScdDer.Draw("same"), YPlusPlusXSCanvas.Write()
YPlusPlusAccCanvas = TCanvas('YPlusPlusAccCanvas','Combined distribution of YPlusPlus variable for events with positive and negative scd derivative (Acc normalisation)')
YPlusPlusAccCanvas.cd(), YPlusPlusAccNegScdDer.Draw(), YPlusPlusAccPosScdDer.Draw("same"), YPlusPlusAccCanvas.Write()
#- Save the YMin variables -#
YMinPosScdDer.SetLineColor(3), YMinXSPosScdDer.SetLineColor(3), YMinAccPosScdDer.SetLineColor(3)
YMinNegScdDer.SetLineColor(2), YMinXSNegScdDer.SetLineColor(2), YMinAccNegScdDer.SetLineColor(2)
YMinCanvas = TCanvas('YMinCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (no normalisation)')
YMinCanvas.cd(), YMinPosScdDer.Draw(), YMinNegScdDer.Draw("same"), YMinCanvas.Write()
YMinXSCanvas = TCanvas('YMinXSCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (XS normalisation)')
YMinXSCanvas.cd(), YMinXSNegScdDer.Draw(), YMinXSPosScdDer.Draw("same"), YMinXSCanvas.Write()
YMinAccCanvas = TCanvas('YMinAccCanvas','Combined distribution of YMin variable for events with positive and negative scd derivative (Acc normalisation)')
YMinAccCanvas.cd(), YMinAccNegScdDer.Draw(), YMinAccPosScdDer.Draw("same"), YMinAccCanvas.Write()
#- Save the YMinMin variables -#
YMinMinPosScdDer.SetLineColor(3), YMinMinXSPosScdDer.SetLineColor(3), YMinMinAccPosScdDer.SetLineColor(3)
YMinMinNegScdDer.SetLineColor(2), YMinMinXSNegScdDer.SetLineColor(2), YMinMinAccNegScdDer.SetLineColor(2)
YMinMinCanvas = TCanvas('YMinMinCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (no normalisation)')
YMinMinCanvas.cd(), YMinMinPosScdDer.Draw(), YMinMinNegScdDer.Draw("same"), YMinMinCanvas.Write()
YMinMinXSCanvas = TCanvas('YMinMinXSCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (XS normalisation)')
YMinMinXSCanvas.cd(), YMinMinXSNegScdDer.Draw(), YMinMinXSPosScdDer.Draw("same"), YMinMinXSCanvas.Write()
YMinMinAccCanvas = TCanvas('YMinMinAccCanvas','Combined distribution of YMinMin variable for events with positive and negative scd derivative (Acc normalisation)')
YMinMinAccCanvas.cd(), YMinMinAccNegScdDer.Draw(), YMinMinAccPosScdDer.Draw("same"), YMinMinAccCanvas.Write()
