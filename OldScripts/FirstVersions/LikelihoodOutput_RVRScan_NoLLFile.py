#! python
import os
import sys
from math import log
#import ROOT
from ROOT import TH1F,TFile,TCanvas,TLegend,gStyle

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Information about the scanned RVR values and the corresponding cross-section
RVRValues =  ["Re(V_{R}) = -0.5","Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1", "Re(V_{R}) = 0.0","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3", "Re(V_{R}) = 0.5"]
RVR =        [-0.5,              -0.3,              -0.2,              -0.1,               0.0,              0.1,              0.2,              0.3,               0.5              ]
MGXS =       [17.9275,           13.3944,           12.06555,          11.25909,           10.90059,         10.97767,         11.49883,         12.49056,          16.1508          ]
MGXSe =      [0.0123100357311,   0.00995808028337,  0.0093464076837,   0.00836607833038,   0.00822214433527, 0.00847293509122, 0.00901976602967, 0.00874682264197,  0.0113081652137  ]
Acceptance = [0.22164,           0.21742,           0.21672,           0.21737,            0.21670,          0.21677,          0.21437,          0.21793,           0.22205          ]

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
print "Will be using file : ",WeightsFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:  title += "Gen"

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"Likelihood_RVRScan_"+title+"_UsingWeights.root"),'recreate')
LLDist        = TH1F('LLRVR',    '-ln(L) for Single Gaus TF (no normalisation -- '+title+' evts)',11,-0.55,0.55)
LLNormDist    = TH1F('LLRVR_XS', '-ln(L) for Single Gaus TF (XS normalisation -- '+title+' evts)',11,-0.55,0.55)
LLNormAccDist = TH1F('LLRVR_Acc','-ln(L) for Single Gaus TF (Acc normalisation -- '+title+' evts)',11,-0.55,0.55)
gStyle.SetOptStat(0)

#Set the bin content originally very low such that these points don't show up in the histograms
#--> Should still try whether TGraph could be used?
for ii in range(LLDist.GetNbinsX()):
  LLDist.SetBinContent(ii+1,-1000)
  LLNormDist.SetBinContent(ii+1,-1000)
  LLNormAccDist.SetBinContent(ii+1,-1000)  

#Initialize the LL and LLNorm which will be created from the weights to 0
LL = []
LLError = []
LLNorm = []
LLNormError = []
LLNormAcc = []
LLNormAccError = []
NrConsEvents = []
for ii in range(len(RVR)):
  LL.append(0)
  LLNorm.append(0)
  LLError.append(0)
  LLNormError.append(0)
  LLNormAcc.append(0)
  LLNormAccError.append(0)
  NrConsEvents.append(0)

#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#" and str(WeightWord[3]) != "0.0" :
    LL[int(WeightWord[1])-1] = LL[int(WeightWord[1])-1] -log(float(WeightWord[3]))
    LLNorm[int(WeightWord[1])-1] = LLNorm[int(WeightWord[1])-1] - log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1])
    LLNormAcc[int(WeightWord[1])-1] = LLNormAcc[int(WeightWord[1])-1] - log(float(WeightWord[3])) + log(MGXS[int(WeightWord[1])-1]) + log(Acceptance[int(WeightWord[1])-1])
    NrConsEvents[int(WeightWord[1])-1] +=1 # NrConsEvents[int(WeightWord[1])-1] +1

print "Obtained LL output = ",LL
print "Obtained LLNorm output =",LLNorm
print "Obtained LLNormAcc output =",LLNormAcc
print "Nr of events considered :",NrConsEvents

#Fill the histograms and calculate the minimums in order to draw the Delta ChiSq
LLMin = float(LL[0])
LLNormMin = float(LLNorm[0])
LLNormAccMin = float(LLNormAcc[0])
for ii in range(len(RVR)):
  LLDist.SetBinContent(LLDist.FindBin(RVR[ii]), float(LL[ii]))
  LLNormDist.SetBinContent(LLNormDist.FindBin(RVR[ii]),float(LLNorm[ii]))
  LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(RVR[ii]),float(LLNormAcc[ii]))
  if float(LL[ii]) < LLMin:
    LLMin = float(LL[ii])
  if float(LLNorm[ii]) < LLNormMin:
    LLNormMin = float(LLNorm[ii])  
  if float(LLNormAcc[ii]) < LLNormAccMin:
    LLNormAccMin = float(LLNormAcc[ii])  
LLDist.SetMinimum(0)
LLDist.SetMarkerStyle(20)
LLDist.SetLineColor(1)
LLDist.SetMarkerColor(1)
LLNormDist.SetMinimum(0)
LLNormDist.SetMarkerStyle(21)
LLNormDist.SetLineColor(2)
LLNormDist.SetMarkerColor(2)
LLNormAccDist.SetMinimum(0)
LLNormAccDist.SetMarkerStyle(22)
LLNormAccDist.SetLineColor(3)
LLNormAccDist.SetMarkerColor(3)
LLDist.Write()
LLNormDist.Write()
LLNormAccDist.Write()
print " Found minimums : ",LLMin,", ",LLNormMin," and ",LLNormAccMin

#Change bin content for combined canvas!
for ii in range(len(RVR)):
  LLDist.SetBinContent(LLDist.FindBin(RVR[ii]), 2*(float(LLDist.GetBinContent(LLDist.FindBin(RVR[ii]))) - LLMin))
  LLNormDist.SetBinContent(LLNormDist.FindBin(RVR[ii]), 2*(float(LLNormDist.GetBinContent(LLNormDist.FindBin(RVR[ii]))) - LLNormMin))
  LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(RVR[ii]), 2*(float(LLNormAccDist.GetBinContent(LLNormAccDist.FindBin(RVR[ii]))) - LLNormAccMin))

#Prepare everything for the Delta ChiSq canvas!
LLDist.GetYaxis().SetTitle("#Delta #chi^{2} = -2ln(L/L_{max})")
LLDist.GetYaxis().SetTitleOffset(1.3)

LLCanvas = TCanvas('LLDistributions','Comparing influence of XS and acceptance normalisation on -ln(L) -- '+title+' evts')
legend = TLegend(0.75,0.7,0.95,0.95)
LLCanvas.cd()
legend.AddEntry(LLDist,'Non-normalized -ln(likelihood)',"p")
legend.AddEntry(LLNormDist,'Cross-section normalized -ln(likelihood)',"p")
legend.AddEntry(LLNormAccDist,'Cross-section and acceptance normalized -ln(likelihood)',"p")
LLDist.SetTitle(LLCanvas.GetTitle())
LLDist.Draw("P")
LLNormDist.Draw("SAMEP")
LLNormAccDist.Draw("SAMEP")
legend.Draw()
LLCanvas.Write()                 
