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

#Information about the scanned top-mass values and the corresponding cross-section
MTopValues = ["m_{top} = 153",  "m_{top} = 163",  "m_{top} = 170",  "m_{top} = 171",  "m_{top} = 172",  "m_{top} = 173",  "m_{top} = 174",  "m_{top} = 175","m_{top} = 183",  "m_{top} = 193"]
MTop =       [153,              163,              170,              171,              172,              173,              174,              175,            183,              193            ]
MGXS =       [8.20916,          9.6299,           10.57123,         10.70485,         10.8257,          10.96469,         11.08428,         11.22448,       12.18068,         13.3046        ]
MGXSe =      [0.00641266950107, 0.00775899974932, 0.00857143063963, 0.00814303595657, 0.00878899028501, 0.00816801717126, 0.00904797742371, 0.008653800078, 0.00931290317946, 0.0103310001752]
Acceptance = [0.16203,          0.19152,          0.21008,          0.21460,          0.21735,          0.21290,          0.21752,          0.22185,        0.23941,          0.26413        ]

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

GenLevel = False
if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:  
  title += "Gen" 
  GenLevel = True

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"Likelihood_TopMass_"+title+"_UsingWeights.root"),'recreate')
LLDist        = TH1F('LL',    '-ln(likelihood) for Single Gaus TF (no normalisation --'+title+' evts)', 41,152.5,193.5)
LLNormDist    = TH1F('LL_XS', '-ln(likelihood) for Single Gaus TF (XS normalisation --'+title+' evts)', 41,152.5,193.5)
LLNormAccDist = TH1F('LL_Acc','-ln(likelihood) for Single Gaus TF (Acc normalisation --'+title+' evts)',41,152.5,193.5)
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
for ii in range(len(MTopValues)):
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
if GenLevel == False: print "Obtained LLNormAcc output =",LLNormAcc
print "Nr of events considered :",NrConsEvents

#Fill the histograms and calculate the minimums in order to draw the Delta ChiSq
LLMin = float(LL[0])
LLNormMin = float(LLNorm[0])
LLNormAccMin = float(LLNormAcc[0])
for ii in range(len(MTop)):
  LLDist.SetBinContent(LLDist.FindBin(MTop[ii]), float(LL[ii]))
  LLNormDist.SetBinContent(LLNormDist.FindBin(MTop[ii]),float(LLNorm[ii]))
  LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(MTop[ii]),float(LLNormAcc[ii]))
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
if GenLevel == False: 
  LLNormAccDist.Write()
  print " Found minimums (LL, LLXS & LLAcc) : ",LLMin,", ",LLNormMin," and ",LLNormAccMin
else:
  print "Found minimums (LL & LLXS) :",LLMin,", ",LLNormMin

#Change bin content for combined canvas!
for ii in range(len(MTop)):
  LLDist.SetBinContent(LLDist.FindBin(MTop[ii]), 2*(float(LLDist.GetBinContent(LLDist.FindBin(MTop[ii]))) - LLMin))
  LLNormDist.SetBinContent(LLNormDist.FindBin(MTop[ii]), 2*(float(LLNormDist.GetBinContent(LLNormDist.FindBin(MTop[ii]))) - LLNormMin))
  LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(MTop[ii]), 2*(float(LLNormAccDist.GetBinContent(LLNormAccDist.FindBin(MTop[ii]))) - LLNormAccMin))
  
#Prepare everything for the Delta ChiSq canvas!
LLDist.GetYaxis().SetTitle("#Delta #chi^{2} = -2ln(L/L_{max})")
LLDist.GetYaxis().SetTitleOffset(1.3)
                
LLCanvas = TCanvas('LLDistributions','Comparing influence of XS and acceptance normalisation on -ln(L) -- '+title+' evts')
legend = TLegend(0.75,0.7,0.95,0.95)
LLCanvas.cd()
legend.AddEntry(LLDist,'Non-normalized -ln(likelihood)',"p")
legend.AddEntry(LLNormDist,'Cross-section normalized -ln(likelihood)',"p")
if GenLevel == False: legend.AddEntry(LLNormAccDist,'Cross-section and acceptance normalized -ln(likelihood)',"p")
LLDist.SetTitle(LLCanvas.GetTitle())
LLDist.Draw("P")
LLNormDist.Draw("SAMEP")
if GenLevel == False: LLNormAccDist.Draw("SAMEP")
legend.Draw()
LLCanvas.Write() 
