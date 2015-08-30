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

NrEvents = 10000 

#File of interest:
LLFile = open(os.path.join(whichDir+'un-normalized_likelihood.out'),'r') 

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"Likelihood_TopMass.root"),'recreate')
LLDist = TH1F('LL_TopMass','Non-normalized -ln(likelihood) for Single Gaus TF',41,152.5,193.5)
LLNormDist = TH1F('NormLL_TopMass','Normalized -ln(likelihood) for Single Gaus TF',41,152.5,193.5)
LLNormAccDist = TH1F('NormLL_Acceptance_TopMass','Normalized -ln(likelihood) for Single Gaus TF (acceptance taken into account)',41,152.5,193.5)
gStyle.SetOptStat(0)

LLMin = []
LLNormMin = []
LLNormAcc = []
for ii in range(LLDist.GetNbinsX()):
  LLDist.SetBinContent(ii+1,-100)
  LLNormDist.SetBinContent(ii+1,-100)
  LLNormAccDist.SetBinContent(ii+1,-1000)  

#Loop over all lines in likelihood file:
for LLLine in LLFile:
  LLWord = LLLine.split()
  #Only interested in files starting with a number
  if str(LLWord[0]) != "#":
    if str(LLWord[1]) != "inf":
      LLDist.SetBinContent(LLDist.FindBin(MTop[int(LLWord[0])-1]),float(LLWord[1]))
      LLNormDist.SetBinContent(LLNormDist.FindBin(MTop[int(LLWord[0])-1]),float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1]))
      LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(MTop[int(LLWord[0])-1]),float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1]) + int(NrEvents)*log(Acceptance[int(WeightWord[1])-1])) 
      if int(LLWord[0]) ==  1:
        LLMin = float(LLWord[1])
        LLNormMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])	
        LLNormAccMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1]) + int(NrEvents)*log(Acceptance[int(LLWord[0])-1)
      if float(LLWord[1]) < LLMin:
        LLMin = float(LLWord[1])
      if (float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])) < LLNormMin:
        LLNormMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])
      if (float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1]) + int(NrEvents)*log(Acceptance[int(WeightWord[1])-1])) < LLNormAccMin:
        LLNormAccMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1] + int(NrEvents)*log(Acceptance[int(WeightWord[1])-1]))

LLDist.Write()
LLDist.SetMinimum(0)
LLDist.SetMarkerStyle(20)
LLDist.SetLineColor(1)
LLDist.SetMarkerColor(1)
LLNormDist.SetMarkerStyle(21)
LLNormDist.SetLineColor(2)
LLNormDist.SetMarkerColor(2)
LLNormAccDist.SetMinimum(0)
LLNormAccDist.SetMarkerStyle(22)
LLNormAccDist.SetLineColor(3)
LLNormAccDist.SetMarkerColor(3)

LLNormDist.Write()
LLNormAccDist.Write()
print " Found minimums : ",LLMin,", ",LLNormMin," and ",LLNormAccMin

#Change bin content for combined canvas!
for ii in range(len(MTop)):
    LLDist.SetBinContent(LLDist.FindBin(MTop[ii]), 2*(float(LLDist.GetBinContent(LLDist.FindBin(MTop[ii]))) - LLMin))
    LLNormDist.SetBinContent(LLNormDist.FindBin(MTop[ii]), 2*(float(LLNormDist.GetBinContent(LLNormDist.FindBin(MTop[ii]))) - LLNormMin))
    LLNormAccDist.SetBinContent(LLNormAccDist.FindBin(MTop[ii]), 2*(float(LLNormAccDist.GetBinContent(LLNormAccDist.FindBin(MTop[ii]))) - LLNormAccMin))

#Scale both histograms to 1000 entries to plot them together!
#LLDist.Scale(1000/LLDist.Integral())
LLDist.GetYaxis().SetTitle("#Delta #chi^{2} = -2ln(L/L_{max})")
LLDist.GetYaxis().SetTitleOffset(1.3)
#LLNormDist.Scale(1000/LLNormDist.Integral())
                
LLCanvas = TCanvas('LLDistributions','Comparing normalized with non-normalized -ln(likelihood)')
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
