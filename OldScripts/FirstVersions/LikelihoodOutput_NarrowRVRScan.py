! python
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
RVRValues = ["Re(V_{R}) = -0.3","Re(V_{R}) = -0.2","Re(V_{R}) = -0.1","Re(V_{R}) = -0.05","Re(V_{R}) = 0.0","Re(V_{R}) = 0.05","Re(V_{R}) = 0.1","Re(V_{R}) = 0.2","Re(V_{R}) = 0.3"]
RVR =       [-0.3,              -0.2,              -0.1,              -0.05,              0.0,              0.05,              0.1,              0.2,              0.3              ]
MGXS =      [13.3944,           12.06555,          11.25909,          11.02784,           10.90059,         10.88228,          10.97767,         11.49883,         12.49056         ]
MGXSe =     [0.00995808028337,  0.0093464076837,   0.00836607833038,  0.0,                0.00822214433527, 0.0                0.00847293509122, 0.00901976602967, 0.00874682264197 ]
#-- Remark --#
# No real need to add acceptance information here since for reco-level events always a couple of events fail the MadWeight computation! -#

NrEvents = 10000 

#File of interest:
LLFile = open(os.path.join(whichDir+'un-normalized_likelihood.out'),'r') 

#ROOT file where all the information will be stored:
Tfile = TFile(os.path.join(whichDir+"Likelihood_RVRNarrowScan.root"),'recreate')
LLDist = TH1F('LL_RVR','Non-normalized -ln(likelihood) for Single Gaus TF',11,-0.55,0.55)
LLNormDist = TH1F('NormLL_RVR','Normalized -ln(likelihood) for Single Gaus TF',11,-0.55,0.55)
gStyle.SetOptStat(0)

LLMin = []
LLNormMin = []
for ii in range(LLDist.GetNbinsX()):
    LLDist.SetBinContent(ii+1,-100)
    LLNormDist.SetBinContent(ii+1,-100)

#Loop over all lines in likelihood file:
for LLLine in LLFile:
    LLWord = LLLine.split()
    #Only interested in files starting with a number
    if str(LLWord[0]) != "#":
        if str(LLWord[1]) != "inf":
            print "Looking at bin ",LLDist.FindBin(RVR[int(LLWord[0])-1])," corresponding to RVR = ",RVR[int(LLWord[0])-1],", which has content ",LLWord[1]
            LLDist.SetBinContent(LLDist.FindBin(RVR[int(LLWord[0])-1]),float(LLWord[1]))
            LLNormDist.SetBinContent(LLDist.FindBin(RVR[int(LLWord[0])-1]),float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1]))
            if int(LLWord[0]) ==  1:
                LLMin = float(LLWord[1])
                LLNormMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])	
            if float(LLWord[1]) < LLMin:
                LLMin = float(LLWord[1])
            if (float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])) < LLNormMin:
                LLNormMin = float(LLWord[1]) + int(NrEvents)*log(MGXS[int(LLWord[0])-1])
LLDist.Write()
LLNormDist.Write()

print " Found minimums : ",LLMin," and ",LLNormMin
print "--> LLDist values :"
for ii in range(LLDist.GetNbinsX()):
  print " ",int(ii)+1," ) Content = ",LLDist.GetBinContent(ii+1)

#Change bin content for combined canvas!
for ii in range(len(RVR)):
    LLDist.SetBinContent(LLDist.FindBin(RVR[ii]), 2*(float(LLDist.GetBinContent(LLDist.FindBin(RVR[ii]))) - LLMin))
    LLNormDist.SetBinContent(LLNormDist.FindBin(RVR[ii]), 2*(float(LLNormDist.GetBinContent(LLNormDist.FindBin(RVR[ii]))) - LLNormMin))

LLDist.SetMinimum(0)
LLDist.SetMarkerStyle(20)
LLDist.SetLineColor(1)
LLDist.SetMarkerColor(1)
LLNormDist.SetMarkerStyle(21)
LLNormDist.SetLineColor(2)
LLNormDist.SetMarkerColor(2)

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
LLDist.SetTitle("Comparing normalized with non-normalized -ln(likelihood)")
LLDist.Draw("P")
LLNormDist.Draw("SAMEP")
legend.Draw()
LLCanvas.Write() 
