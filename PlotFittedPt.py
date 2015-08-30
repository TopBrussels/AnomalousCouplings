#! python
import re
import os
import shutil
import linecache
import math
import sys
#import ROOT
from ROOT import TH1F,TFile,TCanvas

if len(sys.argv) == 1:                               #Python file is seen as 0th argument!
	print "No additional arguments given! --> Need to choose the desired directory"
	os.system('ls -ltr Events')
	print "--------------------------------------------------"
	whichDir = raw_input("Which directory should be used ?  ")
	whichDir = "Events/"+whichDir+"/"
	print " "
else:
	#print "First argument = Directory which will be used:",sys.argv[1]
	whichDir = sys.argv[1]

print "Chosen directory = ",whichDir

#output file used for information:
XMLFile = open(os.path.join(whichDir+'output.xml'),'r')
#XMLFile = open(os.path.join(whichDir+'output_PTOutput.xml'),'r')
#XMLFile = open(os.path.join(whichDir+'ShortFile.xml'),'r')

#ROOT file were all the created histograms will be stored!
Tfile = TFile(os.path.join(whichDir+"FittedPtValuesMW.root"),'recreate')

configCounter=0     #Count the different configurations considered by MadWeight
channelCounter = 1
pTExpUnc = []
pTActualValue = []
ThetaExpUnc = []
ThetaActualValue = []
PhiExpUnc = []
PhiActualValue = []
#InitialPt = [0 for i in range(6)]
#Create 1D histogram:
FittedPtDistribution = TH1F('FittedPt_PermX','Distribution of pT for the different integration points of MadWeight',751,-0.5,750.5)
PtDistCanvas = TCanvas('PtDistCanvas','Comparison between fitted pT and actual pT')
ActualPtPoint = TH1F('ActualPtPoint','Actual Pt point considered by MadWeight',750,0,750)
ActualPtPoint.SetMarkerColor(2)
ActualPtPoint.SetMarkerStyle(34)

for XMLLine in XMLFile:
    XMLWord = XMLLine.split()
    #Get the name of the used Transfer Function
    if re.search( r"Current parametrization", XMLLine):
        print "Name of TF :",XMLWord[4]
        ConsideredTF = XMLWord[4]
    #Get the original Pt information
    #if re.search( r"MG =",XMLLine) and configCounter == 1:
    #    InitialPt[int(XMLWord[6])-3] = str(XMLWord[0])
    #Initialize a counter for the different configurations, and plot histograms for each separately
    if re.search( r"tfset id",XMLLine):
        if configCounter != 0:
            #Set the name correct for each distribution (add configuration number)	
            FittedPtDistribution.SetTitle(str('Transfer Function '+ConsideredTF+' with width = '+str(pTExpUnc[ConsideredParticle-3])+'-- Config '+str(configCounter)))
            FittedPtDistribution.SetName('FittedPt_Config'+str(configCounter))
            FittedPtDistribution.Write()                   #Write the existing distribution
            #Print the actual point and the fitted distribution in 1 canvas
            PtDistCanvas.SetName('PtDistCanvas_Config'+str(configCounter))
            PtDistCanvas.cd()
            #Set the bin-content corresponding to the real value equal to the maximum
            ActualPtPoint.SetTitle(str('Actual Pt point and fitted values using '+ConsideredTF+' TF (width = '+str(pTExpUnc[ConsideredParticle-3])+')'+'-- Config '+str(configCounter)))
            ActualPtPoint.SetBinContent(ActualPtPoint.FindBin(float(pTActualValue[ConsideredParticle-3])),float(FittedPtDistribution.GetMaximum()))
            ActualPtPoint.Draw("P")
            FittedPtDistribution.Draw("same")
            PtDistCanvas.Write()
            PtDistCanvas.Clear()
            PtDistCanvas = TCanvas('PtDistCanvas','Comparison between fitted pT and actual pT')
            #Reset the histogram such that it is empty for the next permutation!
            FittedPtDistribution.Reset()
            ActualPtPoint.Reset()
        configCounter=int(configCounter)+1
    #Initialize a counter for the different channels
    if re.search( r"Current channel of integration",XMLLine):
        channelCounter = XMLWord[5]
    #Get the width of the used Transfer Function (for all particles and kinematic variables)
    if re.search( r"exp. uncertainty",XMLLine):
        if XMLWord[3] == "theta" and configCounter == 1:
#            print "Uncertainty on theta : ",XMLWord[8],"configCounter =",configCounter,"channelCounter = ",channelCounter   #Make sure the first time that uncertainties are identical!
            if channelCounter == 1:
                ThetaExpUnc.append( float(XMLWord[8]))
		ThetaActualValue.append(float(XMLWord[10]))
        if XMLWord[3] == "phi" and configCounter == 1:
#            print "Uncertainty on phi : ",XMLWord[8],"configCounter =",configCounter,"channelCounter = ",channelCounter   #Make sure the first time that uncertainties are identical!
            if channelCounter == 1:
                PhiExpUnc.append( float(XMLWord[8]))
		PhiActualValue.append(float(XMLWord[10]))
        if XMLWord[3] == "pT" and configCounter == 1:
#            print "Uncertainty on p : ",XMLWord[8],"configCounter =",configCounter,"channelCounter = ",channelCounter   #Make sure the first time that uncertainties are identical!
            if channelCounter == 1:
                pTExpUnc.append( float(XMLWord[8]))
		pTActualValue.append(float(XMLWord[10]))
    #Fill the actual distribution with the pT value!
    if re.search( r"transfer_fct", XMLLine):
        for part in xrange(3,9):   #Corresponds with the list [3,4,5,6,7,8]
            if str(part) in str(XMLWord[2]):
                ConsideredParticle = int(part)
        #Change the range of the axis!
        #FittedPtDistribution.GetXaxis().Set(FittedPtDistribution.GetNbinsX(),max(float(pTActualValue[ConsideredParticle-3])-10*float(pTExpUnc[ConsideredParticle-3]),0),min(float(pTActualValue[ConsideredParticle-3])+10*float(pTExpUnc[ConsideredParticle-3]),750))
        #print "Optimal range:",float(pTActualValue[ConsideredParticle-3])-10*float(pTExpUnc[ConsideredParticle-3]),"-",float(pTActualValue[ConsideredParticle-3])+10*float(pTExpUnc[ConsideredParticle-3])
        #print "Bin Content of bin belonging to pT value ",float(XMLWord[4])," :",FittedPtDistribution.GetBinContent(FittedPtDistribution.FindBin(float(XMLWord[4])))
        #if FittedPtDistribution.GetBinContent(FittedPtDistribution.FindBin(float(XMLWord[4]))) == 0:
        FittedPtDistribution.Fill(float(XMLWord[4]))
print "Full list of uncertainties on theta: ",ThetaExpUnc,"versus actual value of ",ThetaActualValue
print "Full list of uncertainties on phi  : ",PhiExpUnc,"versus actual value of ",PhiActualValue
print "Full list of uncertainties on Pt   : ",pTExpUnc,"versus actual value of ",pTActualValue

#Set the name correct for each distribution (add configuration number)	
FittedPtDistribution.SetTitle(str('Transfer Function '+ConsideredTF+' with width = '+str(pTExpUnc[ConsideredParticle-3])+'-- config '+str(configCounter)))
FittedPtDistribution.SetName('FittedPt_Config'+str(configCounter))
FittedPtDistribution.Write()

#Set the bin-content corresponding to the real value equal to the maximum
ActualPtPoint.SetTitle(str('Actual Pt point and fitted values using '+ConsideredTF+' TF (width = '+str(pTExpUnc[ConsideredParticle-3])+')'+'-- Config '+str(configCounter)))
ActualPtPoint.SetBinContent(ActualPtPoint.FindBin(float(pTActualValue[ConsideredParticle-3])),float(FittedPtDistribution.GetMaximum()))

#Print the actual point and the fitted distribution in 1 canvas
PtDistCanvas.SetName('PtDistCanvas_Config'+str(configCounter))
PtDistCanvas.cd()
ActualPtPoint.Draw("P")
FittedPtDistribution.Draw("same")
PtDistCanvas.Write()

shutil.move(Tfile.GetName(),whichDir+'FittedMWPt_TF'+str(ConsideredTF)+'_Width'+str(pTExpUnc[ConsideredParticle-3])+'.root')

#Tfile.Clone(os.path.join(whichDir+"FittedMWPt_TF"+str(ConsideredTF)+"_Width"+str(pTExpUnc[ConsideredParticle-3])+".root"))i
#Tfile.SaveAs(os.path.join(whichDir+"FittedMWPt_TF"+str(ConsideredTF)+"_Width"+str(pTExpUnc[ConsideredParticle-3])+".root"))
