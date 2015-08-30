#! python
import os
import sys
from ROOT import TH1F, TCanvas, TFile

#Directory of interest:
if len(sys.argv) <= 2:
  print "Need to give the directory of interest and the desired pT-cut in command line !"
  sys.exit()
if len(sys.argv) == 1:
  print "Need to give the directory of interest AND the desired pT-cut in command line !"
  sys.exit()
whichDir = sys.argv[1] 
ptCut = sys.argv[2]
print " Interested in directory : ",whichDir

#Input file of interest:
list_dir = []
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".lhco"):
    LHCOFile = open(os.path.join(whichDir+''+file),'r')
print "Will be using lhco file : ",LHCOFile

EventsToDelete = []
h_LeptBPt    = TH1F("PtLeptB","Pt distribution of leptonic b",           250,0,200)
h_LeptonPt   = TH1F("PtLepton","Pt distribution of lepton",              250,0,200)
h_METPt      = TH1F("PtMET","Pt distribution of MET",                    250,0,200)
h_HadrBPt    = TH1F("PtHadrB","Pt distribution of hadronic b",           250,0,200)
h_LightOnePt = TH1F("PtLightOne","Pt distribution of first light quark", 250,0,200)
h_LightTwoPt = TH1F("PtLightTwo","Pt distribution of second light quark",250,0,200)

h_LeptBEta    = TH1F("EtaLeptB","Eta distribution of leptonic b",           250,-6,6)
h_LeptonEta   = TH1F("EtaLepton","Eta distribution of lepton",              250,-6,6)
h_METEta      = TH1F("EtaMET","Eta distribution of MET",                    250,-6,6)
h_HadrBEta    = TH1F("EtaHadrB","Eta distribution of hadronic b",           250,-6,6)
h_LightOneEta = TH1F("EtaLightOne","Eta distribution of first light quark", 250,-6,6)
h_LightTwoEta = TH1F("EtaLightTwo","Eta distribution of second light quark",250,-6,6)

h_LeptBPhi    = TH1F("PhiLeptB","Phi distribution of leptonic b",           250,-0.5,7)
h_LeptonPhi   = TH1F("PhiLepton","Phi distribution of lepton",              250,-0.5,7)
h_METPhi      = TH1F("PhiMET","Phi distribution of MET",                    250,-0.5,7)
h_HadrBPhi    = TH1F("PhiHadrB","Phi distribution of hadronic b",           250,-0.5,7)
h_LightOnePhi = TH1F("PhiLightOne","Phi distribution of first light quark", 250,-0.5,7)
h_LightTwoPhi = TH1F("PhiLightTwo","Phi distribution of second light quark",250,-0.5,7)

#Loop over all lines in LHCO file:
EvtNr = 0
SurvivedCut = True
LeptBPt,  LeptonPt,  METPt,  HadrBPt,  LightOnePt,  LightTwoPt  = 0, 0, 0, 0, 0, 0
LeptBEta, LeptonEta, METEta, HadrBEta, LightOneEta, LightTwoEta = 0, 0, 0, 0, 0, 0
LeptBPhi, LeptonPhi, METPhi, HadrBPhi, LightOnePhi, LightTwoPhi = 0, 0, 0, 0, 0, 0
for lhcoLine in LHCOFile:
  lhcoWord = lhcoLine.split()
  #Only interested in files starting with a number
  if str(lhcoWord[0]) != "#":
    if len(lhcoWord) == 3:
      if str(lhcoWord[1]) != "0":
        if SurvivedCut == True:
          #Only fill the histograms if the entire event survived the cut!
          h_LeptBPt.Fill(float(LeptBPt)),       h_LeptBEta.Fill(float(LeptBEta)),       h_LeptBPhi.Fill(float(LeptBPhi))
          h_LeptonPt.Fill(float(LeptonPt)),     h_LeptonEta.Fill(float(LeptonEta)),     h_LeptonPhi.Fill(float(LeptonPhi))
          h_METPt.Fill(float(METPt)),           h_METEta.Fill(float(METEta)),           h_METPhi.Fill(float(METPhi))
          h_HadrBPt.Fill(float(HadrBPt)),       h_HadrBEta.Fill(float(HadrBEta)),       h_HadrBPhi.Fill(float(HadrBPhi))
          h_LightOnePt.Fill(float(LightOnePt)), h_LightOneEta.Fill(float(LightOneEta)), h_LightOnePhi.Fill(float(LightOnePhi))
          h_LightTwoPt.Fill(float(LightTwoPt)), h_LightTwoEta.Fill(float(LightTwoEta)), h_LightTwoPhi.Fill(float(LightTwoPhi))
      EvtNr = lhcoWord[1]
      SurvivedCut = True
    if len(lhcoWord) == 11:
      if   str(lhcoWord[0]) == "1": LeptBPt,    LeptBEta,    LeptBPhi    = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])
      elif str(lhcoWord[0]) == "2": LeptonPt,   LeptonEta,   LeptonPhi   = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])
      elif str(lhcoWord[0]) == "3": METPt,      METEta,      METPhi      = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])
      elif str(lhcoWord[0]) == "4": HadrBPt,    HadrBEta,    HadrBPhi    = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])
      elif str(lhcoWord[0]) == "5": LightOnePt, LightOneEta, LightOnePhi = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])
      elif str(lhcoWord[0]) == "6": LightTwoPt, LightTwoEta, LightTwoPhi = float(lhcoWord[4]), float(lhcoWord[2]), float(lhcoWord[3])

      if SurvivedCut == True and str(lhcoWord[1]) != "6" and float(lhcoWord[4]) < float(ptCut):
        SurvivedCut = False
        EventsToDelete.append(int(EvtNr))

#canv_PtLeptB, canv_PtLepton, canv_PtMET, canv_PtHadrB, canv_PtLightOne, canv_PtLightTwo = TCanvas("LeptBCanv_Pt","LeptBCanv_Pt"), TCanvas("LeptonCanv_Pt","LeptonCanv_Pt"), TCanvas("METCanv_Pt","METCanv_Pt"), TCanvas("HadrBCanv_Pt","HadrBCanv_Pt"), TCanvas("LightOneCanv_Pt","LightOneCanv_Pt"), TCanvas("LightTwoCanv_Pt","LightTwoCanv_Pt")

#canv_PtLeptB.cd(),    h_LeptBPt.Draw(),    canv_PtLeptB.SaveAs( str(whichDir)+'PtLeptB_PtCut'+str(ptCut)+'.pdf')
#canv_PtLepton.cd(),   h_LeptonPt.Draw(),   canv_PtLepton.SaveAs( str(whichDir)+'PtLepton_PtCut'+str(ptCut)+'.pdf')
#canv_PtMET.cd(),      h_METPt.Draw(),      canv_PtMET.SaveAs( str(whichDir)+'PtMET_PtCut'+str(ptCut)+'.pdf')
#canv_PtHadrB.cd(),    h_HadrBPt.Draw(),    canv_PtHadrB.SaveAs( str(whichDir)+'PtHadrB_PtCut'+str(ptCut)+'.pdf')
#canv_PtLightOne.cd(), h_LightOnePt.Draw(), canv_PtLightOne.SaveAs( str(whichDir)+'PtLightOne_PtCut'+str(ptCut)+'.pdf')
#canv_PtLightTwo.cd(), h_LightTwoPt.Draw(), canv_PtLightTwo.SaveAs( str(whichDir)+'PtLightTwo_PtCut'+str(ptCut)+'.pdf')

#canv_EtaLeptB, canv_EtaLepton, canv_EtaMET, canv_EtaHadrB, canv_EtaLightOne, canv_EtaLightTwo = TCanvas("LeptBCanv_Eta","LeptBCanv_Eta"), TCanvas("LeptonCanv_Eta","LeptonCanv_Eta"), TCanvas("METCanv_Eta","METCanv_Eta"), TCanvas("HadrBCanv_Eta","HadrBCanv_Eta"), TCanvas("LightOneCanv_Eta","LightOneCanv_Eta"), TCanvas("LightTwoCanv_Eta","LightTwoCanv_Eta")

#canv_EtaLeptB.cd(),    h_LeptBEta.Draw(),    canv_EtaLeptB.SaveAs( str(whichDir)+'EtaLeptB_PtCut'+str(ptCut)+'.pdf')
#canv_EtaLepton.cd(),   h_LeptonEta.Draw(),   canv_EtaLepton.SaveAs( str(whichDir)+'EtaLepton_PtCut'+str(ptCut)+'.pdf')
#canv_EtaMET.cd(),      h_METEta.Draw(),      canv_EtaMET.SaveAs( str(whichDir)+'EtaMET_PtCut'+str(ptCut)+'.pdf')
#canv_EtaHadrB.cd(),    h_HadrBEta.Draw(),    canv_EtaHadrB.SaveAs( str(whichDir)+'EtaHadrB_PtCut'+str(ptCut)+'.pdf')
#canv_EtaLightOne.cd(), h_LightOneEta.Draw(), canv_EtaLightOne.SaveAs( str(whichDir)+'EtaLightOne_PtCut'+str(ptCut)+'.pdf')
#canv_EtaLightTwo.cd(), h_LightTwoEta.Draw(), canv_EtaLightTwo.SaveAs( str(whichDir)+'EtaLightTwo_PtCut'+str(ptCut)+'.pdf')

#canv_PhiLeptB, canv_PhiLepton, canv_PhiMET, canv_PhiHadrB, canv_PhiLightOne, canv_PhiLightTwo = TCanvas("LeptBCanv_Phi","LeptBCanv_Phi"), TCanvas("LeptonCanv_Phi","LeptonCanv_Phi"), TCanvas("METCanv_Phi","METCanv_Phi"), TCanvas("HadrBCanv_Phi","HadrBCanv_Phi"), TCanvas("LightOneCanv_Phi","LightOneCanv_Phi"), TCanvas("LightTwoCanv_Phi","LightTwoCanv_Phi")

#canv_PhiLeptB.cd(),    h_LeptBPhi.Draw(),    canv_PhiLeptB.SaveAs( str(whichDir)+'PhiLeptB_PtCut'+str(ptCut)+'.pdf')
#canv_PhiLepton.cd(),   h_LeptonPhi.Draw(),   canv_PhiLepton.SaveAs( str(whichDir)+'PhiLepton_PtCut'+str(ptCut)+'.pdf')
#canv_PhiMET.cd(),      h_METPhi.Draw(),      canv_PhiMET.SaveAs( str(whichDir)+'PhiMET_PtCut'+str(ptCut)+'.pdf')
#canv_PhiHadrB.cd(),    h_HadrBPhi.Draw(),    canv_PhiHadrB.SaveAs( str(whichDir)+'PhiHadrB_PtCut'+str(ptCut)+'.pdf')
#canv_PhiLightOne.cd(), h_LightOnePhi.Draw(), canv_PhiLightOne.SaveAs( str(whichDir)+'PhiLightOne_PtCut'+str(ptCut)+'.pdf')
#canv_PhiLightTwo.cd(), h_LightTwoPhi.Draw(), canv_PhiLightTwo.SaveAs( str(whichDir)+'PhiLightTwo_PtCut'+str(ptCut)+'.pdf')

#Save all histograms in ROOT file
rootFile = TFile(str(whichDir)+'KinematicInformation_PtCut'+str(ptCut)+'.root',"RECREATE")
h_LeptBPt.Write(),  h_LeptonPt.Write(),  h_METPt.Write(),  h_HadrBPt.Write(),  h_LightOnePt.Write(),  h_LightTwoPt.Write()
h_LeptBEta.Write(), h_LeptonEta.Write(), h_METEta.Write(), h_HadrBEta.Write(), h_LightOneEta.Write(), h_LightTwoEta.Write()
h_LeptBPhi.Write(), h_LeptonPhi.Write(), h_METPhi.Write(), h_HadrBPhi.Write(), h_LightOnePhi.Write(), h_LightTwoPhi.Write()

print "Events with a particle with pT < ",str(ptCut)," GeV (excluding MET) :",len(EventsToDelete)
