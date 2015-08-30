#----------------------------------------------#
#   Python script which performs polynomial    #
#   fit on a limited range of the LL TH1F      #
#					       #
#   Needs to be ran using:                     #
#      python FitLikelihood.py Events/blabla/  #
#      MTop (or RVR) y/n (for cut)             #
#                                              #
# ** Update 14/05/2015:                        #
#  - File is replaced by FitDeviationScript.py #
#    (which uses TGraphs, so no empty bins ..) #
#----------------------------------------------#
#! python
import os
import sys
import re           #used for search!
from ROOT import TFile #TH1F,TFile,TCanvas,TLegend,gStyle

#Get all the input from the command line:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest, the considered kinematic variable and whether scdDer cuts should be applied in the command line !"
  print " Correct syntax is : python FitLikelihood.py Events/blabla/ MTop (WithCut)"
  sys.exit()
elif len(sys.argv) == 2:
  print "Need to specify the considered kinematic variable (MTop or RVR)"
  print " Correct syntax is : python FitLikelihood.py Events/blabla/ MTop (WithCut)"
  sys.exit()
elif len(sys.argv) == 3:
  print " --> Considering -ln(L) distribution without any scdDer cut applied !"
  CutApplied = False
else:
  print " --> Considering -ln(L) distribution with scdDer cut applied !"
  CutApplied = True
  whichScdDer = raw_input('Specify which of the second derivatives of the -ln(L) should be positive! \n    1) Only the inner one \n    2) Only the outer one \n    3) Both of them \n    4) No cut desired ... \n --> Choose the corresponding number from the list : ')
  if whichScdDer == "4": CutApplied, AppliedCut = False, ""
  elif whichScdDer == "1": AppliedCut = "ScdDerInner"
  elif whichScdDer == "2": AppliedCut = "ScdDerOuter"
  elif whichScdDer == "3": AppliedCut = "ScdDerBoth"
whichDir = sys.argv[1]
KinVariable = sys.argv[2]

if KinVariable != "MTop" and KinVariable != "RVR":
  print "Need to specify which kinematic variable should be considered (MTop or RVR are the only options!!)"
  KinVariable = raw_input('--> Choose one of the two : ')

#Now change the ROOT analyzer to use the correct directory!!
#--File of interest:
list_dir = os.listdir(whichDir)
RootFileArray, RootFileCounter = [], 0
for file in list_dir:
  if CutApplied == False and file.endswith(".root") and file.startswith("Likelihood"): # eg: '.txt'
    RootFileCounter += 1
    RootFileArray.append(file)
  if CutApplied == True and file.endswith(".root") and file.startswith("FitDeviation"):
    RootFileCounter += 1
    RootFileArray.append(file)

if int(RootFileCounter) == 1:
  RootFile = whichDir+''+RootFileArray[0]
elif int(RootFileCounter) == 0:
  print "No ROOT file found in this directory !"
  sys.exit()
elif int(RootFileCounter) > 1:
  for ii in range(len(RootFileArray)):
    print " ",ii," ) ",RootFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  RootFile = whichDir+''+RootFileArray[int(fileNr)]
print "Will be using file : ",RootFile

#Set title of root file!
title = ""
if whichDir.find("Correct") <= len(whichDir)     and whichDir.find("Correct") > 0:   title = "Correct"
elif whichDir.find("Wrong") <= len(whichDir)     and whichDir.find("Wrong") > 0:     title = "Wrong"
elif whichDir.find("Unmatched") <= len(whichDir) and whichDir.find("Unmatched") > 0: title = "Unmatched"

GenLevel = "false"
if whichDir.find("Reco") <= len(whichDir)  and whichDir.find("Reco") > 0: title += "Reco"
elif whichDir.find("Gen") <= len(whichDir) and whichDir.find("Gen") > 0:
  title += "Gen"
  GenLevel = "true"
title = title+"_"+KinVariable
if CutApplied == True: title = title+"_"+AppliedCut+"CutApplied"
print "GenLevel = ",GenLevel

#-Identify kinematic variable considered (mtop or RVR)
if KinVariable == "RVR":
  KinVar = "$Re(V_R)$"
  FitValues = [-0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5]

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n  3) Many   : [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1] \n --> Choose the correct number : ')
  
  if VarWindow == "1":   FitValues.pop(6), FitValues.pop(4)
  elif VarWindow == "2": FitValues.pop(10), FitValues.pop(0)
  FitValues.pop(8), FitValues.pop(7), FitValues.pop(1), FitValues.pop(0)     #Fit will only be applied on the inner points!
  if VarWindow == "3": FitValues = [-0.1, -0.07, -0.05, -0.02, 0.0, 0.02, 0.05, 0.07, 0.1]

elif KinVariable == "MTop":
  KinVar = "$m_{top}$"
  FitValues = [171,172,173,174,175]

print "FitValues = ",FitValues, "which has length = ",len(FitValues)
fitStep = (float(FitValues[len(FitValues)-1]-FitValues[0]+abs(FitValues[0]-FitValues[1])))/float(len(FitValues))
fitRange = [FitValues[0]-fitStep/2.0,FitValues[len(FitValues)-1]+fitStep/2.0] 
print " Considered fit range is : ", fitRange 

#--Replace ROOT file information in analyzer
RootAnalyzer = open('fitExcludeEmptyBins.C','r')

NewRootAnalyzer = open('output','w')
for RootLine in RootAnalyzer:
  RootWord = RootLine.split()
  if re.search( r".root", RootLine):
    if RootWord[1] == "InputFile":  NewRootAnalyzer.write(RootLine.replace(RootWord[5],RootFile))
    if RootWord[1] == "OutputFile": NewRootAnalyzer.write(RootLine.replace(RootWord[5],whichDir+'LimitedFitResult_'+title+'.root'))
  #Set the correct start and end values for the fit-histogram!
  elif re.search( r"Initialize",RootLine):
    if re.search( r"nrBins", RootLine):     NewRootAnalyzer.write('  int nBins = '+str(len(FitValues))+';         //Initialize nrBins \n')
    elif re.search( r"fitStart", RootLine): NewRootAnalyzer.write('  float fitStart = '+str(fitRange[0])+';  //Initialize fitStart \n')
    elif re.search( r"fitEnd", RootLine):   NewRootAnalyzer.write('  float fitEnd = '+str(fitRange[1])+';    //Initialize fitEnd \n')
  #Get the correct histogram from the InputFile!
  elif re.search( r"TH1F", RootLine) and re.search( r"InputFile",RootLine):
    if CutApplied == False:
      if re.search( r"XS",RootLine):    NewRootAnalyzer.write('  TH1F* LLXS  = (TH1F*) InputFile->Get("LL_XS"); \n')
      elif re.search( r"Acc",RootLine): NewRootAnalyzer.write('  TH1F* LLAcc = (TH1F*) InputFile->Get("LL_Acc"); \n')
      else:                             NewRootAnalyzer.write('  TH1F* LL    = (TH1F*) InputFile->Get("LL"); \n')
    elif CutApplied == True:
      if re.search( r"XS",RootLine):    NewRootAnalyzer.write('  TH1F* LLXS  = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLXSPos'+AppliedCut+'"); \n') 
      elif re.search( r"Acc",RootLine): NewRootAnalyzer.write('  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLAccPos'+AppliedCut+'"); \n')
      else:                             NewRootAnalyzer.write('  TH1F* LL    = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPos'+AppliedCut+'"); \n')
  elif re.search( r" bool GenLevel",RootLine):
    NewRootAnalyzer.write('  bool GenLevel = '+GenLevel+'; \n')
  else:
    NewRootAnalyzer.write(RootLine)
RootAnalyzer.close()
NewRootAnalyzer.close()

os.rename('output','fitExcludeEmptyBins.C'), os.system("root -l -b -q fitExcludeEmptyBins.C")

#-- Now open the created ROOT file and place the minimum (with uncertainty) in a table! --#
FitFile = TFile(os.path.join(whichDir+"LimitedFitResult_"+title+".root"),'r')
TableOutput = open(os.path.join(whichDir+'FitOutput_'+title+'.tex'),'w')

TableOutput.write('\\begin{table}[h!t] \n \\centering \n \\caption{Fit parameters of 2nd degree polynomial ($a_{0} + a_{1}*x + a_{2}*x^{2}$) and corresponding minimum for '+title+' events.} \\label{table::} \n \\begin{tabular}{c|c|c|c|c} \n')
TableOutput.write('  & $a_{0}$ & $a_{1}$ & $a_{2}$ & '+KinVar+' \\\\ \n  \hline \n')

fit = FitFile.Get("fit_LL")
VarPlus1Sig = (fit.GetX(fit.GetMinimum()+0.5, fit.GetMinimumX(), fitRange[1]) - fit.GetMinimumX())
VarMin1Sig =  (fit.GetMinimumX() - fit.GetX(fit.GetMinimum()+0.5, fitRange[0], fit.GetMinimumX()))
if abs(abs(VarPlus1Sig) - abs(VarMin1Sig)) > 0.00001:
  if KinVariable == "RVR": 
    TableOutput.write('  no normalisation & %.3f & %.3f & %.3f & $%.5f^{+%.6f}_{-%.6f}$  \\\\ \n' %( fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2), fit.GetMinimumX(fitRange[0],fitRange[1]), VarPlus1Sig, VarMin1Sig ))
  elif KinVariable == "MTop": 
    TableOutput.write('  no normalisation & %.3f & %.3f & %.3f & $%.3f^{+%.4f}_{-%.4f}$  \\\\ \n' %( fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2), fit.GetMinimumX(fitRange[0],fitRange[1]), VarPlus1Sig, VarMin1Sig ))
else:
  if KinVariable == "RVR": 
    TableOutput.write('  no normalisation & %.3f & %.3f & %.3f & %.5f $\\pm$ %.6f  \\\\ \n' %( fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2), fit.GetMinimumX(fitRange[0],fitRange[1]), abs(VarPlus1Sig) ))
  elif KinVariable == "MTop": 
    TableOutput.write('  no normalisation & %.3f & %.3f & %.3f & %.3f $\\pm$ %.4f  \\\\ \n' %( fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2), fit.GetMinimumX(fitRange[0],fitRange[1]), abs(VarPlus1Sig) ))

fitXS = FitFile.Get("fit_LLXS")
VarXSPlus1Sig = (fitXS.GetX(fitXS.GetMinimum()+0.5, fitXS.GetMinimumX(), fitRange[1]) - fitXS.GetMinimumX())
VarXSMin1Sig =  (fitXS.GetMinimumX() - fitXS.GetX(fitXS.GetMinimum()+0.5, fitRange[0], fitXS.GetMinimumX()))
if abs(abs(VarXSPlus1Sig) - abs(VarXSMin1Sig)) > 0.00001:
  if KinVariable == "RVR": 
    TableOutput.write('  XS normalisation & %.3f & %.3f & %.3f & $%.5f^{+%.6f}_{-%.6f}$ ' %( fitXS.GetParameter(0), fitXS.GetParameter(1), fitXS.GetParameter(2), fitXS.GetMinimumX(fitRange[0],fitRange[1]), VarXSPlus1Sig, VarXSMin1Sig) )
  elif KinVariable == "MTop": 
    TableOutput.write('  XS normalisation & %.3f & %.3f & %.3f & $%.3f^{+%.4f}_{-%.4f}$ ' %( fitXS.GetParameter(0), fitXS.GetParameter(1), fitXS.GetParameter(2), fitXS.GetMinimumX(fitRange[0],fitRange[1]), VarXSPlus1Sig, VarXSMin1Sig) )
else:
  if KinVariable == "RVR": 
    TableOutput.write('  XS normalisation & %.3f & %.3f & %.3f & %.5f $\\pm$ %.6f ' %( fitXS.GetParameter(0), fitXS.GetParameter(1), fitXS.GetParameter(2), fitXS.GetMinimumX(fitRange[0],fitRange[1]), abs(VarXSPlus1Sig) ))
  elif KinVariable == "MTop": 
    TableOutput.write('  XS normalisation & %.3f & %.3f & %.3f & %.3f $\\pm$ %.4f ' %( fitXS.GetParameter(0), fitXS.GetParameter(1), fitXS.GetParameter(2), fitXS.GetMinimumX(fitRange[0],fitRange[1]), abs(VarXSPlus1Sig) ))

if GenLevel == "false":
  fitAcc = FitFile.Get("fit_LLAcc")
  VarAccPlus1Sig = (fitAcc.GetX(fitAcc.GetMinimum()+0.5, fitAcc.GetMinimumX(), fitRange[1]) - fitAcc.GetMinimumX())
  VarAccMin1Sig =  (fitAcc.GetMinimumX() - fitAcc.GetX(fitAcc.GetMinimum()+0.5, fitRange[0], fitAcc.GetMinimumX()))
  if abs(abs(VarAccPlus1Sig) - abs(VarAccMin1Sig)) > 0.00001:
    print "1Sigma interval not similar ! : ", abs(abs(VarAccPlus1Sig) - abs(VarAccMin1Sig))
    if KinVariable == "RVR": 
      TableOutput.write(' \\\\ \n  Acc normalisation & %.3f & %.3f & %.3f & $%.5f^{+%.6f}_{-%.6f}$' %( fitAcc.GetParameter(0), fitAcc.GetParameter(1), fitAcc.GetParameter(2), fitAcc.GetMinimumX(fitRange[0],fitRange[1]), VarAccPlus1Sig, VarAccMin1Sig))
    elif KinVariable == "MTop": 
      TableOutput.write(' \\\\ \n  Acc normalisation & %.3f & %.3f & %.3f & $%.3f^{+%.4f}_{-%.4f}$' %( fitAcc.GetParameter(0), fitAcc.GetParameter(1), fitAcc.GetParameter(2), fitAcc.GetMinimumX(fitRange[0],fitRange[1]), VarAccPlus1Sig, VarAccMin1Sig))
  else:
    if KinVariable == "RVR": 
      TableOutput.write(' \\\\ \n  Acc normalisation & %.3f & %.3f & %.3f & %.5f $\\pm$ %.6f' %( fitAcc.GetParameter(0), fitAcc.GetParameter(1), fitAcc.GetParameter(2), fitAcc.GetMinimumX(fitRange[0],fitRange[1]), abs(VarAccPlus1Sig) ))
    elif KinVariable == "MTop": 
      TableOutput.write(' \\\\ \n  Acc normalisation & %.3f & %.3f & %.3f & %.3f $\\pm$ %.4f' %( fitAcc.GetParameter(0), fitAcc.GetParameter(1), fitAcc.GetParameter(2), fitAcc.GetMinimumX(fitRange[0],fitRange[1]), abs(VarAccPlus1Sig) ))
    
TableOutput.write('\n \\end{tabular} \n\\end{table} \n')
print "Output can be found in the following directory : ",whichDir
