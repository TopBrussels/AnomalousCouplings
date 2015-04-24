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

#-Identify kinematic variable considered (mtop or RVR)
if KinVariable == "RVR":
  KinVar = "$Re(V_R)$"
  ValuesToKeep = [-0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5]

  #Select which window of RVR values was considered!
  VarWindow = raw_input('** Choose the correct RVR-window corresponding to the studied file ** \n** Different options are : \n  1) Wide   : [-0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.5] \n  2) Narrow : [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3] \n --> Choose the correct number : ')
  
  if CutApplied == False:
    if VarWindow == "1":   ValuesToKeep.pop(6), ValuesToKeep.pop(4)
    elif VarWindow == "2": ValuesToKeep.pop(10), ValuesToKeep.pop(0)
  elif CutApplied == True:
    if VarWindow == "1":   ValuesToKeep.pop(10), ValuesToKeep.pop(9), ValuesToKeep.pop(6), ValuesToKeep.pop(4), ValuesToKeep.pop(1), ValuesToKeep.pop(0)
    elif VarWindow == "2": ValuesToKeep.pop(10), ValuesToKeep.pop(9), ValuesToKeep.pop(8), ValuesToKeep.pop(2), ValuesToKeep.pop(1), ValuesToKeep.pop(0)

elif KinVariable == "MTop":
  KinVar = "$m_{top}$"
  ValuesToKeep = [153, 163, 171, 172, 173, 174, 175, 183, 193]
  
  if CutApplied == True: ValuesToKeep.pop(8), ValuesToKeep.pop(7), ValuesToKeep.pop(1), ValuesToKeep.pop(0)

print "Values to keep array = ",ValuesToKeep," which has length = ",len(ValuesToKeep)
fitRange = [ValuesToKeep[0], ValuesToKeep[len(ValuesToKeep)-1]]
print " Considered fit range is : ", fitRange 

if CutApplied == False:
  RejectLine = '  if (reject && ('
  for ii in range(len(ValuesToKeep)):
    if ii == 0:                     RejectLine += ' x[0] != '+str(ValuesToKeep[ii])
    elif ii == len(ValuesToKeep)-1: RejectLine += ' && x[0] != '+str(ValuesToKeep[ii])+') ){ \n'
    else:                           RejectLine += ' && x[0] != '+str(ValuesToKeep[ii])
  print "Correct Reject line = ",RejectLine

#--Replace ROOT file information in analyzer
if CutApplied == False: RootAnalyzer = open('fitExcludeEmptyBins.C','r')
else:                   RootAnalyzer = open('fitExcludeEmptyBins_AfterCut.C','r')

#  TH1F* LL    = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPosScdDerInner");
#  TH1F* LLXS  = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLXSPosScdDerInner");
#  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPosScdDerInner");

NewRootAnalyzer = open('output','w')
for RootLine in RootAnalyzer:
  RootWord = RootLine.split()
  if re.search( r".root", RootLine):
    if RootWord[1] == "InputFile":
      NewRootAnalyzer.write(RootLine.replace(RootWord[5],RootFile))
    if RootWord[1] == "OutputFile":
      NewRootAnalyzer.write(RootLine.replace(RootWord[5],whichDir+'LimitedFitResult_'+title+'.root'))
  elif re.search( r"reject &&",RootLine):
    NewRootAnalyzer.write(RejectLine)
  elif re.search( r"TF1 ",RootLine):
    if CutApplied == False:
      if re.search(   r"XS",RootLine):  NewRootAnalyzer.write('  TF1 *funcXS = new TF1("funcXS",fpol,'+str(fitRange[0])+','+str(fitRange[1])+',3); \n')
      elif re.search( r"Acc",RootLine): NewRootAnalyzer.write('  TF1 *funcAcc = new TF1("funcAcc",fpol,'+str(fitRange[0])+','+str(fitRange[1])+',3); \n')
      else:                             NewRootAnalyzer.write('  TF1 *func = new TF1("func",fpol,'+str(fitRange[0])+','+str(fitRange[1])+',3); \n')
    elif CutApplied == True:
      if re.search(   r"XS",RootLine):  NewRootAnalyzer.write('  TF1 *funcXS = new TF1("funcXS","pol2",'+str(fitRange[0])+','+str(fitRange[1])+'); \n')
      elif re.search( r"Acc",RootLine): NewRootAnalyzer.write('  TF1 *funcAcc = new TF1("funcAcc","pol2",'+str(fitRange[0])+','+str(fitRange[1])+'); \n')
      else:                             NewRootAnalyzer.write('  TF1 *func = new TF1("func","pol2",'+str(fitRange[0])+','+str(fitRange[1])+'); \n')
  elif re.search( r" bool GenLevel",RootLine):
    NewRootAnalyzer.write('  bool GenLevel = '+GenLevel+'; \n')
  elif re.search( r"LikelihoodAfterCuts",RootLine) and CutApplied == True:
      if re.search(   r"XS",RootLine):  NewRootAnalyzer.write('  TH1F* LLXS  = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLXSPos'+AppliedCut+'"); \n')
      elif re.search( r"Acc",RootLine): 
        if GenLevel == "false":  NewRootAnalyzer.write('  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLAccPos'+AppliedCut+'"); \n')
        elif GenLevel == "true": NewRootAnalyzer.write('  TH1F* LLAcc = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPos'+AppliedCut+'"); \n')         #Acc doesn't exist for Gen!
      else:                             NewRootAnalyzer.write('  TH1F* LL    = (TH1F*) InputFile->Get("LikelihoodAfterCuts/SignSecondDerivative/LLPos'+AppliedCut+'"); \n')
  else:
    NewRootAnalyzer.write(RootLine)
RootAnalyzer.close()
NewRootAnalyzer.close()

if CutApplied == False:
  os.rename('output','fitExcludeEmptyBins.C'), os.system("root -l -b fitExcludeEmptyBins.C")
elif CutApplied == True:
  os.rename('output','fitExcludeEmptyBins_AfterCut.C'), os.system("root -l -b fitExcludeEmptyBins_AfterCut.C")

#-- Now open the created ROOT file and place the minimum (with uncertainty) in a table! --#
FitFile = TFile(os.path.join(whichDir+"LimitedFitResult_"+title+".root"),'r')
TableOutput = open(os.path.join(whichDir+'FitOutput_'+title+'.tex'),'w')

TableOutput.write('\\begin{table}[h!t] \n \\centering \n \\caption{Fit parameters of 2nd degree polynomial ($a_{0} + a_{1}*x + a_{2}*x^{2}$) and corresponding minimum for '+title+' events.} \\label{table::} \n \\begin{tabular}{c|c|c|c|c} \n')
TableOutput.write('  & $a_{0}$ & $a_{1}$ & $a_{2}$ & '+KinVar+' \\\\ \n  \hline \n')

fit = FitFile.Get("fit_LL")
TableOutput.write('  no normalisation & '+str(fit.GetParameter(0))+' & '+str(fit.GetParameter(1))+' & '+str(fit.GetParameter(2))+' & '+str(fit.GetMinimumX(fitRange[0],fitRange[1]))+' \\\\ \n')

fitXS = FitFile.Get("fit_LLXS")
TableOutput.write('  XS normalisation & '+str(fitXS.GetParameter(0))+' & '+str(fitXS.GetParameter(1))+' & '+str(fitXS.GetParameter(2))+' & '+str(fitXS.GetMinimumX(fitRange[0],fitRange[1])))

if GenLevel == "false":
  fitAcc = FitFile.Get("fit_LLAcc")
  TableOutput.write(' \\\\ \n  Acc normalisation & '+str(fitAcc.GetParameter(0))+' & '+str(fitAcc.GetParameter(1))+' & '+str(fitAcc.GetParameter(2))+' & '+str(fitAcc.GetMinimumX(fitRange[0],fitRange[1]))+' ')

TableOutput.write('\n \\end{tabular} \n\\end{table} \n')

print "Output can be found in the following directory : ",whichDir
