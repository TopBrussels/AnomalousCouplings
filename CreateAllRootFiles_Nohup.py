##################################################################
#                                                                #
#  Script will do the fitting for the given Events_xx directory  #
#    --> Will loop over all directories existing within this     #
#        given directory                                         #
#    --> Will first even perform the check for uncomplete and    #
#        zero-weight events to be rejected                       #
#    --> Automatically chooses the correct weights_xx.out file   #
#                                                                #
#  !Only important parameter which should be given AND should be #
#  checked in the FitDeviationScript.py file is the range !!     #
#                                                                #
##################################################################

import os
import sys

# Get the Events directory, the considered range and the correction from the input line
if len(sys.argv) < 10:
  print " Need to give the Events directory, the considered range and the correction in the input line"
  print " Correct input : python CreateAllRootFiles_Nohup.py Events_xx RVR/RgR/MTop #evts y/n y/n doublePolFitMacro.C ApplyAccNorm(y/n) ApplyCosTheta(y/n) Range dirString(opt)"
  sys.exit()

whichDir = sys.argv[1]
KinVariable = sys.argv[2]
nEvts = sys.argv[3]
runFitMacro = sys.argv[4]
makeSplitPlots = sys.argv[5]
whichAnalysis = sys.argv[6]
applyAccNorm = sys.argv[7]
applyCosTheta = sys.argv[8]
whichRange = sys.argv[9]

# Can give an optional string which should be retrieved in the considered directories (such as NoCuts/AllCuts/...)
dirString = ""
if len(sys.argv) == 11:
  dirString = sys.argv[10]

list_dir = os.listdir(whichDir)
for VarDir in list_dir:
  print "Can the given string ",dirString," be retrieved in the dirName ",VarDir," .... ? ",VarDir.find(dirString)

  # First make sure whether the uncomplete and zero events are removed!
  # --> Need to make sure that considered VarDir is definitely a directory!
  if not VarDir.endswith(".pdf") and VarDir.find(dirString) > 0:
    os.system('python RemoveEvents.py '+whichDir+''+VarDir+'/')

    # Get the correct weights.out file (_NoZero > _NoUncompleteEvts > weights.out)
    list_dir = os.listdir(whichDir+''+VarDir)
    WeightFileName = ""
    for file_dir in list_dir:
      if os.path.exists(whichDir+''+VarDir+'/weights_NoZero.out'): WeightFileName = "weights_NoZero.out"
      elif os.path.exists(whichDir+''+VarDir+'/weights_NoUncompleteEvts.out'): WeightFileName = "weights_NoUncompleteEvts.out"
      elif os.path.exists(whichDir+''+VarDir+'/weights.out'): WeightFileName = "weights.out"

    if WeightFileName != "":
      # In case the cosTheta reweighting should be applied, create the weights_xx file and change the WeightFileName
      if applyCosTheta == "y" or applyCosTheta == "yes" or applyCosTheta == "Y":
        os.system('python AddCosThetaReweighting.py '+whichDir+''+VarDir+'/ '+KinVariable+' '+VarDir[VarDir.find("MGSample")+8:VarDir.find("_",VarDir.find("MGSample"))]+' '+WeightFileName+' '+nEvts)
        WeightFileName = WeightFileName[:-4]+'_PtCutsApplied_AlsoOnMET_ApplyCosThetaReweighting.out'

      os.system('python FitDeviationScript.py '+str(whichDir)+''+str(VarDir)+'/ '+str(KinVariable)+' '+str(nEvts)+' '+str(runFitMacro)+' '+str(makeSplitPlots)+' '+str(whichAnalysis)+' '+str(applyAccNorm)+' '+str(applyCosTheta)+' '+str(whichDir)+''+str(VarDir)+'/'+str(WeightFileName)+' '+str(whichRange) )
      print "Python command is : ",str('python FitDeviationScript.py '+str(whichDir)+''+str(VarDir)+'/ '+str(KinVariable)+' '+str(nEvts)+' '+str(runFitMacro)+' '+str(makeSplitPlots)+' '+str(whichAnalysis)+' '+str(applyAccNorm)+' '+str(applyCosTheta)+' '+str(whichDir)+''+str(VarDir)+'/'+str(WeightFileName)+' '+str(whichRange))
