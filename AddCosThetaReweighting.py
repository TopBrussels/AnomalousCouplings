############################################
#                                          #
#  RemoveLowPtEvents.py script             #
#  --> Has a two-fold goal:                #
#      1) Removing events below pt-value   #
#      2) Include cos theta weight in file #
#                                          #
############################################
# ! python
import os
import sys

# --------------------------------#
#    Get input from command line  #
# --------------------------------#
if len(sys.argv) <= 4:
  print "Need to give the directory of interest, the anomalous couplings coefficient, the desired pT-cut and the value of the configuration in command line !"
  print "--> Syntax example: python RemoveLowPtEvents.py Events/blabla RVR/RgR Pos05 whichFile(opt) nrEvts(opt)"
  sys.exit()
whichDir, AnomCoef, VarConfig = sys.argv[1], sys.argv[2], sys.argv[3]
print "Stored values are : ", whichDir, ", ", AnomCoef, " and ", VarConfig

# Optional input parameter!
WeightsFile = ""
WeightFileGiven = False
if len(sys.argv) >= 5:
  WeightFileGiven = True
  weightsFileName = sys.argv[4]
  if str(weightsFileName).endswith("NoUncompleteEvts.out") or str(weightsFileName).endswith("NoZero.out") or str(weightsFileName).endswith("weights.out"):
    WeightsFile = open(os.path.join(whichDir+''+weightsFileName), 'r')

NrEvts = 0
NrWeightsGiven = False
if len(sys.argv) >= 6:
  NrWeightsGiven = True
  NrEvts = int(sys.argv[5])

# ----------------------------------------#
#   Get the correct input/output files    #
# ----------------------------------------#
CosThetaCorrFile = open(os.path.join('/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RVRAcceptance/SampleAnalyzer/CosThetaReweighting_'+str(AnomCoef)+'/'+str(VarConfig)+'_'+str(AnomCoef)+'/CosThetaWeight_'+str(VarConfig)+'_'+str(AnomCoef)+'_PtCutsApplied_AlsoOnMET.txt'),'r')
print "Getting cos theta info from : ", CosThetaCorrFile.name

if not WeightFileGiven:
  list_dir = os.listdir(whichDir)
  WeightsFileArray, weightsFileCounter = [], 0
  for file_dir in list_dir:
    if (file_dir.endswith(".out") or file_dir.endswith("NoUncompleteEvts.out") or file_dir.endswith("NoZero.out") ) and file_dir.startswith("weights"):  # eg: '.txt'
      weightsFileCounter += 1
      WeightsFileArray.append(file_dir)

  if int(weightsFileCounter) == 0:
    print "No weights file found in this directory !"
    sys.exit()
  elif int(weightsFileCounter) == 1: WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
  elif int(weightsFileCounter) > 1:
    for ii in range(len(WeightsFileArray)):
      print " ", ii, " ) ",WeightsFileArray[ii]
    fileNr = raw_input('Choose the number of the file of interest! : ')
    WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]), 'r')

print "Will be using file : ", WeightsFile.name, " as starting point!"

WeightsCorr = open(WeightsFile.name[:-4]+'_PtCutsApplied_AlsoOnMET_ApplyCosThetaReweighting.out', 'w')

# -----------------------------------------------------------------------#
#  Get the number of events and configurations from the output.xml file  #
# -----------------------------------------------------------------------#

if not NrWeightsGiven:
  XMLFile = open(os.path.join(whichDir+'output.xml'),'r')
  for xmlLine in XMLFile:
    xmlWord = xmlLine.split()
    if len(xmlWord) > 2 and xmlWord[0] == "nb_exp_events": NrEvts = int(xmlWord[1])
  XMLFile.close()
print "Number of events is : ", NrEvts

# ------------------------------#
#    Save the cos theta weight  #
# ------------------------------#
# First get the info from the output file of MadAnalysis calculation!
CosThetaCorr, CosThetaSum = [], 0
for corrLine in CosThetaCorrFile:
  if int(corrLine.split()[0]) < NrEvts:
    CosThetaCorr.append(float(corrLine.split()[1]))
    CosThetaSum += float(corrLine.split()[1])
CosThetaCorrFile.close()

# Now loop over all the stored cos theta weights and provide normalisation equal to 1 !!
Norm = NrEvts / CosThetaSum
print "Norm is = ", Norm, "with NrEvents and total sum = ", NrEvts, " and ", CosThetaSum
for cosWeight in range(len(CosThetaCorr)):
  CosThetaCorr[cosWeight] *= Norm

for line in WeightsFile:
  word = line.split()
  if str(word[0]) != "#":
    WeightsCorr.write(line.replace(word[4], str(CosThetaCorr[int(word[0])])+'  '+word[4]))
  else: WeightsCorr.write(line)   # Write the first two lines containing text!

WeightsCorr.close()
WeightsFile.close()
