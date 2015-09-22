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
if len(sys.argv) <= 2:
  print "Need to give the directory of interest, the desired pT-cut and the value of the configuration in command line !"
  print "--> Syntax example: python RemoveLowPtEvents.py Events/blabla Pos05 whichFile(opt) nrEvts(opt)"
  sys.exit()
whichDir, VarConfig = sys.argv[1], sys.argv[2]

# Optional input parameter!
NrEvts = 0
NrWeightsGiven = False
if len(sys.argv) >= 4:
  NrWeightsGiven = True
  NrEvts = int(sys.argv[3])

WeightsFile = ""
WeightFileGiven = False
if len(sys.argv) >= 5:
  WeightFileGiven = True
  weightsFileName = sys.argv[4]
  if str(weightsFileName).endswith("NoUncompleteEvts.out") or str(weightsFileName).endswith("NoZero.out") or str(weightsFileName).endswith("weights.out"):
    WeightsFile = open(weightsFileName, 'r')

# ----------------------------------------#
#   Get the correct input/output files    #
# ----------------------------------------#

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

WeightsCorr = open(WeightsFile.name[:-4]+'_CosTheta.out', 'w')

# --------------------------------------------- #
#   Step 1: Store the considered event numbers  #
# --------------------------------------------- #
EventNr = []
highestEvtNr = 0

for weightLine in WeightsFile:
  if weightLine.split()[0] != "#" and len(EventNr) < NrEvts:
    WeightsCorr.write(weightLine)
    if not int(weightLine.split()[0]) in EventNr:
      EventNr.append(int(weightLine.split()[0]))
      if highestEvtNr < int(weightLine.split()[0]):
        highestEvtNr = int(weightLine.split()[0])
  elif weightLine.split()[0] == "#":
    WeightsCorr.write(weightLine)

WeightsCorr.close()
WeightsFile.close()
NrEvts = len(EventNr)
print "--> Highest event number is : ", highestEvtNr
print "--> Length of NrEvts : ", NrEvts

# ------------------------------------------------- #
#   Step 2: Get the corresponding cos theta weight  #
# ------------------------------------------------- #
CosThetaCorrFile = open(os.path.join('/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RgR_SemiMu_CosTheta/SampleAnalyzer/CosThetaReweighting/'+str(VarConfig)+'_SemiMu/CosThetaWeight_'+str(VarConfig)+'_SemiMu_PtCutsApplied_AlsoOnMET.txt'),'r')
print "Getting cos theta info from : ", CosThetaCorrFile.name

CosThetaCorr, CosThetaSum = [], 0
# Initialize the correct number of entries for the CosThetaCorr array:
for ii in range(highestEvtNr+1):
  CosThetaCorr.append(0)

for corrLine in CosThetaCorrFile:
  if int(corrLine.split()[0]) in EventNr:
    CosThetaCorr[int(corrLine.split()[0])] = float(corrLine.split()[1])
    CosThetaSum += float(corrLine.split()[1])
CosThetaCorrFile.close()

# Now loop over all the stored cos theta weights and provide normalisation equal to 1 !!
Norm = NrEvts / CosThetaSum
print "Norm is = ", Norm, "with NrEvents and total sum = ", NrEvts, " and ", CosThetaSum
# for cosWeight in range(highestEvtNr):
#  CosThetaCorr[cosWeight] *= Norm

# ------------------------------------------------------- #
#   Step 3: Add the cos theta weight in the weights file  #
# ------------------------------------------------------- #
CosThetaWeights = open(WeightsFile.name[:-4]+'_CosTheta.out', 'r')
newWeights = open('dummy.out', 'w')
for line in CosThetaWeights:
  word = line.split()
  if str(word[0]) != "#":
    newWeights.write(line.replace(word[4], word[4] + '     ' + str(CosThetaCorr[int(word[0])])))
  else: newWeights.write(line)   # Write the first two lines containing text!

os.system('mv ' + newWeights.name + ' ' + CosThetaWeights.name)
CosThetaWeights.close()
newWeights.close()
