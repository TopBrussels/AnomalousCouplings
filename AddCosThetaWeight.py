############################################
#                                          #
#  AddCosThetaWeight.py script             #
#  --> Include cos theta weight in file    #
#                                          #
############################################
# ! python
import os
import sys

# --------------------------------#
#    Get input from command line  #
# --------------------------------#
if len(sys.argv) < 4:
  print "Need to give the directory of interest, the anomalous couplings coefficient and its value in command line !"
  print "Optionally the weights file can be given! "
  print "--> Syntax example: python AddCosThetaReweighting.py Events/blabla RVR/RgR Pos05 whichFile(opt)"
  sys.exit()
whichDir, AnomCoef, VarConfig = sys.argv[1], sys.argv[2], sys.argv[3]
print "Stored values are : ", whichDir, ", ", AnomCoef, " and ", VarConfig,'\n'

# Optional input parameter!
WeightsFile = ""
WeightFileGiven = False
if len(sys.argv) >= 5:
  WeightFileGiven = True
  weightsFileName = sys.argv[4]
  if str(weightsFileName).endswith(".out") and str(weightsFileName).startswith("weights"):
    WeightsFile = open(os.path.join(whichDir+''+weightsFileName), 'r')

# -----------------------------------------------------------#
#   Save the cos theta weight from the MadAnalysis directory #
#   --> Need to translate the histogram into bin ranges      #
# -----------------------------------------------------------#
MadAnalysisPath = '/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RVRAcceptance/SampleAnalyzer/CosThetaReweighting_'+str(AnomCoef)+'/'+str(VarConfig)+'_'+str(AnomCoef)
CosThetaFile = open(os.path.join(MadAnalysisPath+'/HistoWeight_'+str(VarConfig)+'_'+str(AnomCoef)+'_PtCutsApplied_AlsoOnMET.txt'),'r')
print "Getting cos theta info from : ", CosThetaFile.name

# Now get the info from the output file of MadAnalysis calculation!
CosThetaWeight, CosThetaLowEdge, CosThetaHighEdge = [], [], []
for cosThLine in CosThetaFile:
  cosThWord = cosThLine.split()
  if len(cosThWord) == 3:
    CosThetaLowEdge.append(cosThWord[0])
    CosThetaHighEdge.append(cosThWord[1])
    CosThetaWeight.append(cosThWord[2])
CosThetaFile.close()

# ----------------------------------------#
#   Get the correct input/output files    #
# ----------------------------------------#
if not WeightFileGiven:
  list_dir = os.listdir(whichDir)
  WeightsFileArray, weightsFileCounter = [], 0
  for file_dir in list_dir:
    if file_dir.endswith(".out") and file_dir.startswith("weights"):  # eg: '.txt'
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

# --------------------------------------------#
#   Get the reweighting value for each event  #
# --------------------------------------------#
CosThetaReweighting = []
EvtCounter = []    # Double-check which just counts the event index ==> Should go from 1 to ....
EvtMatchFile = open(os.path.join(whichDir+'EventNrMatching_TTbarJets_SemiLept.txt'),'r')
print "\n Looking at file : ", EvtMatchFile.name
for evtMatchLine in EvtMatchFile:
  evtMatchWord = evtMatchLine.split()
  if len(evtMatchWord) == 7:
    if str(evtMatchWord[5]) == "Correct":    # and int(evtMatchWord[6]) < nEvts:
      EvtCounter.append(int(evtMatchWord[6]))
      for iBin in range(len(CosThetaWeight)):
        if float(evtMatchWord[2]) >= float(CosThetaLowEdge[iBin]) and float(evtMatchWord[2]) <= float(CosThetaHighEdge[iBin]):
          if not (float(CosThetaLowEdge[iBin]) != -1 and float(evtMatchWord[2]) == float(CosThetaLowEdge[iBin])):
            CosThetaReweighting.append(CosThetaWeight[iBin])
          else:
            print "Value is on edge : ", evtMatchWord[2]
print "Saved event number is : ---> Nr of saved is : ", len(EvtCounter)
for it in range(len(EvtCounter)-1):
  if EvtCounter[it]+1 != int(EvtCounter[it+1]):
    print "ERROR: value not +1 ", EvtCounter[it+1]-EvtCounter[it]

print "Will be using file : ", WeightsFile.name, " as starting point! \n"

# -----------------------------------------------------#
#  Now loop over the weights file and store the cosTh  #
# -----------------------------------------------------#
WeightsCorr = open(WeightsFile.name[:-4]+'_PtCutsApplied_AlsoOnMET_CosThetaWeight.out', 'w')

for weightLine in WeightsFile:
  weightWord = weightLine.split()
  if weightWord[0] != "#":
    # print "Looking at event : ", weightLine.split()[0]-1
    WeightsCorr.write(weightWord[0]+'  '+weightWord[1]+'   '+weightWord[2]+'  '+weightWord[3]+'    '+weightWord[4]+'   '+CosThetaReweighting[int(weightLine.split()[0])-1]+'\n')
  else:
    WeightsCorr.write(weightLine)
WeightsCorr.close()