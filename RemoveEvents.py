#########################################
#                                       #
#  Script which will remove bad events  #
#    - Incomplete events                #
#    - Zero-weight events               #
#    - Events not surviving cuts (opt)  #
#  Required inputs are:                 #
#    - whichDir                         #
#                                       #
#  --> If no events have to be removed, #
#      the new weights_xx.out file will #
#      be deleted to avoid overload     #
#                                       #
#  *Files which need to be accessed:    #
#    1) Original weights.out file       #
#    2) output.xml for #evts & configs  #
#                                       #
#########################################

import os
import sys

whichDir = sys.argv[1]

YesOptions = ["y", "yes", "Y", "Yes", "YES"]
# In case extra cuts should be applied it should be mentioned explicitly!
applyExtraCuts = False
if len(sys.argv) >= 3:
  if sys.argv[2] in YesOptions:
    print "\n ** Changing applyExtraCuts variable to True! **"
    applyExtraCuts = True

applyCosTheta = False
if len(sys.argv) >= 4:
  if sys.argv[3] in YesOptions:
    print "\n ** Will be applying cos theta reweighting ! **"
    applyCosTheta = True

# -------------------------------------------------------------------------------#
#  Step 1: Get the number of events and configurations from the output.xml file  #
# -------------------------------------------------------------------------------#
NrConfigs, NrEvts = 0, 0
keepCounting = True
XMLFile = open(os.path.join(whichDir+'output.xml'),'r')
for xmlLine in XMLFile:
  xmlWord = xmlLine.split()
  if len(xmlWord) > 1 and xmlWord[0] == "13" and keepCounting == True: NrConfigs += 1
  # Stop counting as soon as the word 'Permutations' is encountered
  # --> Otherwise the ymm and MM terms are counted as configs as well ...
  if len(xmlWord) > 2 and xmlWord[1] == "Permutations": keepCounting = False
  if len(xmlWord) > 2 and xmlWord[0] == "nb_exp_events": NrEvts = int(xmlWord[1])

print "\nNumber of configs is : ",NrConfigs
print "Number of events is : ",NrEvts

# ---------------------------------------------------------------------------#
#  Step 2: Identify the different types of events which have to be removed  #
# ---------------------------------------------------------------------------#
IncomplEventsToDelete, NrConfigsPerEvent, ZeroEventsToDelete, CutEventsToKeep = [], [], [], []
for ii in range(NrEvts): NrConfigsPerEvent.append(0)

# Loop over all lines in original weights.out file, count the number of configs for each event and identify the zero and non-selected events
weightFile = open(os.path.join(whichDir+'weights.out'),'r')
for weightLine in weightFile:
  weightWord = weightLine.split()
  #Only interested in files starting with a number
  if str(weightWord[0]) != "#":
    NrConfigsPerEvent[int(weightWord[0])-1] += 1
    if str(weightWord[3]) == "0.0" and not int(weightWord[0]) in ZeroEventsToDelete:
      ZeroEventsToDelete.append(int(weightWord[0]))
# Now that the number of configs for each event has been calculated, store the events which don't have the correct number!
for ii in range(NrEvts):
  if NrConfigsPerEvent[ii] != NrConfigs:
    IncomplEventsToDelete.append(ii+1)

# The event numbers not surviving the cut can be obtained from the following file: EventNrMatching_TTbarJets_SemiLept.txt
TypeCutEvts = "empty"
if whichDir.find("Correct") >= 0: 
  TypeCutEvts = "Correct"
elif whichDir.find("Wrong") >= 0:
  TypeCutEvts = "Wrong"
elif whichDir.find("Unmatched") >= 0:
  TypeCutEvts = "Unmatched"

if applyExtraCuts:
  if os.path.exists(os.path.join(whichDir+'EventNrMatching_TTbarJets_SemiLept.txt')):
    evtNrFile = open(os.path.join(whichDir+'EventNrMatching_TTbarJets_SemiLept.txt'),'r')
    for evtNrLine in evtNrFile:
      evtNrWord = evtNrLine.split()
      if len(evtNrWord) == 7:   # For lines with this length, output exists!
        if str(evtNrWord[3]) == "SemiMu" and str(evtNrWord[5]) == TypeCutEvts and str(evtNrWord[1]) == "1":
          CutEventsToKeep.append(int(evtNrWord[6]))
    evtNrFile.close()
  else:
    print "\n   ***   ERROR: Can only remove the non-selected events if an EventNrMatching file is given!!    ***"
    print   "   --> Will skip the event-selection removal but still do the removal of incomplete and zero events!"
    applyExtraCuts = False

print "\n --> Number of events with incomplete weights : ", len(IncomplEventsToDelete)
print " --> Number of events with a weight equal to 0 : ", len(ZeroEventsToDelete)
if applyExtraCuts:
  print " --> Number of events surviving the extra cuts : ", len(CutEventsToKeep)

# ------------------------------------------------------------------#
#  Step 3: Save the cos theta weight from the MadAnalysis directory #
#          --> Need to translate the histogram into bin ranges      #
# ------------------------------------------------------------------#
# Only looking at RECO:
AnomCoef, VarConfig = "RgR", "SM"
CosTh_Norm_CheckedEvts, CosTh_Norm_SelectedEvts, CosTh_Norm_NonSelEvts = 0, 0, 0
CosTheta = []

# Make sure that the EventNrMatching file is present!
if not os.path.exists(os.path.join(whichDir+'EventNrMatching_TTbarJets_SemiLept.txt')):
  print "\n   ***   ERROR: Can only add the cos theta information if an EventNrMatching file is given !!    ***"
  print   "   --> Will skip the cosTheta part, but still continue with the removal of incomplete and zero events!"
  applyCosTheta = False

if applyCosTheta:

  MadAnalysisPath = '/user/aolbrech/AnomalousCouplings/MadAnalysis_v112/RVRAcceptance/SampleAnalyzer/CosThetaReweighting_'+str(AnomCoef)+'/'+str(VarConfig)+'_'+str(AnomCoef)
  CosThetaFile = open(os.path.join(MadAnalysisPath+'/HistoWeight_'+str(VarConfig)+'_'+str(AnomCoef)+'_PtCutsApplied_AlsoOnMET.txt'),'r')

  print "\n CosTheta info taken from : ", CosThetaFile.name

  # Now get the info from the output file of MadAnalysis calculation!
  CosThetaWeight, CosThetaLowEdge, CosThetaHighEdge = [], [], []
  for cosThLine in CosThetaFile:
    cosThWord = cosThLine.split()
    if len(cosThWord) == 3:
      CosThetaLowEdge.append(cosThWord[0])
      CosThetaHighEdge.append(cosThWord[1])
      CosThetaWeight.append(cosThWord[2])
  CosThetaFile.close()

  #   Get the reweighting value for each event

  EvtCounter = []    # Double-check which just counts the event index ==> Should go from 1 to ....
  EvtMatchFile = open(os.path.join(whichDir+'EventNrMatching_TTbarJets_SemiLept.txt'),'r')
  print "\n Looking at file : ", EvtMatchFile.name
  for evtMatchLine in EvtMatchFile:
    evtMatchWord = evtMatchLine.split()
    if len(evtMatchWord) == 7:
      if str(evtMatchWord[5]) == "Correct":
        EvtCounter.append(int(evtMatchWord[6]))
        for iBin in range(len(CosThetaWeight)):
          if float(CosThetaLowEdge[iBin]) <= float(evtMatchWord[2]) <= float(CosThetaHighEdge[iBin]):
            if not (float(CosThetaLowEdge[iBin]) != -1 and float(evtMatchWord[2]) == float(CosThetaLowEdge[iBin])):
              CosTheta.append(CosThetaWeight[iBin])
            else:
              print "Value is on edge : ", evtMatchWord[2]
  for it in range(len(EvtCounter)-1):
    if EvtCounter[it]+1 != int(EvtCounter[it+1]):
      print "ERROR: value not +1 ", EvtCounter[it+1]-EvtCounter[it]

  #   Get the different normalisations
  for iEvt in range(len(CosTheta)):
    # print " Value of evtCounter = ", EvtCounter[iEvt]
    if len(IncomplEventsToDelete) != 0 and len(ZeroEventsToDelete) != 0 and not int(EvtCounter[iEvt]) in ZeroEventsToDelete and not int(EvtCounter[iEvt]) in IncomplEventsToDelete:
      CosTh_Norm_CheckedEvts += float(CosTheta[iEvt])
    if applyExtraCuts and len(CutEventsToKeep) != 0 and int(EvtCounter[iEvt]) in CutEventsToKeep:
      CosTh_Norm_SelectedEvts += float(CosTheta[iEvt])
    if applyExtraCuts and len(CutEventsToKeep) != NrEvts and not int(EvtCounter[iEvt]) in CutEventsToKeep:
      CosTh_Norm_NonSelEvts += float(CosTheta[iEvt])

  CosTh_Norm_CheckedEvts = float((NrEvts - len(ZeroEventsToDelete) - len(IncomplEventsToDelete)))/float(CosTh_Norm_CheckedEvts)
  print "Will be applying normalisation : ", CosTh_Norm_CheckedEvts
  if applyExtraCuts:
    CosTh_Norm_SelectedEvts = float(len(CutEventsToKeep))/float(CosTh_Norm_SelectedEvts)
    CosTh_Norm_NonSelEvts = float(NrEvts - len(CutEventsToKeep))/float(CosTh_Norm_NonSelEvts)

# ----------------------------------------------------------------------#
#  Step 4: Store the complete events and non-zero events in a new file  #
# ----------------------------------------------------------------------#

# In case events have to be deleted open a new file, else skip this step!
if len(IncomplEventsToDelete) != 0 or len(ZeroEventsToDelete) != 0:
  NoIncomplEvtsFile = open(os.path.join(whichDir+'weights_CheckedEvts.out'),'w')
  NoIncomplEvtsFile_CosTh = open(os.path.join(whichDir+'weights_CheckedEvts_CosTheta.out'),'w')
  WeightFile = open(weightFile.name)
  for line in WeightFile:
    word = line.split()
    if str(word[0]) != "#":
      if not int(word[0]) in IncomplEventsToDelete and not int(word[0]) in ZeroEventsToDelete:
        NoIncomplEvtsFile.write(line)
        NoIncomplEvtsFile_CosTh.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' ')
        if applyCosTheta:
          NoIncomplEvtsFile_CosTh.write(str(float(CosTheta[int(word[0])-1])*float(CosTh_Norm_CheckedEvts))+'\n')
    else:
      NoIncomplEvtsFile.write(line)
      NoIncomplEvtsFile_CosTh.write(line)
  NoIncomplEvtsFile.close()
  NoIncomplEvtsFile_CosTh.close()
  if not applyCosTheta:
    # os.system('rm '+NoIncomplEvtsFile_CosTh.name)
    os.system('mv '+NoIncomplEvtsFile_CosTh.name+' '+whichDir+'test.out')
elif applyCosTheta:
  print " In elif case (meaning that there are no events which should be deleted, but cosTheta should still be applied)"
  WeightsFile_CosTheta = open(os.path.join(whichDir+'weights_CosTheta.out'),'w')
  WeightFile = open(weightFile.name)
  for line in WeightFile:
    word = line.split()
    if str(word[0]) != "#":
      WeightsFile_CosTheta.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' ')
      WeightsFile_CosTheta.write(str(float(CosTheta[int(word[0])-1])*float(CosTh_Norm_CheckedEvts))+'\n')
    else:
      WeightsFile_CosTheta.write(line)
  WeightsFile_CosTheta.close()

weightFile.close()

# ------------------------------------------------------------------------------------#
#  Step 5: Now remove from this checked file the events not surviving the extra cuts  #
# ------------------------------------------------------------------------------------#

# The choice of weights_xx.out file depends on whether the weights_CheckedEvts.out exists in the directory!
# If yes, this one should be used!
list_dir = os.listdir(whichDir)
allCheckedWeightFile = open(os.path.join(whichDir+'weights.out'),'r')
for file_dir in list_dir:
  if file_dir.endswith(".out") and file_dir.startswith("weights_CheckedEvts"):
    allCheckedWeightFile = open(os.path.join(whichDir+'weights_CheckedEvts.out'),'r')

if len(CutEventsToKeep) != 0 and applyExtraCuts:
  print "Looking at file : ", allCheckedWeightFile.name
  NoCutsEvtsFile = open(os.path.join(whichDir+'weights_ExtraCuts.out'), 'w')
  NoCutsEvtsFile_CosTh = open(os.path.join(whichDir+'weights_ExtraCuts_CosTheta.out'), 'w')
  DeletedEvtsFile = open(os.path.join(whichDir+'weights_DeletedByExtraCuts.out'), 'w')
  DeletedEvtsFile_CosTh = open(os.path.join(whichDir+'weights_DeletedByExtraCuts_CosTheta.out'), 'w')

  for checkedLine in allCheckedWeightFile:
    checkedWord = checkedLine.split()
    # Only interested in files starting with a number
    if str(checkedWord[0]) != "#":
      if int(checkedWord[0]) in CutEventsToKeep:
        NoCutsEvtsFile.write(checkedLine)
        NoCutsEvtsFile_CosTh.write(checkedWord[0]+' '+checkedWord[1]+' '+checkedWord[2]+' '+checkedWord[3]+' '+checkedWord[4]+' ')
        if applyCosTheta:
          NoCutsEvtsFile_CosTh.write(str(float(CosTheta[int(checkedWord[0])-1])*float(CosTh_Norm_CheckedEvts))+'\n')
      else:
        DeletedEvtsFile.write(checkedLine)
        DeletedEvtsFile_CosTh.write(checkedWord[0]+' '+checkedWord[1]+' '+checkedWord[2]+' '+checkedWord[3]+' '+checkedWord[4]+' ')
        if applyCosTheta:
          DeletedEvtsFile_CosTh.write(str(float(CosTheta[int(checkedWord[0])-1])*float(CosTh_Norm_CheckedEvts))+'\n')
    else:
      NoCutsEvtsFile.write(checkedLine)
      NoCutsEvtsFile_CosTh.write(checkedLine)
      DeletedEvtsFile.write(checkedLine)
      DeletedEvtsFile_CosTh.write(checkedLine)
  print "Done, closing all files ! "
  NoCutsEvtsFile.close(), NoCutsEvtsFile_CosTh.close()
  DeletedEvtsFile.close(), DeletedEvtsFile_CosTh.close()
  if not applyCosTheta:
    os.system('mv '+NoCutsEvtsFile_CosTh.name+' '+whichDir+'test_Keep.out')
    os.system('mv '+DeletedEvtsFile_CosTh.name+' '+whichDir+'test_Delete.out')
allCheckedWeightFile.close()
