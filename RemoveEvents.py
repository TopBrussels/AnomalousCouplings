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

# In case extra cuts should be applied it should be mentioned explicitly!
applyExtraCuts = False
if len(sys.argv) >= 3:
  print "\n ** Changing applyExtraCuts variable to True! **\n"
  applyExtraCuts = True

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

print "Number of configs is : ",NrConfigs
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

# ----------------------------------------------------------------------#
#  Step 3: Store the complete events and non-zero events in a new file  #
# ----------------------------------------------------------------------#

# In case events have to be deleted open a new file, else skip this step!
if len(IncomplEventsToDelete) != 0 or len(ZeroEventsToDelete) != 0:
  NoIncomplEvtsFile = open(os.path.join(whichDir+'weights_CheckedEvts.out'),'w')
  WeightFile = open(weightFile.name)
  for line in WeightFile:
    word = line.split()
    if str(word[0]) != "#":
      if not int(word[0]) in IncomplEventsToDelete and not int(word[0]) in ZeroEventsToDelete:
        NoIncomplEvtsFile.write(line)
    else:
      NoIncomplEvtsFile.write(line)
  NoIncomplEvtsFile.close()
weightFile.close()

# ------------------------------------------------------------------------------------#
#  Step 3: Now remove from this checked file the events not surviving the extra cuts  #
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
  NoCutsEvtsFile = open(os.path.join(whichDir+'weights_ExtraCuts.out'),'w')
  DeletedEvtsFile = open(os.path.join(whichDir+'weights_DeletedByExtraCuts.out'),'w')
  for checkedLine in allCheckedWeightFile:
    checkedWord = checkedLine.split()
    # Only interested in files starting with a number
    if str(checkedWord[0]) != "#":
      if int(checkedWord[0]) in CutEventsToKeep:
        NoCutsEvtsFile.write(checkedLine)
      else:
        DeletedEvtsFile.write(checkedLine)
    else:
      NoCutsEvtsFile.write(checkedLine)
      DeletedEvtsFile.write(checkedLine)
  print "Done, closing all files ! "
  NoCutsEvtsFile.close()
  DeletedEvtsFile.close()
allCheckedWeightFile.close()
