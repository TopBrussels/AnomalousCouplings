#########################################
#                                       #
#  Script which will remove bad events  #
#    - Uncomplete events                #
#    - Zero-weight events               #
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
#  !Add option to perform both removes  #
#   separately (as was the case before  #
#   with two separate files)            #
#                                       #
#########################################

#! python
import os
import sys

whichDir = sys.argv[1]

#--------------------------------------------------------------------------------# 
#  Step 1: Get the number of events and configurations from the output.xml file  #
#--------------------------------------------------------------------------------# 
NrConfigs, NrEvts = 0, 0
keepCounting = True
XMLFile = open(os.path.join(whichDir+'output.xml'),'r')
for xmlLine in XMLFile:
  xmlWord = xmlLine.split()
  if len(xmlWord) > 1 and xmlWord[0] == "13" and keepCounting == True: NrConfigs += 1
  #Stop counting as soon as the word 'Permutations' is encountered
  #--> Otherwise the ymm and MM terms are counted as configs as well ...
  if len(xmlWord) > 2 and xmlWord[1] == "Permutations": keepCounting = False
  if len(xmlWord) > 2 and xmlWord[0] == "nb_exp_events": NrEvts = int(xmlWord[1])

print "Number of configs is : ",NrConfigs
print "Number of events is : ",NrEvts

#----------------------------------------# 
#  Step 2: Remove the uncomplete events  #
#----------------------------------------# 
UncomplEventsToDelete, NrConfigsPerEvent = [], []
for ii in range(NrEvts): NrConfigsPerEvent.append(0)

#Loop over all lines in file and count the number of configs for each event
weightFile = open(os.path.join(whichDir+'weights.out'),'r')
for weightLine in weightFile:
  weightWord = weightLine.split()
  #Only interested in files starting with a number
  if str(weightWord[0]) != "#": NrConfigsPerEvent[int(weightWord[0])-1] += 1

for ii in range(NrEvts):
  if NrConfigsPerEvent[ii] != NrConfigs: UncomplEventsToDelete.append(ii+1)
print "Events with uncomplete weights has length = ",len(UncomplEventsToDelete)
 
#In case events have to be deleted open a new file, else skip this step!
if len(UncomplEventsToDelete) != 0:
  NoUncomplEvtsFile = open(os.path.join(whichDir+'weights_NoUncompleteEvts.out'),'w')
  WeightFile = open(weightFile.name)
  for line in WeightFile:
    word = line.split()
    if str(word[0]) != "#":
      if not int(word[0]) in UncomplEventsToDelete: NoUncomplEvtsFile.write(line)
    else: NoUncomplEvtsFile.write(line)
  NoUncomplEvtsFile.close()
weightFile.close()

#--------------------------------------------------------#
#  Step 3: Check whether events with zero-weights exist  #
#--------------------------------------------------------#
ZeroEventsToDelete = []

#The choice of weights_xx.out file depends on whether the weights_NoUncompleteEvts.out exists in the directory!
#If yes, this one should be used!
list_dir = os.listdir(whichDir)
uncomplWeightFile = open(os.path.join(whichDir+'weights.out'),'r')
for file_dir in list_dir:
  if file_dir.endswith(".out") and file_dir.startswith("weights_NoUncompleteEvts"): uncomplWeightFile = open(os.path.join(whichDir+'weights_NoUncompleteEvts.out'),'r') 

print "Looking at file : ", uncomplWeightFile.name

for uncomplWeightLine in uncomplWeightFile:
  uncomplWeightWord = uncomplWeightLine.split()
  #Only interested in files starting with a number
  if str(uncomplWeightWord[0]) != "#":
    if str(uncomplWeightWord[3]) == "0.0":
      if not int(uncomplWeightWord[0]) in ZeroEventsToDelete: ZeroEventsToDelete.append(int(uncomplWeightWord[0]))
uncomplWeightFile.close()
print "Events with a weight equal to 0 has length = ",len(ZeroEventsToDelete)
print "Events to delete are : ", ZeroEventsToDelete

if len(ZeroEventsToDelete) != 0:
  NoZeroEvtsFile = open(os.path.join(whichDir+'weights_NoZero.out'),'w')
  UncomplWeightFile = open(uncomplWeightFile.name,'r')
  for zeroLine in UncomplWeightFile:
    zeroWord = zeroLine.split()
    if str(zeroWord[0]) != "#":
      if not int(zeroWord[0]) in ZeroEventsToDelete:
        NoZeroEvtsFile.write(zeroLine)
    else: NoZeroEvtsFile.write(zeroLine)
  NoZeroEvtsFile.close()
uncomplWeightFile.close()
