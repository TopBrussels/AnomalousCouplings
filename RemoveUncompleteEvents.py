#! python
import os
import sys

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
if len(sys.argv) == 2:
  print "Need to give the number of configurations and number of events which were considered!"
  sys.exit()
whichDir = sys.argv[1] 
NrConfigs = int(sys.argv[2])
NrEvents = int(sys.argv[3])

print " Interested in directory : ",whichDir
print " Looking at ",NrConfigs," configurations"

#Output file of interest:
outfile = file(os.path.join(whichDir+'weights_NoUncompleteEvts.out'),'w')
print "Using output file : ",outfile

#Input file of interest:
list_dir = []
list_dir = os.listdir(whichDir)
WeightsFileArray = []
weightsFileCounter = 0
for file in list_dir:
  if file.endswith(".out") and file.startswith("weights"): # eg: '.txt'
    weightsFileCounter += 1
    WeightsFileArray.append(file)

if int(weightsFileCounter) == 1:
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
  secondWeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[0]),'r')
elif int(weightsFileCounter) == 0:
  print "No weights file found in this directory !"
  sys.exit()
elif int(weightsFileCounter) > 1:
  for ii in range(len(WeightsFileArray)):
    print " ",ii," ) ",WeightsFileArray[ii]
  fileNr = raw_input('Choose the number of the file of interest! : ')
  WeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')
  secondWeightsFile = open(os.path.join(whichDir+''+WeightsFileArray[int(fileNr)]),'r')

print "Will be using file : ",WeightsFile
print "Will be using second file : ",secondWeightsFile

EventsToDelete = []
NrConfigsPerEvent = []

for ii in range(NrEvents):
  NrConfigsPerEvent.append(0)

print "Start looping over all lines !"
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#":
    NrConfigsPerEvent[int(WeightWord[0])-1] = NrConfigsPerEvent[int(WeightWord[0])-1]+1

for ii in range(NrEvents):
  if NrConfigsPerEvent[ii] != NrConfigs:
    EventsToDelete.append(ii+1)
print "Events with uncomplete weights :",EventsToDelete," which has length = ",len(EventsToDelete)

for line in secondWeightsFile:
  word = line.split()
  if str(word[0]) != "#":
    if not int(word[0]) in EventsToDelete:
      outfile.write(line)
  else:
    outfile.write(line)
