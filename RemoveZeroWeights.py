#! python
import os
import sys

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

#Output file of interest:
outfile = file(os.path.join(whichDir+'weights_NoZero.out'),'w')

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
#Loop over all lines in weights file:
for WeightLine in WeightsFile:
  WeightWord = WeightLine.split()
  #Only interested in files starting with a number
  if str(WeightWord[0]) != "#":
    if str(WeightWord[3]) == "0.0":
      if not int(WeightWord[0]) in EventsToDelete:
       EventsToDelete.append(int(WeightWord[0]))

print "Events with a weight = 0 :",EventsToDelete," which has length = ",len(EventsToDelete)

for line in secondWeightsFile:
  word = line.split()
  if str(word[0]) != "#":
    if not int(word[0]) in EventsToDelete:
      outfile.write(line)
  else:
    outfile.write(line)
#    
#  else:
#    outfile.write(WeightLine) 
