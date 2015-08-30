#! python
import os
import sys

#Directory of interest:
if len(sys.argv) <= 1:
  print "Need to give the directory of interest in command line !"
  sys.exit()
whichDir = sys.argv[1] 
print " Interested in directory : ",whichDir

WeightsFileWide = open(os.path.join(whichDir+'weights_Wide.out'),'r')
WeightsFileNarrow = open(os.path.join(whichDir+'weights_Narrow.out'),'r')
outfile = file(os.path.join(whichDir+'weights_Both.out'),'w')
ZeroRVROutput = file(os.path.join(whichDir+'ZeroRVRWeights.out'),'w')

#Fill the weights of the narrow file in a separate double array
NrNarrowConfigs = 3
NrEvents = 10000
NarrowPart = [["empty" for x in range(NrNarrowConfigs)] for x in range(NrEvents)] 
for narrowLine in WeightsFileNarrow:
  narrowWord = narrowLine.split()
  if str(narrowWord[0]) != "#":
    evtCounter = narrowWord[1]     #Need to store the original event number for filling the array!
    #-- Set the event counter correct --#  
    narrowBuiltLine = "\n"
    if narrowWord[1] == "1": narrowWord[1] = "4"
    if narrowWord[1] == "2": 
      narrowWord[1] = "N"
      narrowBuiltLine = ""
    if narrowWord[1] == "3": narrowWord[1] = "6"
    for ii in range(len(narrowWord)):
      narrowBuiltLine = narrowBuiltLine+narrowWord[ii]+" "
    #-- Now store the line --#
    NarrowPart[int(narrowWord[0])-1][int(evtCounter)-1] = narrowBuiltLine

ConsideredEvent = 0
ConsideredConfig = 0
DifferentZeroRVRValue = 0
RVR =       [-0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
DesiredNr = ["0",  "1",  "2",  "3",  "4",  "5", "6",  "7", "8", "9", "0" ]
for wideLine in WeightsFileWide:
  wideWord = wideLine.split()
  if str(wideWord[0]) == "#":
    #Write away the two start lines
    outfile.write(wideLine)
  else:
#  elif int(wideWord[0]) < 11:
    if int(wideWord[1]) == 1:
      ConsideredConfig = 0
      ConsideredEvent+=1
    if RVR[int(ConsideredConfig)] != -0.5 and RVR[int(ConsideredConfig)] != 0.5: 
      builtLine = "\n"
      if ConsideredEvent == 1 and RVR[int(ConsideredConfig)] == -0.3: builtLine = ""
      wideWord[1] = DesiredNr[int(ConsideredConfig)]
      for ii in range(len(wideWord)):
        builtLine = builtLine+wideWord[ii]+" "
      outfile.write(builtLine)
    #- In case the Re(VR) value of -0.1 is considered, also write away the Re(VR) = -0.05!
    if RVR[int(ConsideredConfig)] == -0.1:
      outfile.write(NarrowPart[int(wideWord[0])-1][0])
      ConsideredConfig+=1
    #- In case the Re(VR) value of 0.0 is considered, also write away the Re(VR) = 0.05!
    elif RVR[int(ConsideredConfig)] == 0.0:
      outfile.write(NarrowPart[int(wideWord[0])-1][2])
      if ConsideredEvent == 1: ZeroRVROutput.write(wideLine)
      else: ZeroRVROutput.write("\n"+wideLine)
      ZeroRVROutput.write(NarrowPart[int(wideWord[0])-1][1])
      NarrowZeroWords = NarrowPart[int(wideWord[0])-1][1].split()
      if NarrowZeroWords[3] != wideWord[3]:
        print "MadWeight output not the same !! "
        DifferentZeroRVRValue+=1
      ZeroRVROutput.write("\n-----------------------------------------")
      ConsideredConfig+=1
    ConsideredConfig+=1

print " ****** Nr of events with different Re(VR) = 0.0 weight : ",DifferentZeroRVRValue," (",ConsideredEvent," considered) ******"
