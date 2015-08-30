#! python
import os

WeightsFileOne = open('weights_FirstPart.out','r')
outfile = file('weights_Combined.out','w')

#Fill the weights of the second file in a separate double array
NrConfigs = 9
NrEvents = 10000
SecondPart = [["empty" for x in range(NrConfigs)] for x in range(NrEvents)] 
WeightsFileTwo = open('weights_SecondPart.out','r')
for secondLine in WeightsFileTwo:
  secondWord = secondLine.split()
  if str(secondWord[0]) != "#":
    SecondPart[int(secondWord[0])-1][int(secondWord[1])-1] = secondLine
  
ConsideredEvent = 0
for weightLine in WeightsFileOne:
  weightWord = weightLine.split()
  if str(weightWord[0]) != "#":
    if int(weightWord[1]) == 1:
      #Stuff about the second file, before the eventNr is changed
      if ConsideredEvent != 0:
        for ii in range(NrConfigs):
          if str(SecondPart[ConsideredEvent-1][ii]) != "empty":
            outfile.write(SecondPart[ConsideredEvent-1][ii])
      #Now update the eventNr and go on with the next one!
      ConsideredEvent = int(weightWord[0])
    #Write away all the lines of the first file!
    if int(weightWord[0]) == ConsideredEvent:
      outfile.write(weightLine)
  #Write away the two start lines
  else:
    outfile.write(weightLine)

#Still need to write away the last event:
for ii in range(NrConfigs):
  if str(SecondPart[ConsideredEvent-1][ii]) != "empty":
    outfile.write(SecondPart[ConsideredEvent-1][ii])
