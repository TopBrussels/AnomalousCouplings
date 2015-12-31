# ----------------------------------------------------------- #
#                                                             #
#  This script will add the scaleFactor to the weights file.  #
#  This will be added as an additional variable such that the #
#  choice can still be made within the analyzer whether the   #
#  scale-factors should be applied or not!                    #
#                                                             #
#  Only a special case should be created for the TTJets       #
#  event since these are splitted in Correct, Wrong or        #
#  Unmatched, but everything should be stored in the file     #
#                                                             #
#  --> Ideally the EventNrMatching file should also hold all  #
#      the systematic information such that it can be chosen  #
#      in this script which systematic should be applied!!    #
# ----------------------------------------------------------- #

import os
import sys
import re

Luminosity, NormFactor, NrEvts = 0, 0, 0
EvtNrMatchingArray = []

EvtNrMatching = open('MWEventNrMatching_TTbarJets_SemiLept.txt')
for evtNr in EvtNrMatching:

    evtNrWord = evtNr.split()
    if re.search(r" * NormFactor = ", evtNr):
        NormFactor = evtNrWord[3]
        print "Norm factor is : ", evtNrWord[3]

    # Only want files with at least 6 input variables, otherwise they have not been processed
    if len(evtNrWord) == 7 and str(evtNrWord[4]) == "Correct":
        EvtNrMatchingArray.append(str(evtNrWord[6]))
EvtNrMatching.close()
print "Stored a total of: ", len(EvtNrMatchingArray)

print " "
# -- Now go into the different directories and change the weights.out file -- #
normFactorAdded = False

origWeights = open('weights_CheckedEvts.out', 'r')
newWeights = open('weights_CheckedEvts_SFAdded.out', 'w')

newWeights.write('TTbarJets_Correct ' + NormFactor + '\n')
for line in origWeights:
    word = line.split()

    if str(word[0]) != "#":
        newWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' ' +
                         EvtNrMatchingArray[int(word[0])-1]+'\n')

origWeights.close()
newWeights.close()