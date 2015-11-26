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

# Values to give as input:
#  1) Which directory
#  2) Which systematic (still to do)

if len(sys.argv) < 3:
    print "--> Should at least give the name of the directory and the sample!!"
    print "    Possible options are: 'TTbarJets_FullHadr', 'TTbarJets_FullLept', 'TTbarJets_SemiLept' "
    print "                          'SingleTop_tWChannel_t', 'SingleTop_tWChannel_tbar' "
    exit()
whichDir = sys.argv[1]
whichSample = sys.argv[2]

# -- First run over the EvtNrMatching file and store the scale factors in an array -- #
EvtNrMatching = open('EventNrMatching_'+whichSample+'.txt', 'r')
EvtNrMatchingArray = []
Luminosity = ""
NormFactor = ""
if str(whichSample) != "TTbarJets_SemiLept":
    for evtNr in EvtNrMatching:
        evtNrWord = evtNr.split()

        #First get the lumi and the normFactor!
        if re.search(r" * Lumi = ", evtNr):
            Luminosity = evtNrWord[3]
            print "Luminosity is : ", Luminosity
        if re.search(r" * NormFactor = ", evtNr):
            NormFactor = evtNrWord[3]
            print "Norm factor is : ", NormFactor

        # Only want files with at least 6 input variables, otherwise they have not been processed
        if len(evtNrWord) == 7:
            EvtNrMatchingArray.append(str(evtNrWord[6]))
else:
    # Special case for TTbarJets_SemiLept since this has been splitted in correct/wrong/unmatched
    # ==> Will need to run over select the correct ones based on the directory name!
    NameOptions = ['Correct','Wrong','Unmatched']
    for iOpt in range(len(NameOptions)):
        if NameOptions[iOpt] in whichDir:
            print "Will be looking at ", NameOptions[iOpt], " events!"
            for evtNrTT in EvtNrMatching:
                evtNrWordTT = evtNrTT.split()

                #First get the lumi and the normFactor!
                if re.search(r" * Lumi = ", evtNrTT):
                    Luminosity = evtNrWordTT[3]
                    print "Luminosity is : ", Luminosity
                if re.search(r" * NormFactor = ", evtNrTT):
                    NormFactor = evtNrWordTT[3]
                    print "Norm factor is : ", NormFactor

                if len(evtNrWordTT) == 7 and str(evtNrWordTT[4]) == NameOptions[iOpt]:
                    EvtNrMatchingArray.append(str(evtNrWordTT[6]))
EvtNrMatching.close()

# -- Check which weights.out files exist in this directory -- #
weightFileArray = []
for weightFile in os.listdir(whichDir):
    if weightFile.startswith('weights.out') and weightFile.endswith('.out'):
        weightFileArray.append(weightFile)
        print "Found file : ", weightFile

if len(weightFileArray) > 1:
    print "Should decide which file should be used!! "
else:
    origWeights = open(whichDir+''+weightFileArray[0], 'r')
    newWeights = open(whichDir+'weights_SFAdded.out', 'w')
    for line in origWeights:
        word = line.split()

        if str(word[0]) != "#" and int(word[0]) <= 117658:
            newWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[int(word[0])-1]+' '+Luminosity+' '+NormFactor+'\n')
        elif str(word[0]) == "#":
            newWeights.write(line)
    origWeights.close()
    newWeights.close()