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

whichDir = "-1"
dirNames = ['TTbarJets', 'WJets', 'ZJets', 'SingleTop']
if len(sys.argv) > 1:
    whichDir = sys.argv[1]
    print "Specified which directory should be considered, so only for ", whichDir, "new weights files will be created!"
    print "This directory contains the following subdirectories that will be considered:"
    for subdir in os.listdir(whichDir):
        if subdir.startswith('Reco_'):
            print "   - ", subdir

sampleNames = ['TTbarJets_FullHadr', 'TTbarJets_FullLept', 'TTbarJets_SemiLept', 'SingleTop_tWChannel_t',
               'SingleTop_tWChannel_tbar', 'SingleTop_tChannel_t', 'SingleTop_tChannel_tbar', 'SingleTop_sChannel_t',
               'SingleTop_sChannel_tbar', 'WJets_1jets', 'WJets_2jets', 'WJets_3jets', 'WJets_4jets', 'ZJets_1jets',
               'ZJets_2jets', 'ZJets_3jets', 'ZJets_4jets']
whichSample = "-1"
if len(sys.argv) > 2:
    whichSample = sys.argv[2]
    print "Specified which sample to considered, so only for ", whichSample, " new weights file will be created!!"

CWUChoice = "-1"
CWUOptions = ['Correct', 'Wrong', 'Unmatched']
if len(sys.argv) > 3:
    CWUChoice = sys.argv[3]
    print "Specified whether Correct/Wrong/Unmatched TTbarJets_SemiLept events should be considered!!"

if str(whichSample) != "-1" and str(whichSample) != "TTbarJets_SemiLept":
    CWUChoice = "0"

Luminosity, NormFactor, NrEvts, consSamples = [], [], [], []
EvtNrMatchingArray = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
for iSample in range(len(sampleNames)):
    print "Looking at : ", sampleNames[iSample]
    if (str(whichSample) == "-1") or (str(whichSample) != "-1" and str(whichSample) in sampleNames[iSample]):
        print " --> Will continue with this one ! \n"

        if str(sampleNames[iSample]) != "TTbarJets_SemiLept":
            EvtNrMatching = open('MWEventNrMatching_'+sampleNames[iSample]+'.txt', 'r')
            consSamples.append(str(sampleNames[iSample]))   # Store the name of the considered sample for weights file!
            sampleSize = len(EvtNrMatchingArray)
            for evtNr in EvtNrMatching:

                # First get the lumi and the normFactor!
                if re.search(r" * Lumi = ", evtNr):
                    Luminosity.append(evtNrWord[3])
                    print "Luminosity is : ", evtNrWord[3]
                if re.search(r" * NormFactor = ", evtNr):
                    NormFactor.append(evtNrWord[3])
                    print "Norm factor is : ", evtNrWord[3]

                # Only want files with at least 6 input variables, otherwise they have not been processed
                if len(evtNrWord) == 7:
                    EvtNrMatchingArray[len(consSamples)-1].append(str(evtNrWord[6]))
            # NrEvts.append(len(EvtNrMatching)-sampleSize)
            print "Stored a total of ", len(EvtNrMatchingArray[len(consSamples)-1]), " events for ", sampleNames[iSample], "!"
            EvtNrMatching.close()
        else:
            # Special case for TTbarJets since this has been splitted in correct/wrong/unmatched!
            print " ... is semiLept thus separate index for correct/wrong/unmatched ..."

            for iCWU in range(len(CWUOptions)):
                print " Looking at CWU index of : ", CWUOptions[iCWU]
                if str(CWUChoice) == "-1" or (str(CWUChoice) != "-1" and str(CWUChoice) in CWUOptions[iCWU]):
                    print " --> Now looking at ", CWUOptions[iCWU], "events ! \n"

                    EvtNrMatching = open('MWEventNrMatching_'+sampleNames[iSample]+'.txt', 'r')
                    consSamples.append(str(sampleNames[iSample])+'_'+CWUOptions[iCWU])
                    print "Array with the considered samples has ", len(consSamples), " entries : ", consSamples
                    for evtNrTT in EvtNrMatching:
                        evtNrWordTT = evtNrTT.split()

                        # First get the lumi and the normFactor!
                        if re.search(r" * Lumi = ", evtNrTT):
                            Luminosity.append(evtNrWordTT[3])
                            print " Luminosity is : ", evtNrWordTT[3]
                        if re.search(r" * NormFactor = ", evtNrTT):
                            NormFactor.append(evtNrWordTT[3])
                            print "Norm factor is : ", evtNrWordTT[3]

                        if len(evtNrWordTT) == 7 and str(evtNrWordTT[4]) == CWUOptions[iCWU]:
                            EvtNrMatchingArray[len(consSamples)-1].append(str(evtNrWordTT[6]))
                    print "Stored a total of ", len(EvtNrMatchingArray[len(consSamples)-1]), " events !"
                    EvtNrMatching.close()

print " "
# -- Now go into the different directories and change the weights.out file -- #
normFactorAdded = False
for mcDir in os.listdir("."):
    if (str(whichDir) == "-1" and not 'Combined' in mcDir and not '.' in mcDir) or \
            (str(whichDir) != "-1" and str(whichDir) in mcDir and not '.' in mcDir):

        print "Looking at directory : ", mcDir
        if str(mcDir) == "TTbarJets":
            combinedSemiMuWeights = open(mcDir+'/weights_CheckedEvts_CombinedSemiLeptEvts_SFAdded.out', 'w')

        for MCsubdir in os.listdir(mcDir):
            if (str(whichSample) == "-1" and MCsubdir.startswith('Reco_')) or \
                    (str(whichSample) != "-1" and str(whichSample) in MCsubdir):
                print " --> Continuing with directory ", MCsubdir

                if not 'TTbarJets_SemiLept' in MCsubdir:

                    origWeights = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts.out', 'r')
                    newWeights = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts_SFAdded.out', 'w')

                    # Get the correct index for the filled arrays!
                    arrayIndex = 999
                    for i, j in enumerate(consSamples):
                        if MCsubdir.startswith('Reco_'+str(j)):
                            arrayIndex = int(i)

                    newWeights.write(consSamples[arrayIndex] + ' ' + NormFactor[arrayIndex] + '\n')
                    for line in origWeights:
                        word = line.split()

                        if str(word[0]) != "#" and int(word[0]) <= 117658:
                            newWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' ' +
                                             EvtNrMatchingArray[arrayIndex][int(word[0])-1]+'\n')
                            if str(whichSample) == "-1":
                                combinedSampleWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '
                                                            + EvtNrMatchingArray[arrayIndex][int(word[0])-1]+'\n')
                        #elif str(word[0]) == "#":
                        #    newWeights.write(line)
                    origWeights.close()
                    newWeights.close()
                else:
                    # Special case for semiLept events because they are splitted in Correct/Wrong/Unmatched!

                    for iCase in range(len(CWUOptions)):
                        print "  * Looping over the different CWU options ", CWUOptions[iCase]
                        if str(CWUChoice) == "-1" or (str(CWUChoice) != "-1" and str(CWUChoice) in CWUOptions[iCase]):
                            if CWUOptions[iCase] in MCsubdir:
                                origWeightsTT = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts.out', 'r')
                                newWeightsTT = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts_SFAdded.out', 'w')

                                # Get the correct index for the filled arrays!
                                arrayIndex = 999
                                for i, j in enumerate(consSamples):
                                    if MCsubdir.startswith('Reco_'+str(j)):
                                        arrayIndex = int(i)

                                newWeightsTT.write(consSamples[arrayIndex] + ' ' + NormFactor[arrayIndex] + '\n')
                                if not normFactorAdded:
                                    combinedSemiMuWeights.write('TTbarJets_SemiLept '+ NormFactor[arrayIndex] + '\n')
                                    normFactorAdded = True

                                for line in origWeightsTT:
                                    word = line.split()

                                    if str(word[0]) != "#":
                                        newWeightsTT.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' ' +
                                                           EvtNrMatchingArray[arrayIndex][int(word[0])-1]+'\n')
                                        if str(CWUChoice) == "-1":
                                            combinedSemiMuWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+
                                                                        word[4] + ' ' + EvtNrMatchingArray[arrayIndex][int(word[0])-1]+'\n')
                                    #elif str(word[0]) == "#":
                                    #    newWeightsTT.write(line)
                                origWeightsTT.close()
                                newWeightsTT.close()

        if str(mcDir) == "TTbarJets":
            combinedSemiMuWeights.close()
            if str(CWUChoice) != "-1":
                os.remove(combinedSemiMuWeights.name)