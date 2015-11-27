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

sampleNames = {'TTbarJets_FullHadr', 'TTbarJets_FullLept', 'TTbarJets_SemiLept', 'SingleTop_tWChannel_t',
               'SingleTop_tWChannel_tbar', 'SingleTop_tChannel_t', 'SingleTop_tChannel_tbar', 'SingleTop_sChannel_t',
               'SingleTop_sChannel_tbar', 'WJets_1jets', 'WJets_2jets', 'WJets_3jets', 'WJets_4jets', 'ZJets_1jets',
               'ZJets_2jets', 'ZJets_3jets', 'ZJets_4jets'}
whichSample = "-1"
if len(sys.argv) > 2:
    whichSample = sys.argv[2]
    print "Specified which sample should be considered, so only for ", whichSample," new weights file will be created!!"

CWUChoice = "-1"
CWUOptions = ['Correct','Wrong','Unmatched']
if len(sys.argv) > 3:
    CWUChoice = sys.argv[3]
    print "Specified whether Correct/Wrong/Unmatched TTbarJets_SemiLept events should be considered!!"

Luminosity, NormFactor, NrEvts, consSamples = [], [], [], []
EvtNrMatchingArray = [[]]
for iSample in range(len(sampleNames)):
    if (str(whichSample) == "-1") or (str(whichSample) != "-1" and str(whichSample) in sampleNames[iSample]):

        EvtNrMatching = open('MWEventNrMatching_'+sampleNames[iSample]+'.txt', 'r')
        if str(sampleNames[iSample]) != "TTbarJets_SemiLept":
            consSamples.append(str(sampleNames[iSample]))   # Store the name of the considered sample for weights file!
            sampleSize = len(EvtNrMatchingArray)
            for evtNr in EvtNrMatching:

                #First get the lumi and the normFactor!
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
        else:
            # Special case for TTbarJets since this has been splitted in correct/wrong/unmatched!

            for iCWU in range(len(CWUOptions)):
                if str(CWUChoice) == "-1" or (str(CWUChoice) != "-1" and str(CWUChoice) in CWUOptions):
                    print "Now looking at ", NameOptions[iCWU], "events !"

                    consSamples.append(str(sampleNames[iSample])+'_'+CWUOptions[iCWU])
                    for evtNrTT in EvtNrMatching:
                        evtNrWordTT = evtNrTT.split()

                        #First get the lumi and the normFactor!
                        if re.search(r" * Lumi = ", evtNrTT):
                            Luminosity.append(evtNrWordTT[3])
                            print "Luminosity is : ", evtNrWordTT[3]
                        if re.search(r" * NormFactor = ", evtNrTT):
                            NormFactor.append(evtNrWordTT[3])
                            print "Norm factor is : ", evtNrWordTT[3]

                        if len(evtNrWordTT) == 7 and str(evtNrWordTT[4]) == CWUOptions[iCWU]:
                            EvtNrMatchingArray[len(consSamples)-1].append(str(evtNrWordTT[6]))
        EvtNrMatching.close()
    else:
        print "Issue with sample name ... --> ", whichSample

# -- Now go into the different directories and change the weights.out file -- #
for mcDir in os.listdir():
    if (str(whichDir) == "-1" and not 'Combined' in mcDir and not '.' in mcDir) or (str(whichDir) != "-1" and str(whichDir) in dirNames):
        combinedMCWeights = open('weights_CheckedEvts_CombinedMC_SFAdded.txt', 'w')  # One output file for all MC!

        print "Looking at directory : ", mcDir
        for MCsubdir in os.listdir(mcDir):
            if (str(whichSample) == "-1" and MCsubdir.startswith('Reco_')) or (str(whichSample) != "-1" and str(whichSample) in sampleNames):
                combinedSampleWeights = open(mcDir+'/weights_CheckedEvts_Combined'+mcDir+'_SFAdded.txt')

                if not 'TTbarJets_SemiLept' in MCsubdir:

                    origWeights = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts.txt')
                    newWeights = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts_SFAdded.txt')

                    # Get the correct index for the filled arrays!
                    arrayIndex = 999
                    for i, j in enumerate(consSamples):
                        if MCsubdir.startswith('Reco_'+str(j)):
                            arrayIndex = int(i)

                    for line in origWeights:
                        word = line.split()

                        if str(word[0]) != "#" and int(word[0]) <= 117658:
                            newWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[arrayIndex][int(word[0])-1]+' '+Luminosity[arrayIndex]+' '+NormFactor[arrayIndex]+' '+consSamples[arrayIndex]+'\n')
                            if str(whichSample) == "-1":
                                combinedSampleWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[arrayIndex][int(word[0])-1]+' '+Luminosity[arrayIndex]+' '+NormFactor[arrayIndex]+' '+consSamples[arrayIndex]+'\n')
                        elif str(word[0]) == "#":
                            newWeights.write(line)
                    origWeights.close()
                    newWeights.close()
                else:
                    # Special case for semiLept events because they are splitted in Correct/Wrong/Unmatched!
                    combinedSemiMuWeights = open(mcDir+'/weights_CheckedEvts_CombinedSemiMuEvts_SFAdded.txt')

                    for iCase in range(len(CWUChoice)):
                        if str(CWUChoice) == "-1" or (str(CWUChoice) != "-1" and str(CWUChoice) in CWUOptions):
                            if CWUChoice[iCase] in MCsubdir:
                                origWeightsTT = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts.txt')
                                newWeightsTT = open(mcDir+'/'+MCsubdir+'/weights_CheckedEvts_SFAdded.txt')

                                # Get the correct index for the filled arrays!
                                arrayIndex = 999
                                for i, j in enumerate(consSamples):
                                    if MCsubdir.startswith('Reco_'+str(j)):
                                        arrayIndex = int(i)

                                for line in origWeights:
                                    word = line.split()

                                    if str(word[0]) != "#" and int(word[0]) <= 117658:
                                        newWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[arrayIndex][int(word[0])-1]+' '+Luminosity[arrayIndex]+' '+NormFactor[arrayIndex]+' '+consSamples[arrayIndex]+'\n')
                                        if str(whichSample) == "-1":
                                            combinedSampleWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[arrayIndex][int(word[0])-1]+' '+Luminosity[arrayIndex]+' '+NormFactor[arrayIndex]+' '+consSamples[arrayIndex]+'\n')
                                        if str(CWUChoice) == "-1":
                                            combinedSemiMuWeights.write(word[0]+' '+word[1]+' '+word[2]+' '+word[3]+' '+word[4]+' '+EvtNrMatchingArray[arrayIndex][int(word[0])-1]+' '+Luminosity[arrayIndex]+' '+NormFactor[arrayIndex]+' '+consSamples[arrayIndex]+'\n')
                                    elif str(word[0]) == "#":
                                        newWeights.write(line)
                            origWeightsTT.close()
                            newWeightsTT.close()

                    # Remove this combined SemiLept weights output file in case only 1 CWU type has been considered
                    combinedSemiMuWeights.close()
                    if str(CWUChoice) != "-1":
                        os.system('rm '+combinedSemiMuWeights.name)

                # Remove this combined MC type weights output file in case only 1 MC sample has been considered
                combinedSampleWeights.close()
                if str(whichSample) != "-1":
                    os.system('rm '+combinedSampleWeights.name)

        # Remove this combined weights output file in case only 1 MC directory has been considered!
        combinedMCWeights.close()
        if str(whichDir) != "-1":
            os.system('rm '+combinedMCWeights.name)