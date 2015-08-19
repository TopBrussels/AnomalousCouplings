######################################################
#                                                    #
#  Python script which will take care of MW-running  #
#    --> List of values given will be considered     #
#    --> Permutations & number of events should be   #
#        set in advance                              #
#    --> Range of interest should be set in advance  #
#                                                    #
######################################################

import os
import re
import sys
import time

VariableList = ["Neg04","Neg03","Neg02","Neg015","Neg01","Neg005","SM","Pos005","Pos01","Pos015","Pos02","Pos03","Pos04"]
DirectoryName = "RgR_CalibrationCurve_AllDeltaTF_RgRMGSample"+VariableList[0]+"_25000Evts_CutsOnMET"
#VariableList = ["Neg03","Neg02","Neg01","SM","Pos01","Pos02","Pos03"]
#DirectoryName = "RVR_ControlCheckForCosTheta_EDepTF_AllDelta_RVRMGSample"+VariableList[0]+"_2500Evts_NarrowShortRange_NoCuts"
FileName = 'RgR'+VariableList[0]+'_RVL1_CutsOnMET.lhco'

MissingLHCOFiles = []

#Loop over each of the elements in this array
for iVar in range(len(VariableList)):

  #Change the DirectoryName for all the folowing variables!
  if not iVar == 0: 
    DirectoryName = DirectoryName.replace(VariableList[iVar-1],VariableList[iVar])
    FileName = FileName.replace(VariableList[iVar-1],VariableList[iVar])
  print "Updated directoryName : ",DirectoryName

  #Now go to the MadWeight_card.dat file and change both the name and the inputfile!
  MWCard = file('Cards/MadWeight_card.dat','r')
  output = file('Cards/output.dat','w')
  for MWLine in MWCard:
    MWWords = MWLine.split()
    if re.search( r"name", MWLine) and MWWords[0] == "name":
      #Now check whether the DirectoryName is already correct or is still using the previous config!       --> Not possible to check because Neg01 is found also in Neg015 string!!
      # ==> Just change it for every loop!
      output.write(MWLine.replace(MWWords[1],DirectoryName))
    elif re.search( r"inputfile",MWLine) and MWWords[0] == "inputfile":
      output.write(MWLine.replace(MWWords[1],FileName))
    else:
      output.write(MWLine)
  output.close()
  os.rename(output.name,MWCard.name)
  MWCard.close()

  #Before running make sure that the inputfile exists, otherwise skip the event and print a big warning at the end of the output file!! 
  if not (os.path.exists(os.path.join("Events/"+FileName))):
    print "  *** Filename ",FileName," does not exist !! "
    MissingLHCOFiles.append(FileName)
  else:
    #Run the MadWeight calculations only for these events !!
    os.system("./bin/madweight.py -1234")
    NrRemainingEvts = os.popen('qstat @cream02 | grep aolbrech | wc -l').read()
    while int(NrRemainingEvts) != 0:
      print " ** Still jobs on the cluster which are running, so should wait to start the -6 step !!"
      print "    - Start time :",time.ctime()
      time.sleep(120)
      NrRemainingEvts = os.popen('qstat @cream02 | grep aolbrech | wc -l').read()
      print "    - Time after waiting : ",time.ctime()," --> Number remaining jobs is now : ",NrRemainingEvts
    os.system("./bin/madweight.py -6 ")
    
    #Before going to the collect or -8 step, make sure no events remain on the cluster
    NrRemainingEvts = os.popen('qstat @cream02 | grep aolbrech | wc -l').read()
    while int(NrRemainingEvts) != 0:
      print " ** Still jobs on the cluster which are running, so should wait to start the -8 step !!"
      print "    - Start time :",time.ctime()
      time.sleep(120)
      NrRemainingEvts = os.popen('qstat @cream02 | grep aolbrech | wc -l').read()
      print "    - Time after waiting : ",time.ctime()," --> Number remaining jobs is now : ",NrRemainingEvts
       
    os.system("./bin/madweight.py -8 ")

print "\n List of missing lhco files is : ",MissingLHCOFiles
