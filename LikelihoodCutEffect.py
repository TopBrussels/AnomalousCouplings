# ----------------------------------------------------------------------------------------------------------- #
#  Script which will compare the different likelihood cut values and help to decide on the ideal cut value!   #
#  Will select the correct ROOT file and then transfer to the .C script to create the information!
# ----------------------------------------------------------------------------------------------------------- #

import sys
import os
import shutil

if len(sys.argv) <= 1:
    print "Need to give the directory of interest !"
    sys.exit()

whichDir = os.path.join(os.getcwd() + '/' + sys.argv[1])
whichDir = whichDir[whichDir.find('user_mnt')+8:]      # Remove the /user_mnt part from the name!

LatexFileName = 'LikelihoodCutEffect.tex'
# Check whether these output files already exist, otherwise delete them !
if os.path.isfile(os.path.join(whichDir + LatexFileName[:-4] + ".pdf")):
    os.system('rm ' + whichDir + LatexFileName[:-4] + '.pdf')
if os.path.isfile(os.path.join(whichDir + LatexFileName)):
    os.system('rm ' + whichDir + LatexFileName)

os.system('g++ -m64 -g -I `root-config --incdir` `root-config --libs` LikelihoodCutEffect.C -o LikelihoodCutEffect')
# --------------------------------------------------------- #
# Check the different ROOT files existing in this directory #
# --> Send each file to the .C script                       #
# --------------------------------------------------------- #
list_dir = os.listdir(sys.argv[1])
fileCounter = 0
consFiles = " "
for file in list_dir:
    if file.endswith(".root") and 'LikelihoodCut' in file:
        consFiles += str(file) + " "
os.system("./LikelihoodCutEffect " + whichDir + ' ' + LatexFileName + ' ' + consFiles)

os.system('pdflatex -interaction=batchmode '+ LatexFileName)
shutil.move(LatexFileName, whichDir)
shutil.move(LatexFileName[:-4] + ".pdf", whichDir)
os.system('rm '+ LatexFileName[:-4] + ".aux")
os.system('rm '+ LatexFileName[:-4] + ".log")
