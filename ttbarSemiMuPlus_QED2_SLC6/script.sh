#!/bin/bash          

##Some general shell commands
STR="Hello World!"
echo $STR    
echo ">> script.sh is checking where it is"
pwd
echo ">> script.sh is checking how much disk space is still available"
df -h
echo ">> script.sh is listing files and directories in the current location"
ls -l
echo ">> script.sh is listing files and directories in userdir on storage element"
ls -l /pnfs/iihe/cms/store/user/$USER

##When accessing files on the storage element it is imporant to execute your code on the /scratch partition of the workernode you are running on. Therefore you need to copy your executable which is accessing/writing root files onto the /scratch partition and execute it there. This is illustrated below.

echo ">> cd into /scratch partition"
cd /scratch
echo ">> ls of /scratch partition"
ls -l

##Create a small root macro

echo "{
  //TFile *MyFile = new TFile(\"testfile.root\",\"RECREATE\"); 
  //MyFile->ls();
  //MyFile->Close(),
  TFile* f=TFile::Open(\"dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/$USER/testfile.root\");
  f->ls();
  f->Close();
} 
" > rootScript.C

cat rootScript.C

echo ">> set root"
##Copied a root version from /user/cmssoft into /localgrid
export ROOTSYS=/localgrid/$USER/cmssoft/root_5.26.00e_iihe_default_dcap/root 
export PATH=$PATH:$ROOTSYS/bin 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib

echo ">> execute root macro"
root -q -l -b -n rootScript.C

echo ">> ls of /scratch partition"
ls -l

echo "copy the file back to the /localgrid sandbox"
#cp testfile.root /localgrid/jmmaes/directsubmissiontest
