import os
import sys

os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C n n Wide 0 &')
os.system('rm nohup.out')
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C n n Wide 15 &')
os.system('rm nohup.out')
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C n n Wide 30 &')
os.system('rm nohup.out')

print " ---------------- Starting with cosTheta but no Acc norm -----------------------"
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C n y Wide 15 &')
os.system('rm nohup.out')
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C n y Wide 30 &')
os.system('rm nohup.out')

print " ---------------- Starting with acc but no costheta norm -----------------------"
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C y n Wide 15 &')
os.system('rm nohup.out')
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C y n Wide 30 &')
os.system('rm nohup.out')

print " ---------------- Starting with both cosTheta and Acc norm -----------------------"
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C y y Wide 15 &')
os.system('rm nohup.out')
os.system('nohup python CreateAllRootFiles_Nohup.py Events_RVRBiasScan_AllPossibleNorms/ RVR 50000 doublePolFitMacro.C y y Wide 30 &')
os.system('rm nohup.out')
