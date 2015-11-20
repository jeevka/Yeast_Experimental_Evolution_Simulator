from __future__ import division
from os import listdir
from os.path import isfile, join
import sys
import re

Files = [ f for f in listdir("Results/")]

N = len(Files)

NFix = 0;TFix = 0
TExt = 0;RFix = 0
RExt = 0;N_Fix = 0
NSeg = 0
Mean_Number_Mutations = 0
K = 0
for i in Files:
    K += 1
    Fname = "Results/" + i
    F1 = open(Fname,"r")
    for j in F1:
        if re.search("TotalMutations:",j):
            MNM = int(j.split()[1])
            Mean_Number_Mutations += int(j.split()[1])
            
        if re.search("NFix:",j):
            NFM = int(j.split()[1])
            NFix += int(j.split()[1])
        
        if re.search("TimeToFix:",j):
            T2F = 0
            try:
                T2F = float(j.split()[1])
                TFix += float(j.split()[1])
                N_Fix += 1
            except:
                pass

        if re.search("RateOfFix:",j):
            RFix += float(j.split()[1])
            ROF = float(j.split()[1])
            
        if re.search("TimetoExt:",j):
            TExt += float(j.split()[1])
            T2E = float(j.split()[1])

        if re.search("RateOfExt:",j):
            RExt += float(j.split()[1])
            ROE = float(j.split()[1])
            
        if re.search("Seggregating:",j):
            NSeg += float(j.split()[1])
            NSM = float(j.split()[1])
    
    print K,"\t",MNM,"\t",NFM,"\t",T2F,"\t",ROF,"\t",ROE,"\t",NSM
            
    F1.close()

print "Total Number of Mutations:",Mean_Number_Mutations/N
print "Number of Fixations:",NFix/N
print "Rate of Fixation:",RFix/N
try:
    print "Time to Fixation:",TFix/N_Fix
except:
    pass
print "Rate of Extinction:",RExt/N
print "Time to Extinction:",TExt/N
print "Seggregating:",NSeg/N