from __future__ import division
import sys
import re

def Analytical_Test(Total_Mutations):
    F1 = open("Mutation_Freq.csv","r")
    
    # Part 1: Store the data
    M = {}; Freq = {};Fix_M = {}
    for i in F1:
        temp = i.split()
        if M.has_key(temp[1]):
            M[temp[1]].append(temp[0])
            Freq[temp[1]].append(float(temp[2]))
            Fix_M[temp[1]] = float(temp[2])
        else:
            M[temp[1]] = [temp[0]]
            Freq[temp[1]] = [float(temp[2])]
            Fix_M[temp[1]] = float(temp[2])
    F1.close()
    
    # To be excluded Mutations from extinction calculation 
    excluded_M = []
    for i in M:
        L = len(M[i])
        if int(M[i][L-1]) == 1000 and Fix_M[i] != 1.0:
            excluded_M.append(i)
    
    # Part II: Time to Fix
    N_fix = 0
    Fixed_Muts = []    
    for i in Fix_M:
        if Fix_M[i] == 1.0:
            N_fix += 1
            Fixed_Muts.append(i)        
    
    print "TotalMutations:",Total_Mutations
    print "NFix:",N_fix
    print "RateOfFix:",N_fix/Total_Mutations
    print "RateOfExt:", 1 - (N_fix/Total_Mutations)
    
    N_G_Fix = 0   
    # Calculate Time to Fix
    if len(Fixed_Muts) > 0:
        for i in Fixed_Muts:
            TN = 0
            for j in Freq[i]:
                TN += 1
                if j == 1.0:
                    break
            N_G_Fix += TN
    
    print "TimeToFix:",N_G_Fix/N_fix
    # Calculate time to extinction
    Time_2_ext = 0
    N_E = 0
    for i in Freq:
        if i not in excluded_M and i not in Fixed_Muts:
            Time_2_ext  += len(Freq[i])
            N_E += 1
    
    Lost_mutations = Total_Mutations - N_E - N_fix - len(excluded_M)
    Time_2_ext += Lost_mutations
    
    print "TimetoExt:",Time_2_ext/(Total_Mutations - N_fix - len(excluded_M))
    
    # Rate of seggregating
    print "Seggregating:",len(excluded_M)/Total_Mutations   
    
    """"    
    if N_fix != 0:
        print "TimeToFix:",N_G_Fix/N_fix
    else:
        print "TimeToFix:","NONE"
    print "RateOfFix:",N_fix/len(M)
    print "TimetoExt:",Time_2_ext/len(M)
    print "RateOfExt:", 1 - (N_fix/len(M))
    """
    
#Analytical_Test(12231)