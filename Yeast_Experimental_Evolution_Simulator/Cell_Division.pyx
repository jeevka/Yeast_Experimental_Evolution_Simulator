from __future__ import division
import sys
import random
from copy import deepcopy
from heapq import *
from stdlib cimport *
from numpy import *
import numpy as np
cimport numpy as np
cimport cython
import gc
import re  
"""
This module memics the mitotic cell division of haploid Yeast cells.
Cells divide and increase in number. During each cell divison, cells
may get mutations based on the given mutation rate. After certain generations,
whole population will go thourgh a bottlenecl where we randomly choose N cells
and restart the whole cycle again.
E.g. Initial popilation (N) = 1000. After five mitotic division, the whole population
will be 32000. Then, 1000 cells wil be chosen randomly the selected cells are allwed to grow
for another 5 generartions. THis whole cycle will repeated number of times.
"""
##############################################################################################
################################ User Defined modules ########################################
##############################################################################################
import Yeast_Simulator
import Sub_Functions
import Mutations
##############################################################################################

# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

##############################################################################################
# Each cell divides at differnt speed. 
##############################################################################################
@cython.boundscheck(False)
def asymmetrical_cell_division(genotypes,cell_groups,PM_fitness,PM,n_generations,n_BN,n_individuals,n_mut,dps,FAMR,mutation_fitness,mutation_rate,BMR,cell_max_age):

    # Output File for Mutation Frequency
    Mut_Freq = open("Mutation_Freq.csv","w")
    
    # Mean Mutation ouput File
    Mean_Mut = open("Mean_Mutations.csv","w")
    
    # Cython variable Declarations
    cdef int i,time,tp
    cdef int small
    cdef int k,h,n_cells, P_cdt, C_cdt

    import heapq
    groups = cell_groups.keys()
    groups.sort()
    
    time = heapq.heappop(groups)
    
    # Number of Bottlenecks
    for i in xrange(n_BN):
        
        #print "Bottleneck :", i
        tp = n_individuals
        
        while tp < dps: 
            # Iterating through the cell groups.
            
            n_cells = len(cell_groups[time])
            for h in xrange(n_cells):
                k = cell_groups[time][h]
                
                # Copying the parents attributes to child.
                genotypes[tp] = genotypes[k]
                genotypes[tp][2] = cell_max_age
                    
                # Copy the mutations from mother to daughter
                PM[tp] = PM[k][:]
                    
                # Introduce mutations
                if random.random() <= mutation_rate:
                   # Sending mother and daughter cells to introduce mutations. 
                   genotypes[k],genotypes[tp],PM, PM_fitness = Mutations.introduce_mutation(genotypes[k],genotypes[tp],BMR,k,tp,n_mut,PM,PM_fitness,FAMR,mutation_fitness)
                   n_mut += 1
                
                # Updating  the next cell division time based on the cell division time for both mother and daughter cells.
                (genotypes[k][0],genotypes[tp][0]) = update_next_CDT(genotypes[k][0],genotypes[k][1],genotypes[tp][0],genotypes[tp][1])                  
                    
                # Group the cells from current time point
                cell_groups, groups = group_the_cells(cell_groups,groups,genotypes[k][0],genotypes[tp][0],k,tp)
                
                # increase the population
                tp += 1
                
                # Come out of the loop once its reached the desired population size.
                if tp == dps:
                   break
            
            del cell_groups[time]
            time = heapq.heappop(groups)
            
        # Sampling - Bottleneck
        if i != n_BN-1: 
            # Random smapling            
            (genotypes,cell_groups,PM) = sampling_individuals(genotypes,n_individuals,i,PM,n_BN)
            dps = n_individuals * (2**5)
            
            # Mutation Frequencies
            calculate_mutation_frequencies(PM,i,n_individuals,Mut_Freq,Mean_Mut)

        else:
            (genotypes1,cell_groups,PM) = sampling_individuals(genotypes,n_individuals,i,PM,n_BN)

            # Fixation analysis: The conditon is to avoid unnecessary fixation calculation where its not needed.
            # Because, Fixation calculation is one of the heavy and time consuming process.
            
            # Mutation Frequencies
            calculate_mutation_frequencies(PM,i,n_individuals,Mut_Freq,Mean_Mut)
            
            # Closing the output files
            Mut_Freq.close()
            Mean_Mut.close()
            
            return genotypes1,PM,PM_fitness,n_mut

        # Its a small mess up. Have to find a better way.
        groups = []
        groups = cell_groups.keys()
        time = heapq.heappop(groups)        
        
#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
cdef sampling_individuals(genotypes,int n_individuals,int i,PM1,n_BN):
    """
    During serial transfers, N number of cells are chosen randomly.
    At the same time mean cell division time of the whole population is calculated
    
    Input:
        Dictionary of lists with cell properties(e.g. cell division time)
        Serial transfer number
        Number of individuals
    
    Output:
        Mean cell division time of the population after random sampling
    """
    # Making the Point mutation and gene duplication Dict of list
    size = len(genotypes)
    size = n_individuals * (2**5)
    PM = {}; 
    for i1 in range(size):
        PM[i1] = []
    
    # Sampling the number of IDs from genotypes dict.
    cdef int n_cells,j
    cdef float total_time = 0
    size = len(genotypes)
    size = n_individuals * (2**5)
    
    sampled_ids = random.random_integers(size-1,size=n_individuals)
    sampled_group = ones((size,4),dtype=int32)
    
    cell_mutations_temp = []
        cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    n_cells = len(sampled_ids)
    n_alive_cells = 0
    for j in xrange(n_cells):
        k = sampled_ids[j]
        
        total_time += genotypes[k][1]
        n_alive_cells += 1
        sampled_group[j] = genotypes[k]
        
        # Copying PM and GD info of cells
        PM[j] = PM1[k][:]

        # Storing the mutations in sampled cells.
        # cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene duplications in sampled cells.
        # cell_gene_duplications_temp.append(cell_gene_duplications[k])
        
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]
    
    mean_CDT = total_time/n_individuals
    LSC = calculate_LSC(mean_CDT)
    # This is to indicate the end of each bottleneck
    print i+1,"\t",mean_CDT,"\t",LSC
    
    return sampled_group,cell_groups,PM


##############################################################################################
# Calculating LSC from Mean Cell divisiobn time after each bottleneck/serial transfer
# In otherwords, LSC is relative growth rate of a mutatnt cell growing in stress condition
# relative to WT cells growing in Normal conditions 
##############################################################################################
cdef calculate_LSC(mean_CDT):
        """
        Input:
            Mean cell division time of a population
        
        Output:
            Growth rate relative to WT cell growing in normal growth media
        """
        
        # For Below formula refer to Jonas Warringer's paper.
        LSC = (math.log(180/60) - math.log(mean_CDT/60))/2
        LSC = LSC/0.3465735902799727        
        
        return LSC
        
##############################################################################################
# This function calculates Frequency of mutations after very bottlenecks
##############################################################################################
cdef calculate_mutation_frequencies(PM,N_BN,N,Mut_Freq,Mean_Mut):
    """
    Input:
        Mutation profiles of each cell in the population
    
    Output:
        Mutation frequency of each mutation appeared in simulation        
    """
    Freq = {}
    N_mut = 0
    
    # Counting Part
    for i in xrange(N):
        for j in PM[i]:
            if Freq.has_key(j):
                Freq[j] += 1
            else:
                Freq[j] = 1
            
            N_mut += 1
            
    # Printing Part. Only the mutations with frequency of >=0.05 are printed
    N_Fixed = 0
    for i in Freq:
        frequency = Freq[i]/N
        if frequency >= 0.05:
            txt = str(N_BN+1) + "\t" + str(i) + "\t" + str(frequency) + "\n"
            Mut_Freq.write(txt)
        if frequency == 1:
            N_Fixed += 1
    
    # Writing Mean Mutations into output file
    txt = str(N_BN+1) + "\t" + str(N_mut/N) + "\t" + "\n"
    Mean_Mut.write(txt)
    
    #if N_BN == 49: 
    #    print "Total Number of Mutations:",N_mut
    #    print "Mean Number of Mutations Per Cell:",N_mut/N
    #    print "Number of Fixed Mutations:",N_Fixed
    
    return 0
###############################################################################################
# To update the next cell division time.
###############################################################################################
cdef inline update_next_CDT(int a1, int a2, int b1, int b2):
     """
     Input:
        Current cell divison time
     Output:
        New cell division time after the mutation
     """
     a1 = a1 + a2
     b1 = b1 + b2
     return a1,b1

###############################################################################################
# To group the cels: Cells which are dividing at "identical" cell division time will be in
# one group. The Cell ids are stored under same id if they have cell divison times are identical
# E.g. {180: [10,15,25,100]} 180 is the cell division time. 10,15,25 and 100 are cell ids.
# When the clock is 180, all these cells will divide at the same time. 
###############################################################################################
import heapq
def group_the_cells(cell_groups,groups,P_cdt,C_cdt,k,tp):
    """
    Input:
        Cell groups
        Parent's cell division time
        Daughter's cell division time
        Cell ids of parent and daughter
    Output:
        New Cell groups
    """
    if P_cdt == C_cdt:               
        # Grouping the mother cells
        if cell_groups.has_key(P_cdt):
            cell_groups[P_cdt].append(k)
            cell_groups[C_cdt].append(tp)
        else:
            cell_groups[P_cdt] = [k]
            cell_groups[C_cdt].append(tp)
            heapq.heappush(groups,P_cdt)
    else:
        # Grouping the mother cells
        if cell_groups.has_key(P_cdt):
            cell_groups[P_cdt].append(k)
        else:
            cell_groups[P_cdt] = [k]
            heapq.heappush(groups,P_cdt)                    
        
        # Grouping the daughter cells.
        if cell_groups.has_key(C_cdt):
            cell_groups[C_cdt].append(tp)
        else:
            cell_groups[C_cdt] = [tp]
            heapq.heappush(groups,C_cdt)
    
    if not P_cdt in groups:
        heapq.heappush(groups,P_cdt)

    if not C_cdt in groups:
        heapq.heappush(groups,C_cdt)        
    
    return cell_groups,groups
