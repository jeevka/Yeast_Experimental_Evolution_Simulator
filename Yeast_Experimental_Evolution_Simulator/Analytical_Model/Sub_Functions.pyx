####################################################################################################
# This module contains the common subroutines for the Yeast Simulator.
####################################################################################################
####################################################################################################
# Importing Python Modules
####################################################################################################
from __future__ import division
import sys
import random
import operator
from copy import copy
from copy import deepcopy
import math
import re
import numpy as np
from heapq import heapify
from numpy import *
#####################################################################################################
# Importing user defined modules.
#####################################################################################################
import Yeast_Simulator
import Mutations

################################## Intialization ####################################################
# Intialization of the individuals before the cell division starts.
# 1. Initial population contains cell at different ages. So the cells divide at different speeds.
# 2. 50% cells are virgin cells. 25% are one time divided. 12.5% are two times divided. and etc ...
#####################################################################################################
def initialize_population(n_individuals,n_generations,cell_division_min,cell_division_max,cell_max_age):
    cdef int i,n
    size = n_individuals * 2**n_generations
    genotypes = ones((size,4),dtype=int32)
    
    import random
    population_ratio = {}
    population_ratio[0] = 0
    population_percentage = n_individuals
    #genotypes = {}
    cell_mutations = []
    cell_gene_deletions = []
    cell_gene_duplications = []
    a = []
    # Getting cell division time proportions based on age.
    #CDT_proportions  = age_cell_divison_time(cell_division_max)

    for n in xrange(n_individuals):
        # Initializing the cell division time and number of generations.
        # Choosing between virgin cells and other cells. Because ~50% cells are virgin cells in the begining.
        # Here n_individuals * 100000 is to get more accurate number upto 5 precision so that the ratio will be better even with small numbers.
            
        genotypes[n] = [cell_division_max,cell_division_max,cell_max_age,0]
        
    # Making the Point mutation and gene duplication Dict of list 
    PM = {}; 
    for i in range(size):
        PM[i] = []
        
    return genotypes,PM

##############################################################################################
# In Reality 10% of the cells wont survive or wont reproduce. This function will do the killings.
##############################################################################################
def inactivate_cells(genotypes,percentage):
    inactive_ids = {}
    n_inactive_cells = int(len(genotypes) * (percentage/100))
    inactive_ids = random.sample(genotypes,n_inactive_cells)    
    for i in range(0,len(inactive_ids)):
        id = inactive_ids[i]
        genotypes[id][2] = 0
    return genotypes

#############################################################################################
# Group the cells based on their cell division time.
#############################################################################################
def group_cells(genotypes,n_individuals):
    cell_group = {}
    n_individuals = int(n_individuals)
    for i in xrange(n_individuals):
        if cell_group.has_key(genotypes[i][0]):
            cell_group[genotypes[i][0]].append(i) 
        else:
            cell_group[genotypes[i][0]] = [i]
            
    return cell_group

#################################################################################################
# Updating the next cell divison time and age
#################################################################################################
"""
def update_CDT_age(parent,daughter):
    # Updating the next cell division time
    parent[0] = int(parent[0]) + int(parent[1])
    daughter[0] = int(daughter[0]) + int(daughter[1])
    # Updating the age
    parent[2] = parent[2] - 1 
    daughter[2] = daughter[2] - 1
    return parent,daughter
"""

#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
"""
def sampling_individuals(genotypes,n_individuals,cell_mutations,cell_gene_duplications):
    # Sampling the number of IDs from genotypes dict.
    sampled_ids = random.sample(genotypes,n_individuals)        
    sampled_group = {}
    cell_mutations_temp = []
    #cell_gene_deletions_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    for j in xrange(len(sampled_ids)):
        k = sampled_ids[j]
       
        sampled_group[j] = genotypes[k][:]
        # Storing the mutations in sampled cells.
        cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene deletions in sampled cells.
        #cell_gene_deletions_temp.append(cell_gene_deletions[k])
        # Storing the gene duplications in sampled cells.
        cell_gene_duplications_temp.append(cell_gene_duplications[k])
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        #print "KLKL:",kl
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]

    #print "Sampling :",len(cell_mutations)
    return sampled_group,cell_groups,cell_mutations_temp,cell_gene_duplications_temp
"""
#################################################################################################
# non-Random : Sampling the desired number of individuals after certain number of generations.
#################################################################################################
"""
def sampling_individuals_test(genotypes,n_individuals,cell_mutations,cell_gene_duplications,small,big):
    # Sampling the number of IDs from genotypes dict.
    sampled_ids = group_the_cells(genotypes,n_individuals,small,big)

    sampled_group = {}
    cell_mutations_temp = []
    #cell_gene_deletions_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    for j in xrange(len(sampled_ids)):
        k = sampled_ids[j]
        sampled_group[j] = genotypes[k][:]
        # Storing the mutations in sampled cells.
        cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene deletions in sampled cells.
        #cell_gene_deletions_temp.append(cell_gene_deletions[k])
        # Storing the gene duplications in sampled cells.
        cell_gene_duplications_temp.append(cell_gene_duplications[k])
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        #print "KLKL:",kl
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]
    #print "Sampling :",len(cell_groups)
    return sampled_group,cell_groups,cell_mutations_temp,cell_gene_duplications_temp

"""
##############################################################################################
# This function will recalculate the cell division time based on its Age
##############################################################################################
"""
def age_based_division_time_change(mother,daughter):
    age = 16 - mother[2]
    mother[1] = mother[1] + (cdt_pd[age-1] * mother[1])

    return mother,daughter
""" 
##############################################################################################
# This function will break the whole population into 3 major groups based on their division time.
##############################################################################################
"""
def group_the_cells(genotypes,n,small,big):
    #fa = 40;me = 20;sl = 40
    break_point = (big - small)/3
    p1 = small + break_point
    p2 = p1 + break_point
    p3 = p2 + break_point
    t = []
    # Faster dividing cells
    faster = {}
    # Medium dividing cells
    medium = {}
    # Slower dividing cells
    slower = {}
    for i in xrange(0,len(genotypes)):
        if genotypes[i][1] < p1:
            faster[i] = genotypes[i][1]
        elif genotypes[i][1] < p2:
            medium[i] = genotypes[i][1]
        else:     
            slower[i] = genotypes[i][1]
    # Decide number of cells in each group based on available population size in each group
    if (float(n)*(fcp/100)) > len(faster) or (float(n)*(mcp/100)) > len(medium) or (float(n)*(scp/100)) > len(slower):
        (n1,n2,n3) = resize_the_cell_groups(n,fcp,mcp,scp,faster,medium,slower)
    else:
        n1 = int(float(n) * (float(fcp)/float(100)))
        n2 = int(float(n) * (float(mcp)/float(100)))
        n3 = int(float(n) * (float(scp)/float(100)))
    #print "Size:", len(faster),len(medium),len(slower)
    # Sampling each group seperately
    #print "NNN:",n
    sampled_ids1 = random.sample(faster,n1)
    sampled_ids2 = random.sample(medium,n2)
    sampled_ids3 = random.sample(slower,n3)    
    
    # merging the arrays
    for i in xrange(0,len(sampled_ids2)):
        sampled_ids1.append(sampled_ids2[i])
    for i in xrange(0,len(sampled_ids3)):
        sampled_ids1.append(sampled_ids3[i])
    
    return sampled_ids1
"""
###############################################################################################    
# This is to re assign the number of cells to be sampled from each group.
###############################################################################################
"""
def resize_the_cell_groups(n,fa,me,sl,faster,medium,slower):
    n1 = int(float(n) * (float(fa)/float(100)))
    n2 = int(float(n) * (float(me)/float(100)))
    n3 = int(float(n) * (float(sl)/float(100)))
    #print "Before:",n1,n2,n3
    #print len(faster),len(medium),len(slower)
    if n1 > len(faster):
       d = n1 - len(faster)
       n1 = n1 - d
       n2 = n2 + d
    if n2 > len(medium):
       d = n2 - len(medium)
       n2 = n2 - d
       # Decide which cells to add
       if n1 > n3:
        n1 = n1 +d
       else:
        n3 = n3 + d
    if n3 > len(slower):
       d = n3 - len(slower)
       n3 = n3 - d
       n2 = n2 + d
    # print "After:",n1,n2,n3   
    return n1,n2,n3
"""
###############################################################################################    
# This is to calculate the mean cell divison time of the cells at each bottlneck.
###############################################################################################
def calculate_mean_cell_division_time(genotypes,output_file,i):
    # "i" is the bottlneck number.
    total_time = 0
    n_cells = 0
    for l in xrange(0,len(genotypes)):
        # Ignore the dead cells.
        if genotypes[l][2] != 0:
            total_time = total_time + genotypes[l][1]
            n_cells = n_cells + 1
    mean_division_time = total_time/n_cells
    print i,"\t",mean_division_time
    
    return mean_division_time
    
###################################################################################################
# Assigning the Haplotype structure.
###################################################################################################
def assign_haplotype(n):
        
    # Chromosome numbers and the length.
    chrom = {
             1:230208,2:813178,3:316617,4:1531918,
             5:576869,6:270148,7:1090946,8:562643,
             9:439885,10:745745,11:666454,12:1078175,
             13:924429,14:784333,15:1091289,16:948062
             }

    haplotype_structure = []
    
    for i in range(n):

        # Choose a chromosome randomly.
        chr_no = random.randint(1,16)

        # Choose a specific position in chromosome randomly.
        bp_position = random.randint(1,chrom[chr_no])
        haplotype_structure.append(str(chr_no) + str(":") + str(bp_position))
        
    return haplotype_structure
#################################################################################################### 