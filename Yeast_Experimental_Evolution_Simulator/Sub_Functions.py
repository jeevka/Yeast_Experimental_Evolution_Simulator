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
import Mutations

################################## Intialization ####################################################
# Intialization of the individuals before the cell division starts.
# 1. Initial population contains cell at different ages. So the cells divide at different speeds.
# 2. 50% cells are virgin cells. 25% are one time divided. 12.5% are two times divided. and etc ...
#####################################################################################################
def initialize_population(n_individuals,n_generations,cell_division_min,cell_division_max,cell_max_age):
    """
    Input:
        cell division time
        number of individuals
        number of generations
        maximum age of the cells
    
    Output:
        Dictionary of lists with basic properties
        e.g. genotypes[n] = [cell_division_max,cell_division_max,cell_max_age,number of mutations]
        e.g. genotype[1] = [180,180,16,0]
    """
    
    size = n_individuals * 2**n_generations
    genotypes = ones((size,4),dtype=int32)
    
    #import random
    population_ratio = {}
    population_ratio[0] = 0
    population_percentage = n_individuals
    #genotypes = {}
    cell_mutations = []
    cell_gene_deletions = []
    cell_gene_duplications = []
    a = []
    
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
    """
    Input:
        dictionary of lists (cells)
        % cells to be dead or wont survive
        
    Output:
        dictionary of lists (cells) with inavtivated cells
    """
    inactive_ids = {}
    n_inactive_cells = int(len(genotypes) * (percentage/100))
    inactive_ids = random.sample(genotypes,n_inactive_cells)    
    for i in range(0,len(inactive_ids)):
        id = inactive_ids[i]
        genotypes[id][2] = 0
    return genotypes

#############################################################################################
# Group the cells based on their cell division time.
# So the cells can divide at the right time
#############################################################################################
def group_cells(genotypes,n_individuals):
    """
    Input:
        dictionary of lists (cells)
        Number of individuals
    
    Output:
        grouped cells
    """
    cell_group = {}
    n_individuals = int(n_individuals)
    for i in xrange(n_individuals):
        if cell_group.has_key(genotypes[i][0]):
            cell_group[genotypes[i][0]].append(i) 
        else:
            cell_group[genotypes[i][0]] = [i]
            
    return cell_group

###############################################################################################    
# This is to calculate the mean cell divison time of the cells at each bottlneck.
###############################################################################################
def calculate_mean_cell_division_time(genotypes,output_file,i):
    """
    Input:
        Dictionary of lists with cell division time
    
    Output:
        Mean cell division time of the population
    """
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
    """
    Input:
        Number to create haplotypes
    
    Output:
        haplotype. E.g. chromsome number and position of the mutation
    """
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
# Convert the fitness values to beneficial and deleterious bases on their proportions
####################################################################################################
def make_fitness_proportions(fitness,MP):
    """
    Input:
        List cell division times
        Proportion of beneficial mutations
    
    Output:
        List of fitness values converted to beneficial or deleterious
    """
    import random as r
    n1 = int(len(fitness) * (100-MP)/100)
    n2 = len(fitness)
    pop = xrange(n2)
    sampled = r.sample(pop,n1)
    for i in sampled:
        fitness[i] *= -1

    return fitness

####################################################################################################
# Convert the fitness values to minutes from random gamma values
####################################################################################################
def covert_to_minutes(fitness):
    """
    Input:
        List of random values from gamma distribution
    
    Output:
        gamma distribution values are converted to mins
    """
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_1(i))
    
    return fitness_mins

####################################################################################################
# Convert the fitness values to minutes
####################################################################################################
def calculate_cell_division_time_1(fitness):
    """
    Input:
        one random value from Gamma distribution
    
    Output:
        converted value back
    """
    cell_division_max = 180
    cell_division_min = 90
    #change_CDT = (cell_division_max-cell_division_min) * alpha
    
    change_CDT = (fitness/(1+fitness)) *  cell_division_max

    if change_CDT > 90:
        change_CDT = 89
        
    return change_CDT

####################################################################################################
# Values bigger than the difference between WT and Mutatant cant be included due to the design of our model
####################################################################################################
def truncate_fitness_effects(fitness,trunc):
    """
    Input:
        List of fitness values
    
    Output:
        Truncated values
        e.g. all the values >89 will be removed
    """
    new_fitness = []
    for i in fitness:
        if i <= trunc:
            new_fitness.append(i)
    
    return new_fitness

####################################################################################################
# Calculating new cell division times based on current cell division time and fitness
####################################################################################################
def calculate_cell_division_time(current_CDT,alpha):
    """
    Input:
        current cell division time of a cell
        Fitness effect of the new mutation
    
    Output:
        New cell division time calculated based on exisiting cell division time
        and fitness effect of the new mutation
    """
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 180
    cell_division_min = 90
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)

    return  new_CDT