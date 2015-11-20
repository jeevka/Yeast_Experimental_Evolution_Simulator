from __future__ import division
import random
import sys
import math
import numpy
import scipy

# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

###################################################################################################
############################### Import user defined Modules. ######################################
################################################################################################### 
import Yeast_Simulator as YS

###################################################################################################
###################### Global variables - Importing from main Modules. ############################
###################################################################################################
#fitness_affecting_mutation_rate = YS.fitness_affecting_mutation_rate
#mutation_fitness = YS.mutation_fitness_random_fitness

###################################################################################################
################### Introducing Mutations Based on Yeast Mutation rate. ###########################
###################################################################################################
def introduce_mutation(mother_cell,daughter_cell,beneficial_mutation_rate,k,tp,n_mut,PM,PM_fitness,FAMR,mutation_fitness):
    beneficial = 0
    deleterious = 0
    fitness_proportion = 0
    
    # Decide whether its fitness altering mutation or neutral.
    if random.randrange(100) < FAMR:
        # Decide Whether its Adavantageous or Disadvantageous or neutral mutation.
        if random.random() <= beneficial_mutation_rate:
            beneficial = 1
        else:
            deleterious = 1
    
    # fix_cal decides whether GDs carry any fitness effects.
    if random.randrange(100) < FAMR:        
        fitness = mutation_fitness[n_mut]
    else:
        fitness = 0
    
    # Thsi is to check neutral model
    fitness = 0 
    ######################################################
    # Choose the cell:Parent or child to have mutation.
    ######################################################
    ran_n = c_libc_random() % 0.99 + 0
    if ran_n <= 0.5:
        # Assigning the Chromosome and the position.
        if fitness != 0:
            mother_cell[1] = YS.calculate_cell_division_time(mother_cell[1],fitness)
            
        PM_fitness[n_mut] = fitness
        mother_cell[3] += 1
        
        # Appending the mutation number
        PM[k].append(n_mut)
        
    else:
        if fitness != 0:
            daughter_cell[1] = YS.calculate_cell_division_time(daughter_cell[1],fitness)
        
        PM_fitness[n_mut] = fitness        
        daughter_cell[3] +=  1
       
        # Appending the mutation number
        PM[tp].append(n_mut)
    
    return mother_cell,daughter_cell,PM,PM_fitness
    
####################################################################################################    
"""
Reference:
Sarah B et al.  Spontaneous Mutations in Diploid Saccharomyces Cerevisiae: More Beneficial Than Expected. 2004 Genetics.
"""
####################################################################################################    
