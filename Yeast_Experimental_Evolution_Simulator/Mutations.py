from __future__ import division
import random
import sys
import math
import numpy
import scipy

###################################################################################################
############################### Import user defined Modules. ######################################
################################################################################################### 
import Sub_Functions as SF
###################################################################################################
################### Introducing Mutations Based on Yeast Mutation rate. ###########################
###################################################################################################
def introduce_mutation(mother_cell,daughter_cell,beneficial_mutation_rate,k,tp,n_mut,PM,PM_fitness,FAMR,mutation_fitness):
    """
    This function decides whether the mutation is fitness affecting or neutral.
    If the mutation is fitness affecting, whether its beneficial or deleterious.
    Whether mother cell or daughter cell will get the mutation
    
    And, t decides whether the mother or daughter cell should get the mutation.
    
    Input:
        Mother and daughter cell
        Proportions of fitness affecting mutations
        Proportions of beneficial mutations.
        Fitness effects
        Number of mutations
    
    Output:
        Mother or daughter cell with extra mutations
        Mother or daughter cell with new cell division time if the mutation is fitness affecting mutation
    """
    
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
    
    ######################################################
    # Choose the cell:Parent or child to have mutation.
    ######################################################
    if random.random() <= 0.5:
        # Assigning the Chromosome and the position.
        if fitness != 0:
            mother_cell[1] = SF.calculate_cell_division_time(mother_cell[1],fitness)
            
        PM_fitness[n_mut] = fitness
        mother_cell[3] += 1
        
        # Appending the mutation number
        PM[k].append(n_mut)
        
    else:
        if fitness != 0:
            daughter_cell[1] = SF.calculate_cell_division_time(daughter_cell[1],fitness)
        
        PM_fitness[n_mut] = fitness        
        daughter_cell[3] +=  1
       
        # Appending the mutation number
        PM[tp].append(n_mut)
    
    return mother_cell,daughter_cell,PM,PM_fitness