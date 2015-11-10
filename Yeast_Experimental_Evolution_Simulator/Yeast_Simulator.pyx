from __future__ import division
import sys
import optparse
import cProfile
import scipy.stats
import random
from copy import deepcopy
from numpy import *

"""
MAIN CONCEPT: Initial population contains n(Ex.n=10) number of yeast cells and each cell divides at 180 minutes.
When the time goes on, each cell divides at its particular time point.
Some cells grows faster(78 mins) and some are slower(102 mins). Cell growth will be stopped when the desird
population size is reached. And then n number of cells will be randomly sampled(bottleneck) and let it grow again.

# IMPORTANT VARIABLES #
genotypes - contains all the information about the each individual
mutations - All the mutation structure, beneficial_mutations - contains the IDs of beneficial mutations
Cell_mutations - cell id and the mutations it carries, gene_duplications - structures of gene duplications.

User defind functions
    import Sub_Functions
    import Cell_Division
"""

import Sub_Functions
import Cell_Division

"""
 ######################## Constant values used in the simualtions ##############################################
 
 n_individuals - Initial population size (e.g. 1000, 10000, 100000)
 number_of_generations - Number of mitotic dicisions before serial transfers
 number_of_bottlenecks - number of serial transfers (bottlenecks) in the whole simulation
 cell_division_min, cell_division_max - mimimum and maximum cell division time/generation time of each cell 
"""

n_individuals = 1000
number_of_generations = 5
number_of_bottlenecks = 50
cell_division_min = 90
cell_division_max = 180
cell_max_age = 16

# Population size before random sampling (bottleneck)
dps = n_individuals * 2**number_of_generations

"""
 How to run the script
 # Compilation
 python setup.py build_ext --inplace
 # Executing
 python Yeast_Simulator.pyx -m 0.004 -s 2 -l 33 -f 50 -b 50 
"""

#################################################################################################################
######################################### INPUT PART ############################################################
#################################################################################################################
parser = optparse.OptionParser("usage: %prog [options] -m 0.04 -s 2 -l 33 -f 50 -b 50")

# Mutation rate
parser.add_option("-m", "--mutation_rate", dest="mutation_rate",default="0",type="float",help="Mutation rate (genome/generation)")

# Shape of the gamma distribution
parser.add_option("-s", "--shape", dest="fitness_shape", default="0",type="int", help="shape value for fitness distribution")

# Scale of the Gamme distribution    
parser.add_option("-l", "--scale", dest="fitness_scale",default="0",type="int",help="Scale value for fitness distribution")

# Proportion of Beneficial Mutation rate in Fitness affecting mutation rate
parser.add_option("-f", "--fitness_affecting_mutations", dest="fitness_affecting_mutations", default="0",type="int", help="proportion of fitness affecting mutations")

# Proportion of Fitness affecting mutations
parser.add_option("-b", "--beneficial_mutations", dest="beneficial_mutations", default="0",type="int", help="proportion of beneficial mutations")

options, args = parser.parse_args()

mutation_rate = options.mutation_rate
mutation_shape = options.fitness_shape
mutation_scale = options.fitness_scale
fitness_affecting_mutation_rate = options.fitness_affecting_mutations
beneficial_mutation_rate = options.beneficial_mutations

##################################################################################################################
###################################### Global variables ##########################################################
##################################################################################################################
# Genotypes is a dict which contains all the details of each cell.
genotypes = {} 

# Contains all the mutations
mutations = []

# Contains beneficial mutation's fitness.
# Index will be the same for beneficial_mutations variable.
mutations_fitness = []

# Cell ID and the beneficial Mutation IDs
cell_beneficial_mutations = []

# Contains cell id and the gene deletion IDs
cell_gene_deletions = []
n_individuals_with_mutation = 0

# Cell group based on their cell division time.
cell_groups = {}

cell_mutations = array([])

#############################################################################################################
#   Initializing the Fitnesss effects and haplotype of genetic variants
#############################################################################################################
mutation_fitness_random_fitness = scipy.stats.gamma.rvs(mutation_shape,loc=0,scale=float(1)/float(mutation_scale),size=dps*0.3* number_of_bottlenecks)
mutation_fitness_random_fitness = Sub_Functions.covert_to_minutes(mutation_fitness_random_fitness)
#mutation_fitness_random_fitness = Sub_Functions.truncate_fitness_effects(mutation_fitness_random_fitness,89)
mutation_fitness  = Sub_Functions.make_fitness_proportions(mutation_fitness_random_fitness,beneficial_mutation_rate)

#############################################################################################################
# Main Program: Contains selection experiment simulation as well as Automated Backcrossings
#############################################################################################################
def Yeast_lab():
    
    # Initialize the cell Population with given cell division time and the total populations
    genotypes,PM = Sub_Functions.initialize_population(n_individuals,number_of_generations,cell_division_min,cell_division_max,cell_max_age)    
    
    # Grouping the cells based on their cell division time.
    cell_groups =  Sub_Functions.group_cells(genotypes,n_individuals)
    
    # Local variables for the subroutine.     
    mutations_fitness = {}; PM_strand = {};
    
    # Total number of mutation events.
    n_mut = 0;

    # Calling the function to allow the cells to divide. Mutations, gene deletions and gene duplications will be introduced during cell division.
    (genotypes,PM,PM_fitness,n_mut) = Cell_Division.asymmetrical_cell_division(genotypes,cell_groups,
                    mutations_fitness,PM,number_of_generations,number_of_bottlenecks,n_individuals,n_mut,dps,fitness_affecting_mutation_rate,mutation_fitness,mutation_rate,beneficial_mutation_rate,cell_max_age)
                    
#############################################################################################################
#############################################################################################################

# This is to run this module normally like "python Yeast_Simulator.py"
if __name__ == "__main__":
    import sys
    Yeast_lab()
