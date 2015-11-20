###################################################################################################################
############################################# MAIN PROGRAM ########################################################
###################################################################################################################
# MAIN CONCEPT: Initial population contains n(Ex.n=10) number of yeast cells and each cell divides at different
# time point (Ex. bw 78-108 mins). When the time goes on, each cell divides at its particular time point.
# Some cells grows faster(78 mins) and some are slower(102 mins). Cell growth will be stopped when the desird
# population size is reached. And then n number of cells will be randomly sampled(bottleneck) and let it grow again.
###################################################################################################################
########################################### IMPORTANT VARIABLES ###################################################
# Genotypes structure goes like this : {id:[next division time, Cell division time, No. of generations,}
# genotypes - contains all the information about the each individual
# mutations - All the mutation structure, beneficial_mutations - contains the IDs of beneficial mutations
# Cell_mutations - cell id and the mutations it carries, gene_duplications - structures of gene duplications.
# cell_gene_duplications - cell id and the gene duplications, gene_deletions - structures of gene deletions
# cell_gene_deletions - cell id and the gene deletions ids, mutations_fitness - mutation id and fitness
# gene_duplication_fitness - gene duplication ID and fitness, gene_deletion_fitness - gene deletion id and fitness.
####################################################################################################################
####################################################################################################################
from __future__ import division
import sys
import optparse
import cProfile
import scipy.stats
import random
from copy import deepcopy
from numpy import *


###################################################################################################################
################################# User defind modules #############################################################
###################################################################################################################
import Sub_Functions
import Cell_Division
import Time_2_Fix_Ext

###################################################################################################################
############################## Global Variables - Constants #######################################################
###################################################################################################################
# Refer the coding documents for all the values use here.
number_of_bottlenecks = 1000

# Number of generations before bottleneck.
# Which means the average number of cell divisions
# for the cells in initial population.
number_of_generations = 5

# Mean cell division time of the WT cells in normal medium
cell_division_min = 90

# Maximum cell divison time of the cells in Stressor medium.
cell_division_max = 180

# Maximum Number of times a cell can divide.
cell_max_age = 16

# Number of chromosomes.
n_chr = 16

# Population size
n_individuals = 10

# Population size before random sampling (bottleneck)
dps = n_individuals * 2**number_of_generations

####################################################################################################################
########################################### INPUT Values ###########################################################
####################################################################################################################

# ID for the simulation (e.g. 1)
#Task_id = sys.argv[1]
# Mutation rate
# mutation_rate =  float(sys.argv[1])
# Shape of the gamma distribution
# mutation_shape = float(sys.argv[2])
# Scale of the Gamme distribution
# mutation_scale = float(sys.argv[3])
# Proportion of Beneficial Mutation rate in Fitness affecting mutation rate
# beneficial_mutation_rate = float(sys.argv[4])
# Proportion of Fitness affecting mutations
# fitness_affecting_mutation_rate = float(sys.argv[5])


parser = optparse.OptionParser("usage: %prog [options] -m 0.04 -s 2 -l 33 -f 50 -b 50")

# How to run
# Compilation
# python setup.py build_ext --inplace
# Executing
# python Yeast_Simulator.pyx -m 0.004 -s 2 -l 33 -f 50 -b 50 

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

# Contains beneficial Mutations
beneficial_mutations = []

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

# Total number of mutation events.
n_mut = 0

cell_mutations = array([])

#############################################################################################################
# Subprograms
#############################################################################################################
def make_fitness_proportions(fitness,MP):
    import random as r
    n1 = int(len(fitness) * (100-MP)/100)
    n2 = len(fitness)
    pop = xrange(n2)
    sampled = r.sample(pop,n1)
    for i in sampled:
        fitness[i] *= -1

    return fitness

def covert_to_minutes(fitness):
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_1(i))
    
    return fitness_mins

def calculate_cell_division_time_1(fitness):
    cell_division_max = 180
    cell_division_min = 90
    #change_CDT = (cell_division_max-cell_division_min) * alpha
    
    change_CDT = (fitness/(1+fitness)) *  cell_division_max

    if change_CDT > 90:
        change_CDT = 89
        
    return change_CDT

def truncate_fitness_effects(fitness,trunc):
    new_fitness = []
    for i in fitness:
        if i <= trunc:
            new_fitness.append(i)
    
    return new_fitness

def calculate_cell_division_time(current_CDT,alpha):
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

#############################################################################################################
#   Initializing the Fitnesss effects and haplotype of genetic variants
#############################################################################################################
mutation_fitness_random_fitness = scipy.stats.gamma.rvs(mutation_shape,loc=0,scale=float(1)/float(mutation_scale),size=dps*0.3* number_of_bottlenecks)
mutation_fitness_random_fitness = covert_to_minutes(mutation_fitness_random_fitness)
#mutation_fitness_random_fitness = truncate_fitness_effects(mutation_fitness_random_fitness,89)
mutation_fitness  = make_fitness_proportions(mutation_fitness_random_fitness,beneficial_mutation_rate)

#############################################################################################################
# Main Program: Contains selection experiment simulation as well as Automated Backcrossings
#############################################################################################################
def Yeast_lab():
    import Sub_Functions
    
    # Initialize the cell Population
    genotypes,PM = Sub_Functions.initialize_population(n_individuals,number_of_generations,cell_division_min,cell_division_max,cell_max_age)    

    # Grouping the cells based on their cell division time.
    cell_groups =  Sub_Functions.group_cells(genotypes,n_individuals)
    
    # de-activate XX% of cells from intial population.
    # genotypes = Sub_Functions.inactivate_cells(genotypes,cell_wont_survive)
    
    # Local variables for the subroutine.     
    mutations_fitness = {}; PM_strand = {};
    
    import Cell_Division
    n_mut = 0;
    
    s_size = n_individuals

    # Calling the function to allow the cells to divide. Mutations, gene deletions and gene duplications will be introduced during cell division.
    (genotypes,PM,PM_fitness,n_mut) = Cell_Division.asymmetrical_cell_division(genotypes,cell_groups,
                    mutations_fitness,PM,number_of_generations,number_of_bottlenecks,n_individuals,n_mut,s_size,dps,fitness_affecting_mutation_rate,mutation_fitness,mutation_rate,beneficial_mutation_rate,cell_max_age)
                    
    Time_2_Fix_Ext.Analytical_Test(n_mut)
    
#############################################################################################################
#############################################################################################################

# This is to run this module normally like "python Yeast_Simulator.py"
if __name__ == "__main__":
    import sys
    Yeast_lab()
