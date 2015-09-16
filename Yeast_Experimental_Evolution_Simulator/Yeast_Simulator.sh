# Compile the Cython codes
python setup.py build_ext --inplace 

echo "\nCompilation is Done\n"

echo "\nSimulating Selection Experiment\n"

# Simulator 
python Yeast_Simulator.pyx -m 0.004 -s 2 -l 33 -f 50 -b 50 >Simulator_Results.csv

echo "\nSelection Experiment is Done \n"

echo "\nPlotting the Results\n"

# Plotting the outputs
R CMD BATCH Simulator_Plots.R

echo "\n ...Finished... \n"
