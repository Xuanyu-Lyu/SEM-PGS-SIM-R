

#READ ME - GE-75
Nov 16, 2023

This version of GeneEvolve will attempt to redo the VT model in the following ways:
* Write a new VT model in GE that allows for cascade-like transmission
* check why i standardized the parental phenotype being passed down and change this if possible
* allow for different maternal and paternal VT paths



Here are what each script does:
GE-75.R - the master script where users change the parameters
GE.Functions_Short.R - 8 different short helper functions




Here are important objects
PAR1 - a list of all the options input by the user
VAR1 - a list of variance parameter values input by the user
ALL - the output from run.gene.evolve

run.gene.evolve - the main function that runs the whole simulation

