# Coarse-grained Simulation Toolkits for Protein Folding
#### Author: Dr. Yang Jiang

### 1. Introduction
This is a package of scripts that are used to create CG models of proteins/ribosomes, optimize CG force field parameters and run MD simulations for protein co- and post-translational folding. All the scripts are ready to use when users have added the directories in `$PATH` and have granted the execution permission (`chmod +x`) for all the scripts. 

### 2. Create CG protein models and tune the force field parameters (*n*<sub>scale</sub>) for a given protein
*  To be able to create a CG Go-like model for a given protein, you need to get the scripts and from Dr. Ed O'Brien lab, including the modified Charmm executable and MMTSB package. Making the CG protein modeling work without Charmm and MMTSB is one of my future development for this toolkit.
*  