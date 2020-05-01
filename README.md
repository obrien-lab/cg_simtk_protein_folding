# Coarse-grained Simulation Toolkits for Protein Folding
#### Author: Dr. Yang Jiang

### 1. Introduction
This is a package of scripts that are used to create CG models of proteins/ribosomes, optimize CG force field parameters and run MD simulations for protein co- and post-translational folding. All the scripts are ready to use when users have added the directories in `$PATH` and have granted the execution permission (`chmod +x`) for all the scripts. 

### 2. Create CG protein models and tune the force field parameters (*n*<sub>scale</sub>) for a given protein
- To be able to create a CG Go-like model for a given protein, you need to get the scripts from Dr. Ed O'Brien lab first, including the modified Charmm executable and MMTSB package. Making the CG protein modeling work without Charmm and MMTSB is one of my future development for this toolkit.

- Scripts will be used in this section:
  - CG_protein_parameterization/**create_cg_protein_model_v34_0.37_nbx3.pl**
    
    Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be only used to build CG model for a single domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more)
  - CG_protein_parameterization/**create_cg_protein_model_v34_nbx3_multidomain.pl**
    
    Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be used to build CG model for a multi-domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more)
  - CG_protein_parameterization/**parallel_temperature_REX.py**
    
    Run parallel temperature replica exchange molecular dynamics (pt-REMD) simulation. This simulation is parallelized using multiple CPU processors. (Learn more)
  - CG_protein_parameterization/**opt_temp.pl**
  
    Optimize the temperature windows for pt-REMD simulation to ensure the good sampling quality around the melting temoerature of the given protein. (Learn more)
  - CG_protein_parameterization/**check_sampling.pl**
  
    Check the sampling quality for pt-REMD simulation. Insufficient sampling will cause problems and inaccuracy in estimating the protein folding stability. (Learn more)
