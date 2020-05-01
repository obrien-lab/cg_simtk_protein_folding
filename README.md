## Coarse-grained Simulation Toolkits for Protein Folding
#### Author: Dr. Yang Jiang

### Table of Contents


### 1. Introduction
This is a package of scripts that are used to create CG models of proteins/ribosomes, optimize CG force field parameters and run MD simulations for protein co- and post-translational folding. All the scripts are ready to use when users have added the directories in `$PATH` and have granted the execution permission (`chmod +x`) for all the scripts. 

### 2. Create CG protein models and tune the force field parameters (*n*<sub>scale</sub>) for a given protein
To be able to create a CG Go-like model for a given protein, you need to get the scripts from Dr. Ed O'Brien lab first, including the modified Charmm executable and MMTSB package. Making the CG protein modeling work without Charmm and MMTSB is one of my future development for this toolkit.
#### 2.1. Tune *n*<sub>scale</sub> for a small single-domain protein that has experimental folding stability reported
- For a small single-domain protein that has experimental folding stability reported, we use an enhanced sampling protocol i.e. replica exchange simulations to extensively sample the protein folding and unfolding and identify the *n*<sub>scale</sub> value to reproduce the experimental protein folding stability.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**create_cg_protein_model_v34_0.37_nbx3.pl** | Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be only used to build CG model for a single domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more) |
| CG_protein_parameterization/**parallel_temperature_REX.py** | Run parallel temperature replica exchange molecular dynamics (pt-REMD) simulation. This simulation is parallelized using multiple CPU processors. ([Learn more](../wikis/help_wiki/parallel_temperature_REX)) |
| CG_protein_parameterization/**opt_temp.pl** | Optimize the temperature windows for pt-REMD simulation to ensure the good sampling quality around the melting temoerature of the given protein. (Learn more) |
| CG_protein_parameterization/**check_sampling.pl** | Check the sampling quality for pt-REMD simulation. Insufficient sampling will cause problems and inaccuracy in estimating the protein folding stability. (Learn more) | 
| CG_protein_parameterization/**analysis_folding_stability.pl** | Estimate the protein folding stability at a given temperature from pt-REMD data using WHAM. (Learn more) |
| CG_protein_parameterization/**scan_nscal_nbx_3_REX.pl** | An automated sript to call `create_cg_protein_model_v34_0.37_nbx3.pl`, `opt_temp.pl`, `parallel_temperature_REX.py` and `analysis_folding_stability.pl` to scan the protein folding stability profile as changing *n*<sub>scale</sub> value. The protein folding stability profile will be used to find the optimized *n*<sub>scale</sub> value for CG model parameterization. (Learn more) | 

- To tune *n*<sub>scale</sub> for a given protein, run `scan_nscal_nbx_3_REX.pl` with a series of *n*<sub>scale</sub> values. Use `scan_nscal_nbx_3_REX.pl` again to analyze the results and get the protein stability profile. The optimized *n*<sub>scale</sub> is the value that reproduces the experimental protein folding stability.

#### 2.2. Tune *n*<sub>scale</sub> for a protein without experimental folding stability reported
- For a protein without experimental folding stability reported, we use a stepwise optimization strategy to tune *n*<sub>scale</sub> according to 5 levels of *n*<sub>scale</sub>. The optimized *n*<sub>scale</sub> of an entire single-domain protein or one domain/interface of a multi-domain protein is identified as the lowest level that maintains the native structure. 
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**create_cg_protein_model_v34_nbx3_multidomain.pl** | Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be used to build CG model for a multi-domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more) | 
| CG_protein_parameterization/**opt_nscal.pl** | An automated sript to find the optimized *n*<sub>scale</sub> for each domain/interface according to 5 levels of *n*<sub>scale</sub> values trained by running the protocol in [Section 2.1](#21-tune-nscale-for-a-small-single-domain-protein-that-has-experimental-folding-stability-reported) for 18 small single-domain proteins. The MD simulator is OpenMM, which is different from that in script `opt_nscal_charmm.pl.pl`. (Learn more) |
| CG_protein_parameterization/**opt_nscal_charmm.pl.pl** | This script has the same function with `opt_nscal.pl` but the levels of *n*<sub>scale</sub> values used in this script was trained by using Charmm for the same protein set. | 

### 3. Temperature quenching simulation
- To estimate the folding rates for a given protein, you need to run temperature quenching simulation where the system is first heated to a very high temperature (usually 800 K) quickly to make protein totally unfolded and then cooled down to the physiological temperature to moniter the refolding of the protein.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**temperature_quenching.py** | Run temperature quenching simulation from the CG native structure of a given protein. This simulation is parallelized using multiple CPU processors. (Learn more) |
| CG_protein_parameterization/**T_quench_nbx_3.pl** | An automated script to build CG model from a pdb file and then run temperature quenching simulations. (Learn more) | 
| CG_protein_parameterization/**analysis_Tq.pl** | Analyze the results of temperature quenching simulations, fit a single- or double- exponential function to the survival probability of the unfolded protein and then estimate the folding rate.  (Learn more) | 

### 4. Create CG ribosome model

### 5. Simulation of co-translational folding

### 6. Simulation of post-translational folding

### 7. Backmapping from coarse-grained model to all-atom model