## Coarse-grained Simulation Toolkits for Protein Folding
#### Author: Dr. Yang Jiang

This is a package of scripts that are used to create coarse-grained (CG) models of proteins/ribosomes, optimize CG force field parameters and run MD simulations for protein co- and post-translational folding. All the scripts are ready to use when users have added the directories in `$PATH` and have granted the execution permission (`chmod +x`) for all the scripts. 

### Table of Contents
  * [1. Create CG protein models and tune the force field parameters (*n*<sub>scale</sub>) for a given protein](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscale-for-a-given-protein)
    + [1.1. Tune *n*<sub>scale</sub> for a small single-domain protein that has experimental folding stability reported](#11-tune-nscale-for-a-small-single-domain-protein-that-has-experimental-folding-stability-reported)
    + [1.2. Tune *n*<sub>scale</sub> for a protein without experimental folding stability reported](#12-tune-nscale-for-a-protein-without-experimental-folding-stability-reported)
  * [2. Temperature quenching simulation](#2-temperature-quenching-simulation)
  * [3. Create CG ribosome model](#3-create-cg-ribosome-model)
  * [4. Simulation of co-translational folding](#4-simulation-of-co-translational-folding)
  * [5. Simulation of post-translational folding](#5-simulation-of-post-translational-folding)
  * [6. Backmapping from coarse-grained model to all-atom model](#6-backmapping-from-coarse-grained-model-to-all-atom-model)
  * [7. Analysis of protein folding trajectories](#7-analysis-of-protein-folding-trajectories)


### 1. Create CG protein models and tune the force field parameters (*n*<sub>scale</sub>) for a given protein
To be able to create a CG Go-like model for a given protein, you need to get the scripts from Dr. Ed O'Brien lab first, including the modified Charmm executable and MMTSB package. Making the CG protein modeling work without Charmm and MMTSB is one of my future development for this toolkit.
#### 1.1. Tune *n*<sub>scale</sub> for a small single-domain protein that has experimental folding stability reported
- For a small single-domain protein that has experimental folding stability reported, we use an enhanced sampling protocol i.e. replica exchange simulations to extensively sample the protein folding and unfolding and identify the *n*<sub>scale</sub> value to reproduce the experimental protein folding stability.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**create_cg_protein_model_v34_0.37_nbx3.pl** | Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be only used to build CG model for a single domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more) |
| CG_protein_parameterization/**parse_cg_prm.py** | Parse the parameters in .prm file and then create a .xml file for OpenMM use. (Learn more) |
| CG_protein_parameterization/**parallel_temperature_REX.py** | Run parallel temperature replica exchange molecular dynamics (pt-REMD) simulation. This simulation is parallelized using multiple CPU processors. ([Learn more](../wikis/help_wiki/parallel_temperature_REX.py)) |
| CG_protein_parameterization/**opt_temp.pl** | Optimize the temperature windows for pt-REMD simulation to ensure the good sampling quality around the melting temoerature of the given protein. (Learn more) |
| CG_protein_parameterization/**check_sampling.pl** | Check the sampling quality for pt-REMD simulation. Insufficient sampling will cause problems and inaccuracy in estimating the protein folding stability. (Learn more) | 
| CG_protein_parameterization/**analysis_folding_stability.pl** | Estimate the protein folding stability at a given temperature from pt-REMD data using WHAM. (Learn more) |
| CG_protein_parameterization/**scan_nscal_nbx_3_REX.pl** | An automated script to call `create_cg_protein_model_v34_0.37_nbx3.pl`, `opt_temp.pl`, `parallel_temperature_REX.py` and `analysis_folding_stability.pl` to scan the protein folding stability profile as changing *n*<sub>scale</sub> value. The protein folding stability profile will be used to find the optimized *n*<sub>scale</sub> value for CG model parameterization. (Learn more) | 

- To tune *n*<sub>scale</sub> for a given protein on PSU ACI cluster, run `scan_nscal_nbx_3_REX.pl` with a series of *n*<sub>scale</sub> values. Use `scan_nscal_nbx_3_REX.pl` again to analyze the results and get the protein stability profile. The optimized *n*<sub>scale</sub> is the value that reproduces the experimental protein folding stability.

#### 1.2. Tune *n*<sub>scale</sub> for a protein without experimental folding stability reported
- For a protein without experimental folding stability reported, we use a stepwise optimization strategy to tune *n*<sub>scale</sub> according to 5 levels of *n*<sub>scale</sub>. The optimized *n*<sub>scale</sub> of an entire single-domain protein or one domain/interface of a multi-domain protein is identified as the lowest level that maintains the native structure. 
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**create_cg_protein_model_v34_nbx3_multidomain.pl** | Create the CG model .psf .top and .prm file that can be used for MD simulations. This script can be used to build CG model for a multi-domain protein. Need to get the modified Charmm and MMTSB installed prior to use. (Learn more) | 
| CG_protein_parameterization/**opt_nscal.pl** | An automated script to find the optimized *n*<sub>scale</sub> for each domain/interface according to 5 levels of *n*<sub>scale</sub> values trained by running the protocol in [Section 1.1](#11-tune-nscale-for-a-small-single-domain-protein-that-has-experimental-folding-stability-reported) for 18 small single-domain proteins. The MD simulator is OpenMM, which is different from that in script `opt_nscal_charmm.pl.pl`. (Learn more) |
| CG_protein_parameterization/**opt_nscal_charmm.pl.pl** | This script has the same function with `opt_nscal.pl` but the levels of *n*<sub>scale</sub> values used in this script was trained by using Charmm for the same protein set. | 

- On PSU ACI cluster, to tune *n*<sub>scale</sub> that is compatible with OpenMM, run `opt_nscal.pl`; To tune *n*<sub>scale</sub> that is compatible with Charmm, run `opt_nscal_charmm.pl`.

### 2. Temperature quenching simulation
- To estimate the folding rates for a given protein, you need to run temperature quenching simulation where the system is first heated to a very high temperature (usually 800 K) quickly to make protein totally unfolded and then cooled down to the physiological temperature to moniter the refolding of the protein.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**temperature_quenching.py** | Run temperature quenching simulation from the CG native structure of a given protein. This simulation is parallelized using multiple CPU processors. ([Learn more](../wikis/help_wiki/temperature_quenching.py)) |
| CG_protein_parameterization/**T_quench_nbx_3.pl** | An automated script to build CG model from a pdb file and then run temperature quenching simulations. (Learn more) | 
| CG_protein_parameterization/**analysis_Tq.pl** | Analyze the results of temperature quenching simulations, fit a single- or double- exponential function to the survival probability of the unfolded protein and then estimate the folding rate.  (Learn more) | 

- To estimate the protein folding rate on PSU ACI cluster, you need to run `T_quench_nbx_3.pl` with optimized *n*<sub>scale</sub> values obtained from [Section 1](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscale-for-a-given-protein) and then run `analysis_Tq.pl` to do the curve fitting.

### 3. Create CG ribosome model
- The CG ribosome model is used to run continuous synthesis of a nascent chain that is parameterized with the CG protein model. The CG ribosome is usually fixed (not allow to move) during the simulation. The force field parameters thus only contain the nonbonding term. To speedup the computation, we usually truncate the ribosome to only contain the tails of P- and A-site tRNA molecules, the entire exit tunnel with a few atoms near the tunnel wall and the surface near the exit that may have contacts with the nascent chain.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_ribosome_parameterization/**gen_50S_pdb.py** | Get he necessary subunits from .cif file of the ribosome and create a pdb file contains those sbuunits. (Learn more) |
| CG_ribosome_parameterization/**fix_orein_50S_pdb.py** | Add missing atoms and rotate/translate the subunits to a desired orientation for E. coli ribosome. (Learn more) |
| CG_ribosome_parameterization/**fix_orein_60S_pdb.py** | Do the same thing with `fix_orein_50S_pdb.py` for S. cerevisiae ribosome. (Learn more) |
| CG_ribosome_parameterization/**create_cg_ribosome_model.py** | Create the CG model for those re-orientated subunits, including .psf, .top, .cor and .prm files. (Learn more) |
| CG_ribosome_parameterization/**truncate_ribosome.py** | Truncate the CG ribosome according to the centroid line of the exit tunnel. (Learn more) |
| CG_ribosome_parameterization/**gen_ribosome_FF_v2.py** | Train collision diameters of CG ribosome beads from given ribosome structures. (Learn more) |

- To create CG ribosome model, you need to download the .cif file of your ribosome from [RCSB PDB](https://www.rcsb.org/). Use `gen_50S_pdb.py` to generate a pdb file that contains the assembly of subunits that you want to model. Note that the 50S in the script name doesn't restrict the scope of usage. You can use it to generate any assembly of subunits you want, no matter what the organism the ribosome belongs to. Then use `fix_orein_50S_pdb.py` for E. coli ribosome and `fix_orein_60S_pdb.py` for S. cerevisiae ribosome to get the re-orientated subunits, followed by using `create_cg_ribosome_model.py` to build the CG model for the assembly of subunits. Finally, use `truncate_ribosome.py` to generate the truncated CG ribosome .cor file. If you are interested in tunning force field parameters for CG ribosome model, use `gen_ribosome_FF_v2.py` to train your parameters using some high-resolution crystal structures or Cryo-EM structures. 

### 4. Simulation of co-translational folding
- The simulation protocol for co-translational folding is called "Continuous Synthesis Protocol" (CSP). We add new amino acid on the A-site tRNA in the truncated CG ribosome and simulate the tRNA translocation and nascent chain elongation at a frequency based on in vivo codon translation time.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Continuous_synthesis_protocol/**continuous_synthesis_v6.py** | Run continuous synthesis of a CG protein on a CG ribosome. Both parallelizations on CPU and GPU are supported. (Learn more) <br>Scripts needed: `Continuous_synthesis_protocol/ribosome_traffic` and `CG_protein_parameterization/parse_cg_prm.py`. |
| Continuous_synthesis_protocol/**ribosome_traffic** | Estimate the real codon translation time by taking into account of the ribosome traffic effects. (Learn more) | 
| Continuous_synthesis_protocol/**visualize_cont_synth.py** | Generate movies of the continuous synthesis process. (Learn more) <br>Scripts needed: `Backmapping/backmap.py`, `Continuous_synthesis_protocol/render_ecoli_RNC.tcl` and `Continuous_synthesis_protocol/render_yeast_RNC.tcl` | 
| Continuous_synthesis_protocol/**render_ecoli_RNC.tcl** | Render the picture of E.coli ribosome-nascent-chain (RNC) complex in VMD. (Learn more) | 
| Continuous_synthesis_protocol/**render_yeast_RNC.tcl** | Render the picture of S. cerevisiae ribosome-nascent-chain (RNC) complex in VMD. (Learn more) |

- To run CSP, you need to prepare the CG protein model for your nascent chain according to [Section 1](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscale-for-a-given-protein) and prepare the CG ribosome model according to [Section 3](#3-create-cg-ribosome-model). All the .psf, .top, .cor and .prm files are required in initialization of CSP. In addition, users have to provide a table of intrinsic codon translation time and the mRNA sequence of the nascent chain.

### 5. Simulation of post-translational folding
- The simulation of post-translational folding is quite simple. The nascent chain structures that are disassociated from the ribosome will be obtained from [CSP](#4-simulation-of-co-translational-folding) and run Langevin dynamics (LD) in implicit water environment at the physiological temperature.
- Scripts will be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Post_translational_folding/**post_trans_single_run_v2.py** | Run a single trajectory of post-translational folding. User can specify a walltime or a threshold to control the termination of the post-translational folding. (Learn more) |
| Post_translational_folding/**PTP_setup_v3.pl** | Automated script to setup post-translation simulations after [CSP](#4-simulation-of-co-translational-folding). (Learn more) <br>Scripts needed: `Post-translational_folding/post_trans_single_run_v2.py` |

- To setup and run post-translation simulations on PSU ACI cluster, use `PTP_setup_v3.pl`.

### 6. Backmapping from coarse-grained model to all-atom model

### 7. Analysis of protein folding trajectories