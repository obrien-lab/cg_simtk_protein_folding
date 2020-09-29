## Coarse-grained Simulation Toolkit for Protein Folding
#### Author: [Dr. Yang Jiang](https://orcid.org/0000-0003-1100-9177)

This is a package of scripts that are used to create coarse-grained (CG) models of proteins/ribosomes, optimize CG force field parameters, run MD simulations for protein co- and post-translational folding and analyze the folding structures. This toolkit is developed based on the original simulation protocols and scripts in [O'Brien Lab](https://obrien.vmhost.psu.edu/). 

An example of the translation and folding process of the *E. coli* D-Ala-D-Ala Ligase B (DDLB, PDBid [4c5c](https://www.rcsb.org/structure/4C5C)) that was simulated by this toolkit is shown below.

<img src="../../wikis/uploads/img/DDLB_folding_example.gif" width="300">

All the scripts are ready to use when users have added the directories (including all the child folders) in `$PATH` (add `export PATH=${PATH}:/path/to/one/child/folder/` in your `~/.bashrc` file) and have granted the execution permission (`chmod -R +x ./cg_simtk_protain_folding/`) for all the scripts. 

 :warning: Some of the scripts need python modules (Python 3.X) and other softerware installed prior to use. Please click `Learn more` in the following script instruction tables to find detailed instruction of usage and basic theory used in the script.

---

### Table of Contents
  * [1. Create CG protein models and tune the force field parameters (*n*<sub>scal</sub>) for a given protein](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscal-for-a-given-protein)
    + [1.1. Tune *n*<sub>scal</sub> for a small single-domain protein that has experimental folding stability reported](#11-tune-nscal-for-a-small-single-domain-protein-that-has-experimental-folding-stability-reported)
    + [1.2. Tune *n*<sub>scal</sub> for a protein without experimental folding stability reported](#12-tune-nscal-for-a-protein-without-experimental-folding-stability-reported)
  * [2. Temperature quenching simulation](#2-temperature-quenching-simulation)
  * [3. Create CG ribosome model](#3-create-cg-ribosome-model)
  * [4. Simulation of co-translational folding](#4-simulation-of-co-translational-folding)
  * [5. Simulation of post-translational folding](#5-simulation-of-post-translational-folding)
  * [6. Backmapping from coarse-grained model to all-atom model](#6-backmapping-from-coarse-grained-model-to-all-atom-model)
  * [7. Analysis of protein folding trajectories](#7-analysis-of-protein-folding-trajectories)


### 1. Create CG protein models and tune the force field parameters (*n*<sub>scal</sub>) for a given protein
- The CG protein model used in the simulations of protein folding is a G&omacr;-like C<sub>&alpha;</sub> model. The force field can be found [here](../../wikis/help_wiki/create_cg_protein_model.py#51-force-field-of-the-c%CE%B1-model). Another CG model ([C<sub>&alpha;</sub>-sidechain model](../../wikis/help_wiki/create_cg_protein_model.py#52-force-field-of-the-c%CE%B1-sidechain-model)) is used to do [backmapping](#6-backmapping-from-coarse-grained-model-to-all-atom-model).
#### 1.1. Tune *n*<sub>scal</sub> for a small single-domain protein that has experimental folding stability reported
- For a small single-domain protein that has experimental folding stability reported, we use an enhanced sampling protocol i.e. replica exchange simulations to extensively sample the protein folding and unfolding and identify the *n*<sub>scal</sub> value to reproduce the experimental protein folding stability.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**create_cg_protein_model.py** | Create the CG model .psf .top .cor and .prm file that can be used for MD simulations. Need to get the Stride software installed prior to use. ([Learn more](../../wikis/help_wiki/create_cg_protein_model.py)) |
| CG_protein_parameterization/**parse_cg_prm.py** | Parse the parameters in .prm file and then create a .xml file for OpenMM use. ([Learn more](../../wikis/help_wiki/parse_cg_prm.py)) |
| CG_protein_parameterization/**parallel_temperature_REX.py** | Run parallel temperature replica exchange molecular dynamics (pt-REMD) simulation. This simulation is parallelized using multiple CPU processors. ([Learn more](../../wikis/help_wiki/parallel_temperature_REX.py)) |
| CG_protein_parameterization/**opt_temp.pl** | Optimize the temperature windows for pt-REMD simulation to ensure the good sampling quality around the melting temperature of the given protein. ([Learn more](../../wikis/help_wiki/opt_temp.pl)) |
| CG_protein_parameterization/**check_sampling.pl** | Check the sampling quality for pt-REMD simulation. Insufficient sampling will cause problems and inaccuracy in estimating the protein folding stability. ([Learn more](../../wikis/help_wiki/check_sampling.pl)) | 
| CG_protein_parameterization/**dGns.pl** | Convert folding propability *P*<sub>N</sub> to folding stability &Delta;*G*<sub>UN</sub> at a given temperature *T*: <br>$`\Delta G_{\text{UN}} (T)=-k_{\text{B}} T \cdot \mathrm{ln}[\frac{P_{\text{N}} (T)}{1-P_{\text{N}} (T)}]`$. |
| CG_protein_parameterization/**analysis_folding_stability.pl** | Estimate the protein folding stability at a given temperature from pt-REMD data using WHAM. ([Learn more](../../wikis/help_wiki/analysis_folding_stability.pl)) |
| CG_protein_parameterization/**run_REX_LD.py** | A simple script to run a Langevin dynamics (LD) simulation. ([Learn more](../../wikis/help_wiki/run_REX_LD.py)) |
| CG_protein_parameterization/**scan_nscal_REX.pl** | An automated script to call `create_cg_protein_model.py`, `opt_temp.pl`, `parallel_temperature_REX.py` and `analysis_folding_stability.pl` to scan the protein folding stability profile as changing *n*<sub>scal</sub> value. The protein folding stability profile will be used to find the optimized *n*<sub>scal</sub> value for CG model parameterization. ([Learn more](../../wikis/help_wiki/scan_nscal_REX.pl)) | 

- To tune *n*<sub>scal</sub> for a given protein on PSU ACI cluster, run `scan_nscal_REX.pl` with a series of *n*<sub>scal</sub> values. Use `scan_nscal_REX.pl` again to analyze the results and get the protein stability profile. The optimized *n*<sub>scal</sub> is the value that reproduces the experimental protein folding stability. [:leftwards_arrow_with_hook:](#table-of-contents)

#### 1.2. Tune *n*<sub>scal</sub> for a protein without experimental folding stability reported
- For a protein without experimental folding stability reported, we use a stepwise optimization strategy to tune *n*<sub>scal</sub> according to 5 levels of *n*<sub>scal</sub>. The optimized *n*<sub>scal</sub> of an entire single-domain protein or one domain/interface of a multi-domain protein is identified as the lowest level that maintains the native structure. 
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**opt_nscal.pl** | An automated script to find the optimized *n*<sub>scal</sub> for each domain/interface according to 5 levels of *n*<sub>scal</sub> values trained by running the protocol in [Section 1.1](#11-tune-nscal-for-a-small-single-domain-protein-that-has-experimental-folding-stability-reported) for 18 small single-domain proteins. The MD simulator is OpenMM, which is different from that in script `opt_nscal_charmm.pl.pl`. ([Learn more](../../wikis/help_wiki/opt_nscal.pl)) |
| CG_protein_parameterization/**opt_nscal_charmm.pl.pl** | This script has the same function with `opt_nscal.pl` but the levels of *n*<sub>scal</sub> values used in this script was trained by using Charmm for the same protein set. | 

- On PSU ACI cluster, to tune *n*<sub>scal</sub> that is compatible with OpenMM, run `opt_nscal.pl`; To tune *n*<sub>scal</sub> that is compatible with Charmm, run `opt_nscal_charmm.pl`. [:leftwards_arrow_with_hook:](#table-of-contents)

### 2. Temperature quenching simulation
- To estimate the folding rates for a given protein, you need to run temperature quenching simulation where the system is first heated to a very high temperature (usually 800 K) quickly to make protein totally unfolded and then cooled down to the physiological temperature to moniter the refolding of the protein.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_protein_parameterization/**temperature_quenching.py** | Run temperature quenching simulation from the CG native structure of a given protein. This simulation is parallelized using multiple CPU processors. ([Learn more](../../wikis/help_wiki/temperature_quenching.py)) |
| CG_protein_parameterization/**T_quench.pl** | An automated script to build CG model from a pdb file and then run temperature quenching simulations. ([Learn more](../../wikis/help_wiki/T_quench.pl)) | 
| CG_protein_parameterization/**analysis_Tq.pl** | Analyze the results of temperature quenching simulations, fit a single- or double- exponential function to the survival probability of the unfolded protein and then estimate the folding rate.  ([Learn more](../../wikis/help_wiki/analysis_Tq.pl)) | 

- To estimate the protein folding rate on PSU ACI cluster, you need to run `T_quench.pl` with optimized *n*<sub>scal</sub> values obtained from [Section 1](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscal-for-a-given-protein) and then run `analysis_Tq.pl` to do the curve fitting. [:leftwards_arrow_with_hook:](#table-of-contents)

### 3. Create CG ribosome model
- The CG ribosome model is used to run continuous synthesis of a nascent chain that is parameterized with the CG protein model. The CG ribosome is usually fixed (not allow to move) during the simulation. The force field parameters thus only contain the nonbonding term. To speedup the computation, we usually truncate the ribosome to only contain the tails of P- and A-site tRNA molecules, the entire exit tunnel with a few atoms near the tunnel wall and the surface near the exit that may have contacts with the nascent chain.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| CG_ribosome_parameterization/**gen_50S_pdb.py** | Get he necessary subunits from .cif file of the ribosome and create a pdb file contains those sbuunits. ([Learn more](../../wikis/help_wiki/gen_50S_pdb.py)) |
| CG_ribosome_parameterization/**fix_orein_50S_pdb.py** | Add missing atoms and rotate/translate the subunits to a desired orientation for *E. coli* ribosome. ([Learn more](../../wikis/help_wiki/fix_orein_50S_pdb.py)) |
| CG_ribosome_parameterization/**fix_orein_60S_pdb.py** | Do the same thing with `fix_orein_50S_pdb.py` for *S. cerevisiae* ribosome. ([Learn more](../../wikis/help_wiki/fix_orein_60S_pdb.py)) |
| CG_ribosome_parameterization/**create_cg_ribosome_model.py** | Create the CG model for those re-orientated subunits, including .psf, .top and .cor files. ([Learn more](../../wikis/help_wiki/create_cg_ribosome_model.py)) |
| CG_ribosome_parameterization/**truncate_ribosome.py** | Truncate the CG ribosome according to the centroid line of the exit tunnel. ([Learn more](../../wikis/help_wiki/truncate_ribosome.py)) |
| CG_ribosome_parameterization/**gen_ribosome_FF.py** | Train collision diameters of CG ribosome beads from given ribosome structures. ([Learn more](../../wikis/help_wiki/gen_ribosome_FF.py)) |

- To create CG ribosome model, you need to download the .cif file of your ribosome from [RCSB PDB](https://www.rcsb.org/). Use `gen_50S_pdb.py` to generate a pdb file that contains the assembly of subunits that you want to model. Note that the 50S in the script name doesn't restrict the scope of usage. You can use it to generate any assembly of subunits you want, no matter what the organism the ribosome belongs to. Then use `fix_orein_50S_pdb.py` for *E. coli* ribosome and `fix_orein_60S_pdb.py` for *S. cerevisiae* ribosome to get the re-orientated subunits, followed by using `create_cg_ribosome_model.py` to build the CG model for the assembly of subunits. Finally, use `truncate_ribosome.py` to generate the truncated CG ribosome .cor file. If you are interested in tunning force field parameters for CG ribosome model, use `gen_ribosome_FF.py` to train your parameters using some high-resolution crystal structures or Cryo-EM structures. [:leftwards_arrow_with_hook:](#table-of-contents)

### 4. Simulation of co-translational folding
- The simulation protocol for co-translational folding is called "Continuous Synthesis Protocol" (CSP). We add new amino acid on the A-site tRNA in the truncated CG ribosome and simulate the tRNA translocation and nascent chain elongation at a frequency based on in vivo codon translation time.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Continuous_synthesis_protocol/**continuous_synthesis_v6.py** | Run continuous synthesis of a CG protein on a CG ribosome. Both parallelizations on CPU and GPU are supported. ([Learn more](../../wikis/help_wiki/continuous_synthesis_v6.py)) <br>Scripts needed: `Continuous_synthesis_protocol/ribosome_traffic` and `CG_protein_parameterization/parse_cg_prm.py`. |
| Continuous_synthesis_protocol/**ribosome_traffic** | Estimate the real codon translation time by taking into account of the ribosome traffic effects. ([Learn more](../../wikis/help_wiki/ribosome_traffic)) | 
| Continuous_synthesis_protocol/**visualize_cont_synth.py** | Generate movies of the continuous synthesis process. ([Learn more](../../wikis/help_wiki/visualize_cont_synth.py)) <br>Scripts needed: `Backmapping/backmap.py`, `Continuous_synthesis_protocol/render_ecoli_RNC.tcl` and `Continuous_synthesis_protocol/render_yeast_RNC.tcl` | 
| Continuous_synthesis_protocol/**render_ecoli_RNC.tcl** | Render the picture of *E. coli* ribosome-nascent-chain (RNC) complex in VMD.  | 
| Continuous_synthesis_protocol/**render_yeast_RNC.tcl** | Render the picture of *S. cerevisiae* ribosome-nascent-chain (RNC) complex in VMD.  |

- To run CSP, you need to prepare the CG protein model for your nascent chain according to [Section 1](#1-create-cg-protein-models-and-tune-the-force-field-parameters-nscal-for-a-given-protein) and prepare the CG ribosome model according to [Section 3](#3-create-cg-ribosome-model). All the .psf, .top, .cor and .prm files are required in initialization of CSP. In addition, users have to provide a table of intrinsic codon translation time and the mRNA sequence of the nascent chain.
- Below is a video of the entire continuous synthesis process using `continuous_synthesis_v6.py` for synthesizing Firfly luciferase (550 residue long) on the *E.coli* ribosome and visualized using `visualize_cont_synth.py` (time is shown in the experimental timescale):

![Example video of synthesizing Firfly luciferase on the E.coli ribosome](../../wikis/uploads/img/Luc_fast_CSP.mp4)

The large subunit of ribosome is shown in silver; The tail of tRNA  is shown in orange; The nascent chain is represented as Cartoon and colored according to the secondary structures; The flexible region of ribosomal protein L24 is shown in green. Only the last frame of the trajectory for each synthesis step at each NC length is shown in this video. [:leftwards_arrow_with_hook:](#table-of-contents)

### 5. Simulation of post-translational folding
- The simulation of post-translational folding is quite simple. The nascent chain structures that are disassociated from the ribosome will be obtained from [CSP](#4-simulation-of-co-translational-folding) and run Langevin dynamics (LD) in implicit water environment at the physiological temperature.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Post_translational_folding/**post_trans_single_run.py** | Run a single trajectory of post-translational folding. User can specify a walltime or a threshold to control the termination of the post-translational folding. ([Learn more](../../wikis/help_wiki/post_trans_single_run.py)) |
| Post_translational_folding/**PTP_setup.pl** | Automated script to setup post-translation simulations after [CSP](#4-simulation-of-co-translational-folding). ([Learn more](../../wikis/help_wiki/PTP_setup.pl)) <br>Scripts needed: `Post-translational_folding/post_trans_single_run.py` |

- To setup and run post-translation simulations on PSU ACI cluster, use `PTP_setup.pl`. [:leftwards_arrow_with_hook:](#table-of-contents)

### 6. Backmapping from coarse-grained model to all-atom model
- CG model has a lot of benifits on saving computational costs and improving sampling efficiency. However, it loses the atomic-level accuracy of the molecular structure. A backmapping strategy presented here can rebuild the all-atom structure from the CG C&alpha; model with high accuracy. The rebuilt all-atom structure can be furthure used for visualization (e.g., `visualize_cont_synth.py`) and simulation (e.g., [Activation Energy Estimation Workflow](https://git.psu.edu/obrien/yang_jiang/activation-energy-estimation-workflow)) at all-atom level.
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Backmapping/**backmap.py** | Backmap the CG C&alpha; structure to its corresponding all-atom structure. ([Learn more](../../wikis/help_wiki/backmap.py)) <br>Scripts needed: `Backmapping/parse_cg_cacb_prm.py` and `CG_protein_parameterization/create_cg_protein_model.py` |
| Backmapping/**parse_cg_cacb_prm.py** | Parse the CG C&alpha;-sidechain model parameters and convert them into OpenMM .xml format. ([Learn more](../../wikis/help_wiki/parse_cg_cacb_prm.py)) |

- To backmap your CG C&alpha; structure, use `backmap.py`. Note that you need to install [PD2](https://github.com/jmacdona/pd2_public) and [Pulchra](http://cssb.biology.gatech.edu/skolnick/files/PULCHRA/index.html) before use this script. [:leftwards_arrow_with_hook:](#table-of-contents)

### 7. Analysis of protein folding trajectories
- To analyze the protein folding process, we usually calculate the order parameters, such as the fraction of native contacts ($`Q`$), fraction of entangment changes ($`G`$) and fraction of chirality changes ($`K`$).
- Some examples of misfolding protein structures with entanglements can be found [**here**](https://sites.google.com/view/vizentanglements/home).
- Scripts to be used in this section:

| Scripts | Instructions |
| ------ | ------ |
| Analysis_protocol/**get_Ep_from_dcd.py** | Get potential energy from a dcd trajectory using OpenMM for a CG model. ([Learn more](../../wikis/help_wiki/get_Ep_from_dcd.py)) |
| Analysis_protocol/**calc_native_contact_fraction.pl** | Calculate $`Q`$ vs. time for a given trajectory. ([Learn more](../../wikis/help_wiki/calc_native_contact_fraction.pl)) |
| Analysis_protocol/**calc_entanglement_number.pl** | Calculate $`G`$ vs. time for a given trajectory. ([Learn more](../../wikis/help_wiki/calc_entanglement_number.pl)) | 
| Analysis_protocol/**calc_chirality_number.pl** | Calculate $`K`$ vs. time for a given trajectory. (Learn more) | 
| Analysis_protocol/**calc_cont_synth_qbb_vs_T.py** | Automated script to calculate $`Q`$ vs. time for [CSP](#4-simulation-of-co-translational-folding) trajectoris. ([Learn more](../../wikis/help_wiki/calc_cont_synth_qbb_vs_T.py)) <br>Scripts needed: `Analysis_protocol/calc_native_contact_fraction.pl` |
| Analysis_protocol/**mrna_silent_mutation.pl** | Do silent mutation for a given mRNA sequence and a mutation scheme, such as fastest translation, slowest translation and random. ([Learn more](../../wikis/help_wiki/mrna_silent_mutation.pl)) |
| Analysis_protocol/**get_co_trans_order_parameters.py** | Collect the pre-calculated order parameters $`Q`$ and $`G`$ and output .npy data files for the co-translation trajectories. ([Learn more](../../wikis/help_wiki/get_co_trans_order_parameters.py)) |
| Analysis_protocol/**build_co_trans_kinetic_model.py** | Assign the metastable states on the co-translational trajectories and generate the visualization of the representative structures. ([Learn more](../../wikis/help_wiki/build_co_trans_kinetic_model.py)) <br>Scripts needed: `Backmapping/backmap.py` |
| Analysis_protocol/**co_trans_JS_divergence.py** | Estimate the Jensen-Shannon divergence (JSD) for the co-translational mirocstates or metastable states in the function of nascent chain length. ([Learn more](../../wikis/help_wiki/co_trans_JS_divergence.py)) |
| Analysis_protocol/**get_post_trans_order_parameters.py** | Collect the pre-calculated order parameters $`Q_{\mathrm{act}}`$ and $`G`$ and output .npy data files for the post-translation trajectories. ([Learn more](../../wikis/help_wiki/get_post_trans_order_parameters.py)) <br>Scripts needed: `Analysis_protocol/calc_native_contact_fraction.pl` and `Analysis_protocol/calc_entanglement_number.pl` |
| Analysis_protocol/**build_post_trans_kinetic_model.py** | Assign the metastable states on the post-translational trajectories and build the Mater equation model for the state probabilities. ([Learn more](../../wikis/help_wiki/build_post_trans_kinetic_model.py)) <br>Scripts needed: `Backmapping/backmap.py` |
| Analysis_protocol/**post_trans_JS_divergence.py** | Estimate the Jensen-Shannon divergence (JSD) for the post-translational mirocstates or metastable states in the function of time. ([Learn more](../../wikis/help_wiki/post_trans_JS_divergence.py)) |
| Analysis_protocol/**get_co_post_folding_pathways.py** | Combine the discrete metastable states trajectories and analyze the folding pathways. ([Learn more](../../wikis/help_wiki/get_co_post_folding_pathways.py)) |
| Analysis_protocol/**get_solubility.py** | Estimate the aggregation, degradation and Hsp70 binding propensities for the post-translational metastable states and estimate the percent of soluble protein in each metastable state. ([Learn more](../../wikis/help_wiki/get_solubility.py)) <br>Scripts needed: `Backmapping/backmap.py` |
| Analysis_protocol/**get_state_probability.py** | Builds the Mater equation model for the state probabilities and runs the bootstrapping resampling for error estimation. ([Learn more](../../wikis/help_wiki/get_state_probability.py)) |
| Analysis_protocol/**calc_specific_activity.py** | Estimates the specific activities, as well as the errors and p-value, by using the $`\Delta G^{\ddagger}`$ obtained from [Activation Energy Estimation Workflow](https://git.psu.edu/obrien/yang_jiang/activation-energy-estimation-workflow) and the state probabilities obtained from `get_state_probability.py`. ([Learn more](../../wikis/help_wiki/calc_specific_activity.py)) |

[:leftwards_arrow_with_hook:](#table-of-contents)
