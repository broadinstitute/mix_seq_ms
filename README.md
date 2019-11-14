# MIX-Seq

This repo contains code associated with our manuscript describing MIX-Seq, our method for multiplexed transcriptional profiling of perturbation responses across mixtures of cancer cell lines.

For SNP-classification code, see [single_cell_classification](https://github.com/acwarren/single_cell_classification)

## Data

The data associated with this analysis will be made available shortly.



## Organization of repo

The code can be organized into config files, helper functions, preprocessing scripts, and analysis/figure generation scripts.


### Configs

- build_figshare_repo.R: script that builds the data repository for Figshare from versioned files stored in our internal data management system (Taiga)

- expt_meta_data.R: Create metadata defining which datasets and comparisons are used in the analysis

- global_params.R: Define global params shared across analysis scripts.

NOTES: 
-In order to run the analyses locally, you should set the global parameter 'local_data' to TRUE in global_params.R. If you have the data from the figshare repo stored locally in a folder called 'data' within your project directory (each data directory unzipped). The code should run. 




### Helper functions

- MixSeq_helpers.R: Define helper functions used throughout the analysis


### Preprocessing

- compute_cell_cycle_change.R: Computes statistics about change in cell cycle phase composition between treat vs control groups for each experiment. Saves these into an intermediate file "module_scores.csv"

- compute_CL_features.R: script used to compile cell line features used in the analysis (e.g. drug sensitivity data, and any genomic features used). The results of this script (as an rds file) are included in the figshare repo.

- compute_collapsed_profiles.R: Precompute sum-collapsed and mean-collapsed expression profiles for each cell line and treatment condition. Saves results as 'summed_counts.rds' and 'avg_profiles_cpm.rds' files

- compute_lfc_pca.R: Computes PCA on matrices of log-fold change values across cell lines and genes for each treatment. 

- compute_viability_sigs.R: Runs differential expression analysis to identify 'viability-related' and 'viability-independent' components of the transcriptional response across cell lines.


## Analysis/fig-gen

There are two scripts that orchestrate running all analyses and generating figs for the manuscript. 

- make_SNP_class_figs.R: Run all scripts used to generate figures associated with SNP-classification analysis (Fig 1, S1, S2)

- make_figs.R: Run all scripts used to generate the remaining figures.