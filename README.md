# COVID: GEMs and MEMs

This repository is supplementing the paper **Genome-scale metabolic modelling of COVID-19: Integration of omics data to generate and analyse COVID-19 specific genome-scale metabolic models**. 


## Data preprocessing (Shell)
* [`make_index.sh`](/Code/Linux/make_index.sh): Downloads the Human transcriptome data and generate the kallisto transcriptome index. 
* [`NHBE.sh`](/Code/Linux/NHBE.sh): Downloads the HBE RNA-seq fastq files and generate the quantification with Kallisto.
* [`293T.sh`](/Code/Linux/293T.sh): Downloads the 293T RNA-seq fastq files and generate the quantification with Kallisto. 
* [`A549.sh`](/Code/Linux/A549.sh): Downloads the A549 RNA-seq fastq files and generate the quantification with Kallisto. 
* [`Calu.sh`](/Code/Linux/Calu.sh): Downloads the Calu-3 RNA-seq fastq files and generate the quantification with Kallisto.  
* [`Lung.sh`](/Code/Linux/Lung.sh): Downloads the Lung RNA-seq fastq files and generate the quantification with Kallisto. 

## Model extraction (Matlab)
* [`Main.m`](/Code/Matlab/Main.m): A script that allows the user to choose a dataset and a MEM used in the extraction process.
* [`Generation_iMAT_alt.m`](/Code/Matlab/Generation_iMAT_alt.m): A function that generate the healthy and the infected model with the iMAT algorithm.
* [`Generation_GIMME_alt.m`](/Code/Matlab/Generation_Gimme_alt.m): A function that generate the healthy and the infected model with the GIMME algorithm.

## Analysis and visualisation (Python)
The analysis of extracted models and the visualisation of the obtained results were performed in Python. The main files are as follows:
* [`model_PCA.ipynb`](model_PCA.ipynb): Performs the PCA analysis on the models obtained with different MEMs.
* [`model_PCA_sampling.ipynb`](model_PCA_sampling.ipynb): Performs the PCA analysis on the model fluxes (the reactions with at least one flux sample different than zero are kept) obtained with different MEMs.
* [`model_statistics.ipynb`](model_statistics.ipynb): Performs the basic anaylses of the extracted models (model sizes, Jaccard index calculation, reaction specificity analysis).
* [`model_statistics_sampling.ipynb`](model_statistics_sampling.ipynb): Similar as [`model_statistics.ipynb`](model_statistics.ipynb) but performs the analysis on the results of flux sampling.
* [`reactions_enrichment.ipynb`](reactions_enrichment.ipynb): Compares pairs (healthy, infected) on a reaction level based on flux samples. Calculates p-values, fold-changes and enrichment. It also extracts active reactions per pairs (reactions that are active - nonzero) in at least one model within a pair.
* [`reactions_subsystems_enrichment.ipynb`](reactions_subsystems_enrichment.ipynb): Produces an XLSX file with enriched reactions per subsystem for all models.
* [`subsystems_enrichment.ipynb`](subsystems_enrichment.ipynb): Performs subsystems enrichment analysis based on the hypergeometric test.
* [`subsystems_enrichment_plot.ipynb`](subsystems_enrichment_plot.ipynb): plots the results of the enrichment analysis.

## Results
* [`results_basic`](results_basic): results of the basic statistical analysis.
* [`results_PCA`](results_PCA): results of the PCA performed on the obtained reactions.
* [`results_PCA_sampling`](results_PCA_sampling): results of the PCA performed on the reactions with nonzero flux values in flux samples.
* [`results_enrichment`](results_enrichment): results of the enrichment analyses.
* [`results_enrichment/reactions_subsystems.xlsx`](results_enrichment/reactions_subsystems.xlsx): full results of the enrichment analysis (for all metabolic reactions).
