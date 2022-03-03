# Genome-scale metabolic modelling of COVID-19

This repository is supplementing the paper **Integration of omics data to generate and analyse COVID-19 specific genome-scale metabolic models**. The main parts of the code in the repository are
* [Data preprocessing (Bash)]: downloading and preprocessing the transcriptomics data used as an input to model extraction methods in the next step.
* **Model extraction (Metlab)**: the integration of the transcriptomics data with the reference model to obtain context-specific models.
* **Analysis and visualisation (Pyhon)**: the analysis of context-specific models and visualisation of the results.

## Data preprocessing (Bash)
* [`make_index.sh`](/Code/Linux/make_index.sh): Downloads the Human transcriptome data and generate the kallisto transcriptome index. 
* [`NHBE.sh`](/Code/Linux/NHBE.sh): Downloads the HBE RNA-seq fastq files and generate the quantification with kallisto.
* [`293T.sh`](/Code/Linux/293T.sh): Downloads the 293T RNA-seq fastq files and generate the quantification with kallisto. 
* [`A549.sh`](/Code/Linux/A549.sh): Downloads the A549 RNA-seq fastq files and generate the quantification with kallisto. 
* [`Calu.sh`](/Code/Linux/Calu.sh): Downloads the Calu-3 RNA-seq fastq files and generate the quantification with kallisto.  
* [`Lung.sh`](/Code/Linux/Lung.sh): Downloads the Lung RNA-seq fastq files and generate the quantification with kallisto. 

## Model extraction (Matlab)
* [`Main.m`](/Code/Matlab/Main.m): A script that allows the user to choose a dataset and a model extraction method for the extraction process.
* [`generate_iMAT.m`](/Code/Matlab/generate_iMAT.m): A function that generates a healthy and an infected model with the iMAT algorithm.
* [`generate_GIMME.m`](/Code/Matlab/generate_GIMME.m): A function that generates a healthy and an infected model with the GIMME algorithm.
* [`generate_INIT.m`](/Code/Matlab/generate_INIT.m): A function that generates a healthy and an infected model with the GIMME algorithm.
* [`generate_tINIT.m`](/Code/Matlab/generate_tINIT.m): A function that generates a healthy and an infected model with the GIMME algorithm.

## Analysis and visualisation (Python)
The analysis of extracted models and the visualisation of the obtained results were performed in Python. The main files are as follows:
* [`model_PCA.ipynb`](model_PCA.ipynb): Performs the PCA analysis on the models obtained with different MEMs.
* [`model_PCA_sampling.ipynb`](model_PCA_sampling.ipynb): Performs the PCA analysis on the model fluxes (the reactions with at least one flux sample different than zero are kept) obtained with different MEMs.
* [`model_statistics.ipynb`](model_statistics.ipynb): Performs the basic analyses of the extracted models (model sizes, Jaccard index calculation, reaction specificity analysis).
* [`model_statistics_sampling.ipynb`](model_statistics_sampling.ipynb): Similar as [`model_statistics.ipynb`](model_statistics.ipynb) but performs the analysis on the results of flux sampling.
* [`reactions_enrichment.ipynb`](reactions_enrichment.ipynb): Compares pairs (healthy, infected) on a reaction level based on flux samples. Calculates p-values, fold-changes and enrichment. It also extracts active reactions per pairs (reactions that are active - nonzero) in at least one model within a pair.
* [`reactions_subsystems_enrichment.ipynb`](reactions_subsystems_enrichment.ipynb): Produces an XLSX file with enriched reactions per subsystem for all models.
* [`subsystems_enrichment.ipynb`](subsystems_enrichment.ipynb): Performs subsystems enrichment analysis based on the hypergeometric test performed on the encriched metabolic reactions.
* [`subsystems_enrichment_plot.ipynb`](subsystems_enrichment_plot.ipynb): plots the results of the metabolic substystems enrichment analysis.

## Results
* [`results_basic`](results_basic): results of the basic statistical analysis.
* [`results_PCA`](results_PCA): results of the PCA performed on the preserved reactions after model extraction process.
* [`results_PCA_sampling`](results_PCA_sampling): results of the PCA performed on the reactions with nonzero flux values after flux sampling on the extracted models.
* [`results_enrichment`](results_enrichment): results of the enrichment analyses.
* [`results_enrichment/reactions_subsystems.xlsx`](results_enrichment/reactions_subsystems.xlsx): full results of the enrichment analysis (for all metabolic reactions).
