These files are the data, scripts and information files needed to perform the analysis described in the article “Forensic Geolocalization of Norwegian Soil Samples by Applying Machine Learning to Microbiomes”. 
-	The folder “fastq_files” contains the demultiplexed fasta files with the sequencing data. Files are named with the combination of named barcodes used for that sample. The files are trimmed and only the primers remain a part of the sequences. 
-	dada2_slurm.sh: the script used to run the dada2 software in an R environment on a HPC. We use Sigma2 - the National Infrastructure for High-Performance Computing and
Data Storage in Norway
-dada2_v1.16_batch.R: The rscript used to do the dada2 processing of the sequence files to produce ASV tables and asign taxonomy using the silva reference databases. 
-	A1_pm.txt: The placement of fall samples on the plate during the library preparation
-	pm_pilot.txt: The placement of summer samples on the plate during the library preparation
-	platemat.txt: The barcode combination used for each well in the plates during library preparation
-	Script_1_Filtering_and_naming.Rmd: The first Rscript to run. In this script the barcode information is used to assign correct names to the samples and the ASV tables are filtered.
-	contamination.R: This script is sourced in the script “Script_1_Filtering_and_naming.Rmd” and uses the package “decontam” to identify contaminant ASVs
-	Script_2.1_Descriptive_plots.Rmd: This script is used to perform several of the analyses. NMDS, PERMANOVA, taxonomic analyses, and Bray-Curtis distances are the most central ones.
-	Script_2.2_Prediction.Rmd: This script is used for the RandomForest classification. 
