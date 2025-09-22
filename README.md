These files are the data, scripts and information files needed to perform the analysis described in the article “Forensic Geolocalization of Norwegian Soil Samples by Applying Machine Learning to Microbiomes”. 

File descriptions:
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

The main steps:

1. dada2_slurm.sh

Input files:
- files in "fastq_files" folder
- dada2_v1.16_batch.R
dada2_slurm.sh is run in the HPC saga provided by Sigma2. This shell script runs the dada2_v1.16_batch.R script in an R environment. The script runs the software dada2 (https://benjjneb.github.io/dada2/tutorial.html).
- To run the dada2_slurm.sh script in saga you need to change the account in line 2 to your account.
- To run the script dada2_v1.16_batch.R you need to provide the path to your fasta files in line 4, and your path to the Silva database files in line 63 and 64.

Output files:
- track_sequence_loss_table.txt
- dada2_tab.txt
- dada2_taxonomy.txt
- ASV.fas
It is also possible to directly run dada2_v1.16_batch.R in a local R environment

2. Script_1_Filtering_and_naming.Rmd

Input files:
- dada2_tab.txt
- dada2_taxonomy.txt
- ASV.fas
- platemap.txt
- pm_pilot.txt
- A1_pm.txt
- contamination.R

Output:
- tab.txt
- taxonomy_filtered.txt
- sequences_filtered.txt
- metaData.txt

3. Script_2.1_Descriptive_plots.Rmd

Input files:
- tab.txt
- taxonomy_filtered.txt
- metaData.txt

Main Output:
- NMDS plots
- PERMaNOVA tables
- Bray_Curtis tables
- Abundance chart

4. Script_2.2_Prediction.Rmd

Input:
- tab.txt
- taxonomy_filtered.txt
- metaData.txt

Output:
- Confusion matrixes for all predictions
- Information table from all predictions
  
