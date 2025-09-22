#!/bin/bash
#SBATCH --account="your_cluster_account"
#SBATCH --time=03:00:00
#SBATCH --job-name=dada2_RD
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --output=dada2.log
#SBATCH --mail-type=END

set -o errexit #make bash exit on any error
set -o nounset #treat unset variable as errors
module --quiet purge
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
Rscript dada2_v1.16_batch.R
