#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_filter_snps_and_make_wide_format.out
#SBATCH -e run_filter_snps_and_make_wide_format.err

cd /scratch.global/haasx092/claudia_gbs_may_2022/220518_snp_calling_results

module load R/3.6.0
Rscript filter_snps_and_make_wide_format.R 220518_snp_calling_normalize.tsv 220518_claudia_gbs_may_2022.csv 220518_claudia_gbs_may_2022_snps.Rdata
