#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_fastqc_trimmed.out
#SBATCH -e run_fastqc_trimmed.err

cd /scratch.global/haasx092/claudia_gbs_may_2022

module load fastqc
module load parallel

cat 220516_claudia_analysis_sample_directory_list.txt | parallel -k --will-cite -j 10 \
        sh ./fastqc_wrapper_script.sh
