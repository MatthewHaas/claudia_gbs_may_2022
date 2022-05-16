#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=36:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_bwa.out
#SBATCH -e run_bwa.err

cd /scratch.global/haasx092/claudia_gbs_may_2022

module load bwa
module load samtools

FASTA='/home/jkimball/shared/WR_Annotation/NCBI_submission/zizania_palustris_13Nov2018_okGsv_renamedNCBI2.fasta'

# FASTA must be indexed first
bwa index $FASTA

for i in $(cat 220516_claudia_analysis_sample_directory_list.txt); do
bwa mem -t 32 $FASTA ${i}/${i}_R1_trimmed.fq  ${i}/${i}_R2_trimmed.fq 2> ${i}/${i}_bwa.err | samtools sort -o ${i}/${i}_sorted.bam 2> ${i}/${i}_samtools_sort.err;
done
