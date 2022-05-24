# Analysis of Claudia's GBS data (May 2022)
Analysis of Claudia's GBS data from May 2022
## Contents
1. [Directory Setup](#Directory-Setup)
2. [Quality Control](#Quality-Control)
3. [Adapter Trimming](#Adapter-Trimming)
4. [Alignment](#Alignment)
5. [Index BAMs](#Index-BAMs)
6. [SNP calling](#SNP-calling)
7. [Filter SNP calls](#Filter-SNP-calls)

## Directory Setup
The data can be found here:
```bash
/home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/
```
We will be working in the [global scratch](https://www.msi.umn.edu/content/scratch-storage) directory so we don't consume our group's storage quota. Just reminder: If you want or need to check the amount of resources the group is using, type `groupquota -g jkimball` on the command line. (Really, `groupquota` is enough, but if `jkimball` isn't your primary group, you'll need to include the second part.<br><br>
The path to our working directory is: `/scratch.global/haasx092/claudia_gbs_may_2022`.

The overall GBS project (Kimball_Project_009) contains projects from multiple group members. Claudia's data begin with either "R" or "S" so we will take advantage of that in order to create a text file containing just her samples.<br><br>
For this step, we will be working in the directory where the data are located.<br>
```bash
cd /home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009
```
We can see all of the files that start with "R" using this code:
```bash
ls R*
```
Similarly, we can see all of the files that start with "S" using this code:
```bash
ls S*
```
**Note:** We are fortunate that these simple patterns don't conflict with samples from other projects in the same submission. For now, things look pretty good. There are only a few potential files that could cause trouble later on because they are different from all the rest, so it's a good idea to keep them in mind. (R109 lacks a hyphen between the R and 109; and a handful of other samples have "rep" after the number.)<br>

Now we can begin the process of writing the sample (file) names to the text file that we can iterate over.<br>
```bash
ls R* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
```
Followed by:
```bash
ls S* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
```
After finishing the entire directoy, I discovered there were a handful of samples that should be included in Claudia's dataset, but do not follow the pattern previously established with the "R" or "S" prefixes. They are: "10-R-2021-pool1", "10-S-2021-pool2", "12A-2021-pool3", and "12B-2021-pool3". So, I need to go back and re-do the setup and FastQC/Adapter Trimming steps so that these are included. Fortunately, this didn't cost too much extra time. I caught the issue early.
```bash
ls 10-R-2021-pool1* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
ls 10-S-2021-pool2* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
ls 12A-2021-pool3* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
ls 12B-2021-pool3* >> /scratch.global/haasx092/claudia_gbs_may_2022/claudia_samples.txt
```

**Note** The `>>` allows output from a command like `ls` to be appended to a file unlke `>` which would overwrite anything in the file that follows the symbol.

Now, we can move back to the global scratch directory where we are working:
```bash
cd /scratch.global/haasx092/claudia_gbs_may_2022
```

This project used paired-end sequencing, so there are actually 2 files (spread out over 2 lines) per sample. They are denoted by "R1" (the forward read) and "R2" (the reverse read). So, as we build the directory structure, we need to keep this in mind because we want to create 1 subdirectory per sample and group both "R1" and "R2" reads from the same sample in each respective subdirectory.

The next step is to move into the `R` statistical environment to go from filenames to a workable CSV file that will be used in the next stage of the working directory setup. I used R version 4.1.0 (`module load R/4.1.0`)
```R
library(data.table)
fread("claudia_samples.txt", header = F) -> x

# Change the column name from V1 to something more informative (filename)
setnames(x, "filename")

# Add a new column called sample_number. It will initially contain the entire filename, but we will work to retain only the sample number
x[, sample_number := filename]

# Strip off first part of filename until sample number begins (S) but do not include it.
x[, sample_number := sub("^.+[S]", "", sample_number)]

# Strip off end of the filename (after the sample number) ... begins with "_R1" or "_R2"
x[, sample_number := sub("_[R1].+$", "", sample_number)]

# Convert sample numbers to numerical and add leading zeros to all samples (to help with sorting).
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]

# Reorder rows in ascending order
x[order(sample_number)] -> x

# Set column order (my personal preference for sample_number to come first)
setcolorder(x, c("sample_number", "filename")) -> x

# Write output to CSV
write.csv(x, file="220516_claudia_analysis_sample_names_and_numbers.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)

# Save table as an R object
save(x, file="220516_claudia_analysis_sample_names_and_numbers.Rdata")
```
After that is done, use the CSV file using bash to create the directory structure.
**Note:** The `echo $i` part is not really necessary. I just included it to watch the progress.

```bash
cat 220516_claudia_analysis_sample_names_and_numbers.csv | cut -f 1 -d , \
	| while read i; do
	d=Sample_$i
	echo $i
	mkdir -p $d
	done
```

Once that is done, you will probably notice that there is a directory called "Sample_sample_number" which is an artefact of the code. I probably could change the code so that the header isn't interpreted as a sample name, but it's also super easy to just remove it after the code finishes. You can easily remove it with a one-liner:
```bash
rm -rf Sample_sample_number
```
Next, you should make a file with the list of directories. This text file will come in handy for future steps of the GBS analysis.
```bash
ls Sample*/ -d | tr -d / > 220516_claudia_analysis_sample_directory_list.txt
```
This next step is necessary because we are working with paired-end reads. We are doing it because the csv file contains 2 lines per sample (one for the forward read and one for the reverse read).
```bash
awk 'FNR%2' 220516_claudia_analysis_sample_names_and_numbers.csv > 220516_claudia_analysis_file_list_every_other.csv
```
_Make sure you open the resulting file using `vi` to manually remove the header line._ Once that is done, we can make symbolic links (symlinks) to point to the data rather than take up disk space by needlessly duplicating the original files. Since this iteration of the analysis only contains the your (Claudia) samples, the first number is 843 because that's where your samples began in the dataset.. The approach I used involves extracting the sample number ("Snumber") from the file name and using that rather than relying on counting iterations through the loop.
```bash
# Make symlinks to GBS data
cat 220516_claudia_analysis_file_list_every_other.csv | cut -f 2 -d , \
	| while read i; do
	STEM=$(echo $i | cut -f 1,2 -d "_")
	Snumber=$(echo $i | rev | cut -f 3 -d "_"| rev | sed 's/^S//g')
	n=$(printf "%04d\n" $Snumber)
	echo $STEM
	ln -s /home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/${STEM}_R1_001.fastq.gz Sample_$n/Sample_${n}_R1.fq.gz
	ln -s /home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/${STEM}_R2_001.fastq.gz Sample_$n/Sample_${n}_R2.fq.gz
	done
```
In the next step, we will move back to the `R` statistical environment to create a sample key.
```R
library(data.table)

# Read in data
x <- fread("220516_claudia_analysis_sample_names_and_numbers.csv")

# Change column names
setnames(x, c("sample_number", "sample_name"))

# Add leading zeros
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]
# Add "Sample_" to each sample number
x[, sample_number := paste0("Sample_", sample_number)]

# Remove beginning the beginning part of the filename to remove the part of the path that is no longer necessary to keep
x[, sample_name := sub("^.+Project_008/", "", sample_name)]

# Remove trailing part of filenames (sample names)---ultimately, we only need one line per sample, not two (a consequence of having 2 files per sample for paired-end reads)
x[, sample_name := sub("_[R1].+$", "", sample_name)]
x[, sample_name := sub("_[R2].+$", "", sample_name)]

# Retain unique values only
x <- unique(x)

# Save to CSV
write.csv(x, file="220516_claudia_analysis_sample_key.csv", row.names = FALSE, sep=",", quote=FALSE)
```
This next bit is just to check if the symlinks are all correct. You can do it manually (by going into each directory and typing `ls -lh`, paying special attention to the first and last to make sure the sample numbers match the "S number" in the filenames provided by UMGC. This way of checking it just quickly puts the same information into a single file so you can view them all at once. Note: I only did the forward reads ("R1") because if they are correct, the reverse reads ("R2") will also be correct.
```bash
for i in $(cat 220516_claudia_analysis_sample_directory_list.txt);
do
ls -lh ${i}/${i}_R1.fq.gz >> check_symlinks_full_paths.txt
done
```

## Quality Control
After we have the directory structure setup, the first thing we want to do is check the FASTQ files with [FastQC](https://github.com/s-andrews/FastQC) which I accomplished using the [run_fastqc.sh](fastqc/run_fastqc.sh) and [fastqc_wrapper_script.sh](fastqc/fastqc_wrapper_script.sh). The results look fine (some flags and warnings are to be expected). The major take-away from this process is that the adapter contamination comes from the Nextera Transposase Sequence. This is helpful information because now we know what sequence to trim in the next step using [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

## Adapter Trimming
To remove the Nextera Transposase Sequence, we used the [run_cutadapt.sh](cutadapt/run_cutadapt.sh) and [cutadapt_wrapper_script.sh](cutadapt/cutadapt_wrapper_script.sh) scripts.

## Alignment
After the adapter sequences have been removed, it was time to perform the alignment using the [run_bwa.sh](alignment/run_bwa.sh) script.

## Index BAMs
Before moving on to the actual SNP calling step, you must first index the bam files which you can accomplish using the [index_bams.sh](index_bams/index_bams.sh) script. Actually, you won't be able to move forward without completing this step. If you try, the SNP calling script [scythe_mpileup.sh](snp_calling/scythe_mpileup.sh) will run, but it will end prematurely if the bam index files are not there.

## SNP calling
Now, we proceed to the actual SNP-calling step using [scythe_mpileup.sh](snp_calling/scythe_mpileup.sh). One parameter to pay particular attention to is the `-q 20` option. This means that _the minimum mapping quality (MQ) for an alignment to be used_ is 20. This is a measurement of the quality of the read being mapped, not the base at the SNP. You can increase the stringency up to `-q 60`, although `-q 20` is acceptable It's a phred-like score and means that the minimum acceptable probability for a read being correct is 99%. Reads with a lower mapping quality are filtered out. Many (if not most) reads will have an even higher probability of being right.**Note:** This script uses [GNU Parallel](https://www.gnu.org/software/parallel/), so make sure you cite the program in any manuscript that uses results from these data. You can get the citation info by running `parallel --citation`. (You'll need to run `module load parallel` first.)

## Filter SNP calls
Once the SNP calling step is done, you will have a list of 2,183 g-zipped `VCF` files (`.vcf.gz`). There is one file per chromosome/scaffold. Most of these don't contain any SNPs at all, so it isn't worth looking at them. They're also quite small (insignificant) in terms of length of genome sequence. Since we renamed the scaffolds, you no longer need to worry about the original scaffold names deliered to us by Dovetail. You will need to make a file like [vcf_file_list.txt](helper_files/vcf_file_list.txt). I made this manually because it's just a list of the files we actually want to look at (instead of all 2,183). The first one is `220518_snp_calling_results_ZPchr0001.vcf.gz`, the second is called `220518_snp_calling_results_ZPchr0002.vcf.gz`, and so forth all the way through `220518_snp_calling_results_ZPchr0016.vcf.gz`. However, the 17th scaffold (which is important because it is greater than 4 Mb in size contains the Northern Wild Rice _sh4_ ortholog) is called `220518_snp_calling_results_ZPchr0458.vcf.gz`. It was originally Scaffold_453, but we didn't include it in the renaming process because it wasn't among the 15 largest scaffolds.  If we had included it, it would have been ZPchr0017.

Anyway, use the script [filter_with_vcftools.sh](filter_vcfs/filter_with_vcftools.sh) to filter the `VCF` files in order to meet your desired parameters. The way the script is currently written, the parameters are:<br>
* Maximum 10% missing data (`--max-missing 0.90`). _I know this is confusing, but it's correct._
* Bi-allelic sites only (`--min-alleles 2 --max-alleles 2`)
* Minor allele frequency is 0.03 (`--maf 0.03`)
* No indels (`--remove-indels`)
* Minimum depth of 8 reads required at a SNP (`--minDP 8`)

**Note:** So far, most of the software programs we have been using so far have already been installed by the Minnesota Supercomputing Institute (MSI). That's why you can use them by calling `module load` and then referring to them in your code simply by calling the name of the program (e.g., `bwa`, `samtools`, or `bcftools`). [`VCFtools`](https://vcftools.github.io/index.html) is different because I had to install it myself and refer to the place where it is installed in my script (`~/vcftools/bin/vcftools`) rather than just using `vcftools`.

Once the `VCF` files have been filtered according to your desired parameters, you can move on to the next step: putting the SNP calls into a `CSV`-formatted SNP matrix. However, I also like working with [plink](https://zzz.bwh.harvard.edu/plink/index.shtml), especially for performing principal component analysis (PCA). As a first step in that analysis, I merge the 17 filtered `VCF` files into a single merged `VCF` file with [concat_filtered_vcfs.sh](filter_vcfs/concat_filtered_vcfs.sh).
