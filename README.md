# Analysis of Claudia's GBS data (May 2022)
Analysis of Claudia's GBS data from May 2022

The data can be found here: `/home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/
`<br><br>
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
rm -rf ample_sample_number
```
Next, you should make a file with the list of directories. This text file will come in handy for future steps of the GBS analysis.
```bash
ls Sample*/ -d | tr -d / > 220516_claudia_analysis_sample_directory_list.txt
```
This next step is necessary because we are working with paired-end reads. We are doing it because the csv file contains 2 lines per sample (one for the forward read and one for the reverse read).
```bash
awk 'FNR%2' 220516_claudia_analysis_sample_names_and_numbers.csv > 220516_claudia_analysis_file_list_every_other.csv
```
_Make sure you open the resulting file using `vi` to manually remove the header line._ Once that is done, we can make symbolic links (symlinks) to point to the data rather than take up disk space by needlessly duplicating the original files. **Note** that when I analyzed the original dataset, I set `n` to start at 73 and then increased that value by 1 with each iteration of the `while` loop. Since this iteration of the analysis only contains the GWAS samples, there are gaps in the sequence of sample numbers, necessitating a different approach. The approach I used involves extracting the sample number ("Snumber") from the file name and using that rather than relying on counting iterations through the loop.
```bash
# Make symlinks to GBS data
cat 220516_claudia_analysis_file_list_every_other.csv | cut -f 9 -d / \
	| while read i; do
	STEM=$(echo $i | rev | cut -f 3,4,5,6,7,8 -d "_" | rev)
	Snumber=$(echo $i | rev | cut -f 3 -d "_"| rev | sed 's/^S//g')
	n=$(printf "%04d\n" $Snumber)
	echo $STEM
	ln -s /home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/${STEM}_R1_001.fastq.gz Sample_$n/Sample_${n}_R1.fq.gz
	ln -s /home/jkimball/data_delivery/umgc/2022-q2/220502_A00223_0826_AHJTHFDSX3/Kimball_Project_009/${STEM}_R2_001.fastq.gz Sample_$n/Sample_${n}_R2.fq.gz
	done
```
