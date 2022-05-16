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
