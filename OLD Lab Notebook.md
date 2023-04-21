Lab Notebook

ls 
# this lists all of the files in my current working directory 
# use this to check if your commands worked when moving files into directory 
# Add a -F which adds / to which directories I can enter
# Add a -lrth lists premissions and other info 

pwd
# "print working directory" shows what folder I am in (Gen 711_Lab) 

CTRL C
# cancel the command line and starts new one

"clear"
# this will clear the terminal 

PS1='$ '
# setting a variable to rid of the text on the command line 

cd
# change directory - write name of folder you want to change into

rm 
# remove - write name of folder you want to delete forever 

ssh kgb1008@ron.sr.unh.edu
# this will connect terminal to the teaching cluster RON - will be in your own folder or workspace 
# New password: Tillywilly0709

****02/10/2023 NOTES:

mkdir NAME
# linex command to make new directory, after space write NAME of new folder you want

cp
# copy file

TAB
# auto completes command line with what you have in directory that matches

.tar
# similar to zip file, compressed data

tar -xf shell_data.tar
# -xf represent options of how to open tar file, will decompress files

man ls
# opens manual listing commands for directory 
# hit q to close out of manual 

head
# shows line of code 

grep 'what you wanna find' what file you want
# searches through files for what you want (like cntl F)
# pike (|) wc

wc
# word count

grep '@SRR097977' SRR097977.fastq | wc -l
# searches file for 'what is in the quotes' in the specific file; the | (pike) wc counts how many times the thing you searched for shows up in your directory

cd ~ 
# brings you back to home directory 

****02/17/2023 NOTES:

relative path
# enetering a path from "step-by-step"

absolute path
# always start with fowards slash, gets you directly to location
# pwd tells you absolute path

cd ../
# to back out of a folder

cd ../../
# back out two directories


cd ../ hit tab WITHOUT CLICKING ENTER
# shows what you can search

cd /home/users/kgb1008/gen711/shell_data/untrimmed_fastq
# absolute way of entering untrimmed_fastq

 Exc 3a - 3 ways to change directories home from untrimmed_fastq:
 1. cd ../../../
 2. cd ~
 3. /home/users/ndb1029/

ls -a
 #  reveal hidden files

/gen711/shell_data/.hidden$ ls
/gen711/shell_data/.hidden$ head youfoundit.txt
# how to open hidden files
# use ls and head command

cat
# outputs whats inside the text file . Used for sticking files together ex: cat youfoundit.txt youfoundit.txt

ls *
# glob / wildcard
# different than ls as it is more specific

/gen711/shell_data/untrimmed_fastq$ ls *.fastq
# will look up terms with just fastq

ls /
# shows all files that server has available

Exc: list applications available on ron that start with "c," "a," and "o"
/gen711/shell_data/untrimmed_fastq$ ls /bin/c* | wc -l
1. starts with c: 90
/gen711/shell_data/untrimmed_fastq$ ls /bin/a* | wc -l
2. starts with a: 45
/gen711/shell_data/untrimmed_fastq$ ls /bin/o* | wc -l
3. starts with o:22

echo 'hi'
# will state the line out

/gen711/shell_data/untrimmed_fastq$ echo *fastq
# will echo what is in folder

>>
echo - gives it back to you, use name and appended 
# echo kim bonanno 
# 
*fastnothing >> newfile
# outputs:
*fastnothing
*fastnothing
# two lines will cause two lines

ctr R
# searches old commands that was run previously
# click ctrl R multiple times to switch results

history
# prints history out to screen

history | grep 'grep'
# finds history that matches grep

grep 'NNNNN' SRR098026.fastq > badreads.fastq

 grep -A -B 'NNNNN' SRR098026.fastq > badreads.fastq

 grep -A2 -B1 'NNNNN' SRR098026.fastq
 # give 2 lines before, 1 after
 cat badreads to see file once saved

 grep -A2 -B1 'N*NNNN' SRR098026.fastq
 # can be anything after

****02/24/2023 NOTES:

cd ../
# one step back in directories 
## gen711/shell_data/untrimmed_fastq$ cd ../ brings you back to (gen711/shell_data$)

realpath gen711/
# helps find absolute path of a specific foler or directory (gen711 folder)

ls -a (folder name)/
# shows hidden files in folders

head and cat before code line (shell_data/.hidden/youfoundit.text)
# shows you what is in a selected file 

How many programs in /bin 1 USE Wildcard (*)
# command = ls/bin/c*

>> double redirect 
# pend file 

tail
#

cp /tmp/*.fastq.gz .
# used to copy all files matching ".fastq.gz" 

gunzip
# unzip compressed files 
## gunzip *fastq.gz - unzips all files with fastq.gz at the end of the file 

Control 'Z'
# stops process 

bg
# puts process in the background to get command line back

top
# shows all processes running on RON - can see who we can yell at for our commands taking long

q
# return to command line from top, less, and other programs but just use for top

SRR2589044_1.fastq  SRR2589044_2.fastq
# _1 indicated paired end reads

less -S filename.fastq
# shows read in a way that makes more sense 

How to find last read in a file 
# tail 
# tail SRR2584863_1.fastq | less -S

How big are these files?
# ls -l = additional information 
# ls -lhS OR ls -l -h = human readable information 

Conda activate genomics 
# activates environment with all the programs downloaded that we need for lab
## (base) kgb1008@ron:~/dc_workshop/data/untrimmed_fastq$ changed to (genomics)

fastqc -h
# list of commands we can do to fastqc files
# written in Java

-o 
# outdir: Creates all output files in the specified output directory. Please note that this directory must exist as the program will not create it.  If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.

cat */summary.txt
# shows all files there 

**** 03/03/23 NOTES:

can highlight lines in notebook and hit control enter, it will run the command in the terminal 

# realpath /tmp/fastqc_output/*fastqc/Images/*png

# get /tmp/fastqc_output/SRR2584863_2_fastqc/Images/per_sequence_quality.png
## this will move file from ron to desktop when you are in desktop directory (sftp kgb1008@ron.sr.unh.edu)

trimmomatic (JAVA language)q
# shows how to use the command 

trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
    SRR_1056_1.trimmed.fastq.gz SRR_1056_1un.trimmed.fastq.gz \
    SRR_1056_2.trimmed.fastq.gz SRR_1056_2un.trimmed.fastq.gz \
    ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20 
# threads allows us to run program on 4 different streams on processor to make things a lot faster 
# SRR_1056_1.fastq name of fastq we will trim 
# The middle lines tells computer where the output should be 
# last line tells computer to look in folder for adapters we just downloaded 
# slidingwindow tells computer how to report quality scores (slide across 4 bp and remove bases with quality score 20 or higher)

ls *gz
# shows zipped files 

cp /opt/anaconda/anaconda/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa .
# goes to folder of adapters with trimmomatic and copied those to terminal so we can remove the adapters from fastq files 

Exercise (when trimmomatic is done)
1. make a folder called ‘trimmed_fastq’ in data
    mkdir trimmed_fastq
2. copy all files that we made (hint: .trim and the move command, and the destination directory will work)
    cp /tmp/trimmed/*fastq.gz . 

**** 03/24 NOTES: Variant Calling - Shell Scripting 
Download Ref. Genome
# cd ~/dc_workshop
# mkdir -p data/ref_genome
## -p = creates pathway for new directory/folder
# curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz 
## -L redirect -o save to specific directory or folder
# gunzip data/ref_genome/ecoli_rel606.fasta.gz
Download more fastq files:
# curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
# tar xvf sub.tar.gz
## tar = untar? tar = archive file, untar = uncompress 
# mv sub/ ~/dc_workshop/data/trimmed_fastq_small
## moves sub directory/folder to new folder called trimmed_fastq_small 
Make Directories:
# mkdir -p results/sam results/bam results/bcf results/vcf
## move files into their own folder 
Activate Environment we will use for Analysis: CONDA
# conda env ls 
## see avalible environments on RON
# conda activate genomics-plus
## activate genomics-plus env 
Index the Reference Genome:
# bwa index data/ref_genome/ecoli_rel606.fasta
## bwa = alignment pathway, speeds up process of aligning reads to reference genome 
Aligning Reads to Ref Genome: 
# bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam
## 1st is ref genome file, 2nd is foward read (has _1 at the end), 3rd is reverse read (has _2 at the end)
## mem = memory saving 
## .aligned.sam is helpful to name files like this to keep track of what you have done to it (aligned and in sam format)
View what is in new folder made (results/sam/SRR2584866.aligned.sam)
# head or cat or less -S (wihtout -S, lines are not wrapped and looks better/easier to read)
To count number of lines in a file *********
# wc -l 
To get an acurate number of lines in a file ******
# grep -v '^@' results/sam/SRR2584866.aligned.sam | wc -l
## inverse grep - shows what doesnt match search text 
## ^@ looks for @ at the beginning of the line 
## | head = shows first few fine not all 35,000
## | wc -l shows the number of lines 
View the SAM/BAM format with samtools view:
# samtools view -s -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
## to manipulate alignnment files 
## -S = 
## -b =
## > = save into 
## converts sam file to bam file 
Sort BAM file by Coordinates with samtools sort:
# samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.sorted.bam 
## -o = output (save location)
Use flagstat to get stats on Alignment:
# samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
Variant Calling:
STEP 1: Calculate the read coverage of positions in the genome 
# bcftools mpileup -Ob -o results/bam/SRR2584866_raw.bcf -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
## -Ob = 
## 1st = file format and name to convert it to, 2nd = reference genome, 3rd = file of interest 

*** 03/31/23 NOTES - Finish Variant Calling and Visualization (SHELL SCRIPTING)
-start in dc_workshop 
STEP 1: Align yout reads to reference genome 
# bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam 
STEP 2: View the SAM/BAM format with samtools view 
# samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
STEP 3: Sort BAM file by coordinates 
# samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 
STEP 4: Activate Conda Environment 
STEP 5: Calculate the read coverage of positions in the genome 
# bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
STEP 6: Detect the single nucleotide variants
# bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf
STEP 7: Filter and report the SNV variants in variant calling format 
# vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf
STEP 8: Index the alignment 
# samtools index results/bam/SRR2584866.aligned.sorted.bam
Viewing with the alignment with tview 
# samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta
OUTPUT:
1         11        21        31        41        51        61        71        
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG
................................................................................
...............................................................................,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,                 ,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  ..............................................................................
                                ................................................
                                                                ,,,,,,,,,,,,,,,,
                                                                  ..............
EXERCISE - Viewing with IGV
1 Download reference genome .bams .bai and .vcf to desktop into a "files_for_igv" folder 
STEP 1: Sign into secure file transfer pathway
NEW TERMINAL 
# sftp kgb1008@ron.sr.unh.edu 
GET for ref genome
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.amb
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.ann
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.bwt
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.fai
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.pac
# /home/users/kgb1008/dc_workshop/data/ref_genome/ecoli_rel606.fasta.sa
STEP 1-4:
# mkdir ~/Desktop/files_for_igv
# cd ~/Desktop/files_for_igv
# sftp YOUR USERNAME@ron.sr.unh.edu
# get /home/unhAW/jtmiller/dc_workshop/data/ref_genome/ecoli_rel606.fasta
IN SSH terminal:
# realpath results/bam/SRR2584866.aligned.sorted.bam*
# realpath realpath results/vcf/SRR2584866_final_variants.vcf*
GET for .bam file 
# /home/users/kgb1008/dc_workshop/results/bam/SRR2584866.aligned.sorted.bam
# /home/users/kgb1008/dc_workshop/results/bam/SRR2584866.aligned.sorted.bam.bai
GET for .vcf
# /home/users/kgb1008/dc_workshop/results/vcf/SRR2584866_final_variants.vcf
IN IGV:
- use genome. load local file. grab both .fasta and .fasta.fai (index)
- use tracks. load local file. grab both .bam and .bam.bai (index)
- use tracks. load local file. grab .vcf file 




