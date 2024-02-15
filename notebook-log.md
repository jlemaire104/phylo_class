A notebook log of steps down for the project for phylogenetics class

# Description of dataset for project:
# Data is from metagenomic shotgun sequencing done on temporal integrated water column samples collected from Lake Mendota Madison, WI

# The original non-redundant set included 3869 genomes but narrowed down to only genomes that were present in 4 or more samples at >2x coverage

# This data was subset to only include the Cyanobacteria phylum - excluding vampirovibrio

# Mapping was done with bowtie2 and parsed with coverm (97%-id and 75% min alignment) using bowtie2 with default settings

# These are .fna files that are filtered and trimmed already - can I still run QC on them? Or use fastqc on the mapped Bam files?


## Running things on GLBRC server Scarcity
# 1. Log-in to WEI VPN
## 2. Log in to scarcity server

ssh jlemaire@scarcity-submit.glbrc.org
# enter PW

## Running quality checks with FastQC - ??Need to determine where/how to use this on glbrc - permission was denied
code:
fastqc seqfile1 seqfile2 .. seqfileN

    # fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN

## Install Alignment Programs
conda install -c bioconda clustalw

# Call the program to see if it works
clustalw

# **************************************************************
 ******** CLUSTAL 2.1 Multiple Sequence Alignments  ********
 **************************************************************


     1. Sequence Input From Disc
     2. Multiple Alignments
     3. Profile / Structure Alignments
     4. Phylogenetic trees

     S. Execute a system command
     H. HELP
     X. EXIT (leave program)
# ***************************************************************     

## run clustalw on test data
# Move into folder with data
cd Documents/phylogenetics-class/data

# Count how many sequences there are
grep ">" primatesAA.fasta | wc -l

# Run clustalW


# Tried to install Muscle but it doesn't seem to be working with the new MacOS M chips
