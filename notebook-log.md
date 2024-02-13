A notebook log of steps down for the project for phylogenetics class

# Description of dataset for project:
# Data is from metagenomic shotgun sequencing done on temporal integrated water column samples collected from Lake Mendota Madison, WI

# The original non-redundant set included 3869 genomes but narrowed down to only genomes that were present in 4 or more samples at >2x coverage

# This data was subset to only include the Cyanobacteria phylum - excluding vampirovibrio

# Mapping was done with bowtie2 and parsed with coverm (97%-id and 75% min alignment) using bowtie2 with default settings

# These are .fna files that are filtered and trimmed already - can I still run QC on them? Or use fastqc on the mapped Bam files?


##Running things on GLBRC server Scarcity
# 1. Log-in to WEI VPN
## 2. Log in to scarcity server

ssh jlemaire@scarcity-submit.glbrc.org
# enter PW

## Running quality checks with FastQC - ??Need to determine where/how to use this on glbrc - permission was denied
# fastqc seqfile1 seqfile2 .. seqfileN

    # fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN




