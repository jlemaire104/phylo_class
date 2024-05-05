# Phylo Class 543 Reproducible Script

## Create multiple sequence alignments

### Prepare sequence files for sample data
1. Compiled sequence data by blasting 16S gene region of Aphanizomenon species and compiling a fasta file of top 10 species plus an outgroup of Syneccochocus found here: /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/16S/Aphani_fasta.fasta file.
    
2. Use the above file as the input to align sequences using both the MUSCLE program and the ClustalW program as specified below. 
   A. The Synnechococcus genus was added as it is a distint relative to Aphanizomenon but still a in the Cyanobacteria phylum.

## 1. Alignment Method: MUSCLE
Download [MUSLCE](https://github.com/rcedgar/muscle/releases/tag/5.1.0) from the Github repository here: [muscle5.1.macos_arm64](https://github.com/rcedgar/muscle/releases).

Note: I have a MacOsx computer with an M chip - may need to adjust the program download version depending on what computer you are currently working with

Troubleshooting: Had issues getting the program to install properly and be executed from the website. Instead used conda to install as follows:

### 1. Install Muscle using Conda
```Bash
conda install muscle
cd
mkdir bin
cp /Users/jacquelinelemaire/Downloads/muscle5.1.macos_arm64 bin/muscle.exe
muscle #check to see if it is successful by calling program
```


### 2. Run MUSCLE on data using input: 
```shell
muscle -align /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/16S/Aphani_fasta.fasta -output ~/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta
```

### Output: 

```
muscle 5.1.osx64 []  8.6Gb RAM, 8 cores
Built Feb 10 2022 07:52:24
(C) Copyright 2004-2021 Robert C. Edgar.
https://drive5.com

Input: 13 seqs, avg length 966, max 1401

00:00 6.0Mb  CPU has 8 cores, running 8 threads
00:01 696Mb   100.0% Calc posteriors
00:01 488Mb   100.0% Consistency (1/2)
00:01 492Mb   100.0% Consistency (2/2)
00:01 492Mb   100.0% UPGMA5
00:02 501Mb   100.0% Refining
```

The MUSCLE multiple sequence alignment stored in the file: ~/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta

## Alignment Method 2: ClustalW

### Install Alignment Programs

```shell
conda install -c bioconda 
clustalw
```

###Call the program to see if it works
```shell
clustalw
```

###Output
```
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
``````

## Running Clustalw Alignment on Aphani 16S sequences

Note: Edited the fasta file to add an underscore between fasta id, geneus and species names because otherwise clustal cuts it off after fasta id.

Troubleshooting: Not necessary for this small dataset but for larger datasets that take longer to run: to stop them from being interrupted the commands were run in a tmux terminal.

## Open a Virtual Terminal or run on a server
```shell
#create a new session
tmux new -s alignment_clustal

#to exit press ctl+b then d on the keyboard

#to reopen tmux terminal session 
tmux attach -t alignment_clustal
```

### Run ClustalW: 
```shell
clustalw -ALIGN -INFILE=/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/16S/Aphani.fasta -OUTFILE=/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/Aphani-aligned-clustal.fasta -OUTPUT=FASTA
```

### Output: 
```shell
There are 12 groups
Start of Multiple Alignment

Aligning...
Group 1: Sequences:   2      Score:13499
Group 2: Sequences:   3      Score:14394
Group 3: Sequences:   4      Score:15092
Group 4: Sequences:   5      Score:15454
Group 5: Sequences:   6      Score:13146
Group 6: Sequences:   2      Score:17736
Group 7: Sequences:   3      Score:17737
Group 8: Sequences:   2      Score:13490
Group 9: Sequences:   3      Score:15584
Group 10: Sequences:   6      Score:15995
Group 11: Sequences:  12      Score:13688
Group 12: Sequences:  13      Score:14568
Alignment Score 439162
firstres = 1 lastres = 1409
FASTA file created!

Fasta-Alignment file created:      
/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/Aphani-aligned-clustal.afa
``````

Open alignment in Rstudio and/or make a tree with the desired method. 

Start with a quick and dirty distance and neighbor joining tree


### Build and Plot Distance-based NJ Trees and Parsimony Trees

### Estimate a distance-based tree using software R package ape

* _ape_ is one of the most widely used phylogenetic software. It is an R package and it has a huge variety of functions
* In particular, it will be used here for distance-based tree estimation methods
  * [Full documentation](https://cran.r-project.org/web/packages/ape/index.html)
  * Procedure here based on this [tutorial](https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf)
* Summary of Software R package _ape_
  * Main distance functions:
    * `nj` (`ape` package): the classical Neighbor-Joining algorithm.
    * [`bionj`](https://www.rdocumentation.org/packages/ape/versions/5.4-1/topics/BIONJ) (`ape`): an improved version of Neighbor-Joining: [Gascuel 1997](https://pubmed.ncbi.nlm.nih.gov/9254330/). It uses information on variances of evolutionary distances
    * [`fastme.bal` and `fastme.ols`](https://www.rdocumentation.org/packages/ape/versions/5.4-1/topics/FastME) (`ape`): minimum evolution algorithms: [Desper and Gascuel, 2002](https://pubmed.ncbi.nlm.nih.gov/12487758/)
    * `hclust` (`stats`): classical hierarchical clustering algorithms including single linkage, complete linkage, UPGMA, and others.


## Within RStudio, install necessary packages:
   A. Install packages: 
```r
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
``` 

B. Load the packages
```r
library(ape)
library(adegenet)
library(phangorn)
```


C. Load in data from alignments: 
   1. Data input must be aligned sequences
   2. Use this command when reading in fasta format sequences
```r
Aphani_muscle <- fasta2DNAbin(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta")

Aphani_clustal <- fasta2DNAbin(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/Aphani-aligned-clustal.fasta")
```


D. Compute the genetic distances
  
```r
#Computing the genetic distances. They choose a Tamura and Nei 1993 model which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate.

#learn more about the distance program

?dist.dna

#Note model options:
#a character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", #"TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", #"paralin", "indel", or "indelblock".
```

```r
# Test the GG95 model 
M <- dist.dna(Aphani_muscle, model="GG95")
C <- dist.dna(Aphani_clustal, model="GG95")

# Test the TN93 model with gamma correction
Mg <- dist.dna(Aphani_muscle, model="TN93", gamma = TRUE)
Cg <- dist.dna(Aphani_clustal, model="TN93", gamma = TRUE)
```



E. Make Neighbor Joining Trees
```r
tree_aphani_muscle <- nj(M)

tree_aphani_clustal <- nj(C)

#Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted

tree_aphani_muscle <- ladderize(tree_aphani_muscle)
tree_aphani_clustal <- ladderize(tree_aphani_clustal)

#view tree info
tree_aphani_muscle
```
Output:
```
#Phylogenetic tree with 13 tips and 11 internal nodes.

#Tip labels:
#  FJ895118.1.Aphanizomenon aphanizomenoides, AF448070.1.Synechococcus_sp._PS845, #FJ895126.1.Aphanizomenon gracile, FJ895128.1.Aphanizomenon gracile, FJ895125.1.Aphanizomenon gracile, #FJ895119.1.Aphanizomenon aphanizomenoides, ...

#Unrooted; includes branch lengths.
```


F. Root all the trees to outgroup
```r
#Try Rooting the Trees

M.rooted <- root(tree_aphani_muscle, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)

C.rooted <- root(tree_aphani_clustal, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)

#check that the tree is now rooted
M.rooted
C.rooted
```

12. Plot the rooted trees
```r

#view the order of tip labels
#coloring the tip labels like this is not very scalable but it works for this small dataset
#Make a vector to color the different species
species_colors_1 <- c("darkblue", "black", "darkgreen", "darkgreen", "darkgreen", "darkblue","darkred", "darkblue", "darkblue", "darkgreen","darkblue", "darkblue", "darkgreen")

plot(M.rooted, cex=.7, edge.color = "black", tip.color = species_colors_1, main="A Neighbor Joining Tree of Aphanizomenon Muscle Alignment")

#view the order of tip labels
C.rooted$tip.label
species_colors_2 <- c("darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkred","darkblue", "darkblue", "darkblue", "darkblue","darkblue", "darkblue", "black")

plot(C.rooted, cex=.7, edge.color = "black", tip.color = species_colors_2, main="A Neighbor Joining Tree of Aphanizomenon Clustal Alignment")
```

### Estimate a parsimony-based tree using software: R package _phangorn_

* _phangorn_ is another widely used phylogenetic software. It is an R package and it has a huge variety of functions.
* Used here for parsimony-based tree estimation methods.
  * Full documentation](https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf)
  * Procedure here based on this [tutorial](https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf)
  * The commands are listed in the [PDF tutorial]((https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf)) that was used as guideline, or as seen in crsl4 class reproducible script [notebook-log.md](https://github.com/crsl4/phylogenetics-class/blob/master/exercises/notebook-log.md). 

A. Install packages
```r
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)
```

B. Load the packages
```r
library(ape)
library(adegenet)
library(phangorn)
```

C. Bring in the  data and convert to phangorn object:
```r
# Step 1: Import aligned sequences and convert them into a phyDat object
Aphani_muscle <- fasta2DNAbin(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta")

Aphani_clustal <- fasta2DNAbin(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/Aphani-aligned-clustal.fasta")

#convert to phangorn object
M2 <- as.phyDat(Aphani_muscle)

C2 <- as.phyDat(Aphani_clustal)
```

D. Create a starting tree and compute the parsimony score. Do this for both the MUSCLE and ClustalW files. 
```r
# Step 2: Construct a starting tree using neighbor-joining

tre.ini.M <- nj(dist.dna(Aphani_muscle,model="TN93", gamma = TRUE))

tre.ini.C <- nj(dist.dna(Aphani_clustal,model="TN93", gamma = TRUE))

#Step 3: Compute the parsimony score of the starting tree
parsimony(tre.ini.M, M2)
#227
parsimony(tre.ini.C, C2)
#225
```

E. Make the maximum parsimony tree:
```r
# Step 4: Search for the tree with maximum parsimony

tre.pars.M <- optim.parsimony(tre.ini.M, M2, method = "fitch", cost = NULL, trace = 1,
  rearrangements = "SPR")

#output:
## Final p-score 225 after 1 nni operations


tre.pars.C <- optim.parsimony(tre.ini.C, C2, method = "fitch", cost = NULL, trace = 1,
  rearrangements = "SPR")

#output:
## Final p-score 222 after  2 nni operations 
``````

F. Root the trees to outgroup 
```r
#Try Rooting the Trees

M.rooted.parsimony <- root(tre.pars.M, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)

C.rooted.parsimony <- root(tre.pars.C, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)

#check that the tree is now rooted
M.rooted.parsimony
C.rooted.parsimony

```
Output:
```
Phylogenetic tree with 13 tips and 11 internal nodes.

Tip labels:
  FJ895118.1_Aphanizomenon_aphanizomenoides, AF448070.1_Synechococcus_sp._PS845, FJ895126.1_Aphanizomenon_gracile, FJ895128.1_Aphanizomenon_gracile, FJ895119.1_Aphanizomenon_aphanizomenoides, EF685373.1_Aphanizomenon_issatschenkoi, ...

Rooted; no branch lengths.
```

G. Plot trees:
```r
#view the order of tip labels
M.rooted.parsimony$tip.label

#Make a vector to color the different species
species_colors <- c("darkblue", "black", "darkgreen", "darkgreen", "darkblue", "darkred","darkblue", "darkblue", "darkblue", "darkgreen","darkgreen", "darkblue", "darkgreen")

# Plot tree and Enhance Tree Visualization
plot(M.rooted.parsimony, cex = 0.7, edge.color = "black", tip.color = species_colors, main="Parsimony Tree of Aphanizomenon Muscle Alignment")