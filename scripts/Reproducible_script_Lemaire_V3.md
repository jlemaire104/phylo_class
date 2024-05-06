
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
```

## Running Clustalw Alignment on Aphani 16S sequences

Note: Edited the fasta file to add an underscore between fasta id, geneus and species names because otherwise clustal cuts it off after fasta id.

Troubleshooting: Not necessary for this small dataset but for larger datasets that take longer to run: to stop them from being interrupted the commands were run in a tmux terminal.

### Open a Virtual Terminal or run on a server
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
```

Open alignment in Rstudio and/or make a tree with the desired method. 

Start with a quick and dirty distance and neighbor joining tree


### Build and Plot Distance-based NJ Trees and Parsimony Trees
Estimate a distance-based tree using software R package ape

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


### Within RStudio, install necessary packages:
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
```

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
```

## Maximum Likelihood Tree Esimation Using RaxML

### Software
```
RAxML-NG

```
### Description
```
Maximum likelihood Next Generation - tree inference tool last updated in Nov 2019
```
### Strengths
```
Fast, good scalability, scores higher when alignments are taxon-rich, good for large datasets due to it reduces memory space and run time
```
### Weaknesses
```
Not as efficient or high likelihoods as IQ-Tree2
```
### Assumptions
```
Iteratively searches by series of Subtree Pruning and Regrafting (SPR) moves
```
### User Choices
```
User must provide outgroup to root data, can select from 3 starting tree options (random topology, tree generated by parsimony-based randomized stepwise addition algorithm, or user-defined provided tree), set the model to run, and can set the seed to create reproducible analysis. User can also choose the number of bootstraps and the boostrap support values used in this new version fbp or tbe or both can be specified
```

###Installing RaxML-NG from Website

Download raxml-ng from here:(https://github.com/amkozlov/raxml-ng). You get a zipped folder: raxml-ng_v1.0.2_macos_x86_64 which I placed in my /Applications folder

We are following HAL 1.3.

Reminders 
1. Branch lengths do not have to do with time
2. need outgroup in order to root tree - otherwise ML will root randomly among the taxa presented

Checking the version
```shell
cd /Applications/raxml-ng_v1.2.1_macos_x86_64

./raxml-ng -v
```
Output:
```
RAxML-NG v. 1.2.1 released on 22.12.2023 by The Exelixis Lab.

Developed by: Alexey M. Kozlov and Alexandros Stamatakis.

Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth, Julia Haag, Anastasis #Togkousidis.

Latest version: https://github.com/amkozlov/raxml-ng

Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml
```
Before we get started, let's first check that the MSA can actually be read and doesn't contain sites with only undetermined characters or sequences with undetermined characters or duplicate taxon names, etc
```shell
./raxml-ng --check --msa /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta --model GTR+G --prefix T1
```

Output:
```
WARNING: Sequences FJ895121.1_Aphanizomenon_aphanizomenoides and FJ895120.1_Aphanizomenon_aphanizomenoides are exactly identical!

WARNING: Duplicate sequences found: 3

NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed)
NOTE: was saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T1.raxml.reduced.phy

Alignment comprises 1 partitions and 1409 sites

Partition 0: noname
Model: GTR+FO+G4m
Alignment sites: 1409
Gaps: 31.43 %
Invariant sites: 87.22 %

Alignment can be successfully read by RAxML-NG.

Execution log saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T1.raxml.log

Analysis started: 27-Apr-2024 17:31:16 / finished: 27-Apr-2024 17:31:16

Elapsed time: 0.013 seconds
```
### Troubleshooting and Optimizing
Run raxml on the new dereplicated file and try running with more starting trees to explore the tree space more

```shell
./raxml-ng  --msa T1.raxml.reduced.phy --model GTR+G --prefix T5 --threads 2 --seed 2 --tree pars{25},rand{25}
```

```shell
#Let us now compare the results of both alternative tree inference runs:
grep "Final LogLikelihood:" T{4,5}.raxml.log

#output:
T4.raxml.log:Final LogLikelihood: -2973.309591
T5.raxml.log:Final LogLikelihood: -2973.308998

#looks good, they are the same!

#inferring bootstrap trees
./raxml-ng --bootstrap --msa T1.raxml.reduced.phy --model GTR+G --prefix T6 --seed 2 --threads 2

#check to see if the trees converged
./raxml-ng --bsconverge --bs-trees T6.raxml.bootstraps --prefix T9 --seed 2 --threads 2 --bs-cutoff 0.01

##Bootstopping test did not converge after 1000 trees

#So maybe I should increase the bootstrapping?
./raxml-ng --bsconverge --bs-trees T6.raxml.bootstraps --prefix T9 --seed 2 --threads 2 --bs-cutoff 0.01

#Bootstopping test did not converge after 1000 trees
#This looks promising, and we
#can expect convergence after few hundred replicates. Luckily, bootstraps are independent,
#and we can thus reuse the 200 BS trees we have already inferred. So let’s add 400 additional
#BS replicate trees.
#IMPORTANT NOTE: It is extremely important to specify a distinct random seed for the second run, otherwise first 200 trees of the second run will be identical to the first run!

./raxml-ng --bootstrap --msa T1.raxml.reduced.phy --model GTR+G --prefix T10 --seed 333
--threads 2 --bs-trees 400

#Now, we can simply concatenate the BS replicate trees from both runs, and re-assess the
#convergence:

cat T6.raxml.bootstraps T10.raxml.bootstraps > allbootstraps
./raxml-ng --bsconverge --bs-trees allbootstraps --prefix T11 --seed 2 --threads 1 --bs-cutoff 0.01

#Bootstopping test did not converge after 1400 trees

#Try doubling this and adding some because we were only about half way to convergence percentage

./raxml-ng --bootstrap --msa T1.raxml.reduced.phy --model GTR+G --prefix T12 --seed 3330 --threads 2 --bs-trees 4000
```
## Troubleshooting and optimizing 
```shell
#an easier way to do all in one:

./raxml-ng --all --msa T1.raxml.reduced.phy --model GTR+G --prefix T14 --seed 2 --threads 2 --bs-trees 4000

#check to see if trees converged this time
./raxml-ng --bsconverge --bs-trees T14.raxml.bootstraps --prefix T15 --seed 2 --threads 2 --bs-cutoff 0.01
```
### Final Version of Command Used
```shell
# Supposedly with --all we can do ML search and boostrap analysis all one step and add --bs-metric to map the bootstrap values on the ML tree

./raxml-ng --all --msa T1.raxml.reduced.phy --model GTR+G --prefix T15 --seed 2 --threads 2 --bs-metric fbp,tbe
```
Output:
```
#This appears to run 1000 bootstraps

Final LogLikelihood: -2973.309591

AIC score: 5998.619182 / AICc score: 5999.635101 / BIC score: 6135.135706
Free parameters (model + branch lengths): 26

WARNING: Best ML tree contains 1 near-zero branches!

Best ML tree with collapsed near-zero branches saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.bestTreeCollapsed
Best ML tree saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.bestTree
All ML trees saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.mlTrees
Best ML tree with Felsenstein bootstrap (FBP) support values saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.supportFBP
Best ML tree with Transfer bootstrap (TBE) support values saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.supportTBE
Optimized model saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.bestModel
Bootstrap trees saved to: /Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.bootstraps
```
Evaluate the Output Tree
```shell
#Was it worth doing ML - evaluate the trees?

./raxml-ng --evaluate --msa T1.raxml.reduced.phy --model GTR+G -prefix T16 --tree T15.raxml.bestTree --threads 2 
```
Output:
```
Bootstopping test did not converge after 4000 trees ?
Is it getting stuck on a hill? Different peaks with same liklihood?
May need to trouble shoot this more in the future
``````

Move files to my results directory
```shell
mv T13.* /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/raxml/aphani
```
Let's look at the tree in R and in FigTree
file = T13.raxml.bestTree
file new tree with bs = T15.raxml.supportFBP
Note: Bootstrap support should be mapped onto the tree already as node label in this T15.raxml.supportFBP

```r
library(ape)
#The output file is a newick tree so use read.tree to open the file

t <- read.tree(file = "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/raxml/aphani/T13.raxml.bestTree")
T2 <- read.tree(file = "/Applications/raxml-ng_v1.2.1_macos_x86_64/T15.raxml.supportFBP") 
```
Root the Tree to the Outgroup
```r
#Try Rooting the Trees

raxml.rooted <- root(T2, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)


#check that the tree is now rooted
raxml.rooted
```
Output:
```
Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
  FJ895127.1_Aphanizomenon_gracile, FJ895128.1_Aphanizomenon_gracile, FJ895126.1_Aphanizomenon_gracile, EF685373.1_Aphanizomenon_issatschenkoi, FJ895119.1_Aphanizomenon_aphanizomenoides, FJ895123.1_Aphanizomenon_aphanizomenoides, ...
Node labels:
  Root, 97, 88, 50, 72, 41, ...

Rooted; includes branch lengths.
```
Plot the Tree
```r
#view the order of tip labels
raxml.rooted$tip.label
raxml.rooted$node.label
species_colors_3 <- c("darkgreen", "darkgreen", "darkgreen", "darkred", "darkblue", "darkblue","darkblue", "darkblue", "darkblue", "black")

plot(raxml.rooted, cex=0.9, edge.color = "black", tip.color = species_colors_3, show.node.label = TRUE, main="A Maximum Likelihood RaxML Tree of Aphanizomenon Muscle Alignment")
```
Note: The BS values are printed on tree but hard to see. Try to troubleshoot this. Maybe try ggtree in the future, but for now view tree in Figtree




# Maximum Likelihood Tree Esimation Using IQ-Tree

1. First run a trial with model finder option to find the best fit model of evolution to use.

   Result: The predicted best fit model for this data was TN+F+G4
2. Run a trial with 5000 bootstraps on the MUSCLE aligned file using the best fit model.
3. Used a high amount of BS trees because RaxML said it was not converging with a low number. This might be overkill, but since this is such a small dataset it is still very fast to make a lot of BS trees. Will need to reduce this number or using automatic bootstopping for larger dataset. 

### Software
```
IQTree2
Download from the website for MacOSX (http://www.iqtree.org/)
```
## Description
```IQ-Tree is a fast and effective stochastic algorithm to infer phylogenetic trees by maximum likelihood
Efficient sampling of local optima in tree space is a core function of IQ-Tree
A combination of hill-climbing and stochastic perturbation optimizes run-time efficiency
Claims higher likelihoods are achieved relative to RAxML and PhyML
Note: all programs used the General Time Reversible (GTR) model of evolution

IQ-Tree combines elements of hill-climbing algorithms, random perturbation of current best trees, and a broad sampling of initial starting trees 
```
### Strengths
```
IQ-TREE finds more higher-likelihood trees than RAxML or PhyML from DNA alignments 
Obtains higher likelihoods than RAxML for 73% of AA alignments 
Flexibility in their user settings (customization) 
Effective tree search algorithm  
Employs small population of candidate trees
Wide collection of models (non-reversible and reversible)
```

### Weaknesses
```
Requires longer CPU times than RAxML for 75% of DNA alignments tested 
Situation is complicated; differences in average CPU times are highly variable
To finish 10 repetitions for 70 tested DNA alignments, IQ-TREE needed 2,020 CPU hours (∼87 CPU days)
VS. RAxML needed 1,870 CPU hours (∼78 CPU days).
```

### Assumptions
```
Some common assumptions include treelikeness (all sites in the alignment have evolved under the same tree), stationarity (nucleotide/amino-acid frequencies remain constant over time), reversibility (substitutions are equally likely in both directions), and homogeneity (substitution rates remain constant over time).
```

### User choices
```
general options:

Option	Usage and meaning
-h or -?	Print help usage.
-s	Specify input alignment file in PHYLIP, FASTA, NEXUS, CLUSTAL or MSF format.
-st	Specify sequence type as either of DNA, AA, BIN, MORPH, CODON or NT2AA for DNA, amino-acid, binary, morphological, codon or DNA-to-AA-translated sequences. This is only necessary if IQ-TREE did not detect the sequence type correctly. Note that -st CODON is always necessary when using codon models (otherwise, IQ-TREE applies DNA models) and you also need to specify a genetic code like this if differed from the standard genetic code. -st NT2AA tells IQ-TREE to translate protein-coding DNA into AA sequences and then subsequent analysis will work on the AA sequences. You can also use a genetic code like -st NT2AA5 for the Invertebrate Mitochondrial Code (see genetic code table).
-t	Specify a file containing starting tree for tree search. The special option -t BIONJ starts tree search from BIONJ tree and -t RANDOM starts tree search from completely random tree. DEFAULT: 100 parsimony trees + BIONJ tree
-te	Like -t but fixing user tree. That means, no tree search is performed and IQ-TREE computes the log-likelihood of the fixed user tree.
-o	Specify an outgroup taxon name to root the tree. The output tree in .treefile will be rooted accordingly. DEFAULT: first taxon in alignment
-pre	Specify a prefix for all output files. DEFAULT: either alignment file name (-s) or partition file name (-q, -spp or -sp)
-nt	Specify the number of CPU cores for the multicore version. A special option -nt AUTO will tell IQ-TREE to automatically determine the best number of cores given the current data and computer.
-ntmax	Specify the maximal number of CPU cores -nt AUTO is allowed to allocate DEFAULT: #CPU cores on the current machine
-seed	Specify a random number seed to reproduce a previous run. This is normally used for debugging purposes. DEFAULT: based on current machine clock
-v	Turn on verbose mode for printing more messages to screen. This is normally used for debugging purposes. DEFAULT: OFF
-quiet	Silent mode, suppress printing to the screen. Note that .log file is still written.
-keep-ident	Keep identical sequences in the alignment. Bu default: IQ-TREE will remove them during the analysis and add them in the end.
-safe	Turn on safe numerical mode to avoid numerical underflow for large data sets with many sequences (typically in the order of thousands). This mode is automatically turned on when having more than 2000 sequences.
-mem	Specify maximal RAM usage, for example, -mem 64G to use at most 64 GB of RAM. By default, IQ-TREE will try to not to exceed the computer RAM size.
```

## Run IQTree Model Finder 
```shell
#Run IQtree with the -m option to use the model finder option
cd /Applications/iqtree-2.2.2.6-MacOSX

bin/iqtree2 -s /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta -m MF
```

Output:
```
Output:
open file to find the best model
more /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta.iqtree
Best-fit model according to BIC: TN+F+G4
```

## Run IQ Tree on Muscle Alignment
```shell
#try running again with the model picked and with bootstraps

bin/iqtree2 -s /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta -m TN+F+G4 --prefix Aph5000 -B 5000
```

Output:
```
Total number of iterations: 103
CPU time used for tree search: 1.041 sec (0h:0m:1s)
Wall-clock time used for tree search: 0.872 sec (0h:0m:0s)
Total CPU time used: 1.147 sec (0h:0m:1s)
Total wall-clock time used: 0.992 sec (0h:0m:0s)

Computing bootstrap consensus tree...
Reading input file T2.splits.nex...
12 taxa and 89 splits.
Consensus tree written to T2.contree
Reading input trees file T2.contree
Log-likelihood of consensus tree: -2977.134

Analysis results written to:
  IQ-TREE report:                T2.iqtree
  Maximum-likelihood tree:       T2.treefile
  Likelihood distances:          T2.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          T2.splits.nex
  Consensus tree:                T2.contree
  Screen log file:               T2.log
```


## Plot the IQ-Tree Output 

1. Within RStudio, install packages:
```r
install.packages("adegenet", dep=TRUE)
```


2. Load the packages:
```r
library(ape)
library(adegenet)
```

3. Read in the best tree from the IQTree output:
```r
iqtree_aphani <- read.tree(file= "/Applications/iqtree-2.2.2.6-MacOSX/Aph5000.treefile")
```

4. Root the trees to the outgroup:
```r
#Try Rooting the Trees

iqtree.rooted <- root(iqtree_aphani, outgroup = "AF448070.1_Synechococcus_sp._PS845", resolve.root = TRUE)


#check that the tree is now rooted
iqtree.rooted
```
Output:
```
Phylogenetic tree with 13 tips and 12 internal nodes.

Tip labels:
  FJ895118.1_Aphanizomenon_aphanizomenoides, AF448070.1_Synechococcus_sp._PS845, FJ895126.1_Aphanizomenon_gracile, FJ895124.1_Aphanizomenon_gracile, FJ895125.1_Aphanizomenon_gracile, FJ895128.1_Aphanizomenon_gracile, ...
Node labels:
  Root, , 50, 45, 98, 89, ...

Rooted; includes branch lengths.
```


5. Plot the tree:

Troubleshooting: Need to figure out what output file to load in and/or what value to call to add bootstrap values to tree

```r
#view the order of tip labels
iqtree.rooted$tip.label
iqtree.rooted$node.label

#coloring the tip labels like this is not very scalable but it works for this small dataset
#Make a vector to color the different species
species_colors_1 <- c("darkblue", "black", "darkgreen", "darkgreen", "darkgreen", "darkgreen","darkgreen", "darkred", "darkblue", "darkblue","darkblue", "darkblue", "darkblue")

plot(iqtree.rooted, cex=.7, edge.color = "black", tip.color = species_colors_1, show.node.label=TRUE, main="A Maximum Likelihood IQTree of Aphanizomenon Muscle Alignment") ##This is not working

#Troubleshooting: Visualize tree in fig tree program instead
```


# Bayesian Inference of Phylogeny Using MrBayes

### Software	
```
Mr. Bayes
Bayesian inference of phylogenetic trees
```

### Description	
```
The program MRBAYES performs Bayesian inference of phylogeny using a variant of Markov chain Monte Carlo. The Markov Chain Monte Carlo (MCMC) then approximates  the posterior probability of trees.In Bayesian analysis, inferences of phylogeny are based upon the posterior probabilities of phylogenetic trees using the Bayes theorem.
```

### Strengths
```
Bayesian inference has several advantages over other phylogenetic inference methods such as easy interpretation of results, ability to incorporate prior information, and computational advantages. The same models of DNA substitution used in Maximum liklihood analyses can be used in a Bayesian analysis of phylogeny. Mr.Bayes not only implements MCMC but also a variant called metropolis-coupled MCMC (MC^3) which are heated: Uses cold and heated chains to explore the space of phylogenetic trees, Can easily explore the space of phylogenetic trees, Mixing these two chains was dramatically improved using MC3.
```

### Weaknesses
```
It can be slow and computationally demanding,
The user impacts are not very clear and it can be confusing as to what to set user options to.
```	
### Assumptions
```
Mr. Bayes uses MCMC to appoximate the posterior probabilites of trees. Uses command line interface. Reads in the standard NEXUS format
```

### User choices
```
User can change assumptions of the substitution model, the prior, and details for the MC3.
```

## Mr Bayes Install

#Run in bash terminal:
```shell
brew install mrbayes --with-open-mpi
# This did not work, got an error:
Error: invalid option: --with-open-mpi

brew reinstall mrbayes
#Error: Cannot install under Rosetta 2 in ARM default prefix (/opt/homebrew)!
#I think this is an error because of the Mac M2 chip

#Try this:
arch -arm64 brew install mrbayes
#It started to install and then this error popped up: Error: gcc: the bottle needs the Apple Command Line Tools to be installed.

#Install command line tools with xcode:
xcode-select --install

#Try to install again:
arch -arm64 brew install mrbayes

##Success!!
which mb
# /opt/homebrew/bin/mb
```

Convert the msa fasta file to nexus file 
```r
msa = read.FASTA("/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.fasta")

write.nexus.data(msa, file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-aligned-muscle.nex", format="dna")
```

2. Create a text file with the mrbayes commands:
```shell
nano #open text editor
#Note that the commands `mcmc;sumt;` must be present at the end so that the mb block is executed.'mcm' runs MCMC and 'sumt' is the command to obtain a summary tree. 
#Copy the mrbayes commands into the text file and adjust parameters as needed:

begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=2 rates=gamma ngammacat=4;
 mcmcp ngen=1000000 samplefreq=10 printfreq=1000 nruns=3 nchains=3 savebrlens=yes;
 outgroup AF448070.1_Synechococcus_sp._PS845;
 mcmc;
 sumt;
end;

#Save as mbblock3.txt
```
   B. Notes on the user choices:
* prset= the priors - which were left ast the default setting for now   
* ngen=1000000 upped this number to tell MrBayes that its robots should each take 1 million steps.
* samplefreq=100 says to only save parameter values and the tree topology every 100 steps.
* printfreq=1000 specifies a progress report every 1000 steps.
* nruns=3 says to just do three independent runs. MrBayes performs two separate analyses by default.
* nchains=3 says that we would like to have 2 heated chains running in addition to the cold chain.
* savebrlens=yes tells MrBayes that we would like it to save branch lengths when it saves the sampled tree topologies.
* outgroup specifies the outgroup to root the tree. 

3. Append the MrBayes block to the end of the nexus files with the data: 
```shell
cat aphani-aligned-muscle.nex /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/scripts/MrBayes/mbblock3.txt > aphani-mb.nex
```

Note: the text of the mbblock3.txt was manually copy and pasted in the bottom of the relevant .nex file. - troubleshoot this issue in the future 

## Run MrBayes:
```shell
mb aphani-mb.nex
```
Output
```
      751000 -- (-3006.036) [-3014.156] (-3001.032) * (-3006.524) [-3001.163] (-3003.739) * (-3004.944) (-3005.384) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      752000 -- (-3002.799) [-3011.061] (-3012.763) * (-2997.738) [-2994.636] (-3011.696) * [-2999.286] (-3006.047) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      753000 -- (-3002.299) [-2997.475] (-3004.641) * [-3001.976] (-3003.403) (-3002.688) * [-3000.645] (-3007.727) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      754000 -- (-3003.695) (-3006.664) [-3001.811] * [-2996.458] (-3005.093) (-3004.084) * [-2994.160] (-3009.436) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      755000 -- (-3003.564) (-3002.840) [-3003.430] * [-2993.963] (-2997.950) (-3002.939) * [-3002.827] (-3002.869) [...1 more local chains...] (...0 remote chains...) -- 0:00:26

      Average standard deviation of split frequencies: 0.003967

      756000 -- (-3005.876) [-2999.516] (-2999.951) * [-2997.319] (-3010.361) (-3005.820) * (-3002.128) (-3012.724) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      757000 -- (-3005.675) (-2998.333) [-2998.192] * (-2998.319) (-3004.156) [-2996.006] * (-2995.393) [-2997.502] [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      758000 -- (-3011.607) [-2997.102] (-3004.772) * (-3005.275) [-2991.495] (-2996.139) * [-3002.989] (-2998.725) [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      759000 -- (-3008.300) [-2994.740] (-3014.189) * (-2997.465) (-3005.659) [-3000.152] * (-3009.448) [-2992.414] [...1 more local chains...] (...0 remote chains...) -- 0:00:26
      760000 -- (-3012.642) [-3006.583] (-2997.996) * (-3005.989) (-3015.443) [-3001.863] * (-3021.554) (-2998.853) [...1 more local chains...] (...0 remote chains...) -- 0:00:25

      Average standard deviation of split frequencies: 0.003848

      761000 -- (-3013.795) [-2995.805] (-3009.188) * (-3006.019) [-3011.588] (-2995.393) * (-3019.149) [-3001.154] [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      762000 -- (-3000.008) [-3005.206] (-3006.563) * (-2998.011) [-3005.765] (-3006.840) * (-2997.536) [-3007.979] [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      763000 -- (-3006.636) [-3004.390] (-2999.439) * [-3000.441] (-3011.530) (-2998.539) * (-2995.937) [-3007.032] [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      764000 -- [-3000.928] (-3006.194) (-3015.477) * (-2998.175) [-2994.207] (-3013.750) * [-3002.285] (-2992.921) [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      765000 -- (-3014.738) (-2998.751) [-3005.426] * [-2999.027] (-2999.291) (-3015.542) * (-3002.394) [-2999.747] [...1 more local chains...] (...0 remote chains...) -- 0:00:25

      Average standard deviation of split frequencies: 0.003866

      766000 -- (-2997.874) [-2994.744] (-3008.997) * (-2996.623) (-3009.607) [-2999.071] * (-3003.885) (-3010.633) [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      767000 -- [-2996.803] (-3000.989) (-2998.836) * [-3011.009] (-3013.220) (-3003.495) * [-2996.788] (-2993.435) [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      768000 -- (-2992.364) [-3005.165] (-3003.279) * (-3002.184) (-3007.303) [-2998.633] * [-3002.698] (-3007.687) [...1 more local chains...] (...0 remote chains...) -- 0:00:25
      769000 -- (-3002.309) [-2999.206] (-3005.418) * (-3005.726) (-3014.196) [-2992.583] * (-3002.971) [-2992.022] [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      770000 -- (-3005.398) [-3000.444] (-3015.744) * (-3001.475) (-3007.205) [-3001.618] * (-3004.464) [-2997.354] [...1 more local chains...] (...0 remote chains...) -- 0:00:24

      Average standard deviation of split frequencies: 0.003880

      771000 -- [-2999.074] (-3010.849) (-3031.298) * (-3006.136) (-3007.047) [-2992.057] * (-3001.486) [-3000.081] [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      772000 -- (-3008.379) (-3005.652) [-3004.978] * [-3001.433] (-2999.124) (-3003.744) * (-2998.985) (-3014.908) [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      773000 -- [-3006.475] (-3003.383) (-2998.432) * (-3013.386) (-3005.266) [-3002.571] * (-3015.314) [-2998.821] [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      774000 -- (-3005.411) (-3009.117) [-3005.593] * [-3006.591] (-3009.386) (-3004.901) * [-2997.950] (-3004.242) [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      775000 -- (-2999.122) (-3003.850) [-2996.601] * [-2995.508] (-3009.125) (-3008.299) * (-2999.682) (-3009.659) [...1 more local chains...] (...0 remote chains...) -- 0:00:24

      Average standard deviation of split frequencies: 0.003874

      776000 -- [-3005.129] (-3007.382) (-3000.989) * (-3001.458) [-3002.088] (-3004.626) * (-3006.541) [-2996.043] [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      777000 -- (-3010.424) [-3001.527] (-3015.963) * (-2998.961) (-3001.910) [-2995.319] * (-3006.470) (-3002.254) [...1 more local chains...] (...0 remote chains...) -- 0:00:24
      778000 -- [-3003.272] (-3017.226) (-3005.565) * [-3009.664] (-3019.646) (-3002.520) * (-3010.423) (-3008.944) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      779000 -- [-2999.125] (-3009.895) (-3011.736) * (-3006.358) [-2999.749] (-3007.862) * (-3006.558) (-3010.788) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      780000 -- [-3005.615] (-3005.245) (-3007.026) * (-3001.963) (-2996.839) [-2999.458] * [-2996.295] (-3004.216) [...1 more local chains...] (...0 remote chains...) -- 0:00:23

      Average standard deviation of split frequencies: 0.003850

      781000 -- (-3003.228) [-2998.007] (-3003.749) * (-3016.569) (-3005.031) [-2996.535] * (-2998.991) (-3002.671) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      782000 -- (-3015.888) (-3010.304) [-2996.896] * (-3001.575) (-3005.470) [-3002.297] * (-2997.050) [-3002.423] [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      783000 -- [-3004.736] (-3001.601) (-2993.851) * [-2998.265] (-2998.763) (-2994.748) * [-2999.996] (-3005.868) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      784000 -- [-2994.317] (-2998.644) (-2997.656) * [-2994.504] (-3004.509) (-3000.971) * (-3015.959) (-2995.547) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      785000 -- (-2993.708) [-2995.423] (-2998.863) * (-3000.482) (-3003.407) [-2999.094] * (-2999.117) [-2992.249] [...1 more local chains...] (...0 remote chains...) -- 0:00:23

      Average standard deviation of split frequencies: 0.003874

      786000 -- (-3002.554) (-3003.972) [-2999.399] * [-3002.363] (-3018.062) (-3002.702) * (-3008.169) [-3007.063] [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      787000 -- (-3000.832) (-2999.616) [-3006.139] * (-3012.823) (-3004.475) [-3005.676] * (-3004.570) (-3001.261) [...1 more local chains...] (...0 remote chains...) -- 0:00:23
      788000 -- (-3002.957) [-2996.861] (-3001.172) * (-3010.309) (-3000.368) [-2999.697] * (-3002.302) [-3004.095] [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      789000 -- (-2999.940) [-3005.503] (-3006.678) * [-3001.808] (-2996.498) (-2998.148) * (-3007.211) [-2998.934] [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      790000 -- [-3004.228] (-3010.115) (-3013.198) * [-3006.668] (-3001.197) (-3007.781) * (-3000.768) [-2999.815] [...1 more local chains...] (...0 remote chains...) -- 0:00:22

      Average standard deviation of split frequencies: 0.003757

      791000 -- [-3006.996] (-3002.513) (-3003.629) * [-3002.019] (-3001.528) (-3009.225) * (-3007.174) [-3007.758] [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      792000 -- (-3020.828) (-3012.334) [-2998.671] * (-3000.168) (-2996.068) [-3002.043] * (-3006.819) [-2999.300] [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      793000 -- (-2999.427) (-3010.678) [-2997.289] * [-2998.951] (-3011.960) (-2999.669) * (-3001.220) [-3006.922] [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      794000 -- (-3016.054) [-3003.777] (-2997.779) * (-3009.428) [-2997.802] (-3009.345) * [-3007.980] (-3002.220) [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      795000 -- (-3005.080) (-3005.614) [-3008.823] * (-3004.026) [-3000.319] (-2999.229) * (-3000.719) (-3011.197) [...1 more local chains...] (...0 remote chains...) -- 0:00:22

      Average standard deviation of split frequencies: 0.003597

      796000 -- (-3000.225) [-2997.888] (-2997.114) * (-2997.914) [-3007.129] (-3009.497) * [-3003.959] (-3008.585) [...1 more local chains...] (...0 remote chains...) -- 0:00:22
      797000 -- [-3003.058] (-2995.672) (-3010.151) * [-2997.919] (-2998.315) (-3007.058) * (-3002.495) [-3005.847] [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      798000 -- [-3000.500] (-3007.571) (-3010.632) * (-2998.748) [-3006.969] (-2998.492) * (-3006.262) (-3010.462) [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      799000 -- [-3003.412] (-2996.826) (-3000.864) * [-3000.460] (-3000.988) (-3000.158) * (-2996.566) [-2996.743] [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      800000 -- (-2994.043) [-3006.368] (-2992.022) * [-3002.016] (-3007.092) (-3000.467) * (-3002.717) (-3006.626) [...1 more local chains...] (...0 remote chains...) -- 0:00:21

      Average standard deviation of split frequencies: 0.003543

      801000 -- (-3008.086) (-3005.614) [-2992.817] * [-3003.263] (-3004.889) (-3004.892) * [-3005.664] (-3018.077) [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      802000 -- (-3001.871) [-3002.283] (-3000.420) * (-2994.076) (-3008.106) [-2999.544] * (-2999.707) (-3000.368) [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      803000 -- (-3003.039) [-2991.976] (-3006.427) * (-3008.482) [-3005.266] (-3001.305) * [-2997.973] (-2997.693) [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      804000 -- (-3006.545) (-3009.524) [-2996.545] * (-3003.393) (-3006.998) [-2990.474] * (-2997.396) (-2995.064) [...1 more local chains...] (...0 remote chains...) -- 0:00:21
      805000 -- [-3000.095] (-3007.644) (-3012.794) * [-2998.427] (-3014.298) (-3001.291) * [-3000.677] (-2999.516) [...1 more local chains...] (...0 remote chains...) -- 0:00:21

      Average standard deviation of split frequencies: 0.003432

      806000 -- (-3005.022) (-3003.216) [-2996.170] * (-3006.222) [-3002.767] (-3001.475) * [-3008.615] (-3005.108) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      807000 -- [-3001.319] (-3008.641) (-3025.112) * (-3017.428) (-2996.237) [-3001.682] * [-3007.601] (-3006.182) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      808000 -- (-3000.487) [-2998.341] (-2997.693) * (-3001.835) (-2998.450) [-2991.652] * [-2994.721] (-3009.355) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      809000 -- (-3005.832) [-3008.114] (-3014.489) * (-3003.807) (-2999.025) [-3004.510] * (-3002.368) (-3007.765) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      810000 -- [-2997.687] (-3001.854) (-3002.746) * (-3013.449) (-3009.107) [-3008.679] * (-3000.845) [-3003.748] [...1 more local chains...] (...0 remote chains...) -- 0:00:20

      Average standard deviation of split frequencies: 0.003364

      811000 -- [-3012.525] (-3001.578) (-2997.748) * (-3004.215) [-2996.252] (-3005.417) * (-3002.363) [-3009.535] [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      812000 -- (-3001.787) [-2998.202] (-3001.657) * [-3002.926] (-3014.927) (-3010.571) * [-2991.829] (-3006.284) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      813000 -- [-3009.650] (-3005.488) (-3000.541) * (-3023.351) [-2996.764] (-3011.678) * [-3003.278] (-3003.977) [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      814000 -- [-2999.983] (-3013.348) (-3011.410) * [-2996.796] (-3004.592) (-3000.660) * (-3009.326) [-3002.762] [...1 more local chains...] (...0 remote chains...) -- 0:00:20
      815000 -- [-3003.172] (-2998.857) (-3002.965) * (-3002.549) (-3002.657) [-3009.391] * (-2999.774) (-3012.975) [...1 more local chains...] (...0 remote chains...) -- 0:00:19

      Average standard deviation of split frequencies: 0.003348

      816000 -- [-3004.507] (-2995.958) (-3001.335) * (-3004.750) (-3003.561) [-2999.622] * (-2999.565) (-3001.908) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      817000 -- (-2997.639) (-3010.258) [-3002.879] * (-3000.049) [-3003.772] (-3008.971) * (-3000.127) (-3002.372) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      818000 -- [-3000.675] (-2996.705) (-3004.241) * (-3003.286) (-3004.152) [-2995.059] * (-3012.567) [-3001.248] [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      819000 -- [-2996.418] (-3002.420) (-3008.045) * [-2992.001] (-3001.452) (-3016.004) * [-2999.483] (-3012.408) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      820000 -- (-3026.365) (-3000.595) [-3001.769] * (-3007.071) (-2999.971) [-3006.369] * [-3005.675] (-2997.526) [...1 more local chains...] (...0 remote chains...) -- 0:00:19

      Average standard deviation of split frequencies: 0.003310

      821000 -- [-2998.492] (-3002.886) (-3000.571) * [-3000.413] (-3001.563) (-3008.019) * (-2997.807) (-3009.464) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      822000 -- [-2996.610] (-3005.834) (-2998.558) * (-3005.276) [-3005.069] (-3012.086) * (-3008.219) (-3007.559) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      823000 -- [-3005.956] (-3008.640) (-2998.274) * [-3007.116] (-2998.698) (-3004.871) * [-2996.179] (-3000.726) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      824000 -- [-2996.979] (-3008.339) (-3001.953) * [-2995.723] (-3008.886) (-3013.939) * (-3011.966) (-3006.456) [...1 more local chains...] (...0 remote chains...) -- 0:00:19
      825000 -- (-2995.410) [-2997.325] (-3008.523) * [-2997.616] (-2996.216) (-2998.469) * (-3006.160) (-3007.144) [...1 more local chains...] (...0 remote chains...) -- 0:00:18

      Average standard deviation of split frequencies: 0.003253

      826000 -- (-3000.643) [-3002.270] (-3005.652) * (-3005.148) (-3002.620) [-2996.159] * (-3015.110) [-3005.261] [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      827000 -- (-3008.482) (-3001.356) [-2994.199] * (-3006.383) (-3002.451) [-3004.253] * (-3016.639) (-3001.823) [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      828000 -- [-3001.335] (-2993.600) (-3007.193) * (-2996.792) (-3009.565) [-3002.506] * (-3013.636) [-2995.459] [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      829000 -- (-3000.213) (-3005.319) [-3001.866] * (-3009.671) (-2992.214) [-2997.659] * [-3003.386] (-2996.800) [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      830000 -- (-2997.532) [-2997.175] (-3014.708) * (-3007.093) (-2997.646) [-3002.344] * (-2992.270) [-2993.842] [...1 more local chains...] (...0 remote chains...) -- 0:00:18

      Average standard deviation of split frequencies: 0.003249

      831000 -- (-3004.941) (-3006.286) [-3002.832] * (-3014.991) [-3000.404] (-3002.078) * [-3003.459] (-3006.322) [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      832000 -- (-2998.791) (-3005.677) [-3001.087] * [-3003.885] (-2997.159) (-3009.407) * (-3000.991) (-3005.052) [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      833000 -- (-3010.650) [-3002.887] (-2997.227) * (-3017.197) (-3004.171) [-3003.368] * (-3000.854) (-3008.672) [...1 more local chains...] (...0 remote chains...) -- 0:00:18
      834000 -- (-3004.771) [-3002.001] (-2996.627) * (-3003.648) (-3000.428) [-2996.380] * [-3007.889] (-3003.627) [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      835000 -- (-3009.929) [-2992.743] (-3007.797) * [-2999.428] (-2998.248) (-3005.741) * (-3004.928) (-3003.986) [...1 more local chains...] (...0 remote chains...) -- 0:00:17

      Average standard deviation of split frequencies: 0.003353

      836000 -- (-3001.740) (-2997.630) [-3002.937] * [-3001.764] (-3004.313) (-3007.061) * (-3007.271) (-3005.841) [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      837000 -- (-2997.278) [-3000.004] (-3002.117) * (-3004.724) (-3002.588) [-3002.956] * [-3003.723] (-3002.053) [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      838000 -- (-2997.150) (-3004.805) [-3008.622] * (-3000.357) (-3012.215) [-2997.984] * (-3005.348) [-2997.099] [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      839000 -- (-2995.333) (-3004.236) [-2998.786] * (-3004.995) [-3007.147] (-3003.234) * (-3004.487) (-3005.549) [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      840000 -- (-2995.358) [-2997.319] (-3006.341) * (-2997.412) (-3019.028) [-3003.346] * (-3006.461) (-3005.199) [...1 more local chains...] (...0 remote chains...) -- 0:00:17

      Average standard deviation of split frequencies: 0.003385

      841000 -- (-3009.676) [-2998.672] (-3005.144) * (-3010.414) (-2995.187) [-3005.623] * (-3005.247) (-3009.930) [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      842000 -- [-2997.612] (-3006.543) (-3013.036) * [-2995.294] (-2999.044) (-3008.038) * (-3007.587) [-2995.465] [...1 more local chains...] (...0 remote chains...) -- 0:00:17
      843000 -- [-3001.300] (-3002.885) (-2998.853) * (-2999.538) [-2998.021] (-3001.802) * (-3021.574) (-3004.483) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      844000 -- [-2994.825] (-3008.998) (-2996.228) * (-2997.208) (-3003.636) [-3006.970] * (-2999.700) (-2997.015) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      845000 -- (-3004.835) (-3003.521) [-2997.331] * [-3001.520] (-3004.718) (-3001.571) * (-3003.190) [-3000.867] [...1 more local chains...] (...0 remote chains...) -- 0:00:16

      Average standard deviation of split frequencies: 0.003412

      846000 -- [-3006.429] (-2998.985) (-2997.884) * (-3000.065) [-3002.965] (-3002.072) * (-3003.401) [-3003.244] [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      847000 -- (-3007.482) [-3010.224] (-2999.306) * [-3004.199] (-3005.631) (-3009.773) * [-2998.154] (-3006.010) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      848000 -- (-3010.881) (-3003.498) [-3002.770] * [-3005.813] (-3001.369) (-3005.604) * (-3012.781) (-3006.242) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      849000 -- (-3001.746) (-2999.651) [-3001.630] * (-3012.019) [-3005.374] (-3005.617) * (-3009.477) (-3001.950) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      850000 -- (-3003.457) (-3004.232) [-2999.644] * (-3003.641) [-2999.057] (-3000.084) * (-2997.901) (-3005.974) [...1 more local chains...] (...0 remote chains...) -- 0:00:16

      Average standard deviation of split frequencies: 0.003406

      851000 -- [-3004.322] (-2998.441) (-3006.069) * (-3011.112) [-2997.065] (-3015.317) * [-3003.030] (-3001.512) [...1 more local chains...] (...0 remote chains...) -- 0:00:16
      852000 -- (-3014.190) (-3007.886) [-2991.935] * (-2999.246) (-3010.025) [-2996.419] * (-3001.190) (-2997.760) [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      853000 -- (-3009.599) [-2991.532] (-3010.101) * (-3006.317) [-3005.404] (-2996.687) * (-2994.730) [-3005.501] [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      854000 -- (-3003.808) [-3001.616] (-3009.722) * [-3001.294] (-2996.923) (-3004.993) * (-2999.910) (-3000.489) [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      855000 -- [-2999.258] (-3004.394) (-3006.133) * (-3001.511) (-3002.059) [-2996.175] * (-3006.733) [-3001.800] [...1 more local chains...] (...0 remote chains...) -- 0:00:15

      Average standard deviation of split frequencies: 0.003340

      856000 -- (-3011.976) (-3000.623) [-3000.807] * (-3002.732) (-3000.963) [-3002.908] * (-3005.320) (-3006.526) [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      857000 -- [-2999.816] (-3003.848) (-3006.330) * [-3005.478] (-2999.192) (-3000.427) * (-3003.715) (-3003.436) [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      858000 -- (-3008.598) (-2998.987) [-3002.016] * (-3007.810) [-2997.991] (-3000.913) * (-2992.887) [-3000.620] [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      859000 -- (-3001.452) (-3001.522) [-3004.144] * (-3000.433) [-3002.359] (-2998.420) * [-3003.883] (-2997.483) [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      860000 -- (-3019.248) (-3001.709) [-2992.394] * [-3002.461] (-2997.921) (-3010.032) * (-3005.850) [-3001.146] [...1 more local chains...] (...0 remote chains...) -- 0:00:15

      Average standard deviation of split frequencies: 0.003328

      861000 -- [-2994.751] (-3005.264) (-3011.353) * (-2999.328) [-3002.444] (-2998.666) * (-3010.250) [-2996.161] [...1 more local chains...] (...0 remote chains...) -- 0:00:15
      862000 -- (-2998.854) [-2999.109] (-3004.427) * (-3011.403) [-2997.870] (-3003.507) * [-3004.953] (-3000.293) [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      863000 -- (-3004.504) [-3002.433] (-3005.469) * [-2998.404] (-3002.879) (-3007.242) * (-3010.949) [-3007.787] [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      864000 -- (-2997.938) (-2998.178) [-3001.347] * [-3008.555] (-3003.060) (-3009.076) * (-3002.608) (-3008.458) [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      865000 -- (-3000.507) [-2998.380] (-3001.449) * (-3010.745) (-3006.143) [-2992.776] * (-3000.227) (-3003.907) [...1 more local chains...] (...0 remote chains...) -- 0:00:14

      Average standard deviation of split frequencies: 0.003339

      866000 -- (-3002.184) [-2996.857] (-2996.220) * (-3001.727) (-3005.960) [-3011.598] * (-3004.509) (-2998.237) [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      867000 -- [-2999.697] (-2999.429) (-3000.718) * (-3016.808) [-3000.600] (-3001.606) * [-3003.847] (-3003.966) [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      868000 -- [-2999.113] (-3005.849) (-2993.960) * [-3006.047] (-3007.830) (-2999.724) * [-3001.343] (-3006.971) [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      869000 -- [-3001.219] (-3008.237) (-3014.522) * [-2997.311] (-3012.250) (-3000.485) * (-3002.829) [-3003.949] [...1 more local chains...] (...0 remote chains...) -- 0:00:14
      870000 -- (-3004.476) [-2998.627] (-3000.048) * (-2991.520) (-2997.268) [-3005.410] * (-3005.354) [-2990.480] [...1 more local chains...] (...0 remote chains...) -- 0:00:14

      Average standard deviation of split frequencies: 0.003397

      871000 -- (-3006.404) (-3018.023) [-3003.007] * (-3003.998) (-3004.728) [-3003.465] * [-3009.239] (-3024.296) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      872000 -- (-3008.528) [-2999.990] (-3006.421) * [-2996.475] (-3006.241) (-2997.470) * (-2996.632) (-3005.999) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      873000 -- (-3000.077) (-2997.349) [-2996.946] * [-2997.651] (-2999.767) (-3025.668) * (-3005.081) [-3012.757] [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      874000 -- [-3007.930] (-3005.743) (-3001.566) * [-2999.752] (-3001.832) (-3001.511) * (-3001.878) (-3001.160) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      875000 -- (-2995.028) [-3000.503] (-3002.803) * (-3008.536) [-3005.815] (-3002.564) * [-3002.869] (-2998.953) [...1 more local chains...] (...0 remote chains...) -- 0:00:13

      Average standard deviation of split frequencies: 0.003300

      876000 -- (-2995.330) [-3005.476] (-3012.756) * [-2998.443] (-3001.386) (-2999.877) * (-3005.268) (-2998.655) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      877000 -- (-2999.608) [-2996.391] (-3010.901) * [-3000.022] (-2996.913) (-3002.897) * (-3009.418) (-3008.946) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      878000 -- [-3002.615] (-3004.280) (-2999.251) * (-2998.861) (-3011.515) [-3008.365] * [-3000.279] (-2998.305) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      879000 -- (-3008.262) (-3004.000) [-3002.883] * (-3005.216) [-2997.265] (-2996.809) * [-2996.732] (-3004.397) [...1 more local chains...] (...0 remote chains...) -- 0:00:13
      880000 -- (-3004.469) [-3001.522] (-3002.466) * (-2998.259) (-2997.381) [-3000.623] * (-3007.137) (-3000.804) [...1 more local chains...] (...0 remote chains...) -- 0:00:12

      Average standard deviation of split frequencies: 0.003383

      881000 -- [-3003.774] (-2998.532) (-2998.675) * [-3001.856] (-3003.996) (-3007.964) * [-2998.018] (-3008.238) [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      882000 -- [-3000.421] (-3009.272) (-2997.301) * [-3001.198] (-3004.249) (-2995.745) * [-3001.022] (-2997.021) [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      883000 -- [-3003.233] (-3008.688) (-3008.412) * (-3007.937) [-3015.756] (-3000.245) * (-3006.920) (-2999.296) [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      884000 -- [-3002.103] (-2991.382) (-3015.421) * [-3001.462] (-3020.687) (-2996.598) * (-3004.922) [-3002.806] [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      885000 -- (-3009.840) (-3005.501) [-3009.763] * (-2999.164) (-3001.740) [-2999.568] * [-3005.900] (-3005.189) [...1 more local chains...] (...0 remote chains...) -- 0:00:12

      Average standard deviation of split frequencies: 0.003308

      886000 -- (-3005.371) (-3000.078) [-2996.915] * (-3004.600) (-3001.014) [-2997.659] * (-2997.461) [-2995.512] [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      887000 -- (-3005.281) (-3003.246) [-3002.556] * [-3006.374] (-3001.045) (-3002.860) * [-3005.725] (-3008.500) [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      888000 -- (-3008.518) [-3001.876] (-3000.208) * [-3003.767] (-3003.429) (-3002.194) * [-3007.081] (-3003.437) [...1 more local chains...] (...0 remote chains...) -- 0:00:12
      889000 -- (-3001.086) [-3002.109] (-3005.520) * (-3008.091) [-2996.587] (-2998.628) * [-2994.949] (-3001.319) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      890000 -- [-2997.445] (-2999.007) (-3015.034) * (-3012.062) (-3006.166) [-3004.334] * [-3000.005] (-3009.478) [...1 more local chains...] (...0 remote chains...) -- 0:00:11

      Average standard deviation of split frequencies: 0.003343

      891000 -- (-3008.085) [-3004.467] (-3007.256) * (-3008.969) (-3003.849) [-3000.584] * [-2996.231] (-3004.792) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      892000 -- (-2996.580) [-2998.809] (-3003.191) * (-3000.852) [-2994.870] (-3002.187) * (-3003.155) (-3007.323) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      893000 -- [-2991.632] (-2997.592) (-3006.915) * [-3000.790] (-3007.153) (-3020.561) * [-2992.319] (-3010.408) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      894000 -- (-3007.691) (-3001.956) [-2997.533] * (-3005.678) (-2996.078) [-3009.842] * [-2993.445] (-3001.495) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      895000 -- [-2994.500] (-3001.646) (-3004.673) * (-3016.798) (-3003.493) [-3005.320] * [-3010.619] (-3003.083) [...1 more local chains...] (...0 remote chains...) -- 0:00:11

      Average standard deviation of split frequencies: 0.003420

      896000 -- (-3007.719) [-3002.549] (-3008.087) * (-3007.189) (-3001.170) [-3007.371] * (-3005.956) [-3007.882] [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      897000 -- (-3007.242) [-3004.628] (-3005.011) * [-2995.129] (-2997.567) (-2996.172) * (-2998.965) (-2999.144) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      898000 -- (-3001.243) [-3018.295] (-3003.239) * [-3000.778] (-3007.674) (-2993.029) * (-3010.436) (-3019.643) [...1 more local chains...] (...0 remote chains...) -- 0:00:11
      899000 -- (-3006.344) [-3001.077] (-2999.119) * (-3013.161) [-2999.569] (-3008.443) * [-3006.617] (-3005.459) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      900000 -- (-3001.068) [-3005.104] (-3004.263) * [-2998.850] (-3003.185) (-3000.484) * (-3002.789) [-3011.077] [...1 more local chains...] (...0 remote chains...) -- 0:00:10

      Average standard deviation of split frequencies: 0.003342

      901000 -- [-2993.841] (-3015.776) (-2994.168) * (-3005.152) (-3004.810) [-2996.502] * (-2998.763) (-3001.248) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      902000 -- (-3010.605) (-3000.744) [-2996.004] * (-3004.352) (-3003.851) [-2998.792] * (-3001.009) [-3003.363] [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      903000 -- (-2995.962) (-3001.171) [-2992.127] * [-3010.212] (-3004.392) (-3011.779) * (-2996.716) (-3003.099) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      904000 -- (-3000.651) (-3006.280) [-2995.131] * (-2999.093) [-3002.546] (-3004.866) * [-2995.906] (-3010.725) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      905000 -- (-3005.122) [-2997.704] (-3004.674) * (-3002.997) (-3005.209) [-3001.394] * (-2998.195) [-2998.614] [...1 more local chains...] (...0 remote chains...) -- 0:00:10

      Average standard deviation of split frequencies: 0.003330

      906000 -- (-2997.363) (-3009.276) [-2996.303] * (-3003.659) (-2997.917) [-3001.363] * [-3002.749] (-3005.398) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      907000 -- (-3012.950) [-3006.808] (-3007.843) * (-2999.147) (-3011.281) [-3004.637] * (-3007.070) (-3009.445) [...1 more local chains...] (...0 remote chains...) -- 0:00:10
      908000 -- [-2992.878] (-3005.291) (-3004.702) * [-3008.837] (-2996.725) (-3002.236) * [-3004.739] (-2999.809) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      909000 -- (-3001.765) (-3010.308) [-3003.566] * (-2997.497) [-3004.461] (-3016.352) * (-2996.230) [-3006.978] [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      910000 -- [-3004.074] (-3008.387) (-3002.731) * (-3007.528) (-2997.439) [-3002.383] * (-3015.067) (-3001.878) [...1 more local chains...] (...0 remote chains...) -- 0:00:09

      Average standard deviation of split frequencies: 0.003366

      911000 -- (-3007.341) [-3001.595] (-3011.760) * [-2999.960] (-3005.840) (-2999.154) * [-3000.404] (-3009.356) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      912000 -- (-3000.056) (-3003.614) [-3000.703] * (-3011.376) (-3001.463) [-2998.107] * [-3000.043] (-3005.809) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      913000 -- [-3003.767] (-2999.566) (-3002.039) * (-3001.955) [-2998.191] (-3002.228) * [-3007.021] (-3013.367) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      914000 -- [-3000.214] (-2999.376) (-3003.197) * (-2997.239) [-3008.356] (-2997.726) * [-3001.227] (-3008.114) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      915000 -- (-3002.056) (-3001.510) [-2999.218] * (-3009.229) (-3009.506) [-3000.647] * [-3002.458] (-2995.668) [...1 more local chains...] (...0 remote chains...) -- 0:00:09

      Average standard deviation of split frequencies: 0.003304

      916000 -- (-3009.915) [-3004.356] (-3017.606) * (-3005.427) (-2999.108) [-3001.189] * [-3003.785] (-3004.592) [...1 more local chains...] (...0 remote chains...) -- 0:00:09
      917000 -- (-3005.787) [-3001.565] (-3015.470) * [-3004.178] (-3018.366) (-3004.952) * [-3009.796] (-3000.339) [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      918000 -- (-2992.501) (-3003.333) [-3005.823] * (-3001.504) (-3006.331) [-3006.854] * (-2999.845) [-2993.635] [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      919000 -- [-3008.471] (-3006.848) (-3000.354) * (-3003.097) [-3002.256] (-3003.530) * [-3004.042] (-3006.604) [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      920000 -- [-2996.917] (-3002.564) (-3016.908) * (-3004.536) [-3002.132] (-2996.707) * (-3009.619) [-3000.798] [...1 more local chains...] (...0 remote chains...) -- 0:00:08

      Average standard deviation of split frequencies: 0.003302

      921000 -- (-3003.825) (-3000.481) [-3002.246] * [-3004.796] (-3007.492) (-2989.924) * (-3008.848) (-3008.528) [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      922000 -- (-3006.838) [-3004.841] (-3015.132) * (-3003.414) (-3000.862) [-2993.935] * (-2996.731) [-2997.953] [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      923000 -- (-2997.776) [-2998.095] (-2997.319) * (-3002.831) [-3005.750] (-3003.799) * (-3003.981) (-3002.226) [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      924000 -- [-3004.048] (-3009.964) (-2996.959) * [-3000.896] (-3011.004) (-3005.371) * [-2999.454] (-3007.844) [...1 more local chains...] (...0 remote chains...) -- 0:00:08
      925000 -- [-3004.882] (-3001.251) (-2995.786) * (-3002.422) [-3001.114] (-3011.094) * [-3001.451] (-3009.708) [...1 more local chains...] (...0 remote chains...) -- 0:00:08

      Average standard deviation of split frequencies: 0.003317

      926000 -- (-3013.953) [-2998.536] (-3002.090) * [-2994.607] (-3003.852) (-3008.129) * (-2996.936) [-3000.776] [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      927000 -- [-3003.680] (-3007.767) (-3001.996) * (-3000.517) (-3003.100) [-2999.336] * (-3006.088) [-3002.952] [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      928000 -- (-2997.316) [-2995.417] (-3006.957) * (-3001.656) [-3012.895] (-3001.657) * [-3002.345] (-3003.815) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      929000 -- (-2997.628) (-3002.205) [-2996.230] * [-3010.080] (-3008.492) (-3008.209) * (-3008.725) (-3006.617) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      930000 -- [-2998.581] (-3015.766) (-3007.451) * (-3005.422) [-2999.473] (-3007.704) * [-2996.353] (-3000.638) [...1 more local chains...] (...0 remote chains...) -- 0:00:07

      Average standard deviation of split frequencies: 0.003158

      931000 -- (-3008.071) [-2998.597] (-3001.229) * (-2997.791) (-3000.843) [-3002.037] * [-2996.632] (-2999.588) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      932000 -- [-2998.374] (-3000.938) (-3005.436) * (-2998.967) [-3002.686] (-3000.569) * (-3003.383) (-3015.743) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      933000 -- (-3003.111) [-2998.300] (-3013.934) * [-3003.417] (-3004.467) (-2998.986) * [-3003.606] (-3020.430) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      934000 -- (-3007.187) [-3001.991] (-3015.882) * (-3014.053) (-3009.362) [-2998.283] * [-3003.966] (-3006.598) [...1 more local chains...] (...0 remote chains...) -- 0:00:07
      935000 -- (-2998.151) (-3007.955) [-3002.303] * [-3002.280] (-3011.113) (-3002.986) * [-2994.832] (-3000.815) [...1 more local chains...] (...0 remote chains...) -- 0:00:07

      Average standard deviation of split frequencies: 0.003071

      936000 -- (-3003.379) [-2997.195] (-3002.551) * (-3010.233) [-3001.981] (-3012.954) * [-3002.354] (-3008.291) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      937000 -- (-3006.758) (-3001.938) [-2999.553] * [-3003.421] (-3000.507) (-3011.113) * [-2998.485] (-3014.291) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      938000 -- (-2995.653) (-3004.034) [-3008.482] * [-3003.279] (-3000.992) (-3005.240) * [-2996.806] (-3000.201) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      939000 -- [-2995.121] (-3007.617) (-3004.328) * [-3003.660] (-2999.566) (-3008.777) * [-2996.687] (-3000.852) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      940000 -- (-3000.915) (-3006.705) [-2999.774] * (-3001.938) (-2992.638) [-3003.441] * (-3005.145) [-2998.402] [...1 more local chains...] (...0 remote chains...) -- 0:00:06

      Average standard deviation of split frequencies: 0.003087

      941000 -- (-3013.314) [-2994.398] (-3010.047) * (-3010.781) (-3001.535) [-3002.933] * [-3001.952] (-3015.220) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      942000 -- [-3000.163] (-2996.189) (-3011.490) * (-2996.381) (-3001.285) [-3000.972] * (-3000.082) (-3002.224) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      943000 -- (-3017.133) [-3002.674] (-3003.484) * (-3000.556) [-2995.766] (-3007.796) * [-3001.779] (-2999.784) [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      944000 -- (-3002.783) (-3000.476) [-2997.531] * (-3008.015) (-3003.823) [-3002.178] * (-2995.281) [-3007.852] [...1 more local chains...] (...0 remote chains...) -- 0:00:06
      945000 -- (-3010.825) (-2997.121) [-3005.819] * [-3002.347] (-3017.705) (-3002.362) * (-3000.891) (-3011.273) [...1 more local chains...] (...0 remote chains...) -- 0:00:05

      Average standard deviation of split frequencies: 0.003046

      946000 -- [-3002.199] (-3012.733) (-3009.376) * [-3000.504] (-3001.132) (-3018.395) * [-3004.194] (-3013.667) [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      947000 -- [-3002.108] (-2996.783) (-3001.904) * [-3018.402] (-3004.619) (-2997.302) * (-3003.083) [-2998.059] [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      948000 -- (-3014.443) (-3004.944) [-2999.832] * (-2999.534) [-3003.623] (-3004.916) * (-3010.670) (-3001.110) [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      949000 -- [-3001.253] (-3003.609) (-3002.944) * (-3010.995) (-3006.065) [-2998.193] * (-3009.189) [-2997.068] [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      950000 -- [-2993.056] (-3002.143) (-3005.031) * (-3003.375) (-2999.898) [-2997.443] * [-2997.954] (-3006.042) [...1 more local chains...] (...0 remote chains...) -- 0:00:05

      Average standard deviation of split frequencies: 0.003117

      951000 -- (-2996.543) (-3003.504) [-2994.189] * (-3013.207) [-3003.995] (-3008.572) * (-3002.488) (-3009.936) [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      952000 -- [-2995.197] (-3001.140) (-3004.514) * (-3002.157) (-3003.933) [-2998.812] * [-3000.087] (-3011.453) [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      953000 -- [-3001.961] (-3006.230) (-3011.757) * (-2996.127) (-3002.757) [-3005.668] * (-3003.742) [-2994.206] [...1 more local chains...] (...0 remote chains...) -- 0:00:05
      954000 -- (-3001.502) [-2997.951] (-3010.886) * (-2995.814) [-3006.257] (-3002.400) * (-3005.735) [-3002.265] [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      955000 -- [-2996.943] (-3010.154) (-3010.020) * (-3005.452) (-3013.422) [-2997.616] * [-2998.258] (-3009.622) [...1 more local chains...] (...0 remote chains...) -- 0:00:04

      Average standard deviation of split frequencies: 0.003188

      956000 -- (-3006.201) (-3010.912) [-2991.140] * [-3000.028] (-3004.816) (-3003.507) * (-2997.922) (-3002.320) [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      957000 -- [-3006.054] (-2996.295) (-2996.966) * [-2995.964] (-2995.440) (-3000.150) * [-3006.877] (-3015.663) [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      958000 -- (-3001.398) [-3005.104] (-3004.895) * [-3009.976] (-3003.803) (-2997.802) * (-3003.741) (-3002.449) [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      959000 -- [-3000.693] (-2999.430) (-2999.119) * [-3010.403] (-3009.068) (-3005.046) * (-3004.649) [-3001.166] [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      960000 -- [-3005.634] (-3000.916) (-3011.551) * (-2997.724) [-2996.109] (-3005.707) * (-2999.355) [-3002.919] [...1 more local chains...] (...0 remote chains...) -- 0:00:04

      Average standard deviation of split frequencies: 0.003162

      961000 -- (-3003.464) [-3010.549] (-3001.713) * (-3003.045) [-3001.428] (-3015.781) * (-3003.348) (-2999.860) [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      962000 -- (-2995.155) (-3008.809) [-3004.062] * [-3001.162] (-3002.273) (-3005.617) * (-3010.865) [-3004.076] [...1 more local chains...] (...0 remote chains...) -- 0:00:04
      963000 -- [-3005.374] (-2995.344) (-3007.605) * (-2998.883) [-2998.030] (-3008.140) * [-3001.101] (-2998.958) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      964000 -- (-3005.757) [-2993.799] (-3007.278) * [-2995.740] (-3012.925) (-3000.810) * (-2999.023) (-3002.041) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      965000 -- (-3001.812) (-3002.370) [-3017.473] * (-3008.055) (-2999.450) [-2996.878] * (-3008.182) (-3003.021) [...1 more local chains...] (...0 remote chains...) -- 0:00:03

      Average standard deviation of split frequencies: 0.003132

      966000 -- (-3017.048) (-3002.907) [-3004.911] * (-3016.196) [-2999.473] (-3007.803) * (-2996.652) (-2999.565) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      967000 -- (-2996.572) (-3009.651) [-3004.959] * (-3014.803) [-2997.284] (-2996.572) * [-2997.801] (-2999.370) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      968000 -- [-2997.579] (-3007.139) (-3001.489) * [-3000.031] (-3000.920) (-2996.488) * [-2996.825] (-3001.556) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      969000 -- [-3002.025] (-3007.018) (-3014.504) * (-3003.768) (-3004.462) [-3022.634] * (-2997.258) [-2994.748] [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      970000 -- (-3001.174) [-3005.265] (-3000.137) * [-3003.686] (-2995.255) (-2997.335) * (-3008.623) (-3000.845) [...1 more local chains...] (...0 remote chains...) -- 0:00:03

      Average standard deviation of split frequencies: 0.003127

      971000 -- (-3005.319) (-3006.011) [-2998.954] * (-3007.583) [-3003.440] (-3002.171) * (-3005.946) (-2996.278) [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      972000 -- (-3001.162) [-2994.815] (-3005.872) * (-3003.756) (-2997.439) [-2996.632] * (-3012.987) [-3001.085] [...1 more local chains...] (...0 remote chains...) -- 0:00:03
      973000 -- (-3012.643) [-3008.867] (-3007.892) * [-3002.388] (-3001.066) (-2997.148) * [-2999.355] (-2998.781) [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      974000 -- [-3001.746] (-2999.849) (-2996.692) * [-2998.152] (-3003.731) (-3004.039) * [-2995.609] (-3005.775) [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      975000 -- [-2994.003] (-3005.525) (-3010.292) * [-3009.616] (-3009.754) (-3000.352) * [-2996.914] (-2993.158) [...1 more local chains...] (...0 remote chains...) -- 0:00:02

      Average standard deviation of split frequencies: 0.003187

      976000 -- (-3006.281) (-3001.640) [-2995.798] * (-3000.814) [-3001.399] (-3000.460) * (-3009.562) [-3003.945] [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      977000 -- [-2994.160] (-3000.868) (-3004.166) * (-3005.854) (-3013.606) [-3016.270] * (-3006.192) [-3002.376] [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      978000 -- (-3006.421) [-2998.591] (-2998.078) * (-3000.786) [-3014.228] (-3004.853) * (-3005.555) [-2998.488] [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      979000 -- (-3002.886) (-2998.553) [-3010.974] * [-2994.998] (-2993.494) (-3005.410) * [-2998.348] (-3005.917) [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      980000 -- (-3012.770) (-3004.347) [-3016.298] * (-3002.781) [-3001.557] (-3015.370) * (-3005.827) (-3003.641) [...1 more local chains...] (...0 remote chains...) -- 0:00:02

      Average standard deviation of split frequencies: 0.003197

      981000 -- (-3007.528) [-2998.383] (-3007.609) * (-3008.448) [-2996.130] (-3001.740) * (-3002.213) (-3004.382) [...1 more local chains...] (...0 remote chains...) -- 0:00:02
      982000 -- (-3009.546) (-2997.036) [-3000.864] * (-3025.182) (-3000.588) [-2997.169] * [-2995.065] (-3005.325) [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      983000 -- (-3000.675) [-2999.783] (-2998.235) * (-3016.983) [-3002.595] (-3014.947) * (-3003.524) [-3004.960] [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      984000 -- [-2999.789] (-3004.288) (-3006.886) * [-2996.572] (-3001.006) (-3008.447) * (-3006.918) [-2999.680] [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      985000 -- (-3004.781) (-3006.127) [-3004.539] * [-2995.230] (-3001.036) (-2996.333) * (-3001.976) (-3003.990) [...1 more local chains...] (...0 remote chains...) -- 0:00:01

      Average standard deviation of split frequencies: 0.003232

      986000 -- (-3001.046) [-3001.075] (-3004.397) * (-3001.821) [-2998.425] (-3004.869) * (-3007.527) (-3004.238) [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      987000 -- [-3004.870] (-3008.690) (-2999.764) * [-3002.175] (-3002.012) (-3006.729) * (-2998.981) (-2999.180) [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      988000 -- (-3003.910) (-3000.662) [-3003.437] * (-3002.411) [-2994.209] (-3000.997) * (-3006.423) [-3003.911] [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      989000 -- [-3021.624] (-3002.562) (-3013.777) * [-3003.006] (-2999.432) (-3004.585) * [-3004.613] (-3002.014) [...1 more local chains...] (...0 remote chains...) -- 0:00:01
      990000 -- (-3013.316) [-2995.438] (-3010.544) * (-3007.007) [-2996.305] (-3004.201) * (-3003.531) [-3002.184] [...1 more local chains...] (...0 remote chains...) -- 0:00:01

      Average standard deviation of split frequencies: 0.003220

      991000 -- (-2999.381) (-3003.064) [-2998.218] * [-2996.233] (-3013.716) (-3007.889) * (-3010.243) [-2999.820] [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      992000 -- [-3001.656] (-3014.868) (-2997.948) * (-2998.511) (-2999.337) [-2999.927] * [-2999.327] (-2995.729) [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      993000 -- (-3005.862) (-3005.679) [-3003.583] * (-3008.093) (-3016.893) [-3001.856] * [-2999.225] (-3002.207) [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      994000 -- [-2999.814] (-2998.351) (-3003.483) * (-3002.650) (-3005.907) [-3005.515] * (-2997.583) [-2999.011] [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      995000 -- [-3004.103] (-3002.026) (-3001.343) * [-3003.479] (-3003.528) (-2994.054) * (-3002.464) (-2998.899) [...1 more local chains...] (...0 remote chains...) -- 0:00:00

      Average standard deviation of split frequencies: 0.003257

      996000 -- (-3007.261) (-3005.973) [-3005.194] * (-3009.331) (-3005.770) [-2996.140] * [-2998.607] (-3004.730) [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      997000 -- (-3001.163) (-2999.166) [-2999.121] * (-3002.955) [-3004.437] (-3008.127) * (-3006.478) [-2998.865] [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      998000 -- (-3005.784) (-3006.007) [-3000.553] * (-3013.003) (-3005.781) [-2998.137] * (-2999.139) (-3006.488) [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      999000 -- (-3003.211) [-2998.294] (-2998.700) * (-3004.551) (-2997.734) [-3010.931] * (-3002.820) [-3004.579] [...1 more local chains...] (...0 remote chains...) -- 0:00:00
      1000000 -- (-3007.385) (-3023.533) [-3002.529] * [-3011.110] (-3011.093) (-3016.529) * [-3002.462] (-3017.585) [...1 more local chains...] (...0 remote chains...) -- 0:00:00

      Average standard deviation of split frequencies: 0.003311

      Analysis completed in 1 mins 48 seconds
      Analysis used 107.09 seconds of CPU time on processor 0
      Likelihood of best state for "cold" chain of run 1 was -2986.26
      Likelihood of best state for "cold" chain of run 2 was -2986.67
      Likelihood of best state for "cold" chain of run 3 was -2987.88

      Acceptance rates for the moves in the "cold" chain of run 1:
         With prob.   (last 100)   chain accepted proposals by move
            32.2 %     ( 24 %)     Dirichlet(Tratio)
            22.2 %     ( 20 %)     Dirichlet(Pi)
            26.2 %     ( 28 %)     Slider(Pi)
            65.7 %     ( 33 %)     Multiplier(Alpha)
            37.9 %     ( 29 %)     ExtSPR(Tau,V)
            23.1 %     ( 25 %)     ExtTBR(Tau,V)
            40.1 %     ( 43 %)     NNI(Tau,V)
            34.0 %     ( 35 %)     ParsSPR(Tau,V)
            27.1 %     ( 30 %)     Multiplier(V)
            58.2 %     ( 53 %)     Nodeslider(V)
            25.6 %     ( 26 %)     TLMultiplier(V)

      Acceptance rates for the moves in the "cold" chain of run 2:
         With prob.   (last 100)   chain accepted proposals by move
            31.5 %     ( 27 %)     Dirichlet(Tratio)
            22.3 %     ( 22 %)     Dirichlet(Pi)
            26.3 %     ( 28 %)     Slider(Pi)
            66.3 %     ( 41 %)     Multiplier(Alpha)
            37.9 %     ( 29 %)     ExtSPR(Tau,V)
            23.4 %     ( 21 %)     ExtTBR(Tau,V)
            40.5 %     ( 42 %)     NNI(Tau,V)
            34.0 %     ( 37 %)     ParsSPR(Tau,V)
            27.1 %     ( 30 %)     Multiplier(V)
            57.9 %     ( 56 %)     Nodeslider(V)
            25.3 %     ( 20 %)     TLMultiplier(V)

      Acceptance rates for the moves in the "cold" chain of run 3:
         With prob.   (last 100)   chain accepted proposals by move
            32.2 %     ( 31 %)     Dirichlet(Tratio)
            22.6 %     ( 20 %)     Dirichlet(Pi)
            26.0 %     ( 21 %)     Slider(Pi)
            65.5 %     ( 34 %)     Multiplier(Alpha)
            38.0 %     ( 39 %)     ExtSPR(Tau,V)
            23.2 %     ( 19 %)     ExtTBR(Tau,V)
            40.7 %     ( 39 %)     NNI(Tau,V)
            34.0 %     ( 30 %)     ParsSPR(Tau,V)
            26.9 %     ( 25 %)     Multiplier(V)
            58.2 %     ( 53 %)     Nodeslider(V)
            25.2 %     ( 23 %)     TLMultiplier(V)

      Chain swap information for run 1:

                   1       2       3
           --------------------------
         1 |            0.75    0.53
         2 |  333050            0.76
         3 |  333781  333169

      Chain swap information for run 2:

                   1       2       3
           --------------------------
         1 |            0.74    0.53
         2 |  332838            0.77
         3 |  333031  334131

      Chain swap information for run 3:

                   1       2       3
           --------------------------
         1 |            0.74    0.53
         2 |  332744            0.77
         3 |  333727  333529

      Upper diagonal: Proportion of successful state exchanges between chains
      Lower diagonal: Number of attempted state exchanges between chains

      Chain information:

        ID -- Heat
       -----------
         1 -- 1.00  (cold chain)
         2 -- 0.91
         3 -- 0.83

      Heat = 1 / (1 + T * (ID - 1))
         (where T = 0.10 is the temperature and ID is the chain number)

      Summarizing trees in files "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.run1.t", "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.run2.t",...,"/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.run3.t"
      Using relative burnin ('relburnin=yes'), discarding the first 25 % of sampled trees
      Writing statistics to files /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.<parts|tstat|vstat|trprobs|con>
      Examining first file ...
      Found one tree block in file "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.run1.t" with 100001 trees in last block
      Expecting the same number of trees in the last tree block of all files

      Tree reading status:

      0      10      20      30      40      50      60      70      80      90     100
      v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
      *********************************************************************************

      Read a total of 300003 trees in 3 files (sampling 225003 of them)
         (Each file contained 100001 trees of which 75001 were sampled)

      General explanation:

      In an unrooted tree, a taxon bipartition (split) is specified by removing a
      branch, thereby dividing the species into those to the left and those to the
      right of the branch. Here, taxa to one side of the removed branch are denoted
      '.' and those to the other side are denoted '*'. Specifically, the '.' symbol
      is used for the taxa on the same side as the outgroup.

      In a rooted or clock tree, the tree is rooted using the model and not by
      reference to an outgroup. Each bipartition therefore corresponds to a clade,
      that is, a group that includes all the descendants of a particular branch in
      the tree.  Taxa that are included in each clade are denoted using '*', and
      taxa that are not included are denoted using the '.' symbol.

      The output first includes a key to all the bipartitions with frequency larger
      or equual to (Minpartfreq) in at least one run. Minpartfreq is a parameter to
      sumt command and currently it is set to 0.10.  This is followed by a table
      with statistics for the informative bipartitions (those including at least
      two taxa), sorted from highest to lowest probability. For each bipartition,
      the table gives the number of times the partition or split was observed in all
      runs (#obs) and the posterior probability of the bipartition (Probab.), which
      is the same as the split frequency. If several runs are summarized, this is
      followed by the minimum split frequency (Min(s)), the maximum frequency
      (Max(s)), and the standard deviation of frequencies (Stddev(s)) across runs.
      The latter value should approach 0 for all bipartitions as MCMC runs converge.

      This is followed by a table summarizing branch lengths, node heights (if a
      clock model was used) and relaxed clock parameters (if a relaxed clock model
      was used). The mean, variance, and 95 % credible interval are given for each
      of these parameters. If several runs are summarized, the potential scale
      reduction factor (PSRF) is also given; it should approach 1 as runs converge.
      Node heights will take calibration points into account, if such points were
      used in the analysis.

      Note that Stddev may be unreliable if the partition is not present in all
      runs (the last column indicates the number of runs that sampled the partition
      if more than one run is summarized). The PSRF is not calculated at all if
      the partition is not present in all runs.The PSRF is also sensitive to small
      sample sizes and it should only be considered a rough guide to convergence
      since some of the assumptions allowing one to interpret it as a true potential
      scale reduction factor are violated in MrBayes.

      List of taxa in bipartitions:

         1 -- FJ895118.1_Aphanizomenon_aphanizomenoides
         2 -- AF448070.1_Synechococcus_sp._PS845
         3 -- FJ895126.1_Aphanizomenon_gracile
         4 -- FJ895128.1_Aphanizomenon_gracile
         5 -- FJ895125.1_Aphanizomenon_gracile
         6 -- FJ895119.1_Aphanizomenon_aphanizomenoides
         7 -- EF685373.1_Aphanizomenon_issatschenkoi
         8 -- FJ895121.1_Aphanizomenon_aphanizomenoides
         9 -- FJ895120.1_Aphanizomenon_aphanizomenoides
        10 -- FJ895124.1_Aphanizomenon_gracile
        11 -- FJ895122.1_Aphanizomenon_aphanizomenoides
        12 -- FJ895123.1_Aphanizomenon_aphanizomenoides
        13 -- FJ895127.1_Aphanizomenon_gracile

      Key to taxon bipartitions (saved to file "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.parts"):

      ID -- Partition
      -------------------
       1 -- *............
       2 -- *.***********
       3 -- ..*..........
       4 -- ...*.........
       5 -- ....*........
       6 -- .....*.......
       7 -- ......*......
       8 -- .......*.....
       9 -- ........*....
      10 -- .........*...
      11 -- ..........*..
      12 -- ...........*.
      13 -- ............*
      14 -- ..***.*..*..*
      15 -- ..***....*..*
      16 -- *....*.**.**.
      17 -- .......**.*..
      18 -- .......**....
      19 -- *......**.*..
      20 -- ..*.*....*..*
      21 -- .....*.....*.
      22 -- *......**.**.
      23 -- ..***....*...
      24 -- *....*.....*.
      25 -- ..*.*....*...
      26 -- ....*....*...
      27 -- ..*......*...
      28 -- ..*.*........
      29 -- .......*..*..
      30 -- ........*.*..
      31 -- *..........*.
      32 -- ..*.........*
      33 -- ....*.......*
      34 -- .........*..*
      35 -- ..*.*.......*
      36 -- ....*....*..*
      37 -- ..*......*..*
      -------------------

      Summary statistics for informative taxon bipartitions
         (saved to file "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.tstat"):

      ID   #obs      Probab.     Sd(s)+      Min(s)      Max(s)   Nruns
      ------------------------------------------------------------------
      14  224875    0.999431    0.000101    0.999360    0.999547    3
      15  224282    0.996796    0.000586    0.996387    0.997467    3
      16  203690    0.905277    0.010046    0.897628    0.916654    3
      17  201514    0.895606    0.002073    0.893748    0.897841    3
      18  142854    0.634898    0.004184    0.630552    0.638898    3
      19  107469    0.477634    0.005003    0.471900    0.481114    3
      20  105722    0.469869    0.006820    0.463314    0.476927    3
      21   94238    0.418830    0.006215    0.413848    0.425794    3
      22   74471    0.330978    0.001716    0.329889    0.332956    3
      23   65000    0.288885    0.008498    0.280143    0.297116    3
      24   47152    0.209562    0.003326    0.206357    0.212997    3
      25   47105    0.209353    0.002234    0.207584    0.211864    3
      26   43148    0.191766    0.001816    0.190064    0.193677    3
      27   42646    0.189535    0.001835    0.188251    0.191637    3
      28   42291    0.187957    0.002006    0.186051    0.190051    3
      29   34936    0.155269    0.000534    0.154745    0.155811    3
      30   34581    0.153691    0.002382    0.151251    0.156011    3
      31   33312    0.148051    0.004020    0.144025    0.152065    3
      32   28053    0.124678    0.002840    0.122305    0.127825    3
      33   27648    0.122878    0.001253    0.122038    0.124318    3
      34   27373    0.121656    0.003487    0.118158    0.125132    3
      35   23889    0.106172    0.003544    0.103425    0.110172    3
      36   23850    0.105999    0.000324    0.105625    0.106199    3
      37   23552    0.104674    0.004612    0.099465    0.108239    3
      ------------------------------------------------------------------
      + Convergence diagnostic (standard deviation of split frequencies)
        should approach 0.0 as runs converge.


      Summary statistics for branch and node parameters
         (saved to file "/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-mb.nex.vstat"):

                                              95% HPD Interval
                                            --------------------
      Parameter      Mean       Variance     Lower       Upper       Median     PSRF+  Nruns
      --------------------------------------------------------------------------------------
      length[1]     0.001582    0.000002    0.000000    0.004461    0.001173    1.000    3
      length[2]     0.256266    0.002491    0.168908    0.357457    0.250852    1.000    3
      length[3]     0.000973    0.000001    0.000000    0.002968    0.000670    1.000    3
      length[4]     0.013027    0.000017    0.005526    0.021293    0.012695    1.000    3
      length[5]     0.000975    0.000001    0.000000    0.002927    0.000672    1.000    3
      length[6]     0.003769    0.000008    0.000000    0.008938    0.003261    1.000    3
      length[7]     0.027155    0.000078    0.011117    0.044491    0.026235    1.000    3
      length[8]     0.001115    0.000001    0.000000    0.003330    0.000780    1.000    3
      length[9]     0.001130    0.000001    0.000000    0.003420    0.000777    1.000    3
      length[10]    0.000980    0.000001    0.000000    0.002939    0.000672    1.000    3
      length[11]    0.001655    0.000002    0.000000    0.004585    0.001252    1.000    3
      length[12]    0.001449    0.000002    0.000000    0.004385    0.000992    1.000    3
      length[13]    0.002622    0.000004    0.000000    0.006664    0.002127    1.000    3
      length[14]    0.058163    0.000429    0.021818    0.100441    0.056170    1.000    3
      length[15]    0.019630    0.000055    0.006240    0.034678    0.018950    1.000    3
      length[16]    0.032525    0.000262    0.000494    0.061170    0.031531    1.000    3
      length[17]    0.003265    0.000004    0.000199    0.007207    0.002883    1.000    3
      length[18]    0.001937    0.000002    0.000000    0.004989    0.001559    1.000    3
      length[19]    0.002940    0.000006    0.000000    0.007649    0.002369    1.000    3
      length[20]    0.003542    0.000006    0.000000    0.008459    0.003026    1.000    3
      length[21]    0.003252    0.000007    0.000000    0.008341    0.002655    1.000    3
      length[22]    0.002991    0.000006    0.000000    0.007831    0.002404    1.000    3
      length[23]    0.002472    0.000004    0.000000    0.006315    0.002003    1.000    3
      length[24]    0.002405    0.000004    0.000001    0.006170    0.001944    1.000    3
      length[25]    0.001292    0.000002    0.000000    0.003836    0.000911    1.000    3
      length[26]    0.000971    0.000001    0.000000    0.002925    0.000677    1.000    3
      length[27]    0.000980    0.000001    0.000000    0.002977    0.000676    1.000    3
      length[28]    0.000979    0.000001    0.000000    0.002937    0.000678    1.000    3
      length[29]    0.001117    0.000001    0.000000    0.003324    0.000777    1.000    3
      length[30]    0.001132    0.000001    0.000000    0.003414    0.000777    1.000    3
      length[31]    0.001491    0.000002    0.000000    0.004291    0.001069    1.000    3
      length[32]    0.000976    0.000001    0.000000    0.002979    0.000675    1.000    3
      length[33]    0.000969    0.000001    0.000000    0.002945    0.000662    1.000    3
      length[34]    0.000975    0.000001    0.000000    0.002956    0.000677    1.000    3
      length[35]    0.000975    0.000001    0.000000    0.002925    0.000680    1.000    3
      length[36]    0.000982    0.000001    0.000000    0.002961    0.000680    1.000    3
      length[37]    0.000953    0.000001    0.000000    0.002867    0.000660    1.000    3
      --------------------------------------------------------------------------------------
      + Convergence diagnostic (PSRF = Potential Scale Reduction Factor; Gelman
        and Rubin, 1992) should approach 1.0 as runs converge. NA is reported when
        deviation of parameter values within all runs is 0 or when a parameter
        value (a branch length, for instance) is not sampled in all runs.


      Summary statistics for partitions with frequency >= 0.10 in at least one run:
          Average standard deviation of split frequencies = 0.003311
          Maximum standard deviation of split frequencies = 0.010046
          Average PSRF for parameter values (excluding NA and >10.0) = 1.000
          Maximum PSRF for parameter values = 1.000


      Clade credibility values:

      /---------------------------------------------------------- AF448070.1_Syne~ (2)
      |
      |                                           /-------------- FJ895126.1_Apha~ (3)
      |                                           |
      |                                           |-------------- FJ895128.1_Apha~ (4)
      |                                           |
      |                            /------100-----+-------------- FJ895125.1_Apha~ (5)
      |                            |              |
      |                            |              |-------------- FJ895124.1_Aph~ (10)
      |-------------100------------+              |
      +                            |              \-------------- FJ895127.1_Aph~ (13)
      |                            |
      |                            \----------------------------- EF685373.1_Apha~ (7)
      |
      |              /------------------------------------------- FJ895118.1_Apha~ (1)
      |              |
      |              |------------------------------------------- FJ895119.1_Apha~ (6)
      |              |
      |              |                            /-------------- FJ895121.1_Apha~ (8)
      \------91------+             /------63------+
                     |             |              \-------------- FJ895120.1_Apha~ (9)
                     |------90-----+
                     |             \----------------------------- FJ895122.1_Aph~ (11)
                     |
                     \------------------------------------------- FJ895123.1_Aph~ (12)


      Phylogram (based on average branch lengths):

      /---------------------------------------------------------- AF448070.1_Syne~ (2)
      |
      |                /- FJ895126.1_Apha~ (3)
      |                |
      |                |--- FJ895128.1_Apha~ (4)
      |                |
      |            /---+- FJ895125.1_Apha~ (5)
      |            |   |
      |            |   |- FJ895124.1_Aph~ (10)
      |------------+   |
      +            |   \- FJ895127.1_Aph~ (13)
      |            |
      |            \------ EF685373.1_Apha~ (7)
      |
      |      /- FJ895118.1_Apha~ (1)
      |      |
      |      |- FJ895119.1_Apha~ (6)
      |      |
      |      |/ FJ895121.1_Apha~ (8)
      \------+|
             || FJ895120.1_Apha~ (9)
             |+
             |\ FJ895122.1_Aph~ (11)
             |
             \- FJ895123.1_Aph~ (12)

      |----------| 0.050 expected changes per site

      Calculating tree probabilities...

      Credible sets of trees (18447 trees sampled):
         50 % credible set contains 402 trees
         90 % credible set contains 4803 trees
         95 % credible set contains 8627 trees
         99 % credible set contains 16197 trees

   Exiting mrbayes block
   Reached end of file

   Tasks completed, exiting program because mode is noninteractive
   To return control to the command line after completion of file processing,
   set mode to interactive with 'mb -i <filename>' (i is for interactive)
   or use 'set mode=interactive'

(base) jacquelinelemaire@Jacquelines-Air mrbayes %
```

Note: If needed convert consensus tree files (con.tre) from Nexus to Newick (.tree) files in [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) by exporting and saving. Or convert in r - but need to troubleshoot how to do this

## Plot the best tree to see how it looks
First read in the consensus tree
```r
mrbayestree_aphani <- read.tree(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/results/aphani-muscle-mrbayes-con.tree")
```
View and reroot the tree
```r
#The tree should already be rooted as we specified the outgroup in the mr bayes file but let's check

mrbayestree_aphani

#Nope, it is not rooted, so let's root the tree

Mrbayes.rooted <- root(mrbayestree_aphani, outgroup = "'AF448070.1_Synechococcus_sp._PS845'", resolve.root = TRUE)

#check that the tree is now rooted
Mrbayes.rooted
#yes
```
Plot the Tree
```r
#view the order of tip labels
Mrbayes.rooted$tip.label
Mrbayes.rooted$

#coloring the tip labels like this is not very scalable but it works for this small dataset
#Make a vector to color the different species
species_colors_B <- c("black", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen","darkred", "darkblue", "darkblue", "darkblue","darkblue", "darkblue", "darkblue")

#Plot the tree
plot(Mrbayes.rooted, cex=.7, edge.color = "black", tip.color = species_colors_B, main="A Bayesian (MrBayes) Tree of Aphanizomenon Muscle Alignment")

#Need to troubleshoot how to add the postier probabilty support values to the tree
```

### FIN


