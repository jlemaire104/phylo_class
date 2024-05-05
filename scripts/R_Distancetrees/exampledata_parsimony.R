## Jackie Lemaire
## Phylogenetics Class 563
## Distance and Parsimony III Class: Distance and Parsimony methods (Part 3: computer lab) 
# 2024-02-27
#Parismony tree

# Software: R package phangorn: Parsimony-based methods
# phangorn is another widely used phylogenetic software
# It is an R package and it has a huge variety of functions
# In particular, today we will use it for parsimony-based tree estimation methods
# Full documentation
# We will follow this great tutorial

#install neccessary packages
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

#loading packages
library(ape)
library(adegenet)
library(phangorn)

#3) Loading the sample data and convert to phangorn object:
dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")
dna2 <- as.phyDat(dna)

#4) We need a starting tree for the search on tree space and compute the parsimony score of this tree (422)
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)

#5) Search for the tree with maximum parsimony:
tre.pars <- optim.parsimony(tre.ini, dna2)
## Final p-score 420 after 2 nni operations

#6) Plot tree:
plot(tre.pars, cex=0.6)




