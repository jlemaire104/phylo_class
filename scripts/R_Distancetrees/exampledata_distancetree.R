## Jackie Lemaire
## Phylogenetics Class 563
## Distance and Parsimony III Class: Distance and Parsimony methods (Part 3: computer lab) 
# 2024-02-27
# Distance Tree

# Software R package ape: Distance-based methods
# ape is one of the most widely used phylogenetic software
# It is an R package and it has a huge variety of functions
# In particular, today we will use it for distance-based tree estimation methods
# Full documentation
# We will follow this great tutorial

#install neccessary packages
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

#load packages
library(ape)
library(adegenet)
library(phangorn)

#loading the sample data
dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")

## Computing the genetic distances. They choose a Tamura and Nei 1993 model which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate (more on Models of Evolution).
D <- dist.dna(dna, model="TN93")

?dist.dna
# Get the NJ tree
tre <- nj(D)

# Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
tre <- ladderize(tre)

# 7) We can plot the tree
plot(tre, cex=.6)
title("A simple NJ tree")

# Summary of Software R package ape
# Main distance functions:

# nj (ape package): the classical Neighbor-Joining algorithm.
# bionj (ape): an improved version of Neighbor-Joining: Gascuel 1997. It uses information on variances of evolutionary distances
# fastme.bal and fastme.ols (ape): minimum evolution algorithms: Desper and Gascuel, 2002
# hclust (stats): classical hierarchical clustering algorithms including single linkage, complete linkage, UPGMA, and others.




