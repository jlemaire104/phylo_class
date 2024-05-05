## Jackie Lemaire
## Phylogenetics Class 563
## Distance and Parsimony III Class: Distance and Parsimony methods (Part 3: computer lab) 
# 2024-02-27
# Distance Tree

## Software
# Software R package ape: Distance-based methods

## Description
# ape is one of the most widely used phylogenetic software
# It is an R package and it has a huge variety of functions
# In particular, today we will use it for distance-based tree estimation methods
## compute a matrix of pairwise genetic distances between the studied taxa, and summarize it using a hierarchical clustering
#algorithm such as UPGMA or Neighbour-Joining (Supplementary Figure S1).

##Strengths
#Advantages: fast (the fastest) and flexible (different genetic distances allow to account for different features of DNA sequence evolution). 

##Weaknesses	
#Limitations: no model comparison (can’t test for the ’best’ tree, or the ’best’ model of evolution); may be inaccurate and highly dependent on the distance and clustering
#algorithm chosen.

##Assumptions	User choices
# We first compute genetic distances using ape’s dist.dna, which proposes no less
# than 15 different genetic distances (see ?dist.dna for details). Here, we use Tamura
# and Nei 1993’s model [8] which allows for different rates of transitions and transversions,

# Summary of Software R package ape
# Main distance functions:
# nj (ape package): the classical Neighbor-Joining algorithm.
# bionj (ape): an improved version of Neighbor-Joining: Gascuel 1997. It uses information on variances of evolutionary distances
# fastme.bal and fastme.ols (ape): minimum evolution algorithms: Desper and Gascuel, 2002
# hclust (stats): classical hierarchical clustering algorithms including single linkage, complete linkage, UPGMA, and others.


#install neccessary packages
install.packages("adegenet", dep=TRUE)
install.packages("phangorn", dep=TRUE)

#load packages
library(ape)
library(adegenet)
library(phangorn)

#loading the sample data
dna <- fasta2DNAbin(file="/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_125_cyanos/combinedtest3.fasta")
dna
class(dna)
as.character(dna)[1:5,1:10]
unclass(dna)[1:5,1:10]
typeof(unclass(dna)[1:5,1:10])


#This results in significant savings in terms of memory required to represent the data:
object.size(as.character(dna))/object.size(dna)
#7.9 bytes

## Computing the genetic distances. They choose a Tamura and Nei 1993 model which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate (more on Models of Evolution).
D <- dist.dna(dna, model="TN93")

# Get the NJ tree
tre <- njs(D)

# Before plotting, we can use the ladderize function which reorganizes the internal structure of the tree to get the ladderized effect when plotted
tre <- ladderize(tre)

# We can plot the tree
plot(tre, cex=.6)
title("A simple NJ tree")

