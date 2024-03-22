## Jackie Lemaire
## Phylogenetics Class 563
## Distance and Parsimony III Class: Distance and Parsimony methods (Part 3: computer lab) 
# 2024-03-14
#RaxML tree

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

#3) Loading the newick tree from raxML output
raxmltree <- read.tree(file = "/Applications/raxml-ng_v1.2.1_macos_x86_64/T3-myseed.raxml.bestTree")

plot(raxmltree)







