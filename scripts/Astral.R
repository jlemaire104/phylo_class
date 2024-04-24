#Astral
# 2024-04-11
#Jackie Lemaire


#Install Java for Mac OS ARM64 
#go to the website and click the download button
#restart terminal
#in terminal
#type
java
#if it is installed properly and working it will call java and show the help info

##Install Astral
cd /Applications/Astral/

For ASTRAL-III
java -jar /Applications/Astral/astral.5.7.8.jar

##Software	
Astral
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees.


##Description	
#The program Astral ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS).
#ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. This predefined set is empirically decided by ASTRAL (but see tutorial on how to expand it). The current code corresponds to **ASTRAL-III**

##Strengths
#Astral uses A leading model of gene evolution is the multi-species coalescent (MSC) (see Degnan and Rosenberg,
#2009). MSC models incomplete lineage sorting (ILS) and the resulting discordance between gene trees and
#the species tree.

##Weaknesses	
##Assumptions
#Mr. Bayes uses MCMC to appoximate the posterior probabilites of trees. Uses command line interface. Reads in the standard NEXUS format

##User choices
#User can change assumptions of the substitution model, the prior, and details for the MC3 




## Mr Bayes Install

#Run in bash terminal:

java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre
java -Djava.library.path=./lib/ -jar astralmp.5.7.8.jar -i test_data/song_mammals.424.gene.tre


#Try this:
arch -arm64 brew install mrbayes
#It started to install and then this error popped up: Error: gcc: the bottle needs the Apple Command Line Tools to be installed.

#Install command line tools with xcode:
xcode-select --install

#Try to install again:
arch -arm64 brew install mrbayes

##Success!!
which mb






#output


