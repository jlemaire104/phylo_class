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
clustalw2 -ALIGN -INFILE=Documents/phylogenetics-class/data/primatesAA.fasta -OUTFILE=Documents/phylo_class/phylo_class/results/Primate_alignment/primatesAA-aligned.fasta -OUTPUT=FASTA

# ClutalW Documentation
DATA (sequences)

-INFILE=file.ext                             :input sequences.
-PROFILE1=file.ext  and  -PROFILE2=file.ext  :profiles (old alignment).


                VERBS (do things)

-OPTIONS            :list the command line parameters
-HELP  or -CHECK    :outline the command line params.
-FULLHELP           :output full help content.
-ALIGN              :do full multiple alignment.
-TREE               :calculate NJ tree.
-PIM                :output percent identity matrix (while calculating the tree)
-BOOTSTRAP(=n)      :bootstrap a NJ tree (n= number of bootstraps; def. = 1000).
-CONVERT            :output the input sequences in a different file format.

                PARAMETERS (set things)

***General settings:****
-INTERACTIVE :read command line, then enter normal interactive menus
-QUICKTREE   :use FAST algorithm for the alignment guide tree
-TYPE=       :PROTEIN or DNA sequences
-NEGATIVE    :protein alignment with negative values in matrix
-OUTFILE=    :sequence alignment file name
-OUTPUT=     :GCG, GDE, PHYLIP, PIR or NEXUS
-OUTORDER=   :INPUT or ALIGNED
-CASE        :LOWER or UPPER (for GDE output only)
-SEQNOS=     :OFF or ON (for Clustal output only)
-SEQNO_RANGE=:OFF or ON (NEW: for all output formats)
-RANGE=m,n   :sequence range to write starting m to m+n
-MAXSEQLEN=n :maximum allowed input sequence length
-QUIET       :Reduce console output to minimum
-STATS=      :Log some alignents statistics to file

***Multiple Alignments:***
-NEWTREE=      :file for new guide tree
-USETREE=      :file for old guide tree
-MATRIX=       :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename
-DNAMATRIX=    :DNA weight matrix=IUB, CLUSTALW or filename
-GAPOPEN=f     :gap opening penalty        
-GAPEXT=f      :gap extension penalty
-ENDGAPS       :no end gap separation pen. 
-GAPDIST=n     :gap separation pen. range
-NOPGAP        :residue-specific gaps off  
-NOHGAP        :hydrophilic gaps off
-HGAPRESIDUES= :list hydrophilic res.    
-MAXDIV=n      :% ident. for delay
-TYPE=         :PROTEIN or DNA
-TRANSWEIGHT=f :transitions weighting
-ITERATION=    :NONE or TREE or ALIGNMENT
-NUMITER=n     :maximum number of iterations to perform
-NOWEIGHTS     :disable sequence weighting

## ==ITERATION==

 A remove first iteration scheme has been added. This can be used to improve the final
 alignment or improve the alignment at each stage of the progressive alignment. During the 
 iteration step each sequence is removed in turn and realigned. If the resulting alignment 
 is better than the  previous alignment it is kept. This process is repeated until the score
 converges (the  score is not improved) or until the maximum number of iterations is 
 reached. The user can  iterate at each step of the progressive alignment by setting the 
 iteration parameter to  TREE or just on the final alignment by seting the iteration 
 parameter to ALIGNMENT. The default is no iteration. The maximum number of  iterations can 
 be set using the numiter parameter. The default number of iterations is 3.
  
 -ITERATION=    :NONE or TREE or ALIGNMENT
 
 -NUMITER=n     :Maximum number of iterations to perform


# Tried to install Muscle but it doesn't seem to be working with the new MacOS M chips
## Installed muscle via bioconda - but need to test that it is working


## Alignment HW - running at least one alignment program on my data 
# Using ClustalW for alignment

# First all the separate fna files need to be combined into one fasta file

cd Documents/phylo_class/phylo_class/results/data/diamond_125_cyanos/
cat *.fna > combined.fasta

clustalw2 -ALIGN -INFILE=Documents/phylo_class/phylo_class/results/data/diamond_125_cyanos/combined.fasta -OUTFILE=Documents/phylo_class/phylo_class/results/cyano-aligned.fasta -TYPE=DNA -QUIET -OUTPUT=FASTA

# Getting an error

# Tried a simple input to test out
clustalw -infile=combined.fasta -seqnos=ON

## -SEQNOS=  :OFF or ON (for Clustal output only)

# The above command worked and began running - need to test out adding this flag to the longer input command