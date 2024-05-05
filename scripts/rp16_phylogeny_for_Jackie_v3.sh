#################################
#################################
# Generate rp16-based phylogeny
# Code from Ben Peterson
#################################
#################################

cd /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/cyano_76/analysis/hmmsearch
mkdir rp16
scripts=/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/scripts # directory where you keep your code

# If necessary, create and then activate a conda environment with hmmer. I installed mine in "hmmer"

conda activate hmmer
PYTHONPATH='' #do we need this?
PERL5LIB='' #do we need this?


#########################
# Use hmmsearch to identify sequences.
#########################

## Work in project directory
cd /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/hmmsearch

#copy over the bins_list.txt file to working directory
cp /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_cyanos_76_fastas/bins_list.txt .

cat /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/databases/rp16_list.txt | while read gene
  do
    if [ ! -d rp16/$gene ]; then
      echo "Working on pulling out" $gene "genes"
      mkdir rp16/$gene
    else
      echo $gene "analysis started previously"
    fi


  cat bins_list.txt | while read bin
  do
    if [ ! -e rp16/$gene/$bin.out ]; then
      echo "Pulling" $gene "out of" $bin
      hmmsearch --tblout rp16/$gene/$bin.out \
                --cpu 2 \
                --cut_tc \
                /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/databases/$gene\_bact.HMM \
                /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/markergenes/ORFs/$bin.faa \
                > rp16/$gene/$bin\HMM_output.txt
    else
      echo "Searched" $bin "for" $gene "already"
    fi
  done
done



#########################
# Pull out amino acid sequences and align them
#########################
# If necessary, create and then activate a conda environment with hmmer. I installed mine in "hmmer"

conda create -n hmmer
conda activate hmmer
conda install biopython

#this worked with python 3.8


cat /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/databases/rp16_list.txt | while read gene
do
  cat bins_list.txt | while read bin
  do
    if [ ! -e rp16/$gene/$bin.faa ]; then
      lineCount=`wc -l < rp16/$gene/$bin.out`
      if [ $lineCount -eq 13 ]; then
        echo "No" $gene "hits in" $bin
      else
        echo "Found" $gene "in" $bin
        python $scripts/extract_protein_hitting_HMM.py \
                                          rp16/$gene/$bin.out \
                                          ../markergenes/ORFs/$bin.faa \
                                          rp16/$gene/$bin.faa \
                                          $bin
      fi
    else
      echo "Pulled" $gene "out of" $bin "already"
    fi
  done
done



#  Prepare files for alignment

cd /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/hmmsearch/rp16
cat /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/databases/rp16_list.txt | while read gene
do
  cat $gene/*.faa > $gene.bins.faa
done


# make the alignment

cat /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/databases/rp16_list.txt | while read gene
do
  muscle -align $gene.bins.faa \
          -output $gene.afa
done

This worked great! And I also aligned the sequences with clustalW

#  concatenate


seqkit concat *.afa > rp16concat.afa

FastTree rp16concat.afa > rp16_fast.tree

# gave error that sequences weren't same length
# This must be because some genomes are missing one or more gene
# How do I fix this without Geneious?


# Ended up using Geneious. See analysis-notes-cyanos.md and /Users/trinamcmahon2017/Research/TrinaComputingProjects/Cyanos/analysis/cyano_phylogeny_rp16.md

#Trina edited the sequences in Geneious to remove genes that were not found in certain samples:
#OK here is the concatenated MSA. I took out some genomes that only one or two genes recovered, and then masked sites that had >30% gaps. The first MCYST is sparse, but figured it was interesting enough to keep in.  I’m a bit worried that the names got garbled, so if you have questions about it LMK.  The letters are fine, it’s just taht GEneious added some numbers to the end which may not match the numbers at the end in the original genomes. If you find something weird, LMK and we can figure it out. But for the purposes of a tree, it will be fine!
#I actually took out the gene that was causing the naming problem. It was the rp22. Maybe there were two copies of it in a MAG??



# Exported an alignment of 199 genomes / 16 ribosomal proteins, 2,498 residues

/Users/trinamcmahon2017/Research/TrinaComputingProjects/Cyanos/analysis/rp16concat-masked-199seqs.fasta

# Looked at JK-distance tree in Geneious and distance matrix

#grep ">" rp_concat_15_curated_mask30.fasta | wc -l
#60

#the rp proteins are all concatenated together now and there seems to be 60 different genomes in the file

#########################
# Generate RAxML tree
#########################

scripts=~/Cyanos/code

source /home/GLBRCORG/trina.mcmahon/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS

cd /Applications/raxml-ng_v1.2.1_macos_x86_64
raxml./ --all \
        --parse \
        --msa /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/hmmsearch/rp16/alignment-files-all/concatenated-alignment/rp_concat_15_curated_mask30.fasta \
        --model PROTGTR+G \
        --bs-trees autoMRE \ &



What each flag means:
The -f option is a very powerful option because, in many cases, it allows you to select what
kind of algorithm RAxML shall execute. If you don't specify ­f at all RAxML will execute the
standard hill climbing algorithm (­f d option)The individual options are listed below.
a (-f a) rapid Bootstrap analysis and search for best­scoring ML tree in one program
run
        -p 283976 = Specify a random number seed for the parsimony inferences. This allows you
to reproduce your results and will help me debug the program.
        -m PROTGAMMAAUTO = The example below will automatically determine which is the best (the one with the highest
likelihood score on the parsimony starting tree) protein substitution model for your
dataset using the base frequencies that come with the models. It will chose among the
following models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62,
MTMAM, LG,  MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, DUMMY, DUMMY2.
These models will not be considered: LG4M, LG4X, PROT_FILE,GTR_UNLINKED, GTR!
        -N autoMRE =
        -x 2381 \
        -T 10 \
        -s rp16concat-masked-199seqs.fasta \
        -n rp16 = Specifies the name of the output file.
        
# oof this takes a LONG TIME

#########################
# Root bacterial tree
#########################

# First we'll need to get our archaeal rp16
# to use for rooting

scripts=~/5M/code/generalUse

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

cd ~/5M/dataEdited/binAnalysis/rp16
mkdir rooting

for gene in $(cat ~/references/rp16/rp16_list.txt)
do

  cd ~/5M/dataEdited/binAnalysis/

  if [ ! -d rp16/rooting/$gene ]; then
    echo "Working on pulling out" $gene "genes"
    mkdir rp16/rooting/$gene
  else
    echo $gene "analysis started previously"
  fi


  for bin in $(echo GCA_000376965.1 GCA_000430905.1)
  do
    binFile=$(awk -F ',' -v bin="$bin" '$3 == bin { print $2 }' ~/references/genomes/bins_metadata.csv \
    | sed 's/.fna/.faa/')

    if [ ! -e rp16/rooting/$gene/$bin.out ]; then
      echo "Pulling" $gene "out of" $bin "using" $binFile
      hmmsearch --tblout rp16/rooting/$gene/$bin.out \
                --cpu 4 \
                --cut_nc \
                ~/references/rp16/$gene\_arch.HMM \
                ~/references/genomes/ORFs/$binFile \
                > rp16/rooting/$gene/$bin\HMM_output.txt
    else
      echo "Searched" $bin "for" $gene "already"
    fi
  done
done


for gene in $(cat ~/references/rp16/rp16_list.txt)
do

  cd ~/5M/dataEdited/binAnalysis/

  for bin in $(echo GCA_000376965.1 GCA_000430905.1)
  do
    if [ ! -e rp16/rooting/$gene/$bin.faa ]; then
      lineCount=`wc -l < rp16/rooting/$gene/$bin.out`
      if [ $lineCount -eq 13 ]; then

        echo "No" $gene "hits in" $bin

      else

        binFile=$(awk -F ',' -v bin="$bin" '$3 == bin { print $2 }' ~/references/genomes/bins_metadata.csv \
                  | sed 's/.fna/.faa/')

        echo "Found" $gene "in" $bin
        python $scripts/extract_protein_hitting_HMM_classification.py \
                                          rp16/rooting/$gene/$bin.out \
                                          ~/references/genomes/ORFs/$binFile \
                                          rp16/rooting/$gene/$bin.faa \
                                          $bin
      fi
    else
      echo "Pulled" $gene "out of" $bin "already"
    fi
  done

  cd rp16/rooting
  rm -f $gene.faa $gene.afa

  cat $gene/*.faa \
      ../$gene.faa \
      > $gene.faa

  muscle -in $gene.faa \
          -out $gene.afa

done


# Generate taxtastic package
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate pipenv
cd ~/5M/dataEdited/binAnalysis/rp16

taxit create -l rp16 \
            -P rp16.refpkg \
            --aln-fasta rp16_masked.afa \
            --tree-stats RAxML_info.rp16 \
            --tree-file RAxML_bipartitions.rp16
conda deactivate

sed 's/AUTO/LG/' phylo_model9hsd9hjn.json > phylo_model9hsd9hjn_edited.json
mv -f phylo_model9hsd9hjn_edited.json phylo_model9hsd9hjn.json

# Then run pplacer
cd ~/5M/dataEdited/binAnalysis/rp16/rooting

conda activate bioinformatics
pplacer -c ../rp16.refpkg \
        rp16_rooting.fasta



#########################
# Generate tree with outgroup
#########################

FastTree rp16_rooting.fasta > rp16_rooting.tree
conda deactivate

cd ~/5M/dataEdited/binAnalysis/rp16/rooting
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 10 \
        -s rp16_rooting.fasta \
        -n rp16_rooting




#################################
#################################
# Generate tree of all bins with Karthik's database
#################################
#################################

# Based this on old workflow (sandbox/bin_analysis_phylogeny.sh)
cp /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16.tree \
    /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree
awk '{ print $1 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/bin_names_to_switch.tsv | while read old_name
do
  new_name=$(awk -v old_name="$old_name" '$1 == old_name { print $2 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/bin_names_to_switch.tsv)
  sed "s/$old_name/$new_name/" /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree \
      > /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree.txt
  mv /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree.txt \
      /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree
done
awk -F ',' '{ print $1 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/binsForAnalysis/newName.csv | while read old_name
do
  new_name=$(awk -F ',' -v old_name="$old_name" '$1 == old_name { print $2 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/binsForAnalysis/newName.csv)
  sed "s/$old_name/$new_name/" /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree \
      > /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree.txt
  mv /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree.txt \
      /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/karthiks_db/bin_phylogeny_rp16_edited.tree
done
