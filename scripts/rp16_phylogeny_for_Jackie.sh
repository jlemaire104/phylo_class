#################################
#################################
# Generate rp16-based phylogeny
# Code from Ben Peterson
#################################
#################################

cd /home/GLBRCORG/trina.mcmahon/Cyanos/analysis/phylogeny
mkdir rp16
scripts=~/Cyanos/code # directory where you keep your code

# If necessary, create and then activate a conda environment with hmmer. I installed mine in "bioinformatics"

source /home/GLBRCORG/trina.mcmahon/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''


#########################
# Use hmmsearch to identify sequences.
#########################

## Work in project directory

cat /home/GLBRCORG/trina.mcmahon/references/rp16/rp16_list.txt | while read gene
  do
    if [ ! -d phylogeny/rp16/$gene ]; then
      echo "Working on pulling out" $gene "genes"
      mkdir phylogeny/rp16/$gene
    else
      echo $gene "analysis started previously"
    fi


  cat bins_list.txt | while read bin
  do
    if [ ! -e phylogeny/rp16/$gene/$bin.out ]; then
      echo "Pulling" $gene "out of" $bin
      hmmsearch --tblout phylogeny/rp16/$gene/$bin.out \
                --cpu 2 \
                --cut_tc \
                ~/references/rp16/$gene\_bact.HMM \
                /home/GLBRCORG/trina.mcmahon/Cyanos/analysis/binAnalysis/ORFs/$bin.faa \
                > phylogeny/rp16/$gene/$bin\HMM_output.txt
    else
      echo "Searched" $bin "for" $gene "already"
    fi
  done
done



#########################
# Pull out amino acid sequences and align them
#########################
scripts=~/Cyanos/code

source /home/GLBRCORG/trina.mcmahon/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

PERL5LIB=''

PYTHONPATH=/opt/bifxapps/python/lib/python3.4/site-packages/
cd ~/Cyanos/analysis/

cat ~/references/rp16/rp16_list.txt | while read gene
do
  cat bins_list.txt | while read bin
  do
    if [ ! -e phylogeny/rp16/$gene/$bin.faa ]; then
      lineCount=`wc -l < phylogeny/rp16/$gene/$bin.out`
      if [ $lineCount -eq 13 ]; then
        echo "No" $gene "hits in" $bin
      else
        echo "Found" $gene "in" $bin
        python $scripts/extract_protein_hitting_HMM_classification.py \
                                          phylogeny/rp16/$gene/$bin.out \
                                          binAnalysis/ORFs/$bin.faa \
                                          phylogeny/rp16/$gene/$bin.faa \
                                          $bin
      fi
    else
      echo "Pulled" $gene "out of" $bin "already"
    fi
  done




# Now pull from the reference genomes (note the file structure is dumb, would need to
# revisit next time to make it better)

  cat /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/ref_list_clean.txt | while read genome
  do
    if [ ! -e /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16/$gene/$genome.faa ]; then
      lineCount=`wc -l < /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16/$gene/$genome.out`
      if [ $lineCount -eq 13 ]; then

        echo "No" $gene "hits in" $genome

      else

        echo "Found" $gene "in" $genome
        python $scripts/extract_protein_hitting_HMM_classification.py \
                                          /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16/$gene/$genome.out \
                                          /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/ORFs/$genome.faa \
                                          /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16/$gene/$genome.faa \
                                          $genome
      fi
    else
      echo "Pulled" $gene "out of" $genome "already"
    fi
  done
done

#  Prepare files for alignment

cd /home/GLBRCORG/trina.mcmahon/Cyanos/analysis/phylogeny/rp16
cat ~/references/rp16/rp16_list.txt | while read gene
do
  cat $gene/*.faa > $gene.bins.faa
done

cd /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16
cat ~/references/rp16/rp16_list.txt | while read gene
do
  cat $gene/*.faa > $gene.reference.faa

done

# copy the references into the main rp16 analysis folder

cd /home/GLBRCORG/trina.mcmahon/Cyanos/data/gtdbtk_refs/phylogeny/rp16
cp *.faa /home/GLBRCORG/trina.mcmahon/Cyanos/analysis/phylogeny/rp16/

# combine the bins and references by gene

cd /home/GLBRCORG/trina.mcmahon/Cyanos/analysis/phylogeny/rp16
cat ~/references/rp16/rp16_list.txt | while read gene
do
  cat $gene.bins.faa $gene.reference.faa > $gene.faa
done

# make the alignment

cat ~/references/rp16/rp16_list.txt | while read gene
do
  muscle -in $gene.faa \
          -out $gene.afa
done


#  concatenate


seqkit concat *.afa > rp16concat.afa

FastTree rp16concat.afa > rp16_fast.tree

# gave error that sequences weren't same length
# This must be because some genomes are missing one or more gene
# How do I fix this without Geneious?


# Ended up using Geneious. See analysis-notes-cyanos.md and /Users/trinamcmahon2017/Research/TrinaComputingProjects/Cyanos/analysis/phylogeny/cyano_phylogeny_rp16.md

# Exported an alignment of 199 genomes / 16 ribosomal proteins, 2,498 residues

/Users/trinamcmahon2017/Research/TrinaComputingProjects/Cyanos/analysis/phylogeny/rp16concat-masked-199seqs.fasta

# Looked at JK-distance tree in Geneious and distance matrix



#########################
# Generate RAxML tree
#########################

scripts=~/Cyanos/code

source /home/GLBRCORG/trina.mcmahon/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 10 \
        -s rp16concat-masked-199seqs.fasta \
        -n rp16 \ &

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

cd ~/5M/dataEdited/binAnalysis/phylogeny/rp16
mkdir rooting

for gene in $(cat ~/references/rp16/rp16_list.txt)
do

  cd ~/5M/dataEdited/binAnalysis/

  if [ ! -d phylogeny/rp16/rooting/$gene ]; then
    echo "Working on pulling out" $gene "genes"
    mkdir phylogeny/rp16/rooting/$gene
  else
    echo $gene "analysis started previously"
  fi


  for bin in $(echo GCA_000376965.1 GCA_000430905.1)
  do
    binFile=$(awk -F ',' -v bin="$bin" '$3 == bin { print $2 }' ~/references/genomes/bins_metadata.csv \
    | sed 's/.fna/.faa/')

    if [ ! -e phylogeny/rp16/rooting/$gene/$bin.out ]; then
      echo "Pulling" $gene "out of" $bin "using" $binFile
      hmmsearch --tblout phylogeny/rp16/rooting/$gene/$bin.out \
                --cpu 4 \
                --cut_nc \
                ~/references/rp16/$gene\_arch.HMM \
                ~/references/genomes/ORFs/$binFile \
                > phylogeny/rp16/rooting/$gene/$bin\HMM_output.txt
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
    if [ ! -e phylogeny/rp16/rooting/$gene/$bin.faa ]; then
      lineCount=`wc -l < phylogeny/rp16/rooting/$gene/$bin.out`
      if [ $lineCount -eq 13 ]; then

        echo "No" $gene "hits in" $bin

      else

        binFile=$(awk -F ',' -v bin="$bin" '$3 == bin { print $2 }' ~/references/genomes/bins_metadata.csv \
                  | sed 's/.fna/.faa/')

        echo "Found" $gene "in" $bin
        python $scripts/extract_protein_hitting_HMM_classification.py \
                                          phylogeny/rp16/rooting/$gene/$bin.out \
                                          ~/references/genomes/ORFs/$binFile \
                                          phylogeny/rp16/rooting/$gene/$bin.faa \
                                          $bin
      fi
    else
      echo "Pulled" $gene "out of" $bin "already"
    fi
  done

  cd phylogeny/rp16/rooting
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
cd ~/5M/dataEdited/binAnalysis/phylogeny/rp16

taxit create -l rp16 \
            -P rp16.refpkg \
            --aln-fasta rp16_masked.afa \
            --tree-stats RAxML_info.rp16 \
            --tree-file RAxML_bipartitions.rp16
conda deactivate

sed 's/AUTO/LG/' phylo_model9hsd9hjn.json > phylo_model9hsd9hjn_edited.json
mv -f phylo_model9hsd9hjn_edited.json phylo_model9hsd9hjn.json

# Then run pplacer
cd ~/5M/dataEdited/binAnalysis/phylogeny/rp16/rooting

conda activate bioinformatics
pplacer -c ../rp16.refpkg \
        rp16_rooting.fasta



#########################
# Generate tree with outgroup
#########################

FastTree rp16_rooting.fasta > rp16_rooting.tree
conda deactivate

cd ~/5M/dataEdited/binAnalysis/phylogeny/rp16/rooting
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
cp /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16.tree \
    /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree
awk '{ print $1 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/bin_names_to_switch.tsv | while read old_name
do
  new_name=$(awk -v old_name="$old_name" '$1 == old_name { print $2 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/bin_names_to_switch.tsv)
  sed "s/$old_name/$new_name/" /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree \
      > /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree.txt
  mv /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree.txt \
      /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree
done
awk -F ',' '{ print $1 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/binsForAnalysis/newName.csv | while read old_name
do
  new_name=$(awk -F ',' -v old_name="$old_name" '$1 == old_name { print $2 }' /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/binsForAnalysis/newName.csv)
  sed "s/$old_name/$new_name/" /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree \
      > /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree.txt
  mv /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree.txt \
      /Users/benjaminpeterson/Documents/research/5M/dataEdited/metagenomes/binAnalysis/phylogeny/karthiks_db/bin_phylogeny_rp16_edited.tree
done
