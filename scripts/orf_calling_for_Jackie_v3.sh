

#####################################
# Call ORFs
# Note this pipeline is modified from Ben Peterson's because his uses
# metapathways right away, and it uses code from Pame Camejo's directory
# (like WTAF?)
#####################################


conda activate prodigal

cd /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/markergenes


# remove the .fna from the bin names and put them in a file

ls *.fna | sed -e 's/\.fna$//' > bins_list.txt

# Run prodigal (v 2.6.3)

cat bins_list.txt | while read bin
do
prodigal -i $bin.fna -o $bin.coords.gbk -a $bin.faa -p meta
done


#####################################
# Pull out ORFs for HMM searches
#####################################

mkdir ORFs
mkdir bins_processed
scripts=/Users/jacquelinelemaire/Documents/phylo_class/phylo_class/scripts    # Use your directory that has scripts in it

cd /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/analysis/cyano_76/markergenes

cat /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_cyanos_76_fastas/bins_list.txt | while read bin
do
  echo "Moving" $bin
  cp /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_cyanos_76_fastas/$bin.faa ORFs/$bin.faa
  cp /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_cyanos_76_fastas/$bin.fna ORFs/$bin.fna
  cp /Users/jacquelinelemaire/Documents/phylo_class/phylo_class/data/diamond_cyanos_76_fastas/$bin.coords.gbk ORFs/$bin.coords.gbk

done

cd ORFs
cat *.faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv
cd ..

