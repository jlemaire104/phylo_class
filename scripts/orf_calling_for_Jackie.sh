

#####################################
# Call ORFs
# Note this pipeline is modified from Ben Peterson's because his uses
# metapathways right away, and it uses code from Pame Camejo's directory
# (like WTAF?)
#####################################


source /home/GLBRCORG/trina.mcmahon/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''

cd /home/GLBRCORG/trina.mcmahon/Cyanos/data/cyanoMAGs_and_refs

#  Moved the compressed files into a separate directory first

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
scripts=~/Cyanos/code    # Use your directory that has scripts in it

cat bins_list.txt | while read bin
do
  echo "Moving" $bin
  cp /home/GLBRCORG/trina.mcmahon/Cyanos/data/cyanoMAGs_and_refs/$bin.faa ORFs/$bin.faa
  cp /home/GLBRCORG/trina.mcmahon/Cyanos/data/cyanoMAGs_and_refs/$bin.fna ORFs/$bin.fna
  cp /home/GLBRCORG/trina.mcmahon/Cyanos/data/cyanoMAGs_and_refs/$bin.coords.gbk ORFs/$bin.coords.gbk

done

cd ORFs
cat *.faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv
cd ..

