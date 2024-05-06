1.	Abstract: 200 words or less
2.	Introduction: convey the background, objective and significance of your analysis
3.	Materials and Methods: describe in details the data and the different methods, software, assumptions or any methodological choices that you did; appropriately cite other papers or software used
4.	Results: clearly and effectively describe the results of your analysis
5.	Discussion: summarize the conclusion of your findings and highlight any weaknesses or limitations of your approach
6.	References
7.	Reproducible script: step-by-step description of the methodology that can be easily followed by others. This can be the link to a github repository or an actual file
Title 
Aphanizomenon Phylogenetic Study and Species Community Composition
study system, variables, result & direction
Phylogenetics Class PLPATH 563
May 2024
Jackie Lemaire
Abstract
 biological rationale, hypothesis, approach, result direction & conclusions
Climate change, agricultural runoff, and species invasions have caused multiple disruptive ecological shocks in Lake Mendota, a freshwater eutrophic lake. The lake suffers from noxious cyanobacterial blooms that produce toxins. Using cyanobacterial genomes recovered from a metagenomic time series spanning 20 years of sampling in Lake Mendota, we have the unique opportunity to investigate the natural community composition. A total of 75 curated cyanobacterial genomes were assembled and used to calculate the coverage of each genome in the 464 metagenome samples. We found that an Aphanizomenon genome had the highest coverage across the dataset by far, and the second highest mean coverage. Microcystis populations had markedly higher microdiversity across the genome, as observed in a shorter time series. We also examined the population-level diversity of the cyanobacteria through time to ask whether new strains have emerged over the 20 years. Analyzing the cyanobacteria phylogenetics throughout this period will help us understand the community composition and relatedness among species and strains. A sample subset of Aphanizomenon species will be used to develop methods to account for species correlations and establish the best pipeline to run on the larger time-series dataset.  

Introduction
relevant & correctly cited background information, question, biological rationale (including biological assumptions), hypothesis (if any), approach
Cyanobacteria are ubiquitous in the environment and an important part of the global ecosystem. They are drivers of biogeochemical cycling (1) as they perform many important functions such as producing oxygen, fixing nitrogen, and exporting organic carbon to the deep ocean. However, some species can form harmful blooms and release toxins into the environment. Cyanotoxins are the potent toxins that are produced by cyanobacteria that include neurotoxins, endotoxins, and hepatotoxins (2). They can cause wide-ranging illness in humans and even death in wild and domestic animals when ingested (3).
Lake Mendota in Madison, Wisconsin is often referred to as one of the most studied lakes in the world (4). Human settlement, urbanization, and farming in the surrounding area of Madison in the mid to late 1800s caused eutrophication in the lake. Scientist E. A. Birge reported Lake Mendota’s first toxic blooms of cyanobacteria, also known as blue-green algae, as early as 1882 (5). Climate and land-use change are two interacting drivers that are altering the microbial communities in Lake Mendota, including the cyanobacteria (8). There is an urgent need to forecast how such drivers will further degrade, or possibly increase, water quality.
I propose using our lab’s unique two-decade long time series metagenomics dataset of the water column in Lake Mendota to analyze the microbial populations and their dynamics over time. Previous researchers in the lab have looked at bacterial community dynamics using 16S rRNA gene amplicon data from the same 472 samples and discovered surprising changes in community phenology and cyanotoxin production (6). The McMahon Lab has the proper infrastructure and expertise in the metagenomics and bioinformatics field of environmental microbiology to fully leverage these large datasets, including access to extensive computational resources. Preliminary results show that a handful of taxa found in lakes all over the world are abundant and persistent in the dataset, and that they have distinct patterns of population-level diversity (i.e. magnitude and distribution of single nucleotide variants). Cyanobacteria are an ancient and diverse phylum of bacteria, and the species community structure shifts seasonally and annually. However, preliminary data shows they have among the lowest levels of nucleotide diversity when they bloom. Prior work in the McMahon lab showed that a diverse array of cyanotoxins are produced throughout the summer and I aim to relate specific populations to toxin production (6). Key bacterial lineages have been associated with cyanobacterial blooms and I will also examine how community- and population-level diversity change over multiple time scales (e.g. days to years). 
This work aims to develop a pipeline to access the community structure and relatedness of species in a natural freshwater system by establishing methods for multiple sequence alignments and evaluating more modern tools for building phylogenetic networks. Water quality is a major concern in Lake Mendota, but its status as a living laboratory presents an opportunity to also learn about fundamental principles driving change at the community level (species dynamics) and population level (evolutionary dynamics). My work will contribute to a broader understanding of how genomes change over time in naturally assembling communities. 

Materials and Methods
Concisely, clearly, & chronologically describes procedure used so that knowledgeable reader could replicate experiment and understand the results. Methods used are appropriate for study. Briefly describes mathematical manipulations or statistical analyses
Multiple Sequence Alignments
A sample sequence dataset of fasta files was compiled by downloading fasta sequences from NCBI GenBank that were the top hits for an Aphanizomenon [Cyanobacteria] 16S gene and an outgroup of Synnecchococcus [Cyanobacteria] 16S gene nucleotide sequence. The sequences consisted of 12 samples with 3 species of Aphanizomenon and one species of Synnechococcus. The sequences were aligned using both clustalW (Clustal 2.1) and Muscle alignment (Muscle 5.1 osx64) software to compare methods. ClustalW is an alignment algorithm that includes sequence weighting, position-specific gap penalties, and can be used with multiple types of input and output file types including fasta files which were used here (Chenna). MUSCLE, is a newer program with high accuracy and high-throughput. The Muscle algorithm increases speed for larger datasets by using a fast distance estimation using kmer counting, makes a progressive alignment using a new profile function referred to as the log-expectation score, and then does refinement using tree-dependent restricted partitioning to output the most accurate alignment (Edgar).

Distance-based Estimation Neighbor Joining Trees 
Each multiple sequence alignment was initially analyzed by making a simple distance-based estimation and neighbor joining tree with the Software R package ape. It is one of the most widely used phylogenetic software packages in R and has a huge variety of functions. To compute the genetic distance between these bacterial species, the Galtier and Gouy 1995 (GG95) model was first used which allows for the G+C content to change through time as well as different rates for transitons and transversions (since cyanobacteria can have variable GC percentage. However, the TN93 model plus gamma (Tamura and Nei 1993) was chosen to move forward with; which allows for different rates of transitions and transversions, heterogeneous base frequencies, and between-site variation of the substitution rate as well as allows for a gamma correction of the inter-site variation in substitution rates. The tree was rooted to the outgroup Synnecchococcus.
Parsimony Trees
Each multiple sequence alignment was next analyzed by making a parsimony-based estimation tree with the Software R package phangorn. This method first constructs a starting tree by using same methods as above to create a neighbor-joining tree using a distance-based estimation with the TN93 model of evolution with gamma correction. The parsimony score is then computed of the starting trees: 227 p-score for the muscle alignment and 225 p-score for the clustal alignment. Then the optim.parsimony function was used to search for the tree with maximum parsimony; the fitch algorithm was used with sub-tree pruning and regrafting for rearrangements. The trees were again rooted to the outgroup Synnecchococcus before plotting.
Maximum Likelihood Trees
The muscle multiple sequence alignment was chosen to move forward with to make maximum likelihood trees with 2 popular programs to compare methods. The first program chosen was the new RAxML-NG v. 1.2.1 (2019) which is a fast and scalable tool for maximum likelihood phylogenetic inference. This method was chosen as it has proven to be faster than the original RAxML algorithim. It shows improved accuracy, flexibility, speed, scalability, and usability compared with other methods (Kozlov) – although this sample dataset is very small and can be easily handled with any of these methods – determining the best method for use when scaling to the larger dataset will prove useful. RAxML-NG has additional features such as the detection of terraces in tree space and the transfer bootstrap support metric (Kozlov). It also now has additional metrics for 22 ‘classical’ GTR-derived models, the GTR+gamma model was chosen as input here, but the program indicted that the GTR+FO+G4m was used once it started running. More running trees were started with (25 parsimony and 25 random) to explore the tree space more and improve convergence. The –all function was called to run the ML search and do bootstrapping at the same time with the boostrap metrics of transfer bootstrap (tbe) and (fbp) to map the support values onto the best-scoring ML tree. TBE appears to be better for very large trees, so for this small dataset the fbp support was used to evaluate branch support as a metric of confidence of the lineage delineation.
  
. 
(Kozlov)
Bayesian Inference Trees


 
1.	Align the sequences
2.	
a.	Compare methods clustalW and Muscle
3.	Make trees to compare and compare tree methods
a.	First do a quick distance method and neighbor joining tree
b.	Then make parsimony trees
4.	Next make better maximum likelihood trees
a.	Compare methods or pick one?
b.	raxML vs IQtree
5.	Make Bayesian Trees
a.	Mr.Bayes and/or RevBayes
6.	If we want to compare multiple gene trees and combine using coalescent methods to make a species tree

Results
The phylogenies are visualized by plotting trees for phylogenetic analysis and methods comparison via multiple different algorithms. The two multiple sequence alignments 

(MSAs) were first analyzed by plotting simple distance-based estimation neighbor joining trees with the muscle alignment (Fig. 1) and the clustalW alignment (Fig. 2). The different alignments showed similar results, with slightly more resolution of divergent clades in the Aphanizomenon aphanizomenoides species lineage as it is split between two clades in the muscle alignment versus creating a monophyletic clade in the clustalW alignment NJ tree. 
The MSAs were then plotted using parsimony-based estimation. This is a more efficient and robust tree-building algorithim which starts with a distance-based NJ tree and then uses a progressive tree searching method with subtree pruning and regrafting (used here) to find the tree with maximum parsimony. This method minimizes the number of genetic transformations in finding the best tree (Goëffon). The muscle alignment parsimony tree (Fig. 3) and the Clustal alignment parsimony tree (Fig. 4) agree with each other and with the NJ tree from the muscle alignment (Fig. 1). although with better resolution. 


The A. issatschenkoi consistently groups with the A. gracile clade in all NJ and parsimony trees. The parsimony trees show a difference among the A. gracile species where 3 of the strains (FJ895124, FJ895124125, and FJ895124126) are shown to be identical replicates, forming their own monophyletic branch that is distinct from the FJ895127 and FJ895128 stains although grouped in the same overarching species clade. There is better separation in minor differences among the A. aphanizomenoides species as well bifurcating into two major clades with subgroupings. 

The methods thus far worked well especially for this small dataset with very closely related species and for a quick view of the results. However, to run a more statistically sound method that is more scalable for running this analysis on a larger dataset in the future and on more divergent species; maximum likelihood estimations are used. This method is better than parsimony because it can perform estimations using a precise (user chosen) evolutionary model that can be adjusted depending on the dataset being used. 




































With a few minor exceptions, contains a concise, well-organized narrative text & tables/figures that highlight key trends/patterns/output. Tables & figures have appropriate legends/ abels & can stand on their own











Discussion

With a few minor exceptions, clearly, concisely, & logically presents all key components: interprets/integrates data; formulates argument for conclusions referring back to biological rationale & by comparing with relevant findings in literature, introduces new literature to discuss or support findings, evaluates methods, evaluates reliability of data, states knowledge generated & implications of results, suggests next investigation steps, includes unique observations, and ends paper with final conclusions
The analyses here combine sequences of 12 cyanobacterial sequences consisting of 3 species of Aphanizomenon and one species of Synnechococcus. Methods were evaluated on this small sample set to optimize and test modern algorithms and programs for phylogenetic analysis of bacterial species composition. 
References
References within body of paper are cited appropriately; references in final citation list are formatted appropriately
References
1.	Sánchez-Baracaldo (2022). Trends in Microbiology, 30(2), 143-157.
2.	Codd, G. A. (2005). Toxicology and applied pharmacology, 203(3), 264-272.
3.	Bell, S. G. (1994). Reviews in Medical Microbiology, 5(4), 256-264.
4.	Brock, T. D. (2012). Springer Science & Business Media.
5.	Van Eyck, M. (2016). College of Letters and Sciences. http://ls. wisc. edu/news/lake-mendota-a-scientific-biography.
6.	Rohwer, R. (2022). Lake iTag measurements over nineteen years, introducing the limony dataset. bioRxiv.
7.	Alvarenga, D. O. (2017). Frontiers in microbiology, 8, 809.
8.	Rohwer, R. R.. (2022). Species invasions shift microbial phenology in a two-decade freshwater time series. bioRxiv.
9.	
Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic acids research, 32(5), 1792-1797.
Chenna, R., Sugawara, H., Koike, T., Lopez, R., Gibson, T. J., Higgins, D. G., & Thompson, J. D. (2003). Multiple sequence alignment with the Clustal series of programs. Nucleic acids research, 31(13), 3497-3500.
Goëffon, A., Richer, J. M., & Hao, J. K. (2008). Progressive tree neighborhood applied to the maximum parsimony problem. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 5(1), 136-145.
Schliep, K. P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics, 27(4), 592-593.
Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), 4453-4455.
Reproducible Script
Fully explained sequence of commands interleaved with comments which make the whole analysis easy to follow and reproduce; details on software installed and versions installed as well as necessary format for the data input files (or links to input data files if data is publicly available)
