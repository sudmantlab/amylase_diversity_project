# Diversity and recent selection on Amylase structural haplotypes in human populations

A repo for code and methods related to analyses as well as for a manuscript



20230125

- Davide
    - TRIED
        - Davide tried clustering
        - didn't help THAT much
        - evaluation with removing bad samples
        - 76% on the first one diploype prediction accuracy, and 100%  a
            - without any clustering!
        - cluster 90% similar
        - can try 80%?
    - will do
        - show prediction as a function of clustering
            - slide DBscan cutoff / clustering cutoff
        - use DBscan but remove the "[].0/[].1"
        - make a reduced copy number of bundles vector to predict
        - show some examples of what's different
- Alessandro
    - UKBIobank HAS the intensities!
    - Spoke w/ Nicola Perastou (sp?) about multiple testing correction
    - downloaded summary stats for 2-3 phenotypes gluten intake and glucose 
    - could we try looking at ALL phenotypes and 
    - will do
        - look into looking at associateion of 20 SNPs w/ ALL phenotypes (or a subset)
        - look into GWAS of 2-3 phenotypes across that locus
        - Compare intensity of 1kg genotypes over that locus with genotypes (intersect the OMNI with UKBB)

- Andrea look into trying to run VERKKO ?
- Peter make the final set of haplotypes, no broken haps! EXCLUDE THOSE
    - Y1 genomes minus 9 fails
    - Y2 genomes 
- Peter Alma Nicolas figgs


20230208
    - Peter
        - 89 or 90 haplotyupes 
        - how does tree change w/ different bunldes

    - Andrea has verkko done
        - will try to extract amylase to look at it

    - Erik
        - Erik will build the graph for Davide
        - odge distance matrix for Davide
        - discovered problem with extracting alignments out of the graph
        - EKG will try to fix

    - Davide
        - remade the tree
        - tried at different cuts of the tree
        - 14 clusters, 11, 3 and 2 works well, can get up to 80 % w/ 14 groups 
        - 80% with first 14 groups
        - Erik suggests using distance matrix based approach
        - Davide tries w/ DBscan, or aglomerative clustering

    - Alessandro
        - really beautiful genotyping
        - I'll give you the right one
        - alessandro will extend
