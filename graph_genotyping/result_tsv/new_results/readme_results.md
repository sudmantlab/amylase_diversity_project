# Results genotyping new graph

graph=`selected_indivs_AMY_region.fa.gz.60ef634.68f91e8.f72a9ab.smooth.final.gfa`

no more bad samples.

- The file [best_geno_tot.tsv](./best_geno/best_geno_tot.tsv) has the following columns: `original_sampleID`, `cosigt_besthaplotypepair (divided by $)`, `cosigt_score`, and `dataset of origin`.

    - `original_sampleID` is the identifier of the sample in the original dataset.
    - `cosigt_besthaplotypepair` is the pair of haplotypes that best match the vector made by short-read mapping of the specific sample.
    - `cosigt_score` is the score assigned by cosigt to the best haplotype pair, ranging from 0 to 1.
    - `dataset` is the name of the dataset that the sample belongs to.


- [Folder](./combos) of combos, it contains the folder for each dataset and in each folder you have:

    - The sorted combos for each sample as `<SAMPLE>combos_sorted.tsv`.
    - The top 10 combos for each sample as `<SAMPLE>combos_sorted_top10.tsv`.
    - The symlink with `<SAMPLE>` is as the `<SAMPLE>combos_sorted.tsv` but not sorted.


Numbers of genotyped samples:

| N of genotyped samples [MODERN] | Dataset            | Number in the source folder  |
|----------|----------------------|---| 
| 2406     | 1K                 |  |
| 698     | 1KGrel       |  698  |
| 838     | GTEX       |  838  |
| 828     | HGDP       |  828 |
| 277     | SGDPv2*       | 279 |


## Samples to check

### SGDP
* LP6005441_DNA_C06.bam does not exist
* LP6005441_DNA_H08.bam cannot create the index  

### 1KGP
* HG02635
* HG03366
* HG03025


| N of genotyped samples [ANCIENT] | Dataset            | Number in the source folder  |
|----------|----------------------|---| 
| 690     | Allentoft                 | 1668 |
| 1673 | Allentofthg19       | 1673   |
| 15     | Excoffier       |  15  |
| 13     | Haak       |  13 |
| 4     | NeanArcaic       | 4 |
| 116     | PosthKrauseHG       | 116 |
| 784     | Reich       | 784 |
| 216     | aDNAhcovtesthg19 (ReichAncients_agdp_subset) | 216 |

