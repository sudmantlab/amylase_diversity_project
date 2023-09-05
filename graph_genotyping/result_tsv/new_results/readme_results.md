# Results genotyping new graph

graph=`selected_indivs_AMY_region.fa.gz.60ef634.68f91e8.f72a9ab.smooth.final.gfa`

no more bad samples.

- The file [best_geno_tot.tsv](./best_geno/best_geno_tot.tsv) has the following columns: `original_sampleID`, `cosigt_besthaplotypepair (divided by $)`, `cosigt_score`, and `dataset of origin`.

    - `original_sampleID` is the identifier of the sample in the original dataset.
    - `cosigt_besthaplotypepair` is the pair of haplotypes that best match the vector made by short-read mapping of the specific sample.
    - `cosigt_score` is the score assigned by cosigt to the best haplotype pair, ranging from 0 to 1.
    - `dataset` is the name of the dataset that the sample belongs to.


- [Folder](./combos) of combos, in each folder you have:

    - The sorted combos for each sample as `<SAMPLE>combos_sorted.tsv`.
    - The top 10 combos for each sample as `<SAMPLE>combos_sorted_top10.tsv`.
    - The symlink with `<SAMPLE>` is as the `<SAMPLE>combos_sorted.tsv` but not sorted.


Numbers of genotyped samples:

| N of genotyped samples [MODERN] | Dataset            | Number in the source folder  |
|----------|----------------------|---| 
| 2504     | 1K                 | 2504 |
| 698     | 1KGrel       |  698  |
| 838     | GTEX       |  838  |
| 828     | HGDP       |  828 |
| 277     | SGDPv2*       | 279 |


## Samples to check

### SGDP
* LP6005441_DNA_C06.bam does not exist
* LP6005441_DNA_H08.bam cannot create the index  

### 1K

There are still some files that are corrupted or missing in `/global/scratch/p2p3/pl1_sudmant/alessandroraveane/re_dwn_1kgp/data/`.

The symlink used for all the files are in `/global/scratch/users/alessandroraveane/input_cram_bam/cram_seq_1000G/all_of_them`.


| N of genotyped samples [ANCIENT] | Dataset            | Number in the source folder  |
|----------|----------------------|---| 
| 1641     | Allentoft  [`stoneage/hg38_bams`]*               | 1668 |
| 1673 | Allentofthg19 [`stoneage/hg19_bams`]      | 1673   |
| 15     | Excoffier       |  15  |
| 13     | Haak       |  13 |
| 4     | NeanArcaic       | 4 |
| 116     | PosthKrauseHG       | 116 |
| 784     | Reich       | 784 |
| 216     | aDNAhcovtesthg19 [`ReichAncients_agdp_subset`] | 216 |

## Samples to check 

### Allentoft [`stoneage/hg38_bams`]

* 27 bams to check `/global/scratch/users/alessandroraveane/bam_seq_ancient/hg38_bams/check_corr/bam_with_issue.tsv`

* RISE509 does not have index (created)