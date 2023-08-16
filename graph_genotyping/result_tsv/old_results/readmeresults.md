# Results cosigt pipeline

This folder contains all the results obtained by running the genotyper on the graph.

The `file cosigt_results.bestgeno.tsv` has four columns that represent: `original_sampleID`, `cosigt_besthaplotypepair`, `cosigt_score`, and `dataset`.

- `original_sampleID` is the identifier of the sample in the original dataset.
- `cosigt_besthaplotypepair` is the pair of haplotypes that best match the vector made by short-read mapping of the specific sample.
- `cosigt_score` is the score assigned by cosigt to the best haplotype pair, ranging from 0 to 1.
- `dataset` is the name of the dataset that the sample belongs to.


This is the table of the datasets:

| N | Dataset            |   |   |   |
|----------|----------------------|---|---|---|
| 3177     | 1KGP                 |   |   |   |
| 1670     | Allentoft_aDNA       |   |   |   |
| 838      | GTEX                 |   |   |   |
| 828      | HGDP                 |   |   |   |
| 13       | Haak_aDNA            |   |   |   |
| 116      | KrausePosth_aDNA     |   |   |   |
| 15       | MarchiExcoffier_aDNA |   |   |   |
| 785      | Reich_aDNA           |   |   |   |
| 176      | SGDP                 |   |   |   |
| 277      | SGDPv2               |   |   |   |
| 216      | highcov              |   |   |   |
| 4        | nea                  |   |   |   |

N.B: There were issues on the analysis of aDNAs - likely because of the coverage - Therefore we re-run the pipeline on a subset of "good" aDNAs which are marked as `highvoc` and `nea`on the table above and in the file. 