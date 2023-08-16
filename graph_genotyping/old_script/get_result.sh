#!/bin/bash

cat graph_genotyper_HGDP/cosigt_smk/results/cosigt_results/HGDP0*/best*.tsv | grep -v "sample" > HGDP_best_geno.tsv

cat graph_geno_separated/graph_genotyper_SGDP/cosigt_smk/results/cosigt_results/LP600*/best*.tsv | grep -v "sample" > SGDP_best_genov2.tsv

cat graph_geno_separated/graph_genotyper_1000G*/cosigt_smk/results/cosigt_results/*/best*.tsv | grep -v "sample" > 1KGP_best_genov2.tsv

for f in *best_geno.tsv; do awk -v var="$f" '{print $0"\t"var}' "$f" | sed -e 's$results/cosigt_results/$$g' -e 's$.x.gaf$$g' -e 's$_best_geno.tsv$$g' -e 's$/$\t$g' >> best_geno_totv2.tmp1.tsv; done


sed "s/haplotype1-/haplotype1_/g" allentoft_best_geno_totv2.tmp1.tsv | sed "s/haplotype2-/haplotype2_/g" | tr '-' '$' > allentoft_best_geno_totv2.tmp2.tsv

cut -f2- allentoft_best_geno_totv2.tmp2.tsv | sed -e 's$ $\t$g' -e 's$_best_genov2.tsv$$g' > allentoft_best_geno_totv2.tsv

paste <(awk '{print $1}' GTEX_best_geno_totv2.tmp1.tsv) <(cut -f3- GTEX_best_geno_totv2.tmp2.tsv) | sed -e 's$ $\t$g' -e 's$_best_genov2.tsv$$g' > GTEX_best_geno_totv2.tsv