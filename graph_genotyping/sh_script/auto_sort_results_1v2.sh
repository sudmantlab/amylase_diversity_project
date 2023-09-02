#!/bin/bash

# This generate a final table with all the best_geno

# Check if the number of arguments is correct
if [ $# -ne 1 ]; then # If not, print an error message and exit
  echo "Usage: ./auto_sort_results_1.sh genotyper_folder"
  exit 1
fi

# process the last sample
genotyper_folder_tmp=$1
genotyper_folder=$(readlink -f ${genotyper_folder_tmp})
name_out_single=$(echo $genotyper_folder | awk -F "_" '{print $NF}')

if [ ! -d best_geno ]; then
  mkdir best_geno
fi

cd best_geno

cat ${genotyper_folder}/cosigt_smk/results/cosigt_results/*/best*.tsv | grep -v "sample" > ${name_out_single}_best_geno.tsv

if [ -e best_geno_tot.tsv ]; then # If temp.txt exists
  rm best_geno_tot.tmp1.tsv
  rm best_geno_tot.tmp2.tsv
  echo "Deleting the previous tmp files"
fi


# add this to a unique table with the other already processed samples
for f in *_best_geno.tsv; \
 do awk -v var="$f" '{print $0"\t"var}' "$f" | \
 sed -e 's$results/cosigt_results/$$g' -e 's$.x.gaf$$g' -e 's$_best_geno.tsv$$g' -e 's$/$\t$g' >> best_geno_tot.tmp1.tsv; \
done

# substitue the haplotype

 
sed -e "s/haplotype1-/haplotype1_/g" -e "s/haplotype2-/haplotype2_/g" -e 's$:\([0-9]\{1,\}\)-\([0-9]\{1,\}\)$:\1_\2$g' best_geno_tot.tmp1.tsv| tr '-' '$' > best_geno_tot.tmp2.tsv

paste <(awk '{print $1}' best_geno_tot.tmp1.tsv) <(cut -f3- best_geno_tot.tmp2.tsv) | sed -e 's$ $\t$g' -e 's$_best_genov2.tsv$$g' > best_geno_tot.tsv