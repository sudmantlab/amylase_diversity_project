#PGR-TK commit 994934bda27da205d1686a0a14a0f10661d7771f, v0.4.0dev brach build
pgr-pbundle-decomp -w 48 -k 56 -r 4 --min-span 28 --min-cov 0 --min-branch-size 8 --bundle-length-cutoff 100 --bundle-merge-distance 10000 AMY1A_region_seq_Y1_Y2_TEST.fa AMY1A_region_seq_Y1_Y2_TEST.
pgr-pbundle-bed2sorted AMY1A_region_seq_Y1_Y2_TEST.bed AMY1A_region_seq_Y1_Y2_TEST.
cat AMY1A_region_seq_Y1_Y2_TEST.ord | awk '{print $1"\t"$1"-"$2}' > AMY1A_region_seq_Y1_Y2_TEST.ord2
pgr-pbundle-bed2svg --track-range 800000 --stroke-width 1.5 AMY1A_region_seq_Y1_Y2_TEST.bed  AMY1A_region_seq_Y1_Y2_TEST. --annotations AMY1A_region_seq_Y1_Y2_TEST.ord2
