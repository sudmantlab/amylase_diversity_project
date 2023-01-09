#PGR-TK commit 75fa20b41592941c9e6eef3f914d97788ee06b86, v0.4.0dev brach build
pgr-pbundle-decomp -w 48 -k 56 -r 4 --min-span 28 --min-cov 0 --min-branch-size 8 --bundle-length-cutoff 100 --bundle-merge-distance 10000 AMY1A_region_seq_Y1_Y2_TEST.fa.gz AMY1A_region_seq_Y1_Y2_TEST.
pgr-pbundle-bed2sorted AMY1A_region_seq_Y1_Y2_TEST.bed AMY1A_region_seq_Y1_Y2_TEST.
cat AMY1A_region_seq_Y1_Y2_TEST.ord | awk '{print $1"\t"$1"-"$2}' > AMY1A_region_seq_Y1_Y2_TEST.ord2
pgr-pbundle-bed2svg --track-range 800000 --stroke-width 1.5 AMY1A_region_seq_Y1_Y2_TEST.bed  AMY1A_region_seq_Y1_Y2_TEST. --annotations AMY1A_region_seq_Y1_Y2_TEST.ord2
pgr-pbundle-bed2dist AMY1A_region_seq_Y1_Y2_TEST.bed AMY1A_region_seq_Y1_Y2_TEST.
pgr-pbundle-bed2svg --track-range 700000 --stroke-width 1.5 AMY1A_region_seq_Y1_Y2_TEST.bed  AMY1A_region_seq_Y1_Y2_TEST. --annotations AMY1A_region_seq_Y1_Y2_TEST.ord2 --ddg-file AMY1A_region_seq_Y1_Y2_TEST.ddg
