# Pangenome Year 1 + Year 2 with `selected_indivs_AMY_region.fa.gz`

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/
PGGB=/home/guarracino/tools/pggb/pggb-736c50d8e32455cc25db19d119141903f2613a63
```

Get the repository:

```shell
cd $DIR_BASE
git clone --recursive https://github.com/sudmantlab/amylase_diversity_project.git
```

Build a pangenome graph:

```shell
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/

$PGGB \
  -i $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/output/pggb_input/selected_indivs_AMY_region.fa.gz \
  -o amy_y1y2.selected_indivs.1 \
  -n 94 -c 2 -t 48 \
  -D /scratch
```

Extract VNTR plus a bit of flanking regions:

```shell
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/amy_y1y2.selected_indivs.1

PREFIX=selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.0ead406.smooth.final
odgi sort -i $PREFIX.og -Y -o $PREFIX.sort.og -P -t 16 --temp-dir /scratch
odgi paths -i $PREFIX.sort.og -Ll > $PREFIX.sort.paths.tsv
bedtools makewindows -g <(cut -f 1,3 $PREFIX.sort.paths.tsv | grep chm13) -w 500 > $PREFIX.sort.paths.chm13.windows.bed
odgi depth -i $PREFIX.sort.og -b $PREFIX.sort.paths.chm13.windows.bed | awk '$4 > 98' | bedtools merge -d 100 | awk '$3-$2 > 500' > $PREFIX.sort.paths.chm13.windows.to_extract.bed 
odgi extract -i $PREFIX.sort.og -b $PREFIX.sort.paths.chm13.windows.to_extract.bed -O -t 16 -P -o - | odgi sort -i - -p gYs -o $PREFIX.sort.extract.og -P -t 16 --temp-dir /scratch
odgi view -i $PREFIX.sort.extract.og -g > $PREFIX.sort.extract.gfa

odgi viz -i $PREFIX.sort.extract.og -o $PREFIX.sort.extract.1D.png -m # To verify that it worked correctly
```

## Leave-One-Out

```shell
mkdir -p $DIR_BASE/amylase_diversity_project/pangenome/pggb/leave-one-out
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/leave-one-out

FASTA=$DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/output/pggb_input/selected_indivs_AMY_region.fa.gz

SEQ_NAMES=($(cut -f1 "${FASTA}.fai"))

for SEQ_TO_REMOVE in "${SEQ_NAMES[@]}"; do
  echo $SEQ_TO_REMOVE

  # Create a new array without the current sequence.
  SEQS=("${SEQ_NAMES[@]/$SEQ_TO_REMOVE}")

  NAME=$(echo $SEQ_TO_REMOVE | tr '#\-:' '_')

  # Use samtools faidx to fetch the remaining sequences and write them to a new FASTA file.
  samtools faidx $FASTA $(echo $SEQS) | bgzip -@ 48 -l 9 > selected_indivs_AMY_region.${NAME}.fa.gz
  samtools faidx selected_indivs_AMY_region.${NAME}.fa.gz
done

# Be careful with the graph made with "selected_indivs_AMY_region.HG02004.verkko_filt_pat.fa_haplotype2_0000813_139386576_140076845_1.fa.gz", because "HG02004.verkko_filt_pat.fa:haplotype2-0000813_139386576_140076845_1" (removed in such a set) is the only H9 haplotype

OUT_DIR=$DIR_BASE/amylase_diversity_project/pangenome/pggb/leave-one-out
ls *fa.gz | while read FASTA; do
  echo $FASTA;

  NAME=$(basename $FASTA .fa.gz)
  
  sbatch -c 48 -p workers --job-name $NAME --wrap "$PGGB -i $OUT_DIR/$FASTA -o $OUT_DIR/pggb.$NAME -n 94 -c 2 -t 48 -D /scratch"
done
```

Layout without labels:

```shell
f=selected_indivs_AMY_region.fa.gz.42c7330.417fcdf.8bc4b72.smooth.final.gfa ; odgi draw -i $f -t 16 -c $f.lay -s $f.lay.svg -b <(sed 's/bundle[0-9]//g' haplotype_all_structures_plotinfo.chm13_colors.bed) -w 30 -B 100
```


# Pangenome Year 1 + Year 2 with `AMY1A_region_seq.fa.gz`

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/
PGGB=/gnu/store/56g6kc0s7gbqflgmk69hw4xp67nxl8ff-pggb-0.4.0+f0db7db-2/bin/pggb
```

Get the repository:

```shell
cd $DIR_BASE
git clone --recursive https://github.com/sudmantlab/amylase_diversity_project.git
```

Build a pangenome graph:

```shell
cd $DIR_BASE/amylase_diversity_project/pangenome/pggb/

$PGGB \
  -i $DIR_BASE/amylase_diversity_project/HPRC_AMY_Sequences/combined_y1_y2_analyses/input/y1_y2_HPRC/AMY1A_region_seq.fa.gz \
  -o amy_y1y2.1 \
  -p 90 -s 5k -k 419 -t 48 -n 200 \
  -D /scratch -Z
```
