## Alignment of short reads using as reference pangenome

- Downloaded three samples from 1KGP with variable lengths in AMY region [HG00438, HG01978, HG02257]
- isolate unaligned reads from .cram files in single .fq (many unaligned were singletons) and map on HPRC multi fasta on AMY region + [chm13v2](https://github.com/marbl/CHM13) - mapper used [mm2-fast](https://github.com/bwa-mem2/mm2-fast)
- we noticed that only ~ 40 unaligned reads were mapped on the AMY region of the HPRC multi-fasta [30, 36 and 40 reads saved]
- we thought, to be consistent to use the same mapper for everything (then we would have re-mapped the unaligned to HPRC multi-fasta only):
    - we converted the whole .cram files into two .fq and we remapped to chm13 and GRCH38 with mm2-fast -> observing differences in unaligned reads (higher number of unaligned on chm13 than GRCH38 - most of the un. reads of GRCH38 are in chm13)
    - we used then bwa-mem2 and we detected the same differences (likely linked to decoy)
    - isolated the reads that map on decoy/HLA but only ~ 40K overlap
    - OT: understand where 1M of unmapped reads map on GRCH38.. (remapping the UL of chm13 on GRCh38 and see where those reads are: [duplicate regions](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1422896709_6Ewf7lC3evtdHGvjsKxZwHBgVEsG&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=fixSeqLiftOverPsl&hgta_table=0&hgta_regionType=genome&position=chr2%3A32%2C915%2C973-32%2C917%2C142&hgta_outputType=primaryTable&hgta_outFileName=), [centromere](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1422896709_6Ewf7lC3evtdHGvjsKxZwHBgVEsG&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=centromeres&hgta_table=0&hgta_regionType=genome&position=chr2%3A32%2C915%2C973-32%2C917%2C142&hgta_outputType=primaryTable&hgta_outFileName=) )
- to do: map all reads on chm13/GRCh38 + haplotype multi-fasta and count how many reads are either on the region of AMY in chm13/GRCh38 or in our haplotypes?  