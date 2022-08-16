## Alignment of short reads using as reference pangenome

- Downloaded three samples from 1KGP with variable lengths in AMY region [HG00438, HG01978, HG02257]
- isolate unaligned reads from .cram files in single .fq (many unaligned were singletons) and map on HPRC multi fasta on AMY region + [chm13v2](https://github.com/marbl/CHM13) - mapper used [mm2-fast](https://github.com/bwa-mem2/mm2-fast)
- we noticed that only ~ 40 unaligned reads were mapped on the AMY region of the HPRC multi-fasta [30, 36 and 40 reads saved] not much improvement.
- we thought, to be consistent to use the same mapper for everything (then we would have re-mapped the unaligned to HPRC multi-fasta only):
    - we converted the whole .cram files into two .fq and we remapped to chm13 and GRCH38 with mm2-fast:  differences were observed in unaligned reads between the two references (higher number of unaligned on chm13 than GRCH38 - however half of the UL reads of chm13 are in GRCh38)
    - we used then bwa-mem2 and we detected the same differences. **OT**: understand where 1M of unmapped reads map on GRCH38:
        - isolated the reads that map on decoy/HLA but only ~ 40K overlap
        - remapping the UL of chm13 on GRCh38 and see where those reads are: [duplicate regions](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1422896709_6Ewf7lC3evtdHGvjsKxZwHBgVEsG&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=fixSeqLiftOverPsl&hgta_table=0&hgta_regionType=genome&position=chr2%3A32%2C915%2C973-32%2C917%2C142&hgta_outputType=primaryTable&hgta_outFileName=), [centromere](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1422896709_6Ewf7lC3evtdHGvjsKxZwHBgVEsG&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=centromeres&hgta_table=0&hgta_regionType=genome&position=chr2%3A32%2C915%2C973-32%2C917%2C142&hgta_outputType=primaryTable&hgta_outFileName=)?)
        - **SOLUTION**: most of these reads goes in chrEBV which is not present on chm13. This is because the LCCs of 1000G are generated with EBV.
- The two references have similar numbers of reads in AMY region. For HG00438 for example we have about 55,044 reads chm13 and 55,041 using GRCh38 aligned regardless of the aligner used. Even if the lengths of the sequences after the liftover are different - chm13 = 438K while GRCh38 = 250K

- Mapping the whole reads on HPRC AMY region + chm13 and extracting how many maps on the haplotypes and on the region of interest give these results: 54,183 and ~ 991 still mapping on chm13 chr1 even if the haplotype of chm13 is present on the multi-fasta extracted by Davide/Erik - however it has a different dimensioin

