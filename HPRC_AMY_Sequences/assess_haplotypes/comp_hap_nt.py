
import pandas as pd
import pdb
import itertools
import numpy as np
from fasta_reader import read_fasta
from Bio.Seq import Seq
from Bio import pairwise2




if __name__=="__main__":
    seqs_by_id = {}
    for item in read_fasta("../AMY1A_region_seq.fa.gz"):
        id = item.defline
        seq = item.sequence
        seqs_by_id[id] = seq

    for i1 in seqs_by_id.keys():
        for i2 in seqs_by_id.keys():
            s1 = Seq(seqs_by_id[i1])
            s2 = Seq(seqs_by_id[i2])
            if s1==s2: continue
            pdb.set_trace()
            alignments = pairwise2.align.globalxx(s1, s2)
            pdb.set_trace()
