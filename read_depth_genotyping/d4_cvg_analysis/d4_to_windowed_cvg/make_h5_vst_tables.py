#!env python3
import pyd4
import numpy as np
import argparse
import pdb
import time
import h5py
import json
import pandas as pd

"""
flat convolution of width n on an array x in O(n) time
padded out to maintain shape
"""

def do_convolve(X, n):
    csum = np.cumsum(X)
    conv = np.r_[np.repeat(csum[n-1],n), csum[n:]-csum[:-n]]
    return conv

def get_windowed_mu(X, w, s):
    mus = do_convolve(X, w)/w
    window_starts = np.arange(0,mus.shape[0],s)
    window_ends = window_starts + s
    return window_starts, window_ends, mus[np.arange(0,mus.shape[0],s)] 


def create(args):
    #contigs = ["chr{c}".format(c=c) for c in range(23)]
    contigs = ["chr1"]

    f = pyd4.D4File(args.d4_file)
    h5_out = h5py.File(args.h5_out, 'w')
    
    for contig in contigs: 

        d4depth = f.load_to_np(contig)
        w_starts, w_ends, windowed_mu = get_windowed_mu(d4depth, args.w, args.s)

        h5_out.create_dataset("{contig}/starts".format(contig=contig),data=w_starts)
        h5_out.create_dataset("{contig}/ends".format(contig=contig),data=w_ends)
        h5_out.create_dataset("{contig}/depth".format(contig=contig), data=windowed_mu)

    h5_out.close()


class h5_genome(object):
    
    def __init__(self,fn):
        self.h5_file = h5py.File(fn,'r')
        self.norm_factor = np.mean(self.h5_file["/chr1/depth"])

    def get_contig_len(self, contig):
        return self.h5_file['{c}/depth'.format(c=contig)].shape[0]
    
    def get_copy(self, contig):
        return 2*self.h5_file['{c}/depth'.format(c=contig)][:] / self.norm_factor
    
    def get_starts(self, contig):
        return self.h5_file['{c}/starts'.format(c=contig)][:]
    
    def get_ends(self, contig):
        return self.h5_file['{c}/ends'.format(c=contig)][:]


def vst(cps1, cps2):
    n1, n2 = cps1.shape[0], cps2.shape[0]
    v1, v2 = np.var(cps1), np.var(cps2)
    vt = np.var(np.r_[cps1, cps2])
    vst = (vt - ((v1*n1 + v2*n2)/(n1+n2))) / vt
    if np.isnan(vst):
        vst=0
    return vst

def calc_vst(args):
    
    #contigs = ["chr{c}".format(c=c) for c in range(23)]
    contigs = ["chr1"]

    p1_gs, p2_gs = [], []
    for f in args.pop1_h5s:
        p1_gs.append(h5_genome(f))
    
    for f in args.pop2_h5s:
        p2_gs.append(h5_genome(f))
    
    df_by_contig = []
    for contig in contigs:
        contig_len = p1_gs[0].get_contig_len(contig)
        p1_array = np.zeros((len(p1_gs),contig_len))
        p2_array = np.zeros((len(p2_gs),contig_len))
        for i, g in enumerate(p1_gs):
            p1_array[i,:] = p1_gs[i].get_copy(contig)
        for i, g in enumerate(p2_gs):
            p2_array[i,:] = p2_gs[i].get_copy(contig)
        
        v1 = np.var(p1_array,0)
        v2 = np.var(p2_array,0)
        vt = np.var(np.r_[p1_array,p2_array],0)
        n1, n2 = len(p1_gs), len(p2_gs)
        vst = (vt - ((v1*n1 + v2*n2)/(n1+n2))) / vt
        vst[np.isnan(vst)] = 0
        
        df = pd.DataFrame({"start":p1_gs[0].get_starts(contig), 
                           "end":p1_gs[0].get_ends(contig),
                           "vst":vst})
        df['contig'] = contig
        df_by_contig.append(df)
    df = pd.concat(df_by_contig)
    df.to_csv(args.fn_out, columns = ["contig", "start", "end", "vst"], sep="\t", index=False) 

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 

    #create h5
    parser_create = subparsers.add_parser("create")
    parser_create.add_argument("--d4_file", required=True)
    parser_create.add_argument("--h5_out", required=True)
    parser_create.add_argument("-w", default=1000, type=int)
    parser_create.add_argument("-s", default=500, type=int)
    parser_create.set_defaults(func=create)

    #calc vst
    parser_create = subparsers.add_parser("calc_vst")
    parser_create.add_argument("--pop1_h5s", required=True, nargs = "+")
    parser_create.add_argument("--pop2_h5s", required=True, nargs = "+")
    parser_create.add_argument("--fn_out", required=True)
    parser_create.set_defaults(func=calc_vst)


    args = parser.parse_args()
    args.func(args)
    
    
