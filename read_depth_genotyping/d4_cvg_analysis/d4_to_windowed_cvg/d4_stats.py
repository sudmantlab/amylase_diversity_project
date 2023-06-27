#!env python3
import pyd4
import numpy as np
import argparse
import pdb
import time
import h5py
import json
import pandas as pd

#import pyarrow as pa
#import pyarrow.csv as csv
import csv
import time as time
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
    window_ends = window_starts + w
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


def get_loci_by_contig(fn_loci):

    loci = json.load(open(fn_loci))
    loci_by_contig = {}

    for locus, inf in loci['loci'].items():
        contig = inf['contig']
        if not contig in loci_by_contig:
            loci_by_contig[contig] = {}
        loci_by_contig[contig][locus] = inf
    
    return loci_by_contig


def genotype(args):
    stime = time.time()
    f = pyd4.D4File(args.d4_file)
    
    t_inf = pd.read_csv(args.fn_stats,sep="\t", header=0)
    alt_corr_factor = t_inf.qmax_filt_mu[t_inf['sample']==args.sample].values[0]
    corr_factor = t_inf.CTRL_region_mean[t_inf['sample']==args.sample].values[0]

    loci_by_contig = get_loci_by_contig(args.fn_loci)
    outrows = []
    for contig, loci in loci_by_contig.items():

        d4depth = f.load_to_np(contig)
        w_starts, w_ends, windowed_mu = get_windowed_mu(d4depth, args.w, args.s)
        windowed_cp = 2*windowed_mu/corr_factor
        windowed_alt_cp = 2*windowed_mu/alt_corr_factor

        for locus, locus_inf in loci.items(): 
            start = locus_inf['start']
            end = locus_inf['end']
            locs = np.where((w_starts>=start)&(w_ends<=end))
            gt = np.mean(windowed_cp[locs])
            gt_alt = np.mean(windowed_alt_cp[locs])
            raw = np.mean(windowed_mu[locs])
            #gt_alt = np.mean(windowed_mu[locs]/corr_factor*alt_corr_factor)
            outrows.append({"locus":locus,
                            "contig":contig,
                            "start":start,
                            "end":end,
                            "cp":gt,
                            "cp_alt":gt_alt,
                            "raw":raw,
                            "sample":args.sample})
            
    t = pd.DataFrame(outrows)
    
    t.to_csv(args.fn_out,
             sep="\t", 
             index=False,
             columns=["locus","contig","start","end","cp","cp_alt","raw","sample"], 
             header=False,
             quoting=csv.QUOTE_NONE)
    print(time.time()-stime)



def coverage(args):
    stime = time.time()
    f = pyd4.D4File(args.d4_file)
    
    t_inf = pd.read_csv(args.fn_stats,sep="\t", header=0)
    corr_factor = t_inf.qmax_filt_mu[t_inf['sample']==args.sample].values[0]
    alt_corr_factor = t_inf.CTRL_region_mean[t_inf['sample']==args.sample].values[0]
    d4depth = f.load_to_np(args.contig)
    w_starts, w_ends, windowed_mu = get_windowed_mu(d4depth, args.w, args.s)
    
    locs = np.where((w_starts>=args.start)&(w_ends<=args.end))
    w_starts = w_starts[locs]
    w_ends = w_ends[locs]
    windowed_mu = windowed_mu[locs]
    windowed_cp = 2*windowed_mu/corr_factor
    windowed_cp_alt = 2*windowed_mu/alt_corr_factor
    
    t = pd.DataFrame({"contig":args.contig,
                      "start":w_starts,
                      "end":w_ends,
                      "cp":windowed_cp,
                      "cp_alt":windowed_cp,
                      "raw":windowed_mu,
                      "sample":args.sample})
    
    #t_pa = pa.Table.from_pandas(t)

    #with pa.CompressedOutputStream(args.fn_out, 'gzip') as out:
    #    opts = csv.WriteOptions(include_header=False)#,delimiter="\t")
    #    csv.write_csv(t_pa, out, opts) 
    #csv.write_csv(t_pa, args.fn_out, csv.WriteOptions(delimiter="\t",include_header=False))
    #csv.write_csv(t_pa, args.fn_out)
    t.to_csv(args.fn_out,
             sep="\t", 
             index=False,
             columns=["contig","start","end","cp","cp_alt","raw","sample"], 
             header=False,
             compression="gzip",
             quoting=csv.QUOTE_NONE)
    print(time.time()-stime)





def stats(args):
    #contigs = ["chr{c}".format(c=c) for c in range(23)]
    contigs = ["chr1"]
    #only rnus on chr 1 for now
    contig = contigs[0]
    f = pyd4.D4File(args.d4_file)
    

    d4depth = f.load_to_np(contig)
    w_starts, w_ends, windowed_mu = get_windowed_mu(d4depth, args.w, args.s)
    quantiles = np.quantile(windowed_mu,q=[0,0.025,0.975,1])
    #quantiles = np.quantile(windowed_mu,q=[0,0.1,0.9,1])
    #quantiles = np.quantile(windowed_mu,q=[0,0.05,0.95,1])
    qmax = quantiles[2]
    qmin = quantiles[1]
    
    qmax_filt_mu = np.mean(windowed_mu[windowed_mu<=qmax])
    qmax_filt_med = np.median(windowed_mu[windowed_mu<=qmax])
    qmax_filt_sd = np.sqrt(np.var(windowed_mu[windowed_mu<=qmax]))
    
    #qrng_filt_mu = np.mean(windowed_mu[(windowed_mu<=qmax)&(windowed_mu>=qmin)])
    #qrng_filt_med = np.median(windowed_mu[(windowed_mu<=qmax)&(windowed_mu>=qmin)])
    #qrng_filt_sd = np.sqrt(np.var((windowed_mu[(windowed_mu<=qmax)&(windowed_mu>=qmin)]))

    ##
    b1 = (w_starts>103800000)&(w_ends<103900001)
    b2 = (w_starts>102500000)&(w_ends<103500001)
    #WRONG COORDS!
    #b1 = (w_starts>103780000)&(w_ends<103500001)
    #b2 = (w_starts>103800000)&(w_ends<103900001)
    x = windowed_mu[np.where(b1|b2)]
    CTRL_region_mean = np.mean(x)
    CTRL_region_med = np.median(x)
    CTRL_region_sd = np.sqrt(np.var(x))
    """ 
    print("qmax_filt_mu:",qmax_filt_mu)
    print("qmax_filt_med:",qmax_filt_med)
    print("qmax_filt_sd:",qmax_filt_sd)
    print("CTRL_region_mean:",CTRL_region_mean)
    print("CTRL_region_med:",CTRL_region_med)
    print("CTRL_region_sd:",CTRL_region_sd)
    """ 
    
    t = pd.DataFrame([{"qmax_filt_mu":qmax_filt_mu,
                       "qmax_filt_med":qmax_filt_med,
                       "qmax_filt_sd":qmax_filt_sd,
                       "CTRL_region_mean":CTRL_region_mean,
                       "CTRL_region_med":CTRL_region_med,
                       "CTRL_region_sd":CTRL_region_sd,
                       "sample":args.sample}])

    t.to_csv(args.fn_out,sep="\t", index=False)


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
    
    #coverage
    parser_create = subparsers.add_parser("coverage")
    parser_create.add_argument("--d4_file", required=True)
    parser_create.add_argument("--fn_stats", required=True)
    parser_create.add_argument("--sample", required=True)
    parser_create.add_argument("--contig", required=True)
    parser_create.add_argument("--start", required=True,type=int)
    parser_create.add_argument("--end", required=True,type=int)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("-w", default=1000, type=int)
    parser_create.add_argument("-s", default=200, type=int)
    parser_create.set_defaults(func=coverage)

    #genotype
    parser_create = subparsers.add_parser("genotype")
    parser_create.add_argument("--d4_file", required=True)
    parser_create.add_argument("--sample", required=True)
    parser_create.add_argument("--fn_stats", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("--fn_loci", required=True)
    parser_create.add_argument("-w", default=1000, type=int)
    parser_create.add_argument("-s", default=200, type=int)
    parser_create.set_defaults(func=genotype)

    #d4 stats
    parser_create = subparsers.add_parser("stats")
    parser_create.add_argument("--d4_file", required=True)
    parser_create.add_argument("--sample", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("-w", default=1000, type=int)
    parser_create.add_argument("-s", default=500, type=int)
    parser_create.set_defaults(func=stats)

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
    
    
