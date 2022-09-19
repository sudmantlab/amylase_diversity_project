import allel
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import pysam
import queue

NODE_COUNT = 0

class haplotype:
    def __init__(self, indiv, hap, rev):
        self.indiv = indiv
        self.hap = hap
        if rev:
            self.hap = self.hap[::-1]
    def is_valid(self):
        if not -1 in self.hap:
            return True
        return False

def get_AF(haps,ix):
    c = 0
    for h in haps:
        c=c+h.hap[ix]
    return float(c)/len(haps)

def get_haplotypes(callset, sample_subset,sample_to_sort_key, min_AF, rev = False):

    #samples = [s for s in callset['samples'] for i in range(2)]
    samples = []
    gts = allel.GenotypeArray(callset['calldata/GT'])
    vcf_afs = callset['variants/AF'][:,0]
    poses = callset['variants/POS']
    
    sites = np.where(callset['variants/AF'][:,0]>min_AF)
    poses = poses[sites]
    gts = gts[sites]
    
    if rev:
        poses = poses[::-1]

    haps = []
    for i,indiv in enumerate(callset["samples"]):
        if (not sample_subset is None) and not indiv in sample_subset:
            continue
        hap1 = haplotype(indiv,gts[:,i,0],rev)
        hap2 = haplotype(indiv,gts[:,i,1],rev)
        if hap1.is_valid():
            haps.append(hap1)
            samples.append(indiv)
        if hap2.is_valid():
            haps.append(hap2)
            samples.append(indiv)
    
    if (not sample_subset is None):
        samples = sorted(samples,key=lambda x: sample_to_sort_key[x])
        haps = sorted(haps, key =lambda x: sample_to_sort_key[x.indiv]) 

    afs = [get_AF(haps,i) for i in range(len(poses))]
    return haps, samples, poses

class node:
    def __init__(self, 
                 pos, 
                 indivs, 
                 level, 
                 snp_idx,
                 f):
        self.f = f
        self.pos = pos
        self.level = level
        self.snp_idx = snp_idx
        self.indivs = indivs
        self.offspring = [None]
        
        self.curr_indiv_to_next_node = {}

        global NODE_COUNT
        self.ix = NODE_COUNT
        NODE_COUNT+=1

    def get_split(self, haps, pos, sep=0, min_AF = 0.01):

        indivs_left = []
        indivs_right = []
        
        for i in self.indivs:
            if haps[i].hap[self.snp_idx+1]==0:
                indivs_left.append(i)
            else:
                indivs_right.append(i)
        
        f_left = float(len(indivs_left))/len(self.indivs) 
        f_right = float(len(indivs_right))/len(self.indivs)
        if f_left > f_right:
            new_indivs_right = indivs_left
            indivs_left = indivs_right
            indivs_right = new_indivs_right
        #print(f_left, f_right) 
        #pdb.set_trace()

        #if (f_left < min_AF) or (f_right <min_AF):
        if (f_left ==0) or (f_right ==0):

            left_node = node(pos, 
                             [i for i in self.indivs],
                             self.level+1, 
                             self.snp_idx+1,
                             self.f)
            for i in left_node.indivs:
                self.curr_indiv_to_next_node[i] = left_node
            self.offspring = [left_node]
            
            return left_node, None 

        left_node = node(pos, 
                         indivs_left, 
                         self.level+1, 
                         self.snp_idx+1,
                         f_left*self.f)
        for i in left_node.indivs:
            self.curr_indiv_to_next_node[i] = left_node
        right_node = node(pos, 
                          indivs_right, 
                          self.level+1,
                          self.snp_idx+1,
                          f_right*self.f)
        for i in right_node.indivs:
            self.curr_indiv_to_next_node[i] = right_node


        self.offspring = [left_node, right_node]
        return left_node, right_node

## assumes the data is phased
## assumes the data is all SNPs (not true for 1kg)
def get_tree(fn_vcf, contig, start, width, sample_subset, sample_to_sort_key, min_AF=0.01, rev=False):
    q = queue.Queue()    
    nodes = []
    
    if rev:
        end = start
        start = start - width
    else:
        end = start + width
    
    loc = "{contig}:{s}-{e}".format(contig=contig,
                                    s=start,
                                    e=end)

    callset = allel.read_vcf(fn_vcf,region=loc,fields="*")
    haps, samples, poses = get_haplotypes(callset, sample_subset, sample_to_sort_key, min_AF, rev=rev)
    
    if rev:
        n = node(end, [i for i in range(len(haps))], 0, -1, 1)
    else:
        n = node(start, [i for i in range(len(haps))], 0, -1, 1)

    nodes_by_ix = {n.ix:n}
    q.put(n)
    
    while q.qsize() >0:
        curr_n = q.get()
        curr_idx = curr_n.snp_idx+1 
        if curr_idx>len(poses)-1:
            continue
        
        left_node, right_node = curr_n.get_split(haps, poses[curr_idx], min_AF=min_AF)
        if not (left_node is None):
            q.put(left_node)
            nodes_by_ix[left_node.ix] = left_node
        if not (right_node is None):
            q.put(right_node)
            nodes_by_ix[right_node.ix] = right_node

    direction = rev and "left" or "right"
    outrows = []
    for ix, curr_n in nodes_by_ix.items():
        for i in curr_n.indivs:
            if curr_n.offspring[0] is None:
                outrows.append({"position":curr_n.pos,
                                "node":curr_n.ix,
                                "next_pos":"NA",
                                "next_node":"NA",
                                "direction":direction,
                                "sample":samples[i]})
            else:
                next_node_ix = curr_n.curr_indiv_to_next_node[i].ix 
                next_node_pos = curr_n.curr_indiv_to_next_node[i].pos 
                outrows.append({"position":curr_n.pos,
                                "node":curr_n.ix,
                                "next_pos":next_node_pos,
                                "next_node":next_node_ix,
                                "sample":samples[i]})
    #outrows = []
    #for ix, curr_n in nodes_by_ix.items():
    #    for i in curr_n.indivs:
    #        if curr_n.offspring[0] is None:
    #            outrows.append({"position":curr_n.pos,
    #                            "node":curr_n.ix,
    #                            "next_pos":"NA",
    #                            "next_node":"NA"})
    #        else:
    #            next_node_ix = curr_n.curr_indiv_to_next_node[i].ix 
    #            next_node_pos = curr_n.curr_indiv_to_next_node[i].pos 
    #            outrows.append({"position":curr_n.pos,
    #                            "node":curr_n.ix,
    #                            "next_pos":next_node_pos,
    #                            "next_node":next_node_ix})

    t = pd.DataFrame(outrows)
    return t
    #callset['variants/AF']
    #callset['variants/POS']
    #callset['variants/is_snp']


if __name__=="__main__":
    
    #also cool https://rdrr.io/rforge/rehh/man/bifurcation.diagram.html
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_vcf",required=True)
    parser.add_argument("--fn_out",required=True)
    parser.add_argument("--contig",required=True)
    parser.add_argument("--label",required=True)
    parser.add_argument("--start_right",required=True,type=int)
    parser.add_argument("--start_left",required=True,type=int)
    parser.add_argument("--width",type=int,default=50000)
    parser.add_argument("--min_AF",type=float,default=0.01)
    #parser.add_argument("--rev",default=False,action='store_true')
    #parser.add_argument('-s','--sample_subset', nargs='+', help='<Required> Set flag', required=False)
    parser.add_argument('-s','--sample_subset', required=False)
    
    args = parser.parse_args()
    sample_subset = None
    sample_to_sort_key = {}
    if not args.sample_subset is None:
        t = pd.read_csv(args.sample_subset,sep="\t", header=0).sort_values('popGrouping')
        sample_subset = [s for s in t['sample'].values]
        sample_to_sort_key = {r['sample']:r["popGrouping"] for i,r in t.iterrows()}


    t_right = get_tree(args.fn_vcf, 
                      args.contig, 
                      args.start_right, 
                      args.width,
                      sample_subset,
                      sample_to_sort_key,
                      min_AF=args.min_AF, 
                      rev = False) 
    t_right['dir'] = "right" 
    t_left = get_tree(args.fn_vcf, 
                      args.contig, 
                      args.start_left, 
                      args.width,
                      sample_subset,
                      sample_to_sort_key,
                      min_AF=args.min_AF, 
                      rev = True) 
    t_left['dir'] = "left" 
    
    t = pd.concat([t_left, t_right])  
    t['label'] = args.label
    t.to_csv(args.fn_out, sep="\t", index=False)
    #get table
    #plot table
    #table to tree

