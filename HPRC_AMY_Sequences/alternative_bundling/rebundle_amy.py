from fasta_reader import read_fasta
import pdb
import pgrtk
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex as to_hex



def output_bed(new_sdb,fout):
    seq_info = new_sdb.seq_info.copy()
    sinfo = list(seq_info.items())
    sinfo.sort(key=lambda x: x[1][0])


    bed_file = open(fout,"w")
    for sid, data in sinfo:
        
        ctg, _, _ = data

        ##print(ctg)
        ctg_items = ctg.split("_")
        ##print(ctg_items[0])
        ctg_bgn = int(ctg_items[-3])
        ctg_end = int(ctg_items[-2])
        ctg_dir = int(ctg_items[-1])
        #assert(ctg_dir==0)

        smps = sid_smps[sid]
        smp_partitions = pgrtk.group_smps_by_principle_bundle_id(smps, 2500, 10000)
        mi = 0
        
        ##plt.figure(figsize=(24,3))
        smp_partitions.reverse()
        for p in smp_partitions:
            b = p[0][0][2]
            e = p[-1][0][3] + shmmrspec["k"]
            bid = p[0][1]
            
            direction = p[0][2]
            
            print(ctg, ctg_bgn+b, ctg_bgn+e, 
                  "{}:{}:{}:{}".format(bid, direction, p[0][3], p[-1][3]), sep="\t", 
                  file=bed_file)
            
    bed_file.close()


def color_graph_by_bundles(principal_bundles,
                           fn_gfa,
                           fn_colors):
    
    
    cmap=plt.get_cmap("nipy_spectral")
    lpb = len(principal_bundles)
    color_theme0 = cmap(np.linspace(0.1, 0.9, lpb))
    #pseudo-randomize the colors
    idx = np.array([_ * (lpb+8011) for _ in range(lpb)]) % lpb
    color_theme = color_theme0[idx]

    v_to_name = {}
    with open(fn_gfa) as f:
        for r in f:
            r = r.strip().split("\t")
            if r[0] != "S":
                continue
            n = r[1]
            v = r[4].split(":")[-1]
            v_to_name[v] = n
    #pdb.set_trace()
    f = open(fn_colors, "w")
    print("Name,Color,Bundle_ID", file=f)
    for bundle in principal_bundles:
        bundle_id = bundle[0]
        for v in bundle[2]:
            vertex = tuple(v[:2])
            color = to_hex(color_theme[bundle_id])
      
            print(v_to_name["{:016x}_{:016x}".format(*vertex)], color, bundle_id, file=f, sep=",")
    f.close()
    #pdb.set_trace()



if __name__=="__main__":
    
    outdir="/global/scratch/users/psudmant/projects/amylase_diversity_project/HPRC_AMY_Sequences/alternative_bundling/output" 
    seq_path = "/global/scratch/users/psudmant/projects/amylase_diversity_project/HPRC_AMY_Sequences/AMY1A_region_seq.fa.gz"    

    shmmrspec = {"w": 48, "k":56, "r":4, "min_span":28 }
    shmmr_id = "{w}_{k}_{r}_{min_span}".format(w = shmmrspec["w"], 
                                               k = shmmrspec["k"], 
                                               r = shmmrspec["r"], 
                                               min_span = shmmrspec["min_span"])
    new_sdb = pgrtk.SeqIndexDB() 
    new_sdb.load_from_fastx(seq_path, 
                           w = shmmrspec["w"], 
                           k = shmmrspec["k"], 
                           r = shmmrspec["r"], 
                           min_span = shmmrspec["min_span"])

    
    print("outputting gfa")
    #new_sdb.generate_mapg_gfa(0, "{out}/AMY1A_region_48_56_4_28.gfa".format(out=outdir))
    new_sdb.generate_principal_mapg_gfa(0, 8,  "{out}/AMY1A_region_PRINCIPAL_{shmmr_id}.gfa".format(out=outdir,
                                                                                                shmmr_id=shmmr_id))
    principal_bundles, sid_smps = new_sdb.get_principal_bundle_decomposition(0,8)
    sid_smps = dict(sid_smps)
    """
     get_principal_bundle_decomposition()
    Returns
    a tuple consist of two lists: (principal_bundles, seqid_smps_with_bundle_id_seg_direction)

        principal_bundles = list of (principal_bundle_id, ave_bundle_position, list_bundle_vertex)

            list_of_bundle_vertex = list of (hash0:u64, hash0:u64, direction:u8)
        
        seqid_smps_with_bundle_id_seg_direction = list of shimmer pairs in the database annotated with principal bundle id and direction
            the elements of the list are ((hash0:u64, hash1:u64, pos0:u32, pos0:u32, direction:0),
                (principal_bundle_is, direction, order_in_the_bundle))
    """
    #Bundles are sets of vertices [bundle id, ave unbdle pos, list of vertices]
    
    
    
    color_graph_by_bundles(principal_bundles, 
                           "{out}/AMY1A_region_PRINCIPAL_{shmmr_id}.gfa".format(out=outdir,
                                                                                shmmr_id=shmmr_id),   
                           "{out}/AMY1A_region_PRINCIPAL_{shmmr_id}_color.csv".format(out=outdir,
                                                                                      shmmr_id=shmmr_id))
 
    fout= "{out}/AMY1A_region_PRINCIPAL_{shmmr_id}_bundles.bed".format(out=outdir,
                                                                       shmmr_id=shmmr_id)
    output_bed(new_sdb,fout)






