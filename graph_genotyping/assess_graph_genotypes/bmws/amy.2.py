#The same as the LCT_data_analysis notebook, but in a script...
import sys
sys.path.append("/global/scratch/users/nicolas931010/amylase_diversity_project/src/bmws/src/bmws")

import estimate, sim
import matplotlib.pyplot as plt
import numpy as np
from betamix import sample_paths, BetaMixture
from sim import sim_and_fit
from scipy.stats import beta
from math import exp, log
import pickle
# genes = ["LCT"]
# obs = {x: np.loadtxt("../../paper/data/data/Britain_" + x + ".txt").astype(int)[::-1] for x in genes}

seed = int(sys.argv[1])

# Load variables from the file
with open('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/amy.environment.pkl', 'rb') as f:
    environment = pickle.load(f)
    
# Access variables from the loaded environment
s_hat = environment['s_hat']
log10_lam = environment['log10_lam']
n = environment['n']
k = environment['k']
Ne = environment['Ne']
paths = environment['paths']

res=sim_and_fit({"s":s_hat}, seed=seed, lam=10 ** log10_lam, n=n, k=k, Ne=Ne, af=paths[seed], em_iterations=10)
np.savetxt("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/bootstrap/amy.s_hat." + str(seed) + ".txt", res["s_hat"])
