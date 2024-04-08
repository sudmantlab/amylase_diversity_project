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

N_BOOSTRAP_REPS=1000
log10_lam=4.5
data = np.loadtxt("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/amy.txt").astype(int)[::-1]
L = len(data)
Ne = np.full(L - 1, 1e4)
# gr = exp(log(100) / (L))
# Ne=np.array([round(10000 * gr ** (10 * int(x / 10))) for x in range(L - 1)])
s_hat, prior = estimate.estimate_em(data, Ne, lam=10 ** log10_lam, em_iterations=10)
paths, _ = sample_paths(s_hat, Ne, data, N_BOOSTRAP_REPS, prior=prior)
nk=[(i, d[0]) for i, d in enumerate(data) if d[0] > 0]
n=[d[1] for d in nk]
k=[d[0] for d in nk]
np.savetxt("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/amy.s_hat.txt", s_hat)
np.savetxt("/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/amy.paths.tsv", np.array(paths), delimiter="\t")

# Save variables to a file
with open('/global/scratch/users/nicolas931010/amylase_diversity_project/graph_genotyping/assess_graph_genotypes/bmws/amy.environment.pkl', 'wb') as f:
    pickle.dump({'s_hat': s_hat, 'log10_lam': log10_lam, 'n': n, 'k': k, 'Ne':Ne, 'paths':paths}, f)
