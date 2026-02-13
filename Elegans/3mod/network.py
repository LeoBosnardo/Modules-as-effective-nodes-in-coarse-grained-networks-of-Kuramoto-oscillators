import pandas as pd
import numpy as np

# ---- input files ----
edges_file = "../connections.xlsx"        # columns: node1, node2 (1..248)
topology_file = "../Elegans.xlsx"  # column: "Topology (3 mod)"

# ---- parameters ----
N = 248

# ---- read files ----
edges = pd.read_excel(edges_file)
topology = pd.read_excel(topology_file)

modules = topology["Topology (3 mod)"].to_numpy().astype(int)

# ---- adjacency matrix (0-based) ----
A = np.zeros((N, N), dtype=int)

for i, j in zip(edges["node1"], edges["node2"]):
    i = int(i) - 1  # 1..248 -> 0..247
    j = int(j) - 1
    A[i, j] = 1
    A[j, i] = 1   # undirected

np.savetxt("A.txt", A, fmt="%d")

# ---- module node lists (0-based) ----
m1 = np.where(modules == 1)[0]
m2 = np.where(modules == 2)[0]
m3 = np.where(modules == 3)[0]

np.savetxt("m1.txt", m1, fmt="%d")
np.savetxt("m2.txt", m2, fmt="%d")
np.savetxt("m3.txt", m3, fmt="%d")
np.savetxt("module.txt", modules, fmt="%d")

# ---- compute k_ij ----
modules_idx = {
    1: np.where(modules == 1)[0],
    2: np.where(modules == 2)[0],
    3: np.where(modules == 3)[0]
}

k = np.zeros((3, 3))

for i in range(1, 4):
    Mi = modules_idx[i]
    for j in range(1, 4):
        Mj = modules_idx[j]
        
        edges_ij = A[np.ix_(Mi, Mj)].sum()
        k[i-1, j-1] = edges_ij / len(Mi)

# ---- save k as 9 values in one column ----
np.savetxt("K.txt", k, fmt="%.6f")
