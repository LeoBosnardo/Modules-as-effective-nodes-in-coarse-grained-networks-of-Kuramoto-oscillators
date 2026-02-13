import pandas as pd
import numpy as np

# ---- input files ----
edges_file = "../connections.xlsx"        # columns: node1, node2 (1..248)
topology_file = "../Elegans.xlsx"  # column: "Topology (10 mod)"

# ---- parameters ----
N = 248

# ---- read files ----
edges = pd.read_excel(edges_file)
topology = pd.read_excel(topology_file)

modules = topology["Topology (10 mod)"].to_numpy().astype(int)

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
m4 = np.where(modules == 4)[0]
m5 = np.where(modules == 5)[0]
m6 = np.where(modules == 6)[0]
m7 = np.where(modules == 7)[0]
m8 = np.where(modules == 8)[0]
m9 = np.where(modules == 9)[0]
m10 = np.where(modules == 10)[0]

np.savetxt("m1.txt", m1, fmt="%d")
np.savetxt("m2.txt", m2, fmt="%d")
np.savetxt("m3.txt", m3, fmt="%d")
np.savetxt("m4.txt", m4, fmt="%d")
np.savetxt("m5.txt", m5, fmt="%d")
np.savetxt("m6.txt", m6, fmt="%d")
np.savetxt("m7.txt", m7, fmt="%d")
np.savetxt("m8.txt", m8, fmt="%d")
np.savetxt("m9.txt", m9, fmt="%d")
np.savetxt("m10.txt", m10, fmt="%d")
np.savetxt("module.txt", modules, fmt="%d")

# ---- compute k_ij ----
modules_idx = {
    1: np.where(modules == 1)[0],
    2: np.where(modules == 2)[0],
    3: np.where(modules == 3)[0],
    4: np.where(modules == 4)[0],
    5: np.where(modules == 5)[0],
    6: np.where(modules == 6)[0],
    7: np.where(modules == 7)[0],
    8: np.where(modules == 8)[0],
    9: np.where(modules == 9)[0],
    10: np.where(modules == 10)[0]
}

k = np.zeros((10, 10))

for i in range(1, 11):
    Mi = modules_idx[i]
    for j in range(1, 11):
        Mj = modules_idx[j]
        
        edges_ij = A[np.ix_(Mi, Mj)].sum()
        k[i-1, j-1] = edges_ij / len(Mi)

# ---- save k as 9 values in one column ----
np.savetxt("K.txt", k, fmt="%.6f")
