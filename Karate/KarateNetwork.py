import networkx as nx
import numpy as np

np.set_printoptions(threshold=np.inf)

G = nx.karate_club_graph()
A = nx.to_numpy_array(G, dtype=int, weight=None)

mr_hi_0based = [n for n, d in G.nodes(data=True) if d['club'] == 'Mr. Hi']
officer_0based = [n for n, d in G.nodes(data=True) if d['club'] == 'Officer']

#mr_hi = [n + 1 for n in mr_hi_0based]
#officer = [n + 1 for n in officer_0based]

khh = sum(1 for u, v in G.edges() if u in mr_hi_0based and v in mr_hi_0based) / 17
kho = sum(1 for u, v in G.edges() if u in mr_hi_0based and v in officer_0based) / 17
koh = kho
koo = sum(1 for u, v in G.edges() if u in officer_0based and v in officer_0based) / 17

print("Mr. Hi (base 0):", mr_hi_0based)
print("Officer (base 0):", officer_0based)

print("khh = ", khh, " kho = ", kho, " koh = ", koh, " koo = ", koo)

np.savetxt("A.txt", A, fmt="%d")
np.savetxt("k.txt", np.array((khh,kho,koh,koo)))
np.savetxt("mr_hi.txt", np.array(mr_hi_0based), fmt="%d")
np.savetxt("officer.txt", np.array(officer_0based), fmt="%d")