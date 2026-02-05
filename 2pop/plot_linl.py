import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.ticker import ScalarFormatter

w = 1.5
p = 0.5
run = 1

folder = f"w{w:.2f}/p{p:.2f}/"

dfl = pd.read_csv(folder + 'linl.csv', header=None)
dfc = pd.read_csv(folder + 'critical.csv', header=None)

l = dfl.iloc[:, 0]
mean, std = dfl.iloc[:, 1], dfl.iloc[:, 2]

lc = dfc.iloc[:, 0]
linc = dfc.iloc[:, 1]

fig, ax = plt.subplots(figsize=(5, 5))
ax.grid(False)

ax.set_xlabel(r'$\lambda$', fontsize=15)
ax.set_ylabel(r'$\lambda_{in}$', fontsize=15)
ax.set_xlim(0, max(l))
ax.set_ylim(0, max(mean+std))

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.fill_between(l, mean - std, mean + std, color='black', alpha=0.2)
plt.plot(l, mean, color='black', linewidth = 1.5)

plt.fill_between(l, linc[0] - linc[1], linc[0] + linc[1], color='turquoise', alpha=0.2)
ax.axhline(y = linc[0], color='turquoise', linestyle='--', linewidth=1.5)

plt.axvspan(xmin = lc[0] - lc[1], xmax = lc[0] + lc[1], color='magenta', alpha=0.2)
ax.axvline(x = lc[0], color='magenta', linestyle='--', linewidth=1.5)

output_image = folder + 'linl.png'
plt.savefig(output_image)

plt.show()
