import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.ticker import ScalarFormatter

df = pd.read_csv('r_lambda.dat', header=None, delim_whitespace=True).values

l = df[:, 0]
lin = df[:, 1]

fig, ax = plt.subplots(figsize=(5, 5))
ax.grid(False)

ax.set_xlabel(r'$\lambda$', fontsize=15)
ax.set_ylabel(r'$\lambda_{in}$', fontsize=15)
ax.set_xlim(0, max(l))

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.plot(l, lin, color = "black")
           
output_image = 'linl.png'
plt.savefig(output_image)

plt.show()
