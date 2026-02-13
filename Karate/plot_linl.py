import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from matplotlib.ticker import ScalarFormatter

w = 1.5
run = 1

folder = f"w{w:.2f}/run_{run}/"

df = pd.read_csv(folder + 'rl.csv', header=None)

l = df.iloc[:, 0]
lin = df.iloc[:, 1]

fig, ax = plt.subplots(figsize=(5, 5))
ax.grid(False)

ax.set_xlabel(r'$\lambda$', fontsize=15)
ax.set_ylabel(r'$\lambda_{in}$', fontsize=15)
ax.set_xlim(0, max(l))

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.plot(l, lin, color = "black")
           
output_image = folder + 'linl.png'
plt.savefig(output_image)

plt.show()
