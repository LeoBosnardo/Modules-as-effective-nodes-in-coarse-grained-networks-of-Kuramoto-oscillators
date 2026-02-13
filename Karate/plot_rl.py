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
mean, std = df.iloc[:, 5], df.iloc[:, 6]
form = df.iloc[:, 7]

triangle = (form == 0)
dot = (form == 1)

fig, ax = plt.subplots(figsize=(5, 5))
ax.grid(False)

ax.set_xlabel(r'$\lambda$', fontsize=15)
ax.set_ylabel(r'$r$', fontsize=15)
ax.set_xlim(0, max(l))
ax.set_ylim(0, 1)

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.fill_between(l, mean - std, mean + std, color='purple', alpha=0.2)

plt.scatter(l[triangle], mean[triangle], marker='^', facecolors='none', edgecolors='purple')

plt.scatter(l[dot], mean[dot], marker='o', color='purple')
           
output_image = folder + 'rl.png'
plt.savefig(output_image)

plt.show()
