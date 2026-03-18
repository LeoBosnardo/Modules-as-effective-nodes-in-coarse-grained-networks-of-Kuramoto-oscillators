import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

w = 1.5
run = 1
l = 3.0

folder = f"w{w:.2f}/run_{run}/l{l:.2f}/"

df = pd.read_csv(folder + 'rt.csv', header=None)

t = df.iloc[:, 0]
lin, r1, r2, r = df.iloc[:, 1], df.iloc[:, 2], df.iloc[:, 3], df.iloc[:, 4]

fig, ax1 = plt.subplots(figsize=(7, 7))
ax1.grid(False)

ax1.set_xlabel('simulation Time', fontsize=15)
ax1.set_ylabel('order parameters', fontsize=15)
ax1.set_xlim(0, max(t))
ax1.set_ylim(0, 1)

ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

ax1.plot(t, r1, color='tab:green', linewidth=1.5)
ax1.plot(t, r2, color='tab:olive', linewidth=1.5)
ax1.plot(t, r, color='purple', linewidth=1.5)

ax1.axhline(y=0.9, color='black', linestyle='--', linewidth=1.5)

ax2 = ax1.twinx() 
ax2.set_ylabel(r'$\lambda_{in}$', fontsize=15)
ax2.plot(t, lin, color='black', linewidth = 1.5)

output_image = folder + 'rt.png'
plt.savefig(output_image)

plt.show()
