import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

df = pd.read_csv('rs_lin_tempo_acoplado.dat', header=None, delim_whitespace=True).values
#df1 = pd.read_csv('rs_lin_tempo_m1.dat', header=None, delim_whitespace=True).values
#df2 = pd.read_csv('rs_lin_tempo_m2.dat', header=None, delim_whitespace=True).values
#df3 = pd.read_csv('rs_lin_tempo_m3.dat', header=None, delim_whitespace=True).values

t = df[:, 0]
#r = df[:, 1]
r1 = df[:, 2]
r2 = df[:, 3]
r3 = df[:, 4]
lin = df[:, 5]

fig, ax1 = plt.subplots(figsize=(5, 5))
ax1.grid(False)

ax1.tick_params(axis='both', which='major', labelsize=13)

ax1.set_xlabel('simulation Time', fontsize=15)
ax1.set_ylabel('order parameters', fontsize=15)
ax1.set_xlim(0, max(t))
ax1.set_ylim(0, 1)

#ax1.plot(t, r, color='purple', linewidth=1.5)
ax1.plot(t, r3, color='darkorange', linewidth=1.5)
ax1.plot(t, r2, color='deepskyblue', linewidth=1.5)
ax1.plot(t, r1, color='deeppink', linewidth=1.5)

ax1.axhline(y=0.9, color='black', linestyle='--', linewidth=1.5)

ax2 = ax1.twinx() 
ax2.set_ylabel(r'$\lambda_{in}$', fontsize=15)
ax2.plot(t, lin, color='black', linewidth = 1.5)
ax2.tick_params(axis='both', which='major', labelsize=13)

plt.tight_layout()

output_image = 'rt.png'
plt.savefig(output_image)

plt.show()
