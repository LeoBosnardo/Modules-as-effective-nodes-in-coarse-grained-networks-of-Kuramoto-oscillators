import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

dfinc = pd.read_csv('epsilon_lambda_m3_inc_v3.dat', header=None, delim_whitespace=True).values

lininc = dfinc[:, 0]
epscinc = dfinc[:, 1]

dfdec = pd.read_csv('epsilon_lambda_m3_dec_v3.dat', header=None, delim_whitespace=True).values

lindec = dfdec[:, 0]
epscdec = dfdec[:, 1]

fig, ax1 = plt.subplots(figsize=(7, 7))
ax1.grid(False)

ax1.set_xlabel(r'$\lambda_{in}$', fontsize=15)
ax1.set_ylabel(r'$\epsilon$', fontsize=15)

ax1.scatter(lininc, epscinc, color='red', marker='^')

ax1.scatter(lindec, epscdec, color='blue', marker='v')

output_image = 'epsc.png'
plt.savefig(output_image)

plt.show()
