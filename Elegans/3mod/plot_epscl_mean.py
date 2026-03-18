import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# read file (3 columns: x, mean, stddev)
data1 = np.loadtxt("epsilon_lambda_m1.dat")
data2 = np.loadtxt("epsilon_lambda_m2.dat")
data3 = np.loadtxt("epsilon_lambda_m3.dat")

x = data1[:,0]
mean1 = data1[:,1]
std1 = data1[:,2]
mean2 = data2[:,1]
std2 = data2[:,2]
mean3 = data3[:,1]
std3 = data3[:,2]

fig, ax1 = plt.subplots(figsize=(5, 5))
ax1.grid(False)

ax1.set_xlabel(r'$\lambda_{in}$', fontsize=15)
ax1.set_ylabel(r'$\epsilon$', fontsize=15)

ax1.tick_params(axis='both', which='major', labelsize=13)

plt.plot(x, mean1, color='deeppink')
plt.fill_between(x, mean1-std1, mean1+std1, alpha=0.3, color='deeppink')

plt.plot(x, mean2, color='deepskyblue')
plt.fill_between(x, mean2-std2, mean2+std2, alpha=0.3, color='deepskyblue')

plt.plot(x, mean3, color='darkorange')
plt.fill_between(x, mean3-std3, mean3+std3, alpha=0.3, color='darkorange')

ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))

output_image = 'epsc.png'
plt.savefig(output_image)

plt.show()
