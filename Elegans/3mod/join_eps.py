import numpy as np

# file lists
files = [
    "epsilon_lambda_m1_dec_v1.dat",
    "epsilon_lambda_m1_dec_v2.dat",
    "epsilon_lambda_m1_dec_v3.dat",
    "epsilon_lambda_m1_inc_v1.dat",
    "epsilon_lambda_m1_inc_v2.dat",
    "epsilon_lambda_m1_inc_v3.dat"
]

# read first file to get x
data0 = np.loadtxt(files[0])
x = data0[:,0]

# store y values
y_all = np.zeros((len(x), len(files)))

for i, f in enumerate(files):
    data = np.loadtxt(f)
    y_all[:, i] = data[:,1]

# statistics
y_mean = np.mean(y_all, axis=1)
y_std  = np.std(y_all, axis=1)

# output
out = np.column_stack((x, y_mean, y_std))
np.savetxt("epsilon_lambda_m1.dat", out)