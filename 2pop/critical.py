import pandas as pd
import glob
import numpy as np

w = 1.5
p = 0.5

folder = f"w{w:.2f}/p{p:.2f}/" 

file_pattern = folder + f"run_*/rl.csv" 
files = sorted(glob.glob(file_pattern))

trigger_data = []

for f in files:
    df = pd.read_csv(f, header=None)
    last_col = df.iloc[:, -1]

    idx = (last_col == 1).idxmax()      
    trigger_data.append([df.iloc[idx, 0], df.iloc[idx, 1]])

temp_df = pd.DataFrame(trigger_data)

final_matrix = pd.DataFrame([
    temp_df.mean().values,
    temp_df.std().values])

final_matrix.to_csv(folder + 'critical.csv', index=False, header=False)