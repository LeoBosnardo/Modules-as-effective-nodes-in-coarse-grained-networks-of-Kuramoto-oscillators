import pandas as pd
import numpy as np

w = 1.5
p = 0.5
run = 5
threshold = 0.01

#file
folder = f"w{w:.2f}/p{p:.2f}/run_{run}/"
file = folder + 'rl.csv'
df = pd.read_csv(file, header=None)

# reorder
df_sorted = df.sort_values(by=0, ascending=True)
df_sorted.to_csv(file, index=False, header=None)

#threshold
last_col_index = df.columns[-1]
condition = df[last_col_index] < threshold
meeting_threshold = condition.idxmax() if condition.any() else len(df)
binary_col = np.where(df.index >= meeting_threshold, 1, 0)
df['Binary_Threshold'] = binary_col
df.to_csv(file, index=False, header=None)