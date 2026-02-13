import pandas as pd
import numpy as np

threshold = 0.01

#file
file = 'rl.csv'
df = pd.read_csv(file, header=None)

# reorder
df = df.sort_values(by=0, ascending=True).reset_index(drop=True)
df.to_csv(file, index=False, header=None)

#threshold
last_col_index = df.columns[-1]
condition = df[last_col_index] < threshold
meeting_threshold = condition.idxmax() if condition.any() else len(df)
binary_col = np.where(df.index >= meeting_threshold, 1, 0)
df['Binary_Threshold'] = binary_col
df.to_csv(file, index=False, header=None)