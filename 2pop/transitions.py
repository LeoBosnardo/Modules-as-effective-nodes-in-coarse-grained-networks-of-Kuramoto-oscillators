import pandas as pd
import glob
import re

file_pattern = "w*/p*/critical.csv"
files = glob.glob(file_pattern)

all_results = []

for f in files:
    parts = re.findall(r"w([-+]?\d*\.\d+|\d+)/p([-+]?\d*\.\d+|\d+)", f)
    
    if parts:
        w_val, p_val = parts[0]
        
        df_critical = pd.read_csv(f, header=None)
        
        stats = df_critical.values.flatten().tolist()
        
        row = [float(w_val), float(p_val)] + stats
        all_results.append(row)

final_df = pd.DataFrame(all_results)

final_df.to_csv('transitions.csv', index=False, header=None)