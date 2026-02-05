import pandas as pd
import glob

w = 1.5
p = 0.5

folder = f"w{w:.2f}/p{p:.2f}/" 

file_pattern = folder + f"run_*/rl.csv" 
files = sorted(glob.glob(file_pattern))

all_data = []
    
for i, f in enumerate(files):

    df = pd.read_csv(f, header=None, usecols=[0, 1])
        
    if i == 0:
        first_column = df.iloc[:, 0]
            
    all_data.append(df.iloc[:, 1])

    combined_df = pd.concat(all_data, axis=1)

    result = pd.DataFrame({
        'Col1': first_column,
        'Mean_Col2': combined_df.mean(axis=1),
        'StdDev_Col2': combined_df.std(axis=1)
    })

output_name = folder + 'linl.csv'
result.to_csv(output_name, index=False, header=None)