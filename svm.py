import pandas as pd

pockets = "/home/r0934354/Documents/overview_not3.2.1.csv"
total = "/home/r0934354/Documents/overview_not3.2.1_faults.csv"

df_pockets = pd.read_csv(pockets, index_col=0)
df_total = pd.read_csv(total, index_col=0)
df = pd.merge(df_pockets, df_total, on='id', how='right')
print(df)