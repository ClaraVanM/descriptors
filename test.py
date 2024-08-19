#third problem
#why still H in exposed atom list?

import pandas as pd

data = pd.read_csv('descriptors3.2.1_newest.csv')
print(data[data['H']!=0])