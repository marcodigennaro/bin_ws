#!/home/mdi0316/anaconda3/bin/python

import sys
import pandas as pd

csv_file = sys.argv[1]
df_file = pd.read_csv(csv_file)
for idx, row in df_file.iterrows():
  print(row)
