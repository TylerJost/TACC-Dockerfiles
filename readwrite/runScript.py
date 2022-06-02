#!/usr/bin/env python3
import pandas as pd
import os
print('Current directory')
print(os.listdir())
df = pd.DataFrame([[1,2,3],[4,5,6],[7,8,9]])
df.to_csv('example.csv')
print('Data Directory:')
print(os.listdir('../data'))
dfNew = pd.read_csv('../data/example2.csv', index_col=0)
print(dfNew.head())
