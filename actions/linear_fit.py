#!/usr/bin/env python3 
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import  pandas as pd
import sys

file_in = sys.argv[1]

df = pd.read_csv(file_in)

x = np.array(df['IS'])
y = np.array(df['TS'])

slope, intercept, r, p, se = stats.linregress(x, y)

y_p = x * slope + intercept 

#MAE  = (np.abs(y_p - y)).mean()
#RMSE = np.sqrt(np.average((y_p-y)**2))
#RMSE = np.sqrt(((y_p-y)**2).mean())
#print(slope, intercept, r**2, MAE, RMSE)

from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
MAE = mean_absolute_error(y,y_p)
RMSE= mean_squared_error(y, y_p, squared=False)
print(slope, intercept, r**2, MAE, RMSE)
