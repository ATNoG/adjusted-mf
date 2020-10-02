# NOTE: to run this code, the data zip file must be previously downloaded (check
# the README for instructions) and must be unzipped into a folder named gas_sensor


import os
import urllib.request
import pandas as pd
import numpy as np

# ---------------------------------------
# GAS SENSOR
# ---------------------------------------
# load data
data = pd.read_csv('gas_sensor/dataset_twosources_downsampled/000_Et_H_CO_n', header=None)

# discard ts (sorted column), temperature and humidity
data=data.drop(columns=[0,1,2])

# get sensors x time matrix
X=data.values.T

# select first  500 first timestamps
X=X[:,:500]

# store dataset
np.save('gas_sensor.npy', X)
