# NOTE: to run this code, file 'tensor.mat' must be previously downloaded
# (check the README for instructions)  and must be renamed 'guangzhou.mat'
# Similarly, file 'mat.csv'  must be downloaded and renamed 'seattle.csv'


import os
import urllib.request
import scipy.io
import numpy as np
import pandas as pd

# load data
tensor_data = scipy.io.loadmat('guangzhou.mat')
X = tensor_data['tensor']

# reshape into matrix
X=X.reshape(X.shape[0],X.shape[1]*X.shape[2])

# select subset of first 5 days
X=X[:,:144*5]

# discard segments with missing data
discard_segments = np.where(np.sum(X==0,1)>1)[0]
X = np.delete(X, discard_segments, axis=0)

# store data
np.save('guangzhou',X)

# ---------------------------------------
# SEATTLE TRAFFIC
# ---------------------------------------
# load data
matrix_data = pd.read_csv('seattle.csv', index_col = 0)
X=matrix_data.values

# select subset of first 500 timestamps
X=X[:,:500]

# store data matrix
np.save('seattle',X)