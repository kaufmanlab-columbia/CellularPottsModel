# Reduce CPM lattice files to a single np binary

import os
from joblib import Parallel, delayed
import numpy as np
import sys

def read_lattice(num):
    d = np.loadtxt('lattice_'+str(num)+'_.txt')
    return d

def reduce_files():
	txt_files = Parallel(n_jobs=n_jobs,verbose=1)(delayed(read_lattice)(i) for i in range(0,total+1,1))
	final_set = np.vstack((txt_files)).astype(np.uint8) # for cell counts > 8-bit range set to 16bit int
	np.save('lattice_condensed',final_set)

total = int(sys.argv[1])
n_jobs = int(sys.argv[2])
reduce_files()