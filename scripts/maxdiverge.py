#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import scipy.stats as ss
from math import sqrt

def contiguous_sets(samples):
    '''return set of contiguous intervals of size 1 to L over samples'''
    sets = []
    L = len(samples)+1
    for i in range(1,L):
        for j in range(L-i):
            s = samples[j:j+i]
            sets.append(s)
    return sets

def col_sums(data):
    '''return vector of column lengths from data frame'''
    sums = []
    for name in data.columns:
        sums.append(sum(data[name]))
    return sums

def load_data(filename):
    '''load CSV into PANDAS dataframe'''
    sys.stderr.write("populating dataframe...")
    data = pd.read_csv(sys.argv[1],index_col=0)
    sys.stderr.write("done\n")
    return data

if len(sys.argv) == 3:
    thresh = float(sys.argv[2])

    data = load_data(sys.argv[1])
    cols = len(data.columns)
    if cols % 2 != 0:
        sys.exit("number of columns must be divisible by two\n")

    ntp  = cols/2
    sets = contiguous_sets(range(ntp))
    sums = col_sums(data)


    for gene in data.iterrows():
        name   = gene[0]
        h_expr = gene[1][0:ntp]
        r_expr = gene[1][ntp:ntp*2]
        h_sums = np.asarray(sums[0:ntp]) # numpy arrays are subscriptable with a list
        r_sums = np.asarray(sums[ntp:ntp*2])

        best_chi2 = 0.0
        best_set  = None       
        expr_diff = 0.0

        # find highest scoring (score = chi2) set 
        for s in sets:
            if min(h_expr) > 5 and min(r_expr) > 5: # chisquare prereq 
                matrix = [[sum(h_expr[s]), sum(r_expr[s])],[sum(h_sums[s]), sum(r_sums[s])]]
                chi2,p,dof,ex = ss.chi2_contingency(matrix, correction=True)
                if chi2 > best_chi2:
                    exp_diff  = np.mean(np.asarray(h_expr[s]) - np.asarray(r_expr[s]))
                    exp_fold  = np.mean(np.asarray(h_expr[s]) / np.asarray(r_expr[s]))
                    best_chi2 = chi2
                    best_set  = ''.join(map(str,s))

        if best_chi2 > thresh: 
            print name,best_set,best_chi2,exp_diff,exp_fold

else:
    print "usage:",sys.argv[0],"<csv> <chi2 cutoff>"
