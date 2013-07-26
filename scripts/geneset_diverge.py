#!/usr/bin/env python

import sys
import numpy as np
from random import shuffle

'''
def contiguous_sets(samples):
    sets = []
    L = len(samples)+1
    for i in range(1,L):
        for j in range(L-i):
            s = samples[j:j+i]
            sets.append(s)
    return sets

def get_norm(n, max_set_size):
    sets = contiguous_sets(range(n))
    norm = np.zeros(n)
    for s in sets:
        if len(s) <= max_set_size:
            norm[s] += 1

    return norm
'''

def calc_accum_diverge(genelist, genediv, genechi, maxdivnum):
    '''calculate acumulated divergence'''
    accum = np.zeros(maxdivnum)
    for name in genelist:
        if name in genediv:
            divs = genediv[name]
            for d in divs:
                accum[d] += genechi[name]
    #accum /= norm # normalize
    return accum

def geneset_diverge(diverge_file, genelist_file):
    '''diverge_file is output from maxdiverge.py'''
    genediv = {}
    genechi = {}
    maxdivlen = 3 # max segment size
    maxdivnum = 6 # number of time points
    with open(diverge_file, 'r') as f:
        for line in f:
            name,div,chi2 = line.strip().split()
            if len(div) <= maxdivlen:
                genediv[name] = map(int,list(div))
                genechi[name] = float(chi2)

    # sample distribution
    genelist = []
    with open(genelist_file, 'r') as f:
        for line in f:
            name = line.strip()
            genelist.append(name)

    sample_accum = calc_accum_diverge(genelist, genediv, genechi, maxdivnum)

    # build null distribution
    rand_accum = np.zeros(maxdivnum)
    sys.stderr.write("building null distribution for " + str(len(genelist)) + " genes...\n")
    numiter = 1000
    for _ in range(numiter):
        randgenes = genediv.keys()
        shuffle(randgenes)
        randgenes = randgenes[0:len(genelist)]
        rand_accum += calc_accum_diverge(randgenes, genediv, genechi, maxdivnum)

    rand_accum /= float(numiter)

    return sample_accum / rand_accum

if len(sys.argv) == 3:
    norm_accum = geneset_diverge(sys.argv[1],sys.argv[2])
    print ','.join(map(str,norm_accum))
else:
    print "usage:",sys.argv[0],"<maxdiverge.py output> <geneset>"
