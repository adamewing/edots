#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import scipy.stats as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from random import sample

'''
def plot_exp(data, randoms, fn, labels):
    means = np.mean(data, axis=0)
    s_num = range(1,len(means)+1)

    #barpngfn = basename(fn) + ".barplot.avg.png"
    linepngfn = basename(fn) + ".lineplot.avg.png"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xticklabels(labels)
    ax.set_ylabel('expression diff: human - rhesus')
    ax.boxplot(data)
    ax.plot(s_num,means)
    plt.savefig(barpngfn)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ylim(min(means)-100,max(means)+100)
    ax.set_xticklabels(labels)
    ax.set_ylabel('avg. expression difference (' + str(len(data)) + ' genes)')
    for random in randoms:
        ax.plot(s_num, np.mean(random, axis=0), c='g', alpha=0.15)

    ax.plot(s_num, means, marker='o', c='k')
    plt.savefig(linepngfn)
'''

def plot_ratio(data, randoms, fn, labels, outdir):
    means = np.mean(data, axis=0)
    print data
    print means
    s_num = range(1,len(means)+1)
    outfn = outdir + '/' + os.path.basename(fn) + ".ratio.avg.png"

    # find interquartile range
    allpoints = []
    for datum in data:
        for exp in datum:
            allpoints.append(datum)
    allpoints = np.asarray(allpoints)
    iles = map(lambda x:x/4.0, range(1,4))

    quantiles = ss.mstats.mquantiles(allpoints,prob=iles)
    bot_quant = quantiles[0]
    top_quant = quantiles[-1]

    if top_quant < max(means)+1.0:
        top_quant = max(means)+1.0
    if bot_quant > min(means)-1.0:
        bot_quant = min(means)-1.0
    print "top_quant:",top_quant,"bot_quant:",bot_quant
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

#    ax.yaxis.set_ticks_position("right")
#    ax.yaxis.set_label_position("right")

    ax.set_xticklabels(labels, size=20)
    ax.set_ylabel('avg. log expr. ratio hu/rh (' + str(len(data)) + ' genes)', size=20)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize('20')

    for random in randoms:
        plt.ylim(min(means)-1.0, max(means)+1.0) # control zoom
        ax.plot(s_num, np.mean(random, axis=0), c='g', alpha=0.15)

    for datum in data:
        plt.ylim(min(means)-1.0, max(means)+1.0) # control zoom
        ax.plot(s_num, datum, c='DodgerBlue', alpha=0.15)

    ax.plot(s_num, means, marker='o', c='k',linewidth=4)
    plt.savefig(outfn, bbox_inches='tight')

'''
def exp_diff(data, ntp, geneset=None):
    if geneset is not None:
        data = data.ix[geneset]

    expdiff = [] 

    for gene in data.iterrows():
        name = gene[0]
        h = np.asarray(gene[1][0:ntp])
        r = np.asarray(gene[1][ntp:ntp*2])
        expdiff.append(h-r)

    return np.asarray(expdiff)
'''

def exp_ratio(data, ntp, geneset=None):
    if geneset is not None:
        data = data.ix[geneset]

    expdiff = [] 

    for gene in data.iterrows():
        name = gene[0]
        h = np.asarray(gene[1][0:ntp])
        r = np.asarray(gene[1][ntp:ntp*2])
        if min(gene[1]) > 0:
            expdiff.append(h/r)

    return np.log(np.asarray(expdiff))

def load_data(csv):
    sys.stderr.write("populating dataframe...")
    data = pd.read_csv(csv,index_col=0)
    sys.stderr.write("done\n")
    return data

def load_geneset(data, infile):
    geneset = []
    with open(infile, 'r') as f:
        for line in f:
            c = line.strip().split()
            if c[0].upper() in data.index:
                geneset.append(c[0].upper())
            else:
                sys.stderr.write("warning: " + c[0] + " not in expression matrix\n")
    return geneset

if __name__ == '__main__':
    if len(sys.argv) == 4:
        data = load_data(sys.argv[1])
        assert len(data.columns) % 2 == 0
        ntp = len(data.columns)/2
        geneset = load_geneset(data, sys.argv[2])
        print "len(geneset):",len(geneset)
        random_exps = []
        random_ratios = []

        outdir = sys.argv[3]
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        for _ in range(200):
            rand_genes = sample(list(data.index),len(geneset))
            #random_exps.append(exp_diff(data, ntp, geneset=rand_genes))
            random_ratios.append(exp_ratio(data, ntp, geneset=rand_genes))

        #genelist_data_exp = exp_diff(data, ntp, geneset=geneset)
        genelist_data_ratio = exp_ratio(data, ntp, geneset=geneset)
#        plot_exp(genelist_data_exp, random_exps, sys.argv[2], list(data.columns))
        plot_ratio(genelist_data_ratio, random_ratios, sys.argv[2], list(data.columns), outdir)
    else:
        print "usage:",sys.argv[0],"<expression csv> <genelist> <outdir>"

