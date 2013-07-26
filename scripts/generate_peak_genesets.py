#!/usr/bin/env python

import sys
import os
from math import log10

'''
builds lists of genes with maximal divergence over width (1..n-1) timepoints
centered at timepoints (0..n)
output files are in peak_genesets/ filename scheme is e.g.:
maxwidth2peak1.txt
'''

class Diverge:
    def __init__(self, diverge):
        self.gene, self.rawpoints, self.chi2, self.expdiff, self.foldexpdiff = diverge.strip().split()
        self.expdiff = float(self.expdiff)
        self.logfoldexpdiff = log10(float(self.foldexpdiff))
        self.width = len(self.rawpoints)
        self.points = map(int, list(self.rawpoints))

if len(sys.argv) == 3:
    diverge = []
    ntp = 0 # number of time points
    with open(sys.argv[1], 'r') as df:
        for line in df:
            dvg = Diverge(line)
            if dvg.width > ntp:
                ntp = dvg.width
            diverge.append(dvg)

    widths = range(1,ntp-1) # 1..ntp-2
    peaks  = range(ntp)   # 0..ntp-1

    outdir = sys.argv[2]
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for w in widths:
        for p in peaks:
            up_filename = outdir + "/maxwidth" + str(w) + "peak" + str(p) + ".up.txt"
            dn_filename = outdir + "/maxwidth" + str(w) + "peak" + str(p) + ".dn.txt"
            up_out = open(up_filename, 'w')
            dn_out = open(dn_filename, 'w')
            for d in diverge:
                if d.width <= w and p in d.points:
                    if d.expdiff > 0:
                        up_out.write(d.gene + "\t" + d.rawpoints +  "\t" + str(d.logfoldexpdiff) + "\t" + str(d.chi2) + "\n")
                    else:
                        dn_out.write(d.gene + "\t" + d.rawpoints +  "\t" + str(d.logfoldexpdiff) + "\t" + str(d.chi2) + "\n")
            up_out.close()
            dn_out.close()
 
else:
    print "usage:",sys.argv[0],"<maxdiverge.py output> <outdir>"
