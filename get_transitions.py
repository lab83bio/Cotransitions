#!/usr/bin/env python3
# coding: utf-8

import sys
import numpy as np
import pandas as pd
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('csv',help='tab-separated file with gene occurrence (genes by rows)')
ap.add_argument('-m','--min_transitions',default=4,
				type=int,help='Minimum number of co-evolutionary transitions in a gene pair')
ap.add_argument('-c','--count_consecutive',action='store_true',
				help='Do not penalize concecutive transitions (e.g. 101)')
args = ap.parse_args()

csv = pd.read_table(args.csv, header=None, index_col=0, comment='#')
ngenes, norgs = len(csv.index), len(csv.columns)

tr = csv.diff(axis=1) #iterative difference
tr = tr.applymap(lambda x: 1 if x>1 else -1 if x<-1 else x) #increasing (1) or decreasing (-1) only

tr_l = tr.values.tolist()

if not args.count_consecutive: #count only once consecutive state transitions (-1,1)
	for r in tr_l:
		for i in range(len(r)-1):
			if (r[i]+r[i+1]==0):
				r[i+1]=0

# sets for fast comparison
t01 = [set(np.nonzero(row > 0)[0]) for row in np.array(tr_l)] # 0->1 transitons
t10 = [set(np.nonzero(row < 0)[0]) for row in np.array(tr_l)] # 1->0 transitons
tt = [len(a | b) for a,b in zip(t01,t10)] #total transitions

sys.stderr.write("done transitions\n")

print('Orthogroup1','Orthogroup2','orgs','t1','t2','c','d','k', sep="\t")
counter = 0
for i in range(ngenes-1):
    for j in range(i+1,ngenes):
        concordant = len(t01[i] & t01[j]) + len(t10[i] & t10[j])
        discordant = len(t10[i] & t01[j]) + len(t01[i] & t10[j])
        k = concordant - discordant
        if abs(k) >= args.min_transitions:
            print(tr.index[i],tr.index[j], norgs, tt[i], tt[j], concordant, discordant, k, sep='\t')
            counter += 1

# All done:
sys.stderr.write("done concordance\n")
summary = "Gene pairs: " + str(int((ngenes*(ngenes-1))/2)) + "; >cutoff: " + str(counter) + "\n" 
sys.stderr.write(summary)
