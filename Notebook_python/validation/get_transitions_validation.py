#!/usr/bin/env python3
# coding: utf-8

import sys
import numpy as np
import pandas as pd
import argparse

Min_t = 4 # minimum number of transitions

ap = argparse.ArgumentParser()
ap.add_argument('csv1',help='tab-separated file with gene occurrence (genes by rows)') 
ap.add_argument('csv2',help='tab-separated validation dataset')
ap.add_argument('-c','--count_consecutive',action='store_true',
                                help='Do not penalize concecutive transitions (e.g. 101)')

args = ap.parse_args()
#args = ap.parse_args("~/Ricerca/OrthoDBv10/Eukaryota.csv.sel.csv.ord.num ~/Ricerca/OrthoDBv10/HogProf/pyprofiler/notebooks/my_humanopt.csv".split())
#args = ap.parse_args("~/Ricerca/OrthoDBv10/Eukaryota.csv.sel.csv.ord.num ~/Ricerca/OrthoDBv10/Validation_ROC/yeastw_opt_trans.csv".split())

csv1 = pd.read_table(args.csv1, header=None, index_col=0, comment='#')

ngenes, norgs = len(csv1.index), len(csv1.columns)

csv2 = pd.read_table(args.csv2,  comment='#', sep='\t')

tr = csv1.diff(axis=1) #iterative difference
tr = tr.applymap(lambda x: 1 if x>1 else -1 if x<-1 else x) #increasing (1) or decreasing (-1) only
tr_l = tr.values.tolist()

if not args.count_consecutive: #count only once consecutive state transitions (-1,1)
        for r in tr_l:
                for i in range(len(r)-1):
                        if (r[i]+r[i+1]==0):
                                r[i+1]=0


t01 = [set(np.nonzero(row > 0)[0]) for row in np.array(tr_l)] # 0->1 transitons
t10 = [set(np.nonzero(row < 0)[0]) for row in np.array(tr_l)] # 1->0 transitons
tt = [len(a | b) for a,b in zip(t01,t10)] #total transitions

sys.stderr.write("done transitions\n")

print('G1','G2','orgs','t1','t2','c','d','k', sep="\t")
counter = 0
for og1,og2 in zip(csv2.og1,csv2.og2):
	try:
		i =  tr.index.get_loc(og1)
	except:
		sys.stderr.write(str(og1)  +" missing\n")
		continue
	else:
		try:
			j =  tr.index.get_loc(og2)
		except:
			sys.stderr.write(str(og2)  +" missing\n")
			continue
		else:
			concordant = len(t01[i] & t01[j]) + len(t10[i] & t10[j])
			discordant = len(t10[i] & t01[j]) + len(t01[i] & t10[j])
			k = concordant - discordant
			if abs(k) >= -1:
				t1 = len(t01[i] | t10[i])
				t2 = len(t01[j] | t10[j])
				print(tr.index[i],tr.index[j], norgs, t1, t2, concordant, discordant, k, sep='\t')
				counter += 1

# All done:
sys.stderr.write("done concordance\n")
#summary = "Gene pairs: " + str(int((ngenes*(ngenes-1))/2)) + "; >cutoff: " + str(counter) + "\n" 
#sys.stderr.write(summary)
