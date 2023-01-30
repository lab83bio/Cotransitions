import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import re
import os
import json
import tqdm
import itertools
from functools import reduce

import requests 
import pickle
from multiprocessing import Pool

import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import textwrap
import traceback
import urllib
import networkx as nx

from Bio import SeqIO

from sklearn.metrics import roc_curve, auc


def transitionstoogs(transitions_file):
    
    transitions = pd.read_table(transitions_file)
    transitions['og1name'] = transitions['ognames'].apply(lambda x: x.split(' --- ')[0])
    transitions['og2name'] = transitions['ognames'].apply(lambda x: x.split(' --- ')[1])
    transitions = transitions.drop('ognames', axis=1)

    ogs = pd.concat([transitions[['Orthogroup1','k_score','p','p.adj','og1name']
                          ].rename(columns={'Orthogroup1':'og','og1name':'ogname'}),
               transitions[['Orthogroup2','k_score','p','p.adj','og2name']
                          ].rename(columns={'Orthogroup2':'og','og2name':'ogname'})]
             )

    ogs['p'] = ogs['p'].astype(float)
    ogs['k_score'] = ogs['k_score'].astype(float)
    ogs['p.adj'] = ogs['p.adj'].astype(float)
    ogs = ogs.sort_values('p.adj').drop_duplicates('og')
    
    ogs = ogs[ogs['k_score']>0]
    
    return ogs

def clustersfromraw(clusters_raw_file, transitions_file):

    df = pd.read_table(clusters_raw_file, header=None)
    df = pd.DataFrame(list(zip(df.index.tolist(), df.values.tolist()))).explode([1]).dropna()
    n = dict(df.groupby(0)[1].count().reset_index().values)
    df['n'] = df[0].apply(lambda x: n.get(x))

    trans = transitionstoogs(transitions_file)
    ognames = dict(trans[['og','ogname']].values)
    df['ogname'] = df[1].apply(lambda x: ognames.get(x))

    # aggiungi le righe per gli score e le transitions dalle transitions
    df['score'] = ''
    df['transition'] = ''

    df.columns = ['cluster','og','n','name','score','transition']
    df = df[['cluster','score','transition','n','og','name']]
    df = df[df['n']>1]

    return df
        
    
        
def scalepval(df, col, minrange, maxrange):
    
    median = df[col].median() 

    if not minrange:
        minrange = 2
    
    if not maxrange:
        maxrange = median-df[col].min()

    df[col] = df[col].apply(
        lambda x: x if x <= maxrange+median else maxrange+median)

    posmin, posmax = minrange, maxrange
    negmin, negmax = -maxrange, -minrange

    finalpos = df[df[col]>=median]
    finalposlog10 = finalpos[col].values
    finalposlog10scaled = ((finalposlog10 - np.min(finalposlog10)) / 
                           (np.max(finalposlog10) - np.min(finalposlog10)) * 
                           (posmax - posmin) + posmin)

    finalneg = df[df[col]<median]
    finalneglog10 = finalneg[col].values
    finalneglog10scaled = ((finalneglog10 - np.min(finalneglog10)) / 
                           (np.max(finalneglog10) - np.min(finalneglog10)) * 
                           (negmax - negmin) + negmin)

    finalpos[f'{col}scaled'] = finalposlog10scaled
    finalneg[f'{col}scaled'] = finalneglog10scaled
    
    final2 = pd.concat([finalpos, finalneg])
    
    return final2

def clustersfromann(clusters_file, level):
    
    clusters = pd.DataFrame(itertools.chain.from_iterable(map(lambda y: list(map(lambda x: [y[0]]+ x.split(f'{level}:')+[y[1].split(':')[0]],':'.join(y[1].split(':')[1:]).strip().strip(' ###').split(' ### '))),
                                                              enumerate(open(clusters_file).readlines())))).dropna(subset=[3])
    clusters.columns = ['cluster','og', 'name', 'feat']
    clusters[['score', 'transition', 'n']] = clusters['feat'].str.split(' ', expand=True)
    clusters['n'] = clusters['n'].apply(lambda x: x.strip('(').strip(')'))
    clusters['og'] = clusters['og']+level
    # clusters[['cluster', 'score', 'transition', 'n', 'og', 'name']].to_csv('Eukaryota_cluster.tsv', sep='\t', index=False)
    clusters = clusters[['cluster', 'score', 'transition', 'n', 'og', 'name']]
    return clusters