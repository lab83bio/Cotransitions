#!/usr/bin/env python3

import pandas as pd
import json, urllib
import tqdm

outfile = '../external/ko_og_pathway.tsv'

kegg_ko_reaction = pd.read_table('http://rest.kegg.jp/link/ko/pathway', 
                                 header=None, names=['reaction', 'ko'])

ko_set = {x.split(':')[1] for x in kegg_ko_reaction['ko'].to_list()}

def get_og(ko):
    try:
        j = pd.read_table(f"https://www.orthodb.org/tab?query={ko}&level=2759")
                        
        j = j.groupby('pub_og_id')['og_name'].count().reset_index()
        j.columns = ['og', 'count']
        j['ko'] = ko
        return j[['ko', 'og', 'count']]
    except:
        j = pd.DataFrame([[ko, None, None]])
        j. columns = ['ko', 'og', 'count']
        return j
    
for x in tqdm.tqdm(list(ko_set)):
    if list(ko_set).index(x) == 0:
        koog = get_og(x)
        koog.to_csv(outfile, index=False, sep='\t')
    else:
        koog = get_og(x)
        koog.to_csv(outfile, index=False, sep='\t', mode='a', header=False)
