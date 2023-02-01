#!/usr/bin/env python3
# coding: utf-8

import warnings
warnings.simplefilter("ignore")

import pandas as pd
import numpy as np
from Bio import Entrez

import argparse
import itertools
import requests
import gzip
import shutil
import tqdm

import omadb
import pyham

import sys
import os

import sys

ap = argparse.ArgumentParser()
ap.add_argument('-d', '--taxonomy_database', default='ncbi',
               help='Taxonomy database ("oma" for oma level, "ncbi" to filter species based on taxid)') 
ap.add_argument('-t','--taxonomy_level',default='root',
               help='Taxonomy level (e.g. root, Eukaryota, Mammalia)')
ap.add_argument('-n','--newick_path',default='speciestree.nwk')
ap.add_argument('-o','--orthoxml_path',default='oma-hogs.orthoXML')
args = ap.parse_args()

def HOG_data(hog):
    
    return hog.hog_id, [(gene.unique_id, gene.prot_id, 
                         gene.genome.taxid) 
             for gene in hog.get_all_descendant_genes()]

def HOGs_at_root_level(hogs):
    
    return [HOG_data(hog) for hog in hogs]

def HOGs_at_level(hogs, taxonomy_level):
    
    hogs_data = []
    for hog in tqdm.tqdm(hogs):
        
        descendant_levels = [i.name for i in 
            hog.get_all_descendant_hog_levels()]

        if taxonomy_level not in descendant_levels:
            continue
            
        if descendant_levels.index(taxonomy_level) == 0:
            hogs_data.append(HOG_data(hog))

        else:
            genome = hamutil.get_ancestral_genome_by_name(taxonomy_level)
            for hog2 in hog.get_at_level(genome):
                hogs_data.append((hog.hog_id, HOG_data(hog2)[1]))
    
    return hogs_data

def getlineage(taxidlist):
    
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.efetch(db="Taxonomy", id=taxidlist, retmode="xml")
    taxDf = pd.json_normalize(Entrez.read(handle))
    return [(row[0], row[5].split('; ')) 
            for row in taxDf.values.tolist()]

def taxidlevel(taxidlist):
    
    return Entrez.read(Entrez.esearch(db="Taxonomy", 
            term=taxidlist, retmode="xml"))['IdList'][0]

taxonomy_level = args.taxonomy_level
taxonomy_database = args.taxonomy_database 
newick_path = args.newick_path
orthoxml_path = args.orthoxml_path
    
taxonomy_level_oma = 'root' if taxonomy_database == 'ncbi' else taxonomy_level

print(f'OMA TRANSLATOR', file=sys.stderr)
print(f'Translation at level {taxonomy_level}, database {taxonomy_database}', file=sys.stderr)
print(f'Download {newick_path} and {orthoxml_path}', file=sys.stderr)

r = requests.get(f'https://omabrowser.org/All/{orthoxml_path}.gz')
open(f'{orthoxml_path}.gz' , 'wb').write(r.content)

with gzip.open(f'{orthoxml_path}.gz', 'rb') as f_in:
    with open(orthoxml_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

r = requests.get(f'https://omabrowser.org/All/{newick_path}')
open(newick_path , 'wb').write(r.content)

hamutil = pyham.Ham(newick_path, orthoxml_path, tree_format="newick", 
                    use_internal_name=True, species_resolve_mode="OMA")

os.remove(orthoxml_path)

print(f'Extracting HOGs data from {orthoxml_path}', file=sys.stderr)

hogs = hamutil.get_list_top_level_hogs()

allhogs = (HOGs_at_root_level(hogs) if taxonomy_level_oma == 'root' 
           else HOGs_at_level(hogs, taxonomy_level_oma))
allhogs_df = pd.DataFrame([(hog[0], [f"{i[2]}_x:{i[1]}" 
                for i in hog[1]]) for hog in allhogs])

print(f'Downloading HOGs descriptions', file=sys.stderr)

HOGS_descriptions = omadb.OMARestAPI.HOGs(omadb.Client()).list().as_dataframe()
HOGS_descriptions.to_csv('HOGs_list.tsv', sep='\t', index=None)
HOGS_descriptions = dict(HOGS_descriptions[['roothog_id',
                    'description']].values.tolist())

allhogs_df[2] = allhogs_df[0].apply(lambda x: HOGS_descriptions.get(int(x)))
allhogs_df[0] = allhogs_df[0].apply(lambda x: f'HOG:C{str(x).zfill(7)}')
allhogs_df[3] = taxidlevel(taxonomy_level)

ogs = allhogs_df[[0,3,2]]
og2genes = allhogs_df.explode(1)[[0,1]]

print(f"\t- Orthogroups found: {len(ogs)}", file=sys.stderr)
print(f"\t- Genes found: {len(og2genes)}", file=sys.stderr)
print(f"\t- Genomes found: {len(set(og2genes[1].apply(lambda x: x.split('_')[0])))}", file=sys.stderr)

if taxonomy_database == 'ncbi':
    
    print(f'Filtering taxids', file=sys.stderr)
    
    taxis = list(set([int(i.split('_')[0]) for i in og2genes[1].tolist()]))
    taxis_lineage = [getlineage(t.tolist()) 
                     for t in np.array_split(taxis, int(len(taxis)/400)+1)]
    taxis_lineage = list(itertools.chain(*taxis_lineage))

    taxonomy_df = pd.DataFrame([(t[0],';'.join(t[1])) for t in taxis_lineage])
    taxonomy_df.to_csv('HOGs_lineage.tsv', sep='\t', index=False, header=None)
    
    taxids = list(set(taxonomy_df[taxonomy_df[1].str.contains(
        taxonomy_level)][0].tolist()))

    og2genes[2] = og2genes[1].apply(lambda x: int(x.split('_')[0]))
    og2genes = og2genes[og2genes[2].astype(str).isin(taxids)]
    
    print(f"\t\t- Orthogroups found: {len(ogs[ogs[0].isin(set(og2genes[0].tolist()))])}", file=sys.stderr)
    print(f"\t\t- Genes found: {len(og2genes)}", file=sys.stderr)
    print(f"\t\t- Genomes found: {len(set(og2genes[1].apply(lambda x: x.split('_')[0])))}", file=sys.stderr)
    
outfile = f'HOGs_{taxonomy_level}_{taxonomy_database}'
print(f'Saving files at {outfile}_*', file=sys.stderr)
og2genes[[0,1]].to_csv(f'{outfile}_OG2genes.tab.gz', 
                    sep='\t', index=False, header=None, compression='gzip')
ogs[ogs[0].isin(set(og2genes[0].tolist()))].to_csv(f'{outfile}_OGs.tab.gz', 
                    sep='\t', index=False, header=None, compression='gzip')