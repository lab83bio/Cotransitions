#!/usr/bin/env python3

import pandas as pd
import requests
import tqdm
import pickle
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Orthogroup data download (Orthodb v101)')
parser.add_argument('-i', '--input_file', help='orthogroups list in a table')
parser.add_argument('-o', '--output_file', help='output (saved in pickle format)')
parser.add_argument('-p', '--pool_threads', type=int, help='n. of parallel threads')
args = parser.parse_args()

allogs_file = args.input_file # 'external/Eukaryota.ogs.list'
odb_file = args.output_file #'external/odbinfos.pickle'

def odbinfo(og):
    try:
        r = requests.get(f"https://v101.orthodb.org/group?id={og}").json()
        return [og, r]
    except:
        return [og, []]
    
allogs = set(pd.read_table(allogs_file)['og'])
odinfolist = list(tqdm.tqdm(Pool(args.pool_threads).imap(odbinfo, allogs), 
                            total=len(list(allogs))))

with open(odb_file, 'wb') as handle:
    pickle.dump(odinfolist, handle, protocol=pickle.HIGHEST_PROTOCOL)