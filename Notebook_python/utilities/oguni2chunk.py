import os
import sys
import glob
import tqdm
import warnings
import pandas as pd
from multiprocessing import Pool

warnings.filterwarnings('ignore')

oguninfile = sys.argv[1] # 'goa_uniprot_all.gaf.gz'
oguni = sys.argv[2] # 'og_uni_Eukaryota.tsv'
ogunigo = sys.argv[3] # Eukaryota.raxml.RL.csv.num.transitions.annotated.unigo
proc = sys.argv[4] # processes

print(sys.argv[1:])

abspath = os.path.dirname(os.path.abspath(oguninfile))

os.chdir(abspath)
os.mkdir('uni_go')

chunk_size=500000
os.system(f'zcat {oguninfile} | split -dl {chunk_size} - uni_go/chunk')
og_uni = pd.read_table(os.path.basename(oguni))

def merge_chunk_unigo(x):
    df = pd.read_table(x, 
                       comment='!', 
                       header=None,
                       low_memory=False)
    ch = og_uni.merge(df, left_on='id', right_on=1)
    ch['part'] = x
    return ch

total = (pd.concat(list(tqdm.tqdm(Pool(int(proc)).imap(merge_chunk_unigo, 
                                 glob.glob('uni_go/chunk*')) , 
                                  total=len(glob.glob('uni_go/chunk*'))))
                  ).drop_duplicates())


ch_gr = pd.merge(total.groupby(['og', 8, 4])[6].apply(list).reset_index(),
                  total.groupby(['og', 8, 4])['id'].apply(list).reset_index())
ch_gr['len'] = ch_gr[6].apply(len)
ch_gr['len_set'] = ch_gr['id'].apply(lambda x: len(set(x)))

def listone(og):
    d = ch_gr[ch_gr['og']==og]
    maxs = dict(d.groupby(['og', 8]).max().reset_index()[[8, 'len_set']].values)
    d['ratio'] = d.apply(lambda x: x['len_set']/maxs.get(x[8]), axis=1)
    return d
    
totalistone = list(tqdm.tqdm(Pool(int(proc)).imap(listone, set(ch_gr['og'])), total=len(set(ch_gr['og']))))


pd.concat(totalistone).to_csv(f"{os.path.basename(ogunigo)}", 
                              index=False, 
                              sep='\t', 
                              header=['og', 'type', 'go', 'evidence', 'uni', 'len', 'len_set','ratio']
                             )
