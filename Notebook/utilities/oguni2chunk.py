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
# os.mkdir('uni_go')

# chunk_size=500000
# batch_no=1
# dfog = pd.read_csv(os.path.basename(oguninfile), 
#                    chunksize=chunk_size, 
#                    comment='!', 
#                    compression='gzip', 
#                    sep ='\t', header=None)
# for chunk in tqdm.tqdm(dfog):
#     chunk.to_csv(f'uni_go/chunk{str(batch_no)}.csv', index=False, sep ='\t')
#     batch_no+=1

og_uni = pd.read_table(os.path.basename(oguni))

def merge_chunk_unigo(x):
    df = pd.read_table(x, low_memory=False)
    ch = og_uni.merge(df, left_on='id', right_on='1'
                     ).groupby(['og', '4', '8'])['6'].apply(
        lambda x:','.join(list(set(x)))).reset_index()
    return ch

total = (pd.concat(list(tqdm.tqdm(Pool(int(proc)).imap(merge_chunk_unigo, 
                                 glob.glob('uni_go/chunk*.csv')) , total=len(glob.glob('uni_go/chunk*.csv'))))
                  ).drop_duplicates())

total.to_csv(ogunigo, 
                 sep ='\t', 
                 index=False, 
                 header=None,
                 names=['og', 'go', 'type', 'evidence'])