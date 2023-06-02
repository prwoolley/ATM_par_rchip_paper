import numpy as np
import pandas as pd

chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

# Inputs are extended_signal numpy arrays
beds = [
'sh.august2022.dom_transcript.promoter.',
'sh.august2022.dom_transcript.tss.',
'sh.august2022.dom_transcript.terminator.',
'sh.august2022.dom_transcript.tts.',
'sh.august2022.gene.promoter.',
'sh.august2022.gene.tss.',
'sh.august2022.gene.terminator.',
'sh.august2022.gene.tts.'
]


def extend(a, b, c):
    return np.ones(b-a)*c

for i in chrs:
    dna_vec = np.load('/home1/05664/pwoolley/work/scripts/chip.smoothing/dna_vec/gc_vec/'+i+'_vec.npy')
    max_len = len(dna_vec)
    del dna_vec

    # Adding zeros to the end of each array if needed
    for bed in beds:
        result = np.zeros(max_len)
        loc = pd.read_csv('/home1/05664/pwoolley/work/scripts/chip.smoothing/inputs/SHSY5Y/rnaseq/'+bed+'sorted_merged.bed',sep = '\t',names=['chr','start','end'],header=None,index_col=None) # Reading gene/transcript bedgraphs
        loc = loc[loc['chr']==i] # saving rows with chromosomes that match the current one
        loc = loc[['start','end']].to_numpy() # saving info as numpy matrix for multi column list comprehension
        for row in loc:
            result[int(row[0])-1:int(row[1])-1] = 1        

        np.save('/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/'+bed+i+'.npy',result)

