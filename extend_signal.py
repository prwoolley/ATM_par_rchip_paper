import pandas as pd
import numpy as np
import scipy.signal as signal
print('libraries imported')

chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

chips = [
#'shsy5y.par.atmi.ppois.',
#'shsy5y.par.dmso.ppois.',
#'shsy5y.rchip.atmi.ppois.',
#'shsy5y.rchip.dmso.ppois.',
#'shsy5y.atac.SRR5819663_SRR5819664_treat_pileup.',
#'shsy5y.ctcf.SRR6334830.ppois.',
#'shsy5y.h3k27ac.SRR3363257.ppois.',
#'shsy5y.h3k27me3.SRR3363258.ppois.',
#'shsy5y.h3k4me1.SRR3363256.ppois.',
#'shsy5y.h3k4me3.SRR3363255.ppois.',
#'shsy5y.igg.SRR4291521_SRR4291522.ppois.',
#'shsy5y.pol2.SRR6322542_SRR6322543.ppois.',
#'shsy5y.s96.SRR4291525_SRR4291526.ppois.',
'shsy5y.rchip.atmi_nac.ppois.',
'shsy5y.rchip.dmso_nac.ppois.'
]

def extend(a, b, c):
    return np.ones(b-a)*c

for chip in chips:
    for chrom in chrs:
        signal = pd.read_csv('./inputs/SHSY5Y/prepped/'+chip+chrom+'.bdg',sep = '\t',names=['chr','start','end','signal'],header=None,index_col=None)
        signal = signal[['start','end','signal']].to_numpy()
        result = [extend(int(row[0]),int(row[1]),row[2]) for row in signal]
        result = np.hstack(result)
        np.save('/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/extended_signal/'+chip+'extended_signal.'+chrom+'.npy',result)


