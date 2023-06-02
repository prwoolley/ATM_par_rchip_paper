# Library imports
import numpy as np
import pandas as pd


chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

# Inputs are from signal numpy arrays
chips = [
'shsy5y.par.atmi.ppois.',
'shsy5y.par.dmso.ppois.',
'shsy5y.rchip.atmi.ppois.',
'shsy5y.rchip.dmso.ppois.',
'shsy5y.atac.SRR5819663_SRR5819664_treat_pileup.',
'shsy5y.ctcf.SRR6334830.ppois.',
'shsy5y.h3k27ac.SRR3363257.ppois.',
'shsy5y.h3k27me3.SRR3363258.ppois.',
'shsy5y.h3k4me1.SRR3363256.ppois.',
'shsy5y.h3k4me3.SRR3363255.ppois.',
'shsy5y.igg.SRR4291521_SRR4291522.ppois.',
'shsy5y.pol2.SRR6322542_SRR6322543.ppois.',
'shsy5y.s96.SRR4291525_SRR4291526.ppois.',
'shsy5y.rchip.atmi_nac.ppois.',
'shsy5y.rchip.dmso_nac.ppois.'
]

# Rolling mean over chip signal
for chip in chips:
    for chrom in chrs:
        dna_vec = np.load('/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/trimmed/'+chip+'trimmed.'+chrom+'.npy')
        dna_vec = pd.DataFrame(dna_vec, columns = ['vec'])
        dna_vec['vec'] = dna_vec['vec'].rolling(100).mean()
        dna_vec = dna_vec.fillna(0)
        dna_vec = dna_vec['vec'].to_numpy()
        np.save('/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/'+chip+'rolling_mean_100.'+chrom+'.npy',dna_vec)

