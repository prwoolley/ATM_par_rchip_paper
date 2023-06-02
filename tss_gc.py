# Script with several functions that collects statistics for x (and y) arrays based on bin values for z (and w)

import numpy as np
import pandas as pd

chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

constants = pd.read_csv('/home1/05664/pwoolley/work/scripts/bedgraph.standardize/standardize_constants.txt',sep='\t')
var = [
'/home1/05664/pwoolley/work/scripts/chip.smoothing/dna_vec/rolling_mean_100/gc_vec.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/sh.august2022.gene.abundance.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.par.atmi.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.par.dmso.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.rchip.atmi.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.rchip.dmso.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.rchip.atmi_nac.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rolling_mean_100/shsy5y.rchip.dmso_nac.ppois.rolling_mean_100.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/sh.august2022.gene.tss.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/sh.august2022.gene.tts.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/sh.august2022.gene.promoter.',
'/scratch/05664/pwoolley/SHSY5Y_ChIP/chip.smoothing/rna_seq_vec/sh.august2022.gene.terminator.'
]

def z_bin_x_stat(array,bin_no,range_bin): # array[:,(z,x)]
    if bin_no == False:
        array = array[array[:, 0].argsort()] # sort array by z vector
        bin_no = len(array)//250_000 # Determining bin no
        bins = np.linspace(0,bin_no*250_000,bin_no+1) # creating equally SIZED bins
        bins = np.digitize(np.linspace(0,len(array),len(array)),bins) # assigning bins based on index value
        u = np.array([np.mean(array[bins==i,1]) for i in range(1,bin_no+1)]) # mean x value in each bin
        s = np.array([np.std(array[bins==i,1]) for i in range(1,bin_no+1)]) # standard deviation of x in each bin
        m = np.array([np.median(array[bins==i,1]) for i in range(1,bin_no+1)])
        X = np.array([len(array[bins==i,1]) for i in range(1,bin_no+1)])
        u_z = np.array([np.mean(array[bins==i,0]) for i in range(1,bin_no+1)]) # mean z value in each bin
        return(np.vstack((u_z, u,s,m,X)).T)
    elif bin_no:
        bins = np.linspace(range_bin[0],range_bin[1],bin_no+1) # creating equally spaced bins of a set number from min(x) to max(x)
        array[:,0] = np.digitize(array[:,0],bins) # assigning bins to x array
        u = np.array([np.mean(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        s = np.array([np.std(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        m = np.array([np.median(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        X = np.array([len(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        return(np.vstack((u,s,m,X)).T)


for i in range(2,8):
    genome = []
    for chrom in chrs:
        z_chr = np.load(var[0]+chrom+'.npy').astype('float32') # GC
        x_chr = np.load(var[i]+chrom+'.npy').astype('float32') # ChIP
        tss = np.load(var[8]+chrom+'.npy').astype('intc') # TSS
        locs = np.where((tss>0) & (x_chr>0))
        del tss
        z_chr = z_chr[locs]
        x_chr = x_chr[locs]
        del locs
        D = constants[constants['signal_type']==var[i].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
        x_chr = (x_chr*D).astype('float32') # operating on signal
        genome.append(np.vstack((z_chr,x_chr))) # appending the chromosome signal to the array
        del z_chr,x_chr
    genome = np.hstack(genome).T
    bins = z_bin_x_stat(genome,False,False) # Bins of equal size, not spacing
    x_name = var[i].split('/')[-1] # ChIP name
    gc_abun = pd.DataFrame()
    gc_abun[x_name+'tss_gc_bin.u'] = bins[:,0]
    gc_abun[x_name+'tss_gc.u'] = bins[:,1]
    gc_abun[x_name+'tss_gc.s'] = bins[:,2]
    gc_abun[x_name+'tss_gc.m'] = bins[:,3]
    gc_abun[x_name+'tss_gc.X'] = bins[:,4]
    gc_abun.to_csv('./outputs/tss_gc/'+x_name+'tss_gc.250K.NAs.txt',sep='\t',index=False)
    gc_abun.dropna().to_csv('./outputs/tss_gc/'+x_name+'tss_gc.250K.txt',sep='\t',index=False)

