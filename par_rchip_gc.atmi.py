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


def zw_bin_xy_stat_MODIFIED(array): # array[:,(z,x,w_y)]
    array = array[array[:,0].argsort()] # sorting by GC value
    zbin_no = len(array)//500_000
    zbins = np.linspace(0,zbin_no*500_000,zbin_no+1)
    zbins = np.digitize(np.linspace(0,len(array),len(array)),zbins)
    out_array = []
    for i in range(1,zbin_no+1):
        tmp = array[zbins==i] # tmp[zbin==i,(z,x,w_y)]
        tmp = tmp[tmp[:,2].argsort()] # sorting by RChIP value
        wbin_no = len(tmp)//25_000 # number of 1k sized bins
        wbins = np.linspace(0,wbin_no*25_000,wbin_no+1) # creating equally SIZED bins
        wbins = np.digitize(np.linspace(0,len(tmp),len(tmp)),wbins) # assigning bins based on index value
        u_z = np.array([np.mean(tmp[:,0])]*wbin_no)
        u_x = np.array([np.mean(tmp[(wbins==j),1]) for j in range(1,wbin_no+1)])
        #s_x = np.array([np.std(tmp[(wbins==j),1]) for j in range(1,wbin_no+1)])
        #m_x = np.array([np.median(tmp[(wbins==j),1]) for j in range(1,wbin_no+1)])
        #X_x = np.array([len(tmp[(wbins==j),1]) for j in range(1,wbin_no+1)])
        u_y = np.array([np.mean(tmp[(wbins==j),2]) for j in range(1,wbin_no+1)])
        #s_y = np.array([np.std(tmp[(wbins==j),2]) for j in range(1,wbin_no+1)])
        #m_y = np.array([np.median(tmp[(wbins==j),2]) for j in range(1,wbin_no+1)])
        #X_y = np.array([len(tmp[(wbins==j),2]) for j in range(1,wbin_no+1)])
        out_array.append(np.vstack((u_z,u_x,u_y)).T)
        #out_array.append(np.vstack((u_z,u_x,s_x,m_x,X_x,u_y,s_y,m_y,X_y)).T)
    out_array = np.vstack(out_array)
    return(out_array)



gc_abun = pd.DataFrame()
genome = []
for chrom in chrs:
    z_chr = np.load(var[0]+chrom+'.npy').astype('float32') # GC
    x_chr = np.load(var[2]+chrom+'.npy').astype('float32') # PAR
    y_chr = np.load(var[4]+chrom+'.npy').astype('float32') # RChIP
    locs = np.where((x_chr>0) & (y_chr>0) & (z_chr>0.45))
    z_chr = z_chr[locs]
    x_chr = x_chr[locs]
    y_chr = y_chr[locs]
    del locs
    D = constants[constants['signal_type']==var[2].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
    x_chr = (x_chr*D).astype('float32') # operating on signal
    D = constants[constants['signal_type']==var[4].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
    y_chr = (y_chr*D).astype('float32') # operating on signal
    genome.append(np.vstack((z_chr,x_chr,y_chr))) # appending the chromosome signal to the array
    del z_chr,x_chr,y_chr
genome = np.hstack(genome).T
bins = zw_bin_xy_stat_MODIFIED(genome) #
x_name = var[2].split('/')[-1] # ChIP name
y_name = var[4].split('/')[-1] # ChIP name
gc_abun['bin_gc_u'] = bins[:,0]
gc_abun[x_name+'global_gc.u'] = bins[:,1]
#gc_abun[x_name+'global_gc.s'] = bins[:,2]
#gc_abun[x_name+'global_gc.m'] = bins[:,3]
#gc_abun[x_name+'global_gc.X'] = bins[:,4]
gc_abun[y_name+'global_gc.u'] = bins[:,2]
#gc_abun[y_name+'global_gc.s'] = bins[:,6]
#gc_abun[y_name+'global_gc.m'] = bins[:,7]
#gc_abun[y_name+'global_gc.X'] = bins[:,8]
gc_abun.to_csv('par_x.rchip_wy.gc_z.atmi.figures.500K_25K_0.45GC.NAs.txt',sep='\t',index=False)
gc_abun.dropna().to_csv('par_x.rchip_wy.gc_z.atmi.figures.500K_25K_0.45GC.txt',sep='\t',index=False)




