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

def z_bin_x_stat(array,bin_no,bin_range): # array[:,(z,x)], bin_range:[bin_min,bin_max]
    if bin_no == False: # Useful for Abundance plots, unnecessary for GC plots
        array = array[array[:, 0].argsort()] # sort array by z vector
        bin_no = len(array)//1_000_000 # Determining bin no
        bins = np.linspace(0,bin_no*1_000_000,bin_no+1) # creating equally SIZED bins
        bins = np.digitize(np.linspace(0,len(array),len(array)),bins) # assigning bins based on index value
        u = np.array([np.mean(array[bins==i,1]) for i in range(1,bin_no+1)]) # mean x value in each bin
#        s = np.array([np.std(array[bins==i,1]) for i in range(1,bin_no+1)]) # standard deviation of x in each bin
#        m = np.array([np.median(array[bins==i,1]) for i in range(1,bin_no+1)])
#        X = np.array([len(array[bins==i,1]) for i in range(1,bin_no+1)])
        u_z = np.array([np.mean(array[bins==i,0]) for i in range(1,bin_no+1)]) # mean z value in each bin
        return(np.vstack((u_z, u)).T)
        #return(np.vstack((u_z, u,s,m,X)).T)
    elif bin_no:
        bins = np.linspace(bin_range[0],bin_range[1],bin_no+1) # creating equally spaced bins of a set number from min(x) to max(x)
        array[:,0] = np.digitize(array[:,0],bins) # assigning bins to x array
        u = np.array([np.mean(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        s = np.array([np.std(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        m = np.array([np.median(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        X = np.array([len(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        return(np.vstack((u,s,m,X)).T)

def z_bin_xy_stat(array,bin_no,bin_range): # array[:,(z,x,y)]
    if bin_no == False: # Useful for Abundance plots, unnecessary for GC plots
        array = array[array[:, 0].argsort()] # sort array by z vector
        bin_no = len(array)//100_000 # Determining bin no
        bins = np.linspace(0,bin_no*100_000,bin_no+1) # creating equally SIZED bins
        bins = np.digitize(np.linspace(0,len(array),len(array)),bins) # assigning bins based on index value
        u_z = np.array([np.mean(array[bins==i,0]) for i in range(1,bin_no+1)]) # mean z value in each bin
        u_x = np.array([np.mean(array[bins==i,1]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        #s_x = np.array([np.std(array[bins==i,1]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        #m_x = np.array([np.median(array[bins==i,1]) for i in range(1,bin_no+1)])
        #X_x = np.array([len(array[bins==i,1]) for i in range(1,bin_no+1)])
        u_y = np.array([np.mean(array[bins==i,2]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        #s_y = np.array([np.std(array[bins==i,2]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        #m_y = np.array([np.median(array[bins==i,2]) for i in range(1,bin_no+1)])
        #X_y = np.array([len(array[bins==i,2]) for i in range(1,bin_no+1)])
        #return(np.vstack((u_z,u_x,s_x,m_x,X_x,u_y,s_y,m_y,X_y)).T)
        return(np.vstack((u_z,u_x,u_y)).T)

    elif bin_no:
        bins = np.linspace(bin_range[0],bin_range[1],bin_no+1) # creating equally spaced bins of a set number from min(x) to max(x)
        array[:,0] = np.digitize(array[:,0],bins) # assigning bins to x array
        u_x = np.array([np.mean(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        s_x = np.array([np.std(array[array[:,0]==i,1]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        m_x = np.array([np.median(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        X_x = np.array([len(array[array[:,0]==i,1]) for i in range(1,bin_no+1)])
        u_y = np.array([np.mean(array[array[:,0]==i,2]) for i in range(1,bin_no+1)]) # mean y value in each x bin
        s_y = np.array([np.std(array[array[:,0]==i,2]) for i in range(1,bin_no+1)]) # standard deviation of y in each x bin
        m_y = np.array([np.median(array[array[:,0]==i,2]) for i in range(1,bin_no+1)])
        X_y = np.array([len(array[array[:,0]==i,2]) for i in range(1,bin_no+1)])
        return(np.vstack((u_x,s_x,m_x,X_x,u_y,s_y,m_y,X_y)).T)


#gc_abun = pd.DataFrame()
## PAR and RChIP where both are present
#genome = []
#for chrom in chrs:
#    z_chr = np.load(var[0]+chrom+'.npy').astype('float32') # GC
#    x_chr = np.load(var[3]+chrom+'.npy').astype('float32') # PAR ChIP
#    y_chr = np.load(var[5]+chrom+'.npy').astype('intc') # RChIP
#    locs = np.where((y_chr>0) & (x_chr>0))
#    z_chr = z_chr[locs]
#    x_chr = x_chr[locs]
#    y_chr = y_chr[locs]
#    del locs
#    D = constants[constants['signal_type']==var[3].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
#    x_chr = (x_chr*D).astype('float32') # operating on signal
#    D = constants[constants['signal_type']==var[5].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
#    y_chr = (y_chr*D).astype('float32') # operating on signal
#    genome.append(np.vstack((z_chr,x_chr,y_chr))) # appending the chromosome signal to the array
#    del z_chr,x_chr,y_chr
#genome = np.hstack(genome).T
#bins = z_bin_xy_stat(genome,False,False) # Bins of equal size, not spacing
#gc_abun['par_n_rchip.gc_bin.global_gc.u'] = bins[:,0]
#gc_abun['par_n_rchip.par_dmso.global_gc.u'] = bins[:,1]
##gc_abun['par_n_rchip.par_dmso.global_gc.s'] = bins[:,2]
##gc_abun['par_n_rchip.par_dmso.global_gc.m'] = bins[:,3]
##gc_abun['par_n_rchip.par_dmso.global_gc.X'] = bins[:,4]
#gc_abun['par_n_rchip.rchip_dmso.global_gc.u'] = bins[:,2]
##gc_abun['par_n_rchip.rchip_dmso.global_gc.s'] = bins[:,6]
##gc_abun['par_n_rchip.rchip_dmso.global_gc.m'] = bins[:,7]
##gc_abun['par_n_rchip.rchip_dmso.global_gc.X'] = bins[:,8]
#gc_abun.to_csv('./outputs/rchip_n_par/par_n_rchip.dmso.gc.figures.100K.NAs.txt',sep='\t',index=False)
#gc_abun.dropna().to_csv('./outputs/rchip_n_par/par_n_rchip.dmso.gc.figures.100K.txt',sep='\t',index=False)

# PAR where no RChIP is present
gc_abun = pd.DataFrame()
genome = []
for chrom in chrs:
    z_chr = np.load(var[0]+chrom+'.npy').astype('float32') # GC
    x_chr = np.load(var[3]+chrom+'.npy').astype('float32') # PAR ChIP
    y_chr = np.load(var[5]+chrom+'.npy').astype('intc') # RChIP
    locs = np.where((y_chr==0) & (x_chr>0))
    z_chr = z_chr[locs]
    x_chr = x_chr[locs]
    del locs,y_chr
    D = constants[constants['signal_type']==var[3].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
    x_chr = (x_chr*D).astype('float32') # operating on signal
    genome.append(np.vstack((z_chr,x_chr))) # appending the chromosome signal to the array
    del z_chr,x_chr
genome = np.hstack(genome).T
bins = z_bin_x_stat(genome,False,False) # Bins of equal size, not spacing
gc_abun['par_n_rchipC.gc_bin.global_gc.u'] = bins[:,0]
gc_abun['par_n_rchipC.par_dmso.global_gc.u'] = bins[:,1]
#gc_abun['par_n_rchipC.par_dmso.global_gc.s'] = bins[:,2]
#gc_abun['par_n_rchipC.par_dmso.global_gc.m'] = bins[:,3]
#gc_abun['par_n_rchipC.par_dmso.global_gc.X'] = bins[:,4]
gc_abun.to_csv('./outputs/rchip_n_par/par_n_rchipC.dmso.gc.figures.1mil.NAs.txt',sep='\t',index=False)
gc_abun.dropna().to_csv('./outputs/rchip_n_par/par_n_rchipC.dmso.gc.figures.1mil.txt',sep='\t',index=False)

## RChIP where no PAR is present
#gc_abun = pd.DataFrame()
#genome = []
#for chrom in chrs:
#    z_chr = np.load(var[0]+chrom+'.npy').astype('float32') # GC
#    x_chr = np.load(var[5]+chrom+'.npy').astype('float32') # RChIP
#    y_chr = np.load(var[3]+chrom+'.npy').astype('intc') # PAR
#    locs = np.where((y_chr==0) & (x_chr>0))
#    z_chr = z_chr[locs]
#    x_chr = x_chr[locs]
#    del locs,y_chr
#    D = constants[constants['signal_type']==var[5].split('/')[-1].split('ppois')[0]+'ppois']['D_depth_normalization_factor'].item() # getting operation factor
#    x_chr = (x_chr*D).astype('float32') # operating on signal
#    genome.append(np.vstack((z_chr,x_chr))) # appending the chromosome signal to the array
#    del z_chr,x_chr
#genome = np.hstack(genome).T
#bins = z_bin_x_stat(genome,False,False) # Bins of equal size, not spacing
#gc_abun['rchip_n_parC.gc_bin.global_gc.u'] = bins[:,0]
#gc_abun['rchip_n_parC.rchip_dmso.global_gc.u'] = bins[:,1]
##gc_abun['rchip_n_parC.rchip_dmso.global_gc.s'] = bins[:,2]
##gc_abun['rchip_n_parC.rchip_dmso.global_gc.m'] = bins[:,3]
##gc_abun['rchip_n_parC.rchip_dmso.global_gc.X'] = bins[:,4]
#gc_abun.to_csv('./outputs/rchip_n_par/rchip_n_parC.dmso.gc.figures.1mil.NAs.txt',sep='\t',index=False)
#gc_abun.dropna().to_csv('./outputs/rchip_n_par/rchip_n_parC.dmso.gc.figures.1mil.txt',sep='\t',index=False)

