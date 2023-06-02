# Script calculates the "effective GC content" of a DNA string using the method described in DOI: 10.1101/gr.220673.117

# General workflow is to load the DNA sequence, convert the AT/CG to 0/1 in a numpy array, and multiply each position by a weight vector.
# The first portion converts chromosome 1 into a np.array(), which is time consuming since it uses a generic for loop.
# Once it is run, the np.array() is saved and can be reloaded, skipping the step.

# Library imports
import numpy as np
import pandas as pd
#from os import listdir
#import swifter


chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
for chrom in chrs:
    # Loading the DNA sequence:
    with open('/home1/05664/pwoolley/work/FASTA/Ensembl/Homo_sapiens.GRCh38.dna.chromosome.'+chrom+'.fa','r') as chrom_file:
        chrom_file.readline() # readline since the first line isn't the sequence
        dna = chrom_file.read().replace('\n','') # save every line after the first as the chrom_seq variable, removing all new line characters

    # Turning sequence string into a vector
    gc = ['G','g','C','c']
    # Running on the entire chromosome
    dna_list = []
    for j in dna:
        dna_list.append(1 if j in gc else 0)
    dna_vec = np.array(dna_list)
    np.save('./dna_vec/gc_vec.chr'+chrom+'.npy', dna_vec)

