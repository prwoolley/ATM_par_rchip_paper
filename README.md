# ATM_par_rchip_paper
Repository containing Python scripts and RMarkdown document for ATM_par_rchip paper.

Figures 1 and 2 are found on "AT_paper.Swish_analysis.rmd"

Figure 6 uses the python scripts to process Bedgraph files.

"sample_prep.sh" separates the bedgraph files by chromosome.

"extend_signal.py" extends the chromosome bedgraph files so that each numpy row represents a nucleotide position on the genome.

"trim.py" adds extra 0's to the end of each numpy array so each array matches the length of the GC numpy arrays, which cover the entire genome.

"rolling_mean.py" takes trimmed numpy arrays and creates a rolling mean through each of them of window size 100.

"make_rna_seq_vec.py" converts bedgraph files of transcript abundance into numpy arrays that match the length of the genome.

"make_annotation_vec.py" marks features on the genome that are linked to genes with RNA seq transcript abundance.

"make_gc_vec.py" converts the Hg38 fasta sequence into numpy arrays (one per chromosome) with 0 for A/T and 1 for G/C.

"tss_gc.py" prepares data for plotting ChIP against promoter/TSS GC content.

"tss_abun.py" prepares data for plotting ChIP against promoter/TSS transcript abundance.

"rchip_n_par.dmso.py" prepares data for plotting locations where PAR and RChIP intersect and locations where they don't.

"rchip_n_par.atmi.py" same as above, but for ATMi treated cells.

"par_rchip_gc.dmso.py" prepares data for plotting PAR signal against RChIP and GC content

"par_rchip_gc.atmi.py" same as above, but for ATMi treated cells.
