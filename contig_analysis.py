#!/usr/bin/env python3.5
# Plots the distribution of contig sizes from a fasta file
# and analyzes percentages
# JP Tomkins July 6, 2018

import matplotlib.pyplot as plt
import seaborn as sns

def get_fastadata(fastafile):
    '''Makes a list of tuples, each with a list containing fasta
    header data and a string containing DNA seq ([info], 'seq')'''
    with open(fastafile, 'r') as fi:
        return [(part[0],
        part[2].replace('\n', ''))
        for part in
        [entry.partition('\n')
        for entry in fi.read().split('>')[1:]]]

seqdata = get_fastadata("contigs_all.fa")

#  Array of seq lengths for seqdata
seq_lens = [len(x[1]) for x in seqdata]

seqs_under_50k = [x for x in seq_lens if x < 50000]
seqs_under_250k = [x for x in seq_lens if x < 250000]
seqs_over_250k = [x for x in seq_lens if x > 250000]
seqs_over_300k = [x for x in seq_lens if x > 300000]
seqs_over_400k = [x for x in seq_lens if x > 400000]

print("seqs_under_50k: ", round(len(seqs_under_50k)/18000*100, 2),"%")
print("seqs_under_250k: ", round(len(seqs_under_250k)/18000*100, 2),"%")
print("seqs_over_250k: ", round(len(seqs_over_250k)/18000*100, 2),"%")
print("seqs_over_300k: ", round(len(seqs_over_300k)/18000*100, 2),"%")
print("seqs_over_300k: ", len(seqs_over_300k))
print("seqs_over_400k: ", len(seqs_over_400k))

# Make plot
sns.kdeplot(data=seq_lens, shade=True, legend=False)
plt.show()

