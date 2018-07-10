#!/usr/bin/env python3.5

'''Summarizes the results of a blastn output and can be run
as a batch script to analyze each of an entire set of blastn
output files in csv format. Also generates a summary of the
entire set of blastn output files.'''

from os import listdir

import pandas as pd


def get_dataframe(csvfile):
    '''Puts a csv file into a pandas dataframe'''
    try:
        with open(csvfile, 'r') as fi:
            # Open *csv as a pandas data frame
            return pd.read_csv(fi)
    except IOError:
        print('Cannot open', csvfile)


def set_csv_header():
    '''Set the csv data file header line for blastn results.'''
    try:
        with open(csv_outfilename, 'w') as fo:
            header = ['file_id', 'aln_ident', 'qseq_ident', 'aln_len',
                      'qseqret', 'qseqall', 'num_qseqs', 'num_hits',
                      'hitfreq', 'overall_ident']
            fo.write(",".join(header) + "\n")
    except IOError:
            print('Cannot open', csv_outfilename)


def get_blast_data(csv_filename, fasta_filename):
    '''Takes a blastn output file in csv format along with it's
    corresponding fasta file and generates a set of basic stats
    for the blast results and the sequences used to generate them.'''

    def _get_fastadata(fastafile):
        '''Makes a list of tuples, each with a list containing fasta
        header data and a string containing DNA seq ([info], 'seq')'''
        try:
            with open(fastafile, 'r') as fi:
                return [(part[0],
                        part[2].replace('\n', ''))
                        for part in
                        [entry.partition('\n')
                        for entry in fi.read().split('>')[1:]]]
        except IOError:
            print('Cannot open', fasta_file)

    # Get and setup the pandas dataframe
    df = get_dataframe(csv_filename)
    seqdata = _get_fastadata(fasta_filename)
    file_id = "pan" + csv_filename[4:7]

    num_hits = df.count()[0]
    num_qseqs = len(seqdata)

    # Assign column headings - replaces default [0,1,2...etc]
    df.columns = \
        'qseqid qstart qend mismatch gapopen pident nident length qlen'.split()

    # Speeds calc and removes this col as an obj dtype
    del df['qseqid']

    # Get overall qseq aln identity and put in new column called 'ident'
    df['ident'] = df.apply(lambda row: (row['nident']/row['length']
                           if row['length'] > row['qlen']
                           else row['nident']/row['qlen']), axis=1)

    # Create variables for outputs
    total_aligned = (df['ident'].mean() * df['qlen'].mean()) * num_hits
    total_seq_len = sum(len(s[1]) for s in (seqdata))
    ave_seqlen = ((total_seq_len/num_qseqs) * 100, 2)
    ave_aln_ident = df['pident'].mean()
    med_aln_ident = df['pident'].median()
    min_aln_ident = df['pident'].min()
    max_aln_ident = df['pident'].max()
    ave_aln_len = df['length'].mean()
    med_aln_len = df['length'].median()
    min_aln_len = df['length'].min()
    max_aln_len = df['length'].max()
    ave_qseqret = df['qlen'].mean()
    perc_aln = (ave_aln_len/ave_qseqret) * 100
    ave_qseqall = total_seq_len/num_qseqs
    ave_hitfreq = num_hits/num_qseqs * 100
    ave_qseq_ident = df['ident'].mean() * 100

    # Send a human readable report to standard out
    print("Ave aln ident   : ", round(ave_aln_ident, 2))
    print("Median aln ident: ", round(med_aln_ident, 2))
    print("Min aln ident   : ", round(min_aln_ident, 2))
    print("Max aln ident   : ", round(max_aln_ident, 2))
    print("Ave aln len     : ", round(ave_aln_len, 2))
    print("Min aln len     : ", round(min_aln_len, 2))
    print("Max aln len     : ", round(max_aln_len, 2))
    print("Ave Perc aln    : ", round(perc_aln, 2))
    print("Median aln len  : ", med_aln_len)
    print("Ave qseqret len : ", round(ave_qseqret, 2))
    print("Ave qseqall len : ", round(ave_qseqall, 2))
    print("Num queryseqs   : ", num_qseqs)
    print("Num queryhits   : ", num_hits)
    print("Ave hit freq    : ", round(ave_hitfreq, 2))
    print("Ave qseq ident  : ", round(ave_qseq_ident, 2))

if __name__ == '__main__':
    '''Run on all blastn csv and fasta
    files in the working directory.'''

    csv_files = sorted([file for file in listdir('.')
                        if file.endswith('csv')])
    if not csv_files:
        print("No files with *.csv found!")

    fasta_files = sorted([file for file in listdir('.')
                          if file.endswith('fa')])
    if not fasta_files:
        print("No files with *.fa found!")

    for csv_filename, fa_filename in zip(csv_files, fasta_files):
        print("<{}>".format(csv_filename))
        print("<{}>".format(fa_filename))
        get_blast_data(csv_filename, fa_filename)
        print()
