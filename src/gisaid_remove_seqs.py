#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gisaid_remove_seqs.py - remove short ands sequences with too man N's from fasta file.
author: Bill Thompson
license: GPL 3
copyright: 2024-02-01
"""
import argparse
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import TextIOWrapper

def GetArgs() ->  argparse.Namespace:
    """Get commandline arguments

    Returns
    -------
    argparse.Namespace
        An argparse record.
    """
 
    def ParseArgs(parser: argparse.ArgumentParser) -> argparse.Namespace:
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)
                
        parser = Parser(description='Remove short sequences and sequences with too many N\'s in a row.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input file (required). Fasta file from GISAID.',
                            type = str)
        parser.add_argument('-m', '--meta_file',
                            required = True,
                            help = ' '.join(['Meta file (required). GISAID meta file corresponding to FASTA file.',
                                              'Column names must be reformatted with R .name_repair = "universal"']),
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = False,
                            help = 'Output fasta file. Default: write to screen.',
                            type = str)
        parser.add_argument('-r', '--replace',
                            required = False,
                            default = False,
                            help = "Replace U's with T's in sequence.",
                            action = 'store_true')
        parser.add_argument('-n', '--max_N_pct',
                            required = False,
                            default = 20,
                            help = "Maximum percent of allowed number of N's. (default = 0.1, 0 to ignore)",
                            type = float)
        parser.add_argument('-l', '--min_seq_length',
                            required = False,
                            default = 0,
                            help = "Minimum sequence length. Sequences shorter than this are removed. (default = 0, 0 to ignore)",
                            type = int)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def process_seqs(seq_recs: SeqIO.FastaIO.FastaIterator,
                 df: pd.DataFrame,
                 out_handle: TextIOWrapper,
                 max_n_content: float,
                 min_seq_length: int) -> None:
    """Remove sequences with too many N's and short sequences.

    Parameters
    ----------
    seq_recs : SeqIO.FastaIO.FastaIterator
        An iterator for a fasta file.
    df : pd.DataFrame
        A dataframe cretaed by R, read_tsv_chunked.
    out_handle : TextIOWrapper
        Output fasta file handle
    max_n_content : float
        Max N content fraction.
    min_seq_length : int
        Minimum sequence length.
    """
    for seq_rec in seq_recs:
        items1 = seq_rec.description.split('|')
        virus_name = items1[0]

        df2 = df[df['Virus.name'] == virus_name]
        if df2.shape[0] == 0:
            print('Missing id:', virus_name, file = sys.stderr)
            continue
        elif np.isnan(df2['Is.complete.'].values[0]):
            print('Incomplete id:', virus_name, file = sys.stderr)
            continue
        elif np.isnan(df2['N.Content'].values[0]):
            print('Missing N content:', virus_name, file = sys.stderr)
            continue
        elif df2['N.Content'].values[0] > max_n_content:
            print('High N content:', virus_name, df2['N.Content'].values[0], file = sys.stderr)
            continue
        elif df2['Sequence.length'].values[0] < min_seq_length:
            print('Short sequence:', virus_name, df2['Sequence.length'].values[0], file = sys.stderr)
            continue
        else:
            location_str = '|'.join(df2['Location'].values[0].split(' / '))
            new_desc = '|'.join([virus_name,
                                 str(df2['Accession.ID'].values[0]),
                                 str(df2['Collection.date'].values[0]),
                                 str(df2['Pango.lineage'].values[0]),
                                 location_str])
            new_id = str(df2['Accession.ID'].values[0])
            new_rec = SeqRecord(id = new_id,
                                description = new_desc,
                                seq = seq_rec.seq)

            SeqIO.write(new_rec, out_handle, 'fasta')       
        
def main():
    args = GetArgs()
    
    if args.output_file is not None:
        f = open(args.output_file, 'w')
    else:
        f = sys.stdout
    
    rec_iterator = SeqIO.parse(args.input_file, 'fasta')
    
    df = pd.read_csv(args.meta_file, 
                     sep = '\t')
    
    process_seqs(rec_iterator, df, f, args.max_N_pct, args.min_seq_length)
    f.close()
    
    print()

if __name__ == "__main__":
    main()
