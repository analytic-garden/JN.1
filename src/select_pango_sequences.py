#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
select_pango_sequences.py - select sequnecs from GIID fSAasta file
author: Bill Thompson
license: GPL 3
copyright: 2024-01-24
"""
import argparse
import sys
import gzip
from Bio import SeqIO
import pandas as pd
from bisect import bisect_left

def GetArgs() -> argparse.Namespace:
    """Get commandline arguments

    Returns
    -------
    argparse.Namespace
        An argparse record.
    """
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description = 'Select sequences from GISAID fasta file.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input fasta file (required). A gzipped fasta file from GISAID',
                            type = str)
        parser.add_argument('-m', '--meta_file',
                            required = True,
                            help = 'Input meta file (required). A tab delimited meta file for the pango variant.',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = True,
                            help = 'Output fasta file (required). Selected fasta sequences.',
                            type = str)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def select_sequences(fasta_file: str,
                     output_file: str,
                     fasta_ids: list[str]) -> None:
    """Seqrch a gzipped fasta file for records contained in fasta_ids. Matching records are written to output_file.

    Parameters
    ----------
    fasta_file : str
        A gzipped fasta file downloaded from GISAID.
    output_file : str
        An output fasta file containing matching records.
    fasta_ids : list[str]
        A list of fasta ids.
        
    Note
    ----
    GISAID fasta headers may have spaces in the country names. SeqIO parse breaks the id on whitespace. 
    Use description istead to match
    """
    
    def BinarySearch(target: str) -> int:
        """A helper function. Binary seach of the fasta id list.

        Parameters
        ----------
        target : str
            The ID to be found.

        Returns
        -------
        int
            The position in the fasta_ids list or -1 if not found.
        """
        i = bisect_left(fasta_ids, target)
        if i != len(fasta_ids) and fasta_ids[i] == target:
            return i
        else:
            return -1
    
    f_out = open(output_file, 'w')
            
    with gzip.open(fasta_file, 'rt') as f:
        count = 0
        for seq_rec in SeqIO.parse(f, 'fasta'):
            if BinarySearch(seq_rec.description) != -1:   # SeqIO uses the whole header as descripton.
                SeqIO.write(seq_rec, f_out, 'fasta')
            count += 1
            if count % 10000 == 0:
                sys.stdout.write("\r %d" %  count)
                sys.stdout.flush()
            
    f_out.close()
            
def main():
    args = GetArgs()
    
    # construct bfasta ID like GISAID.
    df = pd.read_csv(args.meta_file, sep = '\t')
    pango_ids = sorted(['|'.join([x, y, z]) 
                    for x,y,z in zip(list(df['Virus.name']), list(df['Collection.date']), list(df['Submission.date']))
                        if isinstance(y, str) and isinstance(z, str)])
    
    select_sequences(args.input_file, args.output_file, pango_ids)
    print()
    
if __name__ == "__main__":
    main()
