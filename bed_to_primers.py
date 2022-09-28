#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse

def get_primers(reference, bed_file):
    df_bed = pd.read_csv(bed_file, header = None, names = ["chrom", "chrom_start", "chrom_end", "name", "score", "strand", "primer"], sep = "\t")
    df_bed["primer_seq"] = ""
    df_bed["size"] = ""
    df_bed["%gc"] = "UNK"
    df_bed["tm"] = "UNK"
    sequence = []
    fasta_sequences = SeqIO.parse(open(reference),'fasta')
    for seq in fasta_sequences:
        sequence.append(str(seq.seq))

    for i in range(0, df_bed.shape[0]):
        if "LEFT" in df_bed["name"][i]:
            primer_seq = sequence[0][df_bed["chrom_start"][i]:df_bed["chrom_end"][i]]
            size = df_bed["chrom_end"][i] - df_bed["chrom_start"][i]
            df_bed["primer"][i] = primer_seq
            df_bed["size"][i] = size
        if "RIGHT" in df_bed["name"][i]:
            primer_seq = sequence[0][df_bed["chrom_start"][i]:df_bed["chrom_end"][i]]
            primer_seq = Seq(primer_seq).reverse_complement()
            size = df_bed["chrom_end"][i] - df_bed["chrom_start"][i]
            df_bed["primer"][i] = primer_seq
            df_bed["size"][i] = size

    df_bed = df_bed[["name", "chrom", "primer", "size", "%gc", "tm"]]
    df_bed.columns = [["name", "pool", "seq", "size", "%gc", "tm"]]
    df_bed.to_csv(bed_file.replace("bed", "tsv"), index = None, sep = "\t", na_rep = "UNK")

    ### Write out fasta file of primers
    with open(bed_file.replace("bed", "fasta"), "w") as out_fasta:
        for i in range(0,df_bed.shape[0]):
            out_fasta.write(">" + df_bed.loc[i][0] + "\n" + str(df_bed.loc[i][2]) + "\n")

def parse_arguments(system_args):
    """
    Parse command line arguments
    """
    usage = """Takes a bed file and associated reference and creates primer sequence files.

    Example usage:
    bed_to_primers.py -i neb_vss1a_primer.bed -r NC_045512.2_sequence.fasta """


    parser = argparse.ArgumentParser(description= usage, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", "--bed_file", metavar="STR", dest="bed_file", default=None, required = True, help="bed file.")
    parser.add_argument("-r", "--reference", metavar="STR", dest="reference", default=None, required = True, help="Name of reference fasta file")

    args = parser.parse_args(system_args)
    return args

def main(args):
    """
    Generate results
    """
    # All args
    bed_file = args.bed_file
    reference = args.reference

    get_primers(reference, bed_file)

if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)
