#!/usr/bin/env python
# coding: utf-8
"""
Created on Tue Dec 14 19:32:12 2021

@author: James.Pettengill
"""

from Bio import SeqIO
import pandas as pd
import subprocess
from io import BytesIO
import argparse
import sys
import os

def get_results(in_file, primer_amplicon, id_file, number_reads, read_length, fragment_length):
    ### Get reference genomes out of in_file
    print("Getting reference ids and proportions")
    ref_ids = {}
    with open(id_file, "r") as infile:
        for line in infile:
            ref_ids[line.split("\t")[0]] = line.split("\t")[1].replace("\n", "")
     
    for ref_id, proportion in ref_ids.items():
        with open(ref_id + ".fasta", "w") as f:
            fasta_sequences = SeqIO.parse(open(in_file),'fasta')
            for seq in fasta_sequences:
                if ref_id in seq.id:
                    SeqIO.write([seq], f, "fasta")
    
    
    ### Get primers in in_silico_pcr.pl format
    print("Reading primer tsv file")
    neb = pd.read_csv(primer_amplicon, header = 0, sep = "\t")
    neb['output_prefix'] = neb.name.replace("_LEFT.*|_RIGHT.*", "", regex = True)
    neb_left = neb[neb["name"].str.contains("_LEFT")]
    neb_left = neb_left[["name", "seq", "output_prefix"]]
    neb_right = neb[neb["name"].str.contains("_RIGHT")]
    neb_right = neb_right[["name", "seq", "output_prefix"]]
    neb_formatted = pd.merge(neb_left, neb_right, how = "inner", on = "output_prefix")
    neb_formatted = neb_formatted[["seq_x", "seq_y", "output_prefix"]]
    neb_formatted.to_csv(os.path.basename(primer_amplicon).split(".")[0] + "_for_in_silico_pcr.tsv", index = None, header = False, sep = "\t")
    
    ### Run in_silico_pcr.pl
    print("Running in_silico_pcr.pl")
    for ref_id, proportion in ref_ids.items():
        cmd = "in_silico_PCR.pl -p " + primer_amplicon.split(".")[0] + "_for_in_silico_pcr.tsv -s " + ref_id + ".fasta"
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        results = proc.stdout.read()
        try:
            tmp = pd.read_csv(BytesIO(results), sep='\t', header=0)
            tmp.to_csv(ref_id + "_" + os.path.basename(primer_amplicon).split(".")[0] + ".bed", index = None, sep = "\t")
        except:
            print("\n***ERROR****\nProblems reading in_silico_pcr.pl results - make sure it is installed")
            sys.exit(1)
    ### Generate amplicon fasta file from bed file
    bed_info = {}
    for ref_id, proportion in ref_ids.items():
        bed_info[ref_id] = []
        with open(ref_id + "_" + os.path.basename(primer_amplicon).split(".")[0] + ".bed", "r") as bed_file:
            for line in bed_file:
                if "AmpId" not in line and "No amplification" not in line:
                    name = line.split("\t")[0]
                    begin = line.split("\t")[2]
                    end = line.split("\t")[3].replace("\n", "")
                    bed_info[ref_id].append([name, begin, end])
    
    ### Write out a single file of amplicons adjusting for the proportion of each variant
    output_prefix = '_'.join('{}_{}'.format(key, value) for key, value in ref_ids.items())
    number_amplicons = []
    with open(output_prefix + "_" + "amplicons.fasta", "w") as outfile:
        for ref_id, values in bed_info.items():
            fasta_sequences = SeqIO.parse(open(ref_id + ".fasta"),'fasta')
            proportion = int(ref_ids[ref_id])
            for seq in fasta_sequences:
                if ref_id in seq.id:
                    for i in range(1, proportion):
                        for value in values:
                                ## Get count of amplicons to scale number of reads
                                number_amplicons.append(value)
                                amplicon = seq.seq[int(float(value[1])-1):int((float(value[1]) + float(value[2]))-1)]
                                outfile.write(">" + ref_id + "_" + value[0] + "\n" + str(amplicon) + "\n")
    number_amplicons = len(number_amplicons)
    ### Generate amplicon files for each genome
    for ref_id, values in bed_info.items():
        fasta_sequences = SeqIO.parse(open(ref_id + ".fasta"),'fasta')
        for seq in fasta_sequences:
            if ref_id in seq.id:
                with open(ref_id + "_" + primer_amplicon.split(".")[0] + "amplicons.fasta", "w") as outfile:
                    for value in values:
                        amplicon = seq.seq[int(float(value[1])-1):int((float(value[1]) + float(value[2]))-1)]
                        outfile.write(">" +  ref_id + "_" + value[0] + "\n" + str(amplicon) + "\n")
    
    
    # ### Run art and adjust coverage to account for total number of reads desired
    # input_id = output_prefix + "_" + "amplicons.fasta"
    # output_id = output_prefix + "art_out"
    # print("Running art for ", input_id)
    # try:
    #     subprocess.call(["art_illumina", "-ss", "MSv3", "-p", "-na", "-i", input_id, "-l", str(read_length), "-s", "10", "-m", str(fragment_length), "-f", str(round(int(number_reads)/int(number_amplicons))), "-o", output_id])
    # except: 
    #     print("\n***ERROR***\nProblems running art - make sure it is installed and the fasta input file isn't empty.")
    #     sys.exit(1)
    # ### Generate single art files per reference genome
    # for ref_id, proportion in ref_ids.items():
    #     input_id = ref_id + "_" + primer_amplicon.split(".")[0] + "amplicons.fasta"
    #     output_id = ref_id + "_" + primer_amplicon.split(".")[0] + "art_out"
    #     subprocess.call(["art_illumina", "-ss", "MSv3", "-p", "-na", "-i", input_id, "-l", str(read_length), "-s", "10", "-m", str(fragment_length),   "-f", str(1),  "-o", output_id])

def parse_arguments(system_args):
    """
    Parse command line arguments
    """
    usage = """Takes a SC2 reference fasta, tsv file with amplicon id, start position, and length (output from in_silico_pcr.pl); \
        tsv file with amplicon id and number of copies used to simulate differential coverage.  \
        Must have in_silico_pcr.pl and art installed and in in your path.

    Example usage:
    generate_reads.py -f relevantGISAID.fa -i cdc_epi_isl.txt -p neb_vss1a.primer.tsv """


    parser = argparse.ArgumentParser(description= usage, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-f", "--in-file", metavar="FILE", dest="in_file", default=None, required = True, help="Path to reference SC2 reference fastas.")
    parser.add_argument("-i", "--id-file", metavar="FILE", dest="id_file", default=None, required = True, help="Path to new line delimited list of epi isls in SC2 reference fastas.")
    parser.add_argument("-p", "--primer-amplicon", metavar="FILE", dest="primer_amplicon", default=None,  help="Path to file with of format https://raw.githubusercontent.com/nebiolabs/VarSkip/main/neb_vss1a.primer.tsv")
    parser.add_argument("-n", "--number-reads", metavar="INT", dest="number_reads", default=10000,  help="Total number of reads to generate via art.  Default is 10,000")
    parser.add_argument("-l", "--read-length", metavar="INT", dest="read_length", default=150,  help="Read length for art simulations.")
    parser.add_argument("-r", "--fragment-length", metavar="INT", dest="fragment_length", default=400,  help="Expected fragment length.")


    # parser.add_argument("-f", "--frequency-amplicons", metavar="FILE", dest="frequency_amplicons", default=None, help="Path to file with amplicon id and frequency")

    args = parser.parse_args(system_args)
    return args

def main(args):
    """
    Generate results
    """
    # All args
    in_file = args.in_file
    id_file = args.id_file
    primer_amplicon = args.primer_amplicon
    number_reads = args.number_reads
    fragment_length = args.fragment_length
    read_length = args.read_length
    # frequency_amplicons = args.frequency_amplicons

    get_results(in_file, primer_amplicon, id_file, number_reads, read_length, fragment_length)

if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)






