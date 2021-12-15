#!/usr/bin/env python
# coding: utf-8
"""
Created on Tue Dec 14 19:32:12 2021

@author: James.Pettengill
"""
import random
import argparse
import sys

def generate_replicates(in_file, num_reps, out_file):
    variants = []
    with open(in_file, "r") as infile:
        for line in infile:
            variants.append(line.replace("\n", ""))
    max_variants = len(variants)
    res = {}
    for i in range(0, int(num_reps)):
        res[i] = []
        samples = [0]
        for j in range(0, max_variants):
            if j < (max_variants + 1):
                # print(sum(samples))
                a = random.randint(0, (100 - sum(samples)))
                samples.append(a)
                res[i].append([variants[j], a])
        remainder = 100 - sum(samples)
        list_index = random.randint(0, max_variants - 1)
        res[i][list_index][1] = res[i][list_index][1] + remainder
    
    for key, value in res.items():
       with open(out_file + "_replicate_" + str(key) + ".tsv", "w") as output_file:
           # output_file.write("Variant\tProportion\n")
           for val in value:
                output_file.write(str(val[0]) + "\t" + str(val[1]) + "\n")
    
def parse_arguments(system_args):
    """
    Parse command line arguments
    """
    usage = """Takes an integer denoting the number of variants to model and another integer denoting
    the number of replicates.

    Example usage:
    generate_replicates.py -i epi_isls.txt -n 10 -o replicates_5_variants_10_reps.tsv """


    parser = argparse.ArgumentParser(description= usage, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", "--in-file", metavar="INT", dest="in_file", default=None, required = True, help="File with new line delimited list of epi isl ids.")
    parser.add_argument("-n", "--num-reps", metavar="INT", dest="num_reps", default=None, required = True, help="Number of replicates to generate.")
    parser.add_argument("-o", "--out-file", metavar="STR", dest="out_file", default=None,  help="Prefix for output file - replicate number will be appended.")

    args = parser.parse_args(system_args)
    return args

def main(args):
    """
    Generate results
    """
    # All args
    in_file = args.in_file
    num_reps = args.num_reps
    out_file = args.out_file

    generate_replicates(in_file, num_reps, out_file)

if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)

    