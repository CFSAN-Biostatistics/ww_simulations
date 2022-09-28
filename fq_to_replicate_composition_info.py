#!/usr/bin/env python
# coding: utf-8

"""
Created on Feb. 4 2022

@author: James.Pettengill
"""

import os
import sys

# print(sys.argv[0])
# print(len(sys.argv))
infile = sys.argv[1]
print(infile)

var_1 = infile.split("_")[0]
prop_1 = infile.split("_")[1]
var_2 = infile.split("_")[2]
prop_2 = infile.split("_")[3]
var_3 = infile.split("_")[4]
prop_3 = infile.split("_")[5]
var_4 = infile.split("_")[6]
prop_4 = infile.split("_")[7]
var_5 = infile.split("_")[8]
prop_5 = infile.split("_")[9].replace("art", "")

with open(infile.replace("art_out1.fq", "composition_info.tsv"), "w") as outfile:
    outfile.write(var_1 + "\t" + prop_1 + "\n" + var_2 + "\t" + prop_2 + "\n" + var_3 + "\t" + prop_3 + "\n" + var_4 + "\t" + prop_4 + "\n" + var_5 + "\t" + prop_5 )

