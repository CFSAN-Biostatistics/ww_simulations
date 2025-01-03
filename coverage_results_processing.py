#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import pandas as pd


# In[38]:


# os.chdir("../")
os.getcwd()


# In[ ]:


results_directories = [
    "./art_L75_F150_ARTICv4/",
    "./art_L75_F150_nebvss1a/",
    "./art_L75_F150_QIAseq/",
    "./art_L150_F300_QIAseq/",
    "./art_L150_F300_ARTICv4/",
    "./art_L150_F300_nebvss1a/",
]
var_names_df = pd.read_csv("who_GISAID_mapping.txt", header=0, sep="\t", dtype="str")
var_names_df["who"] = var_names_df["who"].str.upper()

for results_directory in results_directories:
    all_res_df = pd.DataFrame()
    print(results_directory)
    read_length = results_directory.split("_")[1].replace("L", "")
    fragment_length = results_directory.split("_")[2].replace("F", "").replace("/", "")
    platform = "illumina"
    amplicon_panel = results_directory.split("_")[3].replace("/", "")
    rep_info = pd.DataFrame()
    dirs = [
        i
        for i in os.listdir(results_directory)
        if "_results" in i and "primer" not in i
    ]
    dirs = [i for i in dirs if "tsv" not in i]
    num_vars = len(dirs[0].split("_"))
    for res_dir in dirs:
        """Get replicate variation composition info"""
        # var_prop = []
        # for j in range(0, (num_vars - 1)):
        #     if (j % 2) == 0:
        #         var_prop.append([res_dir.split("_")[j], res_dir.split("_")[j+1]])
        #         var_prop_df = pd.DataFrame(var_prop)
        #         var_prop_df.columns = ["GISAID", "KnownProportion"]
        #         var_prop_df["Rep"] = res_dir
        #         var_prop_df["GISAID"] = var_prop_df["GISAID"].str.upper()
        #         rep_info = rep_info.append(var_prop_df)
        #         rep_info_df = pd.merge(rep_info, var_names_df, on = "GISAID", how = "inner")
        #         rep_info_df = rep_info_df[["who", "GISAID", "Rep", "KnownProportion"]]
        """ Get coverage files"""
        pos_info_files = []
        contents = os.listdir(results_directory + "/" + res_dir)
        for f in contents:
            if "_pos-coverage-quality.tsv" in f:
                pos_info_files.append(results_directory + res_dir + "/" + f)
        """ Get coverage info"""
        for idx, f in enumerate(pos_info_files):
            rep = f.split("/")[2]
            print(rep)
            res_df = pd.read_csv(
                f, header=None, sep="\t", names=["Position", "Coverage", "Quality"]
            )
            res_df["Rep"] = rep
            res_df = res_df[res_df["Coverage"] == 0]
            # print(res_df.shape)
            # pos_info_df = pd.merge(pos_info_df, rep_info_df, on = ["Rep"], how = "left")
            all_res_df = all_res_df.append(res_df)
            # pos_info_df["Platform"] = platform
            # pos_info_df["ReadLength"] = read_length
            # pos_info_df["FragmentLength"] = fragment_length
            # pos_info_df["AmpliconPanel"] = amplicon_panel
            # pos_info_df['KnownProportion'] = pos_info_df['KnownProportion'].fillna(0)
    all_res_df.to_csv(
        results_directory.replace("./", "").replace("/", "") + "_coverage_results.tsv",
        sep="\t",
        index=None,
        na_rep="UNK",
    )

# print(var_names_df)


# In[ ]:


# all_res_df.to_csv("coverage_results.tsv", sep = "\t", index = None, na_rep = "UNK")
