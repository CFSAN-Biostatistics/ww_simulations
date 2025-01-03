#!/usr/bin/env python
# coding: utf-8

# In[37]:


import os
import pandas as pd
import sys


results_directory = sys.argv[1]
read_length = results_directory.split("_")[1].replace("L", "")
fragment_length = results_directory.split("_")[2].replace("F", "").replace("/", "")
platform = "illumina"
amplicon_panel = results_directory.split("_")[3].replace("/", "")
# print(read_length, fragment_length, platform, amplicon_panel)


# Read in pango/who/GISAID mapping file

# In[40]:


var_names_df = pd.read_csv("who_GISAID_mapping.txt", header=0, sep="\t", dtype="str")
var_names_df["who"] = var_names_df["who"].str.upper()
# print(var_names_df)


# Parse replicate names

# In[41]:


rep_info = pd.DataFrame()
dirs = [
    i for i in os.listdir(results_directory) if "_results" in i and "primer" not in i
]
num_vars = len(dirs[0].split("_"))
for res_dir in dirs:
    var_prop = []
    for j in range(0, (num_vars - 1)):
        if (j % 2) == 0:
            var_prop.append([res_dir.split("_")[j], res_dir.split("_")[j + 1]])
    var_prop_df = pd.DataFrame(var_prop)
    var_prop_df.columns = ["GISAID", "KnownProportion"]
    var_prop_df["Rep"] = res_dir
    var_prop_df["GISAID"] = var_prop_df["GISAID"].str.upper()
    rep_info = rep_info.append(var_prop_df)
# print(rep_info)


# Merge replicate names and mapping file

# In[42]:


rep_info_df = pd.merge(rep_info, var_names_df, on="GISAID", how="inner")
rep_info_df = rep_info_df[["who", "GISAID", "Rep", "KnownProportion"]]
# print(rep_info_df.head(10))


# Get freyja results

# In[43]:


freyja_files = []
for res_dir in dirs:
    contents = os.listdir(results_directory + "/" + res_dir + "/")
    for f in contents:
        if "freyja.demix" in f:
            freyja_files.append(results_directory + res_dir + "/" + f)

            freyja_res_df = pd.DataFrame()
for idx, f in enumerate(freyja_files):
    freyja_res_dict = {}
    rep = f.split("/")[2]
    freyja_rep_res = []
    with open(f, "r") as freyja_res:
        for line in freyja_res:
            if "summarized" in line:
                to_parse = line.replace("summarized\t", "").replace("\n", "").split(")")
                num_vars = len(to_parse)
                for i in range(0, num_vars - 1):
                    var_info = (
                        to_parse[i]
                        .replace("[('", "")
                        .replace(", ('", "")
                        .replace("'", "")
                    )
                    variant = var_info.split(", ")[0]
                    abund = var_info.split(", ")[1]
                    if variant != "Other":
                        freyja_rep_res.append([variant, abund])
                    else:
                        pass
    freyja_res_dict[rep] = freyja_rep_res
    tmp_df = pd.DataFrame.from_dict(freyja_rep_res)
    if tmp_df.empty:
        pass
    else:
        tmp_df.columns = ["who", "Abundance"]
        tmp_df["who"] = tmp_df["who"].str.upper()
        tmp_df["Abundance"] = tmp_df["Abundance"].astype(float) * 100
        tmp_total = sum(tmp_df["Abundance"])
        new_row = {"who": "OTHER", "Abundance": 100 - tmp_total}
        tmp_df = tmp_df.append(new_row, ignore_index=True)
        tmp_df["RepNumber"] = idx
        tmp_df["Rep"] = rep
        tmp_df["Method"] = "freyja"
        freyja_res_df = freyja_res_df.append(tmp_df)
freyja_res_df = pd.merge(freyja_res_df, rep_info_df, on=["who", "Rep"], how="left")
freyja_res_df.drop_duplicates(subset=["who", "Rep"], inplace=True, keep="first")
# freyja_res_df.to_csv("freyja_results.tsv", index = None, sep = "\t", na_rep = "UNK")


# In[44]:


# print(freyja_res_df.head(10))


# all_Covid_bracken results

# In[45]:


allCovid_bracken_files = []
for res_dir in dirs:
    contents = os.listdir(results_directory + "/" + res_dir + "/")
    for f in contents:
        if "allCovid_bracken" in f:
            allCovid_bracken_files.append(results_directory + res_dir + "/" + f)

allCovid_bracken_df = pd.DataFrame()
for idx, f in enumerate(allCovid_bracken_files):
    rep = f.split("/")[2]
    try:
        res_df = pd.read_csv(f, header=None, sep="\t")
        res_df = res_df[res_df[3] == "P"]
        res_df = res_df[[0, 5]]
        res_df.columns = ["Abundance", "who"]
        res_df["who"] = res_df["who"].str.upper()
        res_df["who"] = res_df["who"].str.replace(" ", "")
        res_df["RepNumber"] = idx
        res_df["Rep"] = rep
        res_df["Method"] = "allCovid_bracken"
        res_df.drop_duplicates(inplace=True)
        allCovid_bracken_df = allCovid_bracken_df.append(res_df)
    except Exception as e:
        print(f"Error processing file {f}: {e}", file=sys.stderr)
allCovid_bracken_df = pd.merge(
    allCovid_bracken_df, rep_info_df, on=["who", "Rep"], how="left"
)
# allCovid_bracken_df.drop_duplicates(subset=['who', "Rep"], inplace = True, keep='first')
# allCovid_bracken_df.to_csv("allCovid_bracken_df_results.tsv", index = None, sep = "\t", na_rep = "UNK")


# majorCovid_bracken results

# In[46]:


majorCovid_bracken_files = []
for res_dir in dirs:
    contents = os.listdir(results_directory + "/" + res_dir + "/")
    for f in contents:
        if "majorCovid_bracken" in f:
            majorCovid_bracken_files.append(results_directory + res_dir + "/" + f)

majorCovid_bracken_df = pd.DataFrame()
for idx, f in enumerate(majorCovid_bracken_files):
    rep = f.split("/")[2]
    try:
        res_df = pd.read_csv(f, header=None, sep="\t")
        res_df = res_df[res_df[3] == "P"]
        res_df = res_df[[0, 5]]
        res_df.columns = ["Abundance", "who"]
        res_df["who"] = res_df["who"].str.upper()
        res_df["who"] = res_df["who"].str.replace(" ", "")
        res_df["RepNumber"] = idx
        res_df["Rep"] = rep
        res_df["Method"] = "majorCovid_bracken"
        majorCovid_bracken_df = majorCovid_bracken_df.append(res_df)
    except Exception as e:
        print(f"Error processing file {f}: {e}", file=sys.stderr)
majorCovid_bracken_df = pd.merge(
    majorCovid_bracken_df, rep_info_df, on=["who", "Rep"], how="left"
)
# majorCovid_bracken_df.drop_duplicates(subset = ["who", "Rep"], inplace = True, keep='first')
# majorCovid_bracken_df.to_csv("majorCovid_bracken_df_results.tsv", index = None, sep = "\t", na_rep = "UNK")


# kallisto results

# In[47]:


kallisto_abundance_files = []
for res_dir in dirs:
    contents = os.listdir(results_directory + "/" + res_dir + "/")
    for f in contents:
        if "kallisto_abundance" in f:
            kallisto_abundance_files.append(results_directory + res_dir + "/" + f)

kallisto_abundance_df = pd.DataFrame()
kallisto_mapping = pd.read_csv("kallisto_who_mapping.tsv", sep="\t", header=0)
kallisto_mapping["who"] = kallisto_mapping["who"].str.upper()
kallisto_mapping["kallisto"] = kallisto_mapping["kallisto"].str.upper()
for idx, f in enumerate(kallisto_abundance_files):
    rep = f.split("/")[2]
    try:
        res_df = pd.read_csv(f, header=0, sep="\t")
        total = sum(res_df["tpm"])
        res_df["target_id_simplified"] = res_df["target_id"].str.replace("_.*", "")
        res_df_new = res_df.pivot_table(
            index="target_id_simplified", values="tpm", aggfunc="sum"
        )
        res_df_new["kallisto"] = res_df_new.index
        res_df_new.reset_index(drop=True, inplace=True)
        res_df_new["Abundance"] = res_df_new["tpm"] / total
        res_df = res_df_new[["kallisto", "Abundance"]]
        res_df["kallisto"] = res_df["kallisto"].str.upper()
        res_df = pd.merge(res_df, kallisto_mapping, on="kallisto", how="left")
        res_df = res_df.pivot_table(index="who", values="Abundance", aggfunc="sum")
        res_df["Abundance"] = res_df["Abundance"] * 100
        res_df["who"] = res_df.index
        res_df.reset_index(drop=True, inplace=True)
        res_df["RepNumber"] = idx
        res_df["Rep"] = rep
        res_df["Method"] = "kallisto"
        kallisto_abundance_df = kallisto_abundance_df.append(res_df)
    except Exception as e:
        print(f"Error processing file {f}: {e}", file=sys.stderr)
kallisto_abundance_df = pd.merge(
    kallisto_abundance_df, rep_info_df, on=["who", "Rep"], how="left"
)
# kallisto_abundance_df.drop_duplicates(subset = ["who", "Rep"], inplace = True, keep='first')
# kallisto_abundance_df.to_csv("kallisto_abundance_df_results.tsv", index = None, sep = "\t")


# linear deconvolution results

# In[48]:


linearDeconvolution_abundance_files = []
for res_dir in dirs:
    contents = os.listdir(results_directory + "/" + res_dir + "/")
    for f in contents:
        if "linearDeconvolution_abundance" in f:
            linearDeconvolution_abundance_files.append(
                results_directory + res_dir + "/" + f
            )

linearDeconvolution_abundance_df = pd.DataFrame()
linear_deconvolution_mapping = pd.read_csv(
    "lineardeconvolution_who_mapping.txt", sep="\t", header=0
)
linear_deconvolution_mapping["who"] = linear_deconvolution_mapping["who"].str.upper()
linear_deconvolution_mapping["LinearDeconvolution"] = linear_deconvolution_mapping[
    "LinearDeconvolution"
].str.upper()
for idx, f in enumerate(linearDeconvolution_abundance_files):
    rep = f.split("/")[2]
    try:
        res_df = pd.read_csv(f, header=None, sep=" ")
        res_df.columns = ["LinearDeconvolution", "Abundance"]
        res_df = res_df[res_df["Abundance"] > 0]
        res_df["LinearDeconvolution"] = res_df["LinearDeconvolution"].str.upper()
        res_df = pd.merge(
            res_df, linear_deconvolution_mapping, on="LinearDeconvolution", how="left"
        )
        res_df_new = res_df.pivot_table(index="who", values="Abundance", aggfunc="sum")
        res_df_new["who"] = res_df_new.index
        res_df_new.reset_index(drop=True, inplace=True)
        res_df_new["RepNumber"] = idx
        res_df_new["Rep"] = rep
        res_df_new["Method"] = "LinearDeconvolution"
        linearDeconvolution_abundance_df = linearDeconvolution_abundance_df.append(
            res_df_new
        )
    except Exception as e:
        print(f"Error processing file {f}: {e}", file=sys.stderr)
linearDeconvolution_abundance_df = pd.merge(
    linearDeconvolution_abundance_df, rep_info_df, on=["who", "Rep"], how="left"
)
# linearDeconvolution_abundance_df.drop_duplicates(subset = ["who", "Rep"], inplace = True, keep='first')
# linearDeconvolution_abundance_df.to_csv("linearDeconvolution_abundance_df_results.tsv", index = None, sep = "\t", na_rep = "UNK")


# Merge results

# In[49]:


df = pd.concat(
    [
        linearDeconvolution_abundance_df,
        allCovid_bracken_df,
        majorCovid_bracken_df,
        kallisto_abundance_df,
        freyja_res_df,
    ],
    ignore_index=True,
)
df["Platform"] = platform
df["ReadLength"] = read_length
df["FragmentLength"] = fragment_length
df["AmpliconPanel"] = amplicon_panel
df["KnownProportion"] = df["KnownProportion"].fillna(0)


# In[50]:


df.to_csv(
    platform
    + "_L"
    + read_length
    + "_F"
    + fragment_length
    + "_"
    + amplicon_panel
    + "_results.tsv",
    sep="\t",
    index=None,
    na_rep="UNK",
)


# In[55]:


print(len(df["RepNumber"].unique()))


# In[52]:


# print(allCovid_bracken_df.head(10))


# In[53]:


# print(linearDeconvolution_abundance_df.head(10))b


# In[54]:


print()
