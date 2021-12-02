
from os import path
import pandas as pd
import numpy as np
import argparse
import glob


def get_full_gene_list(mutsig_file):

    ms_df = pd.read_csv(mutsig_file, sep="\t")
    genes = ms_df["gene"].unique()

    return sorted(list(genes))


def get_all_maf_files(manifest_path):
    maf_dirname = path.dirname(manifest_path)
    maf_files = glob.glob(path.join(maf_dirname, "TCGA-*.txt"))
    return maf_files


def get_barcode(file_path):
    with open(file_path, "r") as f:
        # Get the third line
        for _ in range(4):
            ln = f.readline()
        # Search the third line for the barcode column
        idx = 0
        for colname in ln.split("\t"):
            if colname.strip() == "Tumor_Sample_Barcode":
                break
            idx += 1
        # Get the fourth line
        ln = f.readline().split("\t")
        
    barcode = ln[idx]
    return barcode


def sanitize(score_str):
    res = max(score_str.split("|"))
    if res.strip() == ".":
        res = "0.0"
    return res

def get_mutation_scores(maf_file):

    print(maf_file)
    # load the file
    df = pd.read_csv(maf_file, sep="\t", skiprows=3, encoding="ISO-8859-1")
    
    # find the _rankscore columns
    rankscore_cols = [col for col in df.columns if "rankscore" in col]

    # index by genes
    df.set_index(["Hugo_Symbol"], inplace=True)

    unscored = (df[rankscore_cols].isnull().sum(axis=1) == len(rankscore_cols))  
    df.loc[ (unscored & (df["Variant_Classification"] == "Silent")), rankscore_cols ] = "0.0" 
    df.loc[ (unscored & (df["Variant_Classification"] != "Silent")), rankscore_cols ] = "0.5" 


    result = df[rankscore_cols].astype(str)
    result = result.applymap(sanitize)
    result = result.astype(float)

    result["agg_rankscore"] = np.nanmean(result[rankscore_cols], axis=1) 

    return result[["agg_rankscore"]]


def tabulate_annotations(maf_files, gene_ls):

    barcodes = sorted(list(set([get_barcode(mf) for mf in maf_files])))
    result = pd.DataFrame(index=gene_ls, columns=barcodes)
    result.loc[:,:] = 0.0
 
    for mf in maf_files:
 
        barcode = get_barcode(mf)
        score_df = get_mutation_scores(mf)

        # in case of redundancies: groupby gene; take the max
        score_df = score_df.groupby(["Hugo_Symbol"]).max()
        score_df.columns = [barcode]

        result[barcode] = score_df[barcode]


    result[result.isnull()] = 0.0
    return result


if __name__=="__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("maf_manifest")
    parser.add_argument("mutsig_file")
    parser.add_argument("output_tsv")

    args = parser.parse_args()

    all_genes = get_full_gene_list(args.mutsig_file)

    maf_files = get_all_maf_files(args.maf_manifest)

    df = tabulate_annotations(maf_files, all_genes)

    df.to_csv(args.output_tsv, sep="\t") 

 
