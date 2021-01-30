
from os import path
import pandas as pd
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


def get_patient_id(file_path):

    fname = path.splitext(path.split(file_path)[0])[0]

    return fname[:12]


def tabulate_annotations(maf_files, gene_ls):

    patient_ids = sorted(list(set([get_patient_id(mf) for mf in maf_files])))
    result = pd.DataFrame(index=gene_ls, columns=patient_ids)
    
    for mf in maf_files:
        # load the file
        # find the _rankscore columns
        # assign heuristic rankscores to rows that have all-nan rankscores
        # compute mean rankscore for each mutation
        # index by gene
            # just in case: groupby gene; take the max

        # update the `result` dataframe with the values from this patient-specific dataframe
        # (use pandas.DataFrame.update)

    return

if __name__=="__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("maf_manifest")
    parser.add_argument("mutsig_file")
    parser.add_argument("output_tsv")

    args = parser.parse_args()

    all_genes = get_full_gene_list(args.mutsig_file)

    maf_files = get_all_maf_files(args.maf_manifest)

    df = tabulate_annotations(maf_files, all_genes)

    df.to_csv(args.output_tsv)   

 
