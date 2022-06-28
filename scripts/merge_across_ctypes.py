
from os import path
import pandas as pd
import numpy as np
import argparse
import h5py


PRIMARY_TUMOR_TYPES = set(["01","03","09"])


def standardize_gene_id(gene_id, suffix):

    std_id = gene_id.split("|")[0]
    std_id = std_id.upper()
    std_id += suffix

    return std_id


def standardize_antibody_id(ab_id, suffix):
    gene_ab = ab_id.split("|")
    std_gene = gene_ab[0].upper()
    std_ab = gene_ab[1].upper()
    std_ab = std_ab.replace("_","-")

    return "_".join([std_gene, suffix])


def get_patient_id(barcode):
    std_id = barcode.upper()
    std_id = std_id.replace("_","-")
    std_id = "-".join(std_id.split("-")[:3])
    return std_id


def get_sample_type(barcode):
    return barcode.replace("_","-").split("-")[3][:2]


"""
Discard IDs that don't seem to be valid
"""
def is_valid_gene_id(gene_id):
    return ("?" not in gene_id)
        

def read_mRNAseq_data(data_filename):

    gene_col = "gene"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_gene_id(x, "_mrnaseq")
    df[gene_col] = df[gene_col].map(sgi_func)


    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).agg(np.nanmean)   
    
    # transpose the df
    df = df.transpose()

    df.index.rename("barcode", inplace=True)

    df["patient"] = df.index.map(get_patient_id)
    df["sample_type"] = df.index.map(get_sample_type)

    primary_tumor_idx = df["sample_type"].map(lambda x: x in PRIMARY_TUMOR_TYPES)
    df = df.loc[primary_tumor_idx, :]
    barcodes = list(df.index)

    assert len(df["patient"].unique()) == df.shape[0]
    df.set_index(["patient"], inplace=True)
    df.drop("sample_type", axis=1, inplace=True)

    # transpose the df (again...)
    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    
    return df, barcodes


def read_RPPA_data(data_filename):

    gene_col = "Composite.Element.REF"

    # standardize the gene/antibody IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_antibody_id(x, "rppa")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard any unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).agg(np.nanmean)   

    # transpose the df
    df = df.transpose()
    
    df.index.rename("barcode", inplace=True)

    df["patient"] = df.index.map(get_patient_id)
    df["sample_type"] = df.index.map(get_sample_type)

    primary_tumor_idx = df["sample_type"].map(lambda x: x in PRIMARY_TUMOR_TYPES)
    df = df.loc[primary_tumor_idx, :]
    barcodes = list(df.index)

    assert len(df["patient"].unique()) == df.shape[0]
    df.set_index(["patient"], inplace=True)
    df.drop("sample_type", axis=1, inplace=True)

    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]

    return df, barcodes


def read_Methylation_data(data_filename):
    gene_col = "Hybridization REF"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t", skiprows=[1])

    sgi_func = lambda x: standardize_gene_id(x, "_methylation")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).agg(np.nanmean) 

    # transpose the df
    df = df.transpose()

    # Parse barcodes
    df.index.rename("barcode", inplace=True)    
    df["patient"] = df.index.map(get_patient_id)
    df["sample_type"] = df.index.map(get_sample_type)

    # Keep only records for primary tumor samples
    primary_tumor_idx = df["sample_type"].map(lambda x: x in PRIMARY_TUMOR_TYPES)
    df = df.loc[primary_tumor_idx, :]
    barcodes = list(df.index)

    assert len(df["patient"].unique()) == df.shape[0]
    df.set_index(["patient"], inplace=True)
    df.drop("sample_type", axis=1, inplace=True)

    df = df.transpose() 
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    return df, barcodes


def read_CopyNumber_data(data_filename):
    gene_col = "Gene Symbol"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_gene_id(x, "_cna")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard some extra columns
    df.drop(labels=["Locus ID", "Cytoband"], axis='columns', inplace=True)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]

    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).mean()   

    df = df.transpose()

    # Parse barcodes
    df.index.rename("barcode", inplace=True)
    df["patient"] = df.index.map(get_patient_id)
    df["sample_type"] = df.index.map(get_sample_type)

    # Keep only records for primary tumor samples
    primary_tumor_idx = df["sample_type"].map(lambda x: x in PRIMARY_TUMOR_TYPES)
    df = df.loc[primary_tumor_idx, :]
    barcodes = list(df.index)

    assert len(df["patient"].unique()) == df.shape[0]
    df.set_index(["patient"], inplace=True)
    df.drop("sample_type", axis=1, inplace=True)

    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    return df, barcodes


def read_maf_data(data_filename):

    df = pd.read_csv(data_filename, sep="\t", index_col=0)

    sgi_func = lambda x: standardize_gene_id(x, "_mutation")

    sample_types = df.columns.map(get_sample_type)
    patient_ids = df.columns.map(get_patient_id)

    primary_tumor_cols = sample_types.map(lambda x: x in PRIMARY_TUMOR_TYPES)

    df = df.loc[:,primary_tumor_cols]
    patient_ids = patient_ids[primary_tumor_cols]
    barcodes = list(df.columns)
    df.columns = patient_ids

    df.index = df.index.map(sgi_func)

    return df, barcodes



def read_data(data_filename, data_type_str):

    if data_type_str=="mRNAseq_Preprocess":
        df, barcodes = read_mRNAseq_data(data_filename)
    elif data_type_str == "RPPA_AnnotateWithGene":
        df, barcodes = read_RPPA_data(data_filename)
    elif data_type_str == "Methylation_Preprocess":
        df, barcodes = read_Methylation_data(data_filename)
    elif data_type_str == "CopyNumber_Gistic2":
        df, barcodes = read_CopyNumber_data(data_filename)
    elif data_type_str == "Mutation_Packager_Oncotated_Calls":
        df, barcodes = read_maf_data(data_filename)
    else:
        raise ValueError

    return df, barcodes


def read_all_data(all_data_files, all_cancer_types, data_type):

    ctype_ls = []
    combined = pd.DataFrame()
    barcode_ls = []

    for data_file, cancer_type in zip(all_data_files, all_cancer_types):
        new_df, new_barcodes = read_data(data_file, data_type)
        ctype_ls += [cancer_type]*new_df.shape[1]
        barcode_ls += new_barcodes
        combined = pd.concat((combined, new_df), axis=1)

    
    return combined, ctype_ls, barcode_ls



def write_hdf(dataframe, ctype_ls, barcode_ls, filename):

    with h5py.File(filename, "w") as f:

        dset = f.create_dataset("data", dataframe.shape, dtype=float)
        dset[:,:] = dataframe.values

        columns = f.create_dataset("columns", dataframe.columns.shape,
                                              dtype=h5py.string_dtype('utf-8'))
        columns[:] = dataframe.columns.values

        idx = f.create_dataset("index", dataframe.index.shape,
                                        dtype=h5py.string_dtype('utf-8'))
        idx[:] = dataframe.index.values

        ctypes = f.create_dataset("cancer_types", len(ctype_ls),
                                  dtype=h5py.string_dtype('utf-8'))
        ctypes[:] = ctype_ls

        barcodes = f.create_dataset("barcodes", len(barcode_ls),
                                    dtype=h5py.string_dtype('utf-8'))
        barcodes[:] = barcode_ls

    return


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("data_type", help="Specifies the type of omic data",
                                     choices=["mRNAseq_Preprocess",
                                              "RPPA_AnnotateWithGene",
                                              "Methylation_Preprocess",
                                              "CopyNumber_Gistic2",
                                              "Mutation_Packager_Oncotated_Calls"])
    parser.add_argument("output_hdf", help="path to the output HDF file")
    parser.add_argument("--csv-files", help="A list of CSV file paths; the data to merge",
                                       nargs="+")
    parser.add_argument("--cancer-types", help="A list of cancer types corresponding to the CSV files",
                                          nargs="+")

    args = parser.parse_args()


    df, ctype_ls, barcode_ls = read_all_data(args.csv_files, 
                                             args.cancer_types, 
                                             args.data_type)

    write_hdf(df, ctype_ls, barcode_ls, args.output_hdf)


