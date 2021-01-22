
import h5py
import pandas as pd
import argparse


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

    return std_gene + suffix + "_" + std_ab


def standardize_patient_id(patient_id):
    std_id = patient_id.upper()
    std_id = std_id.replace("_","-")
    std_id = "-".join(std_id.split("-")[:3])
    return std_id

"""
Discard IDs that don't seem to be valid
"""
def is_valid_gene_id(gene_id):
    return ("?" not in gene_id)
        

def read_mRNAseq_data(data_filename):

    gene_col = "gene"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_gene_id(x, "_MRNA_SEQ")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    df.set_index(gene_col, inplace=True)
    
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).mean()   

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()

    df.index.rename("patient", inplace=True)

    # group by patient; take the mean over replicates
    df = df.groupby("patient").mean() 

    # transpose the df (again...)
    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    return df


def read_RPPA_data(data_filename):

    gene_col = "Composite.Element.REF"

    # standardize the gene/antibody IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_antibody_id(x, "_PROT_RPPA")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard any unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    df.set_index(gene_col, inplace=True)
    
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).mean()   

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()
    df.index.rename("patient", inplace=True)

    df = df.groupby("patient").mean() 

    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]

    return df


def read_Methylation_data(data_filename):
    gene_col = "Hybridization REF"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t", skiprows=[1])

    sgi_func = lambda x: standardize_gene_id(x, "_METH")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    df.set_index(gene_col, inplace=True)
    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).mean()   

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()
    df.index.rename("patient", inplace=True)

    df = df.groupby("patient").mean()

    df = df.transpose() 
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    return df


def read_CopyNumber_data(data_filename):
    gene_col = "Gene Symbol"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    sgi_func = lambda x: standardize_gene_id(x, "_CNV")
    df.loc[:,gene_col] = df[gene_col].map(sgi_func)

    # Discard some extra columns
    df.drop(labels=["Locus ID", "Cytoband"], axis='columns', inplace=True)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]

    # Eliminate some redundant gene measurements by taking averages
    df = df.groupby(gene_col).mean()   

    # standardize the patient IDs
    df.columns = df.columns.map(standardize_patient_id)

    # If a patient has multiple records, take their average
    df = df.transpose()
    df.index.rename("patient", inplace=True)
    df = df.groupby("patient").mean() 

    df = df.transpose()
    srt_cols = sorted(df.columns)
    df = df[srt_cols]
    return df



def read_data(data_filename, data_type_str):

    if data_type_str=="mRNAseq_Preprocess":
        df = read_mRNAseq_data(data_filename)
    elif data_type_str == "RPPA_AnnotateWithGene":
        df = read_RPPA_data(data_filename)
    elif data_type_str == "Methylation_Preprocess":
        df = read_Methylation_data(data_filename)
    elif data_type_str == "CopyNumber_Gistic2":
        df = read_CopyNumber_data(data_filename)
    else:
        raise ValueError

    return df


def read_all_data(all_data_files, all_cancer_types, data_type):

    col_dict = {}
    combined = pd.DataFrame()

    for data_file, cancer_type in zip(all_data_files, all_cancer_types):
        new_df = read_data(data_file, data_type)
        col_dict[cancer_type] = new_df.columns.tolist()
        combined = pd.concat((combined, new_df), axis=1)

    
    return combined, col_dict


def write_hdf(dataframe, col_dict, filename):

    with h5py.File(filename, "w") as f:

        for ctype, cols in col_dict.items():
               
            dset = f.create_dataset(ctype+"/data", dataframe.loc[:,cols].values.shape,
                                                   dtype=float)
            dset[:,:] = dataframe.loc[:,cols].values

            columns = f.create_dataset(ctype+"/columns", (len(cols),),
                                       dtype=h5py.string_dtype('utf-8'))
            columns[:] = cols


        idx = f.create_dataset("index", dataframe.index.shape,
                                        dtype=h5py.string_dtype('utf-8'))

        idx[:] = dataframe.index.values
    return


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("data_type", help="Specifies the type of omic data",
                                     choices=["mRNAseq_Preprocess",
                                              "RPPA_AnnotateWithGene",
                                              "Methylation_Preprocess",
                                              "CopyNumber_Gistic2"])
    parser.add_argument("output_hdf", help="path to the output HDF file")
    parser.add_argument("--csv-files", help="A list of CSV file paths; the data to merge",
                                       nargs="+")
    parser.add_argument("--cancer-types", help="A list of cancer types corresponding to the CSV files",
                                          nargs="+")

    args = parser.parse_args()


    df, col_dict = read_all_data(args.csv_files, 
                                 args.cancer_types, 
                                 args.data_type)

    write_hdf(df, col_dict, args.output_hdf)


