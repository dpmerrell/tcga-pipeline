# takes all of the CSVs for a single cancer type;
# creates one big dataframe from them;
# and dumps them into a HDF file

import pandas as pd
import argparse
import h5py


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

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()

    df.index.rename("patient", inplace=True)

    # group by patient; take the mean over replicates
    df = df.groupby("patient").mean() 

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

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()
    df.index.rename("patient", inplace=True)

    df = df.groupby("patient").mean() 
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

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()
    df.index.rename("patient", inplace=True)

    df = df.groupby("patient").mean() 
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
    print(df) 

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()
    df.index.rename("patient", inplace=True)
    df = df.groupby("patient").mean() 
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


def read_all_data(all_data_files, all_data_types):

    combined = pd.DataFrame()

    for i, data_file in enumerate(all_data_files):
        data_type = all_data_types[i]
        new_df = read_data(data_file, data_type)
        combined = pd.concat((combined, new_df), axis=1)

    print(combined)
    return combined



def write_hdf(dataframe, filename):

    with h5py.File(filename, "w") as f:
        
        dset = f.create_dataset("data", dataframe.shape)
        dset[:,:] = dataframe.values

        columns = f.create_dataset("columns", dataframe.columns.shape,
                                              dtype=h5py.string_dtype('utf-8'))
        columns[:] = dataframe.columns.values

        idx = f.create_dataset("index", dataframe.index.shape,
                                        dtype=h5py.string_dtype('utf-8'))
        idx[:] = dataframe.index.values

    return


if __name__=="__main__":

    parser = argparse.ArgumentParser("Get data from text files; merge it; and put it in an HDF file")
    parser.add_argument("--csv-files", nargs="+", help="a list of path to the sqlite database")
    parser.add_argument("--data-types", nargs="+", help="path to the text file")
    parser.add_argument("hdf_file", help="path to the output HDF file")

    args = parser.parse_args()

    df = read_all_data(args.csv_files, args.data_types) 

    write_hdf(df, args.hdf_file)


