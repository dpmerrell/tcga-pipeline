# takes all of the CSVs for a single cancer type;
# creates one big dataframe from them;
# and dumps them into a HDF file

import pandas
import argparse
import h5py

def standardize_gene_id(gene_id):

    std_id = gene_id.split("|")[0]
    std_id = std_id.upper()

    return std_id


def standardize_patient_id(patient_id):
    std_id = patient_id.upper()
    std_id = std_id.replace("_","-")
    std_id = "-".join(std_id.split("-")[:3])
    return std_id


def is_valid_gene_id(gene_id):
    return gene_id.isalnum()


def read_mRNAseq_data(data_filename):

    gene_col = "gene"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    df.loc[:,gene_col] = df[gene_col].map(standardize_gene_id)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    df.set_index(gene_col, inplace=True)

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()

    return df


def read_data(data_filename, data_type_str):

    if data_type_str=="mRNAseq_Preprocess":
        df = read_mRNAseq_data(data_filename)

    return df


def read_all_data(all_data_files, all_data_types):

    combined = pd.DataFrame()

    for i, data_file in enumerate(data_files):
        data_type = data_types[i]
        
        new_df = read_data(data_file, data_type)
        combined = pd.concat((combined, new_df), axis=1)

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


