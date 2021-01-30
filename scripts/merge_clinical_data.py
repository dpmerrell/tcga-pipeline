
import argparse
import pandas as pd
import h5py


def standardize_patient_id(patient_id):
    std_id = patient_id.upper()
    std_id = std_id.replace("_","-")
    std_id = "-".join(std_id.split("-")[:3])
    return std_id

def read_clinical_data(data_filename):

    index_col = "bcr_patient_barcode"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t", skiprows=[1,2])
    
    # index the df by gene
    df.set_index(index_col, inplace=True)

    # standardize the patient IDs and put the columns in sorted order
    df.columns = df.columns.map(standardize_patient_id)
    srt_cols = sorted(df.columns)
    df = df[srt_cols]

    return df


def read_all_data(data_filenames):

    df_list = [read_clinical_data(filename) for filename in data_filenames]
    
    return df_list


def store_all_data(df_list, cancer_type_list, output_filename):

    with h5py.File(output_filename, "w") as out_f:

        for df, ctype in zip(df_list, cancer_type_list):
           
            df = df.astype(str)
 
            dset = out_f.create_dataset(ctype+"/data", df.shape,
                                        dtype=h5py.string_dtype('utf-8'))
            dset[:,:] = df.values

            rowset = out_f.create_dataset(ctype+"/index", df.index.shape,
                                          dtype=h5py.string_dtype('utf-8'))
            rowset[:] = df.index.values
            
            colset = out_f.create_dataset(ctype+"/columns", df.columns.shape,
                                          dtype=h5py.string_dtype('utf-8'))
            colset[:] = df.columns.values

    return 


if __name__=="__main__":

    parser = argparse.ArgumentParser("Get data from text files; merge it; and put it in an HDF file")
    parser.add_argument("--csv-files", nargs="+", help="a list of path to the sqlite database")
    parser.add_argument("--cancer-types", nargs="+", help="a list of cancer type names")
    parser.add_argument("hdf_file", help="path to the output HDF file")

    args = parser.parse_args()

    df_list = read_all_data(args.csv_files) 

    store_all_data(df_list, args.cancer_types, args.hdf_file)


