
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


def read_all_data(data_filenames, all_cancer_types):

    ctype_ls = []
    combined = pd.DataFrame()

    for data_file, cancer_type in zip(data_filenames, all_cancer_types):
        new_df = read_clinical_data(data_file)
        ctype_ls += [cancer_type]*new_df.shape[1]
        combined = pd.concat((combined, new_df), axis=1)

    return combined, ctype_ls


def store_all_data(combined_data, cancer_type_list, output_filename):

    combined_data = combined_data.astype(str)
    
    with h5py.File(output_filename, "w") as out_f:
 
        dset = out_f.create_dataset("data", combined_data.shape,
                                    dtype=h5py.string_dtype('utf-8'))
        dset[:,:] = combined_data.values

        rowset = out_f.create_dataset("index", combined_data.index.shape,
                                      dtype=h5py.string_dtype('utf-8'))
        rowset[:] = combined_data.index.values
        
        colset = out_f.create_dataset("columns", combined_data.columns.shape,
                                      dtype=h5py.string_dtype('utf-8'))
        colset[:] = combined_data.columns.values

        ctype_set = out_f.create_dataset("cancer_types", (len(cancer_type_list),),
                                         dtype=h5py.string_dtype('utf-8'))
        ctype_set[:] = cancer_type_list

    return 


if __name__=="__main__":

    parser = argparse.ArgumentParser("Get data from text files; merge it; and put it in an HDF file")
    parser.add_argument("--csv-files", nargs="+", help="a list of paths to clinical CSV files")
    parser.add_argument("--cancer-types", nargs="+", help="a list of cancer type names")
    parser.add_argument("hdf_file", help="path to the output HDF file")

    args = parser.parse_args()

    combined_df, ctype_ls = read_all_data(args.csv_files, args.cancer_types) 

    store_all_data(combined_df, ctype_ls, args.hdf_file)


