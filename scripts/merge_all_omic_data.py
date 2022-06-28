
import argparse
import numpy as np
import h5py


def get_patients_and_features(hdf_path_list):

    print("Getting patients and features")
    
    features = []
    ctype_to_patients = {}
    for hdf_path in hdf_path_list:

        with h5py.File(hdf_path, "r") as f:
            
            features += list(f["index"][:].astype(str))
            new_patients = f["columns"][:]
            new_ctypes = f["cancer_types"][:]

            for pat, ct in zip(new_patients, new_ctypes):
                if ct not in ctype_to_patients.keys():
                    ctype_to_patients[ct] = set([pat])
                else:
                    ctype_to_patients[ct].add(pat)

    return ctype_to_patients, features



def initialize_output(output_path, all_patients, all_ctypes, features):

    print("Initializing output")

    f_out = h5py.File(output_path, 'w')

    # OMIC DATA GROUP
    index = f_out.create_dataset("omic_data/features", shape=(len(features),),
                                 dtype=h5py.string_dtype('utf-8'))
    index[:] = features

    dataset = f_out.create_dataset("omic_data/data", shape=(len(features), len(all_patients)))
    dataset[:,:] = np.nan

    columns = f_out.create_dataset("omic_data/instances", shape=(len(all_patients),),
                                   dtype=h5py.string_dtype('utf-8'))
    columns[:] = all_patients

    ctypes = f_out.create_dataset("omic_data/cancer_types", shape=(len(all_ctypes),),
                                  dtype=h5py.string_dtype('utf-8'))
    ctypes[:] = all_ctypes

    # BARCODE GROUP
    omic_types = [str(feat).split("_")[-1] for feat in features]
    unique_omic_types = sorted(set(omic_types))

    barcodes = f_out.create_dataset("barcodes/data", shape=(len(unique_omic_types), len(all_patients)),
                                    dtype=h5py.string_dtype('utf-8'))
    barcodes[:,:] = ""
    
    barcode_columns = f_out.create_dataset("barcodes/instances", shape=(len(all_patients),),
                                           dtype=h5py.string_dtype('utf-8'))
    barcode_columns[:] = all_patients
    
    barcode_index = f_out.create_dataset("barcodes/features", shape=(len(unique_omic_types),),
                                    dtype=h5py.string_dtype('utf-8'))
    barcode_index[:] = unique_omic_types

    return f_out


def compute_chunks(idx_ls):

    idx_ls_srt = sorted(idx_ls)
    chunks = []

    cur_start = idx_ls_srt[0]
    cur_end = idx_ls_srt[0]
    for idx in idx_ls_srt[1:]:
        if idx == cur_end + 1:
            cur_end += 1
        else:
            chunks.append([cur_start, cur_end])
            cur_start = idx
            cur_end = idx
    chunks.append([cur_start, cur_end])

    return chunks


def add_dataset(input_path, patient_to_col, f_out, leading, lagging):

    print("Adding data from ", input_path)
    with h5py.File(input_path, "r") as f_in:
        
        # Update the leading row index
        leading += len(f_in["index"]) 

        # Need to build maps between
        # columns of f_in and columns of f_out.
        # The columns of f_in are a subset of the 
        # columns of f_out.
        new_data = f_in["data"][:,:]
        in_out_idx = [patient_to_col[pat] for pat in f_in["columns"][:]]
        out_in_idx = {out_idx: in_idx for in_idx, out_idx in enumerate(in_out_idx)}

        # Writing to HDF is much more efficient if we
        # write "chunks" of columns at a time
        # (rather than individual columns)
        out_column_chunks = compute_chunks(in_out_idx)
        for chunk in out_column_chunks:
            mapped_chunk = [out_in_idx[out_idx] for out_idx in range(chunk[0],chunk[1]+1)]
            f_out["omic_data/data"][lagging:leading, chunk[0]:chunk[1]+1] = new_data[:,mapped_chunk]
   
        barcode_omics = f_out["barcodes/features"][:].astype(str)
        omic_type = f_in["index"][:1].astype(str)[0].split("_")[-1]
        barcode_omic_row = np.argwhere(barcode_omics == omic_type)[0][0]
        f_out["barcodes/data"][barcode_omic_row, in_out_idx] = f_in["barcodes"][:]

    # Update the lagging row index
    lagging = leading

    return leading, lagging


def add_all_datasets(input_path_list, ctype_to_patients, features, output_path):

    # Some bookkeeping data structures
    ctypes = sorted(list(ctype_to_patients.keys()))
    all_patients = sum([sorted(list(ctype_to_patients[k])) for k in ctypes], start=[])
    patient_to_ctype = {pat: k for k in ctypes for pat in ctype_to_patients[k]}
    all_ctypes = [patient_to_ctype[pat] for pat in all_patients] 
    patient_to_col = {pat: idx for idx, pat in enumerate(all_patients)}

    f_out = initialize_output(output_path, all_patients, all_ctypes, features)
    
    leading = 0
    lagging = 0
    for i, input_path in enumerate(input_path_list):
        leading, lagging = add_dataset(input_path, patient_to_col, f_out, 
                                       leading, lagging)

    f_out.close()

    return
        

if __name__=="__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--hdf-path-list", nargs="+", help="list of HDF files to merge")
    argparser.add_argument("--output", help="path to the output: a combined HDF file")
   
    args = argparser.parse_args()

    ctype_to_patients, features = get_patients_and_features(args.hdf_path_list) 

    add_all_datasets(args.hdf_path_list, ctype_to_patients, features, args.output)


