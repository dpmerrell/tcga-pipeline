
import argparse
import numpy as np
import h5py


def get_patients_and_features(hdf_path_list):

    print("Getting patients and features")
    patients = {}
    features = []
    for hdf_path in hdf_path_list:

        with h5py.File(hdf_path, "r") as f:
            
            features += list(f["index"][:])
            
            for k in f.keys():
                if k != "index":
                    if k not in patients.keys():
                        patients[k] = set([])

                    patients[k] |= set(f[k]['columns'])

    return patients, features



def initialize_output(output_path, patients, features):

    print("Initializing output")

    f_out = h5py.File(output_path, 'w')

    # Create the 'index' group in the HDF file
    index = f_out.create_dataset("index", shape=(len(features),),
                                 dtype=h5py.string_dtype('utf-8'))
    index[:] = features

    # Create the cancer type-specific datasets
    for ctype, p_set in patients.items():

        dataset = f_out.create_dataset(ctype+"/data", shape=(len(features), len(p_set)))
        dataset[:] = np.nan

        columns = f_out.create_dataset(ctype+"/columns", shape=(len(p_set),),
                                       dtype=h5py.string_dtype('utf-8'))
        columns[:] = sorted(list(p_set))

    return f_out


def add_dataset(input_path, f_out, leading, lagging):

    print("Adding data from ", input_path)
    with h5py.File(input_path, "r") as f_in:
        ctypes = [k for k in f_in.keys() if k != "index"]
        
        # Update the leading row index
        leading += len(f_in["index"]) 

        # Add data for each cancer type
        for ct in ctypes:

            # Get the correct column indices
            patient_encoder = {pat: i for i, pat in enumerate(f_out[ct]["columns"])}
            patient_idx = [patient_encoder[pat] for pat in f_in[ct]["columns"]]
            
            # Insert data from f_in to the correct parts of f_out
            for i, pat in enumerate(patient_idx): 
                f_out[ct+"/data"][lagging:leading, pat] = f_in[ct+"/data"][:,i]

    # Update the lagging row index
    lagging = leading

    return leading, lagging


def add_all_datasets(input_path_list, patients, features, output_path):

    f_out = initialize_output(output_path, patients, features)
    
    leading = 0
    lagging = 0
    for i, input_path in enumerate(input_path_list):
        leading, lagging = add_dataset(input_path, f_out, 
                                       leading, lagging)

    f_out.close()

    return
        

if __name__=="__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--hdf-path-list", nargs="+", help="list of HDF files to merge")
    argparser.add_argument("--output", help="path to the output: a combined HDF file")
   
    args = argparser.parse_args()

    patients, features = get_patients_and_features(args.hdf_path_list) 

    add_all_datasets(args.hdf_path_list, patients, features, args.output)


