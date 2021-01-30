
import argparse
import h5py
import numpy as np
import grsn

def load_all_data(input_hdf):

    print("Loading data: ", args.input_hdf)

    with h5py.File(args.input_hdf, "r") as f_in:

        ctypes = [k for k in f_in.keys() if k != "index"]
        patients_ls = [f_in[ct]['columns'][:] for ct in ctypes]
        
        n_patients = sum([len(p) for p in patients_ls])
        idx = f_in['index'][:]
        
        combined = np.empty((idx.size, n_patients)) 

        leading = 0
        lagging = 0
        for pat, ct in zip(patients_ls, ctypes):
            leading += len(pat)
            combined[:, lagging:leading] = f_in[ct]['data'][:]
            lagging = leading
            
    return combined, ctypes, patients_ls, idx


def dump_data(arr, ctypes, patients, idx, out_hdf):
    
    print("Saving results: ", out_hdf)

    with h5py.File(out_hdf, "w") as f_out:
        
        # Store the feature set
        rows = f_out.create_dataset("index", shape=idx.shape,
                             dtype=h5py.string_dtype('utf-8'))
        rows[:] = idx

        leading = 0
        lagging = 0

        # For each cancer type:
        for ct, pat in zip(ctypes, patients):

            leading += len(pat)

            # Store the data values
            dset = f_out.create_dataset(ct+"/data",
                                shape=(arr.shape[0], len(pat)))

            dset[:] = arr[:, lagging:leading]
 
            # Store the columns
            columns = f_out.create_dataset(ct+"/columns",
                                   shape=(len(pat),),
                                   dtype=h5py.string_dtype('utf-8'))
            columns[:] = pat

            lagging = leading

    return 

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_hdf", help="path to the input HDF file")
    parser.add_argument("output_hdf", help="path to the output HDF file")
    parser.add_argument("--set-size", help="size of the rank-invariant set", type=int, default=100)
    parser.add_argument("--iterations", help="number of iterations", type=int, default=4)

    args = parser.parse_args()

    arr, ctypes, patients, idx = load_all_data(args.input_hdf)

    print("Performing GSRN")
    result = grsn.grsn(arr, args.set_size, args.iterations)

    dump_data(result, ctypes, patients, idx, args.output_hdf)
    
    
