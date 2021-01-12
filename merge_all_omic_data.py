
import argparse
import pandas as pd
import h5py


def build_dataframe(cols, rows, data, all_features):

        df = pd.DataFrame(index=all_features, columns=cols, dtype=float)

        chunk_size = 1000
        lagging = 0
        leading = chunk_size
        while leading < len(rows):
           df.loc[rows[lagging:leading], cols] = data[lagging:leading,:]
           lagging = leading
           leading = min(len(rows), leading + chunk_size)

        # put the columns (i.e., patients) in sorted order
        srt_cols = sorted(cols)
        df = df[srt_cols]

        return df


def get_all_features(hdf_path_list):

    all_features = set([])
    for hdf_path in hdf_path_list:
        with h5py.File(hdf_path, "r") as f:
            new_features = set(f['index'])
            all_features |= new_features

    return sorted(list(all_features))


def add_dataset(input_path, group_name, all_features, out_f):
    
    with h5py.File(input_path, "r") as in_f:
        cols = in_f["columns"]
        rows = in_f["index"]
        data = in_f["data"]

        df = build_dataframe(cols, rows, data, all_features)

        dset = out_f.create_dataset(group_name + "/data", df.shape)
        dset[:,:] = df.values

        colset = out_f.create_dataset(group_name+"/columns", df.columns.shape,
                                      dtype=h5py.string_dtype('utf-8'))
        colset[:] = df.columns.values

    return


def add_all_datasets(input_path_list, group_name_list, all_features, output_path):
    
    with h5py.File(output_path, 'w') as out_f:
        for i, input_path in enumerate(input_path_list):
            add_dataset(input_path, group_name_list[i], all_features, out_f)
    
        rowset = out_f.create_dataset("index", (len(all_features),),
                                      dtype=h5py.string_dtype('utf-8'))
        rowset[:] = all_features

    return
        

if __name__=="__main__":

    argparser = argparse.ArgumentParser("Take all of the cancer-type-specific HDF files and combine them into one!")
    argparser.add_argument("--hdf-path-list", nargs="+", help="list of HDF files to merge")
    argparser.add_argument("--group-name-list", nargs="+", help="list of group names") 
    argparser.add_argument("--output", help="path to the output: a combined HDF file")
   
    args = argparser.parse_args()
    all_features = get_all_features(args.hdf_path_list) 
    print(len(all_features), "UNIQUE FEATURES")    
    add_all_datasets(args.hdf_path_list, args.group_name_list, all_features, args.output)


