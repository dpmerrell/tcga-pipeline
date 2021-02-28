

from sklearn.cluster import AgglomerativeClustering
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib import pyplot as plt
from copy import copy
import numpy as np
import argparse
import json
import h5py
import sys

sys.setrecursionlimit(10000)

def cluster_to_ordering(ch, n_s, data_rows):

    def rec_cto(children, n_samples):
        left = children[-1,0]
        right = children[-1,1]

        if left < n_samples:
            left_result = [left]
            left_var = np.var(data_rows[left,:])
            left_n = 1
        else:
            left_result, left_var, left_n = rec_cto(children[:(left-n_samples+1),:], n_samples)

        if right < n_samples:
            right_result = [right]
            right_var = np.var(data_rows[right,:])
            right_n = 1
        else:
            right_result, right_var, right_n = rec_cto(children[:(right-n_samples+1),:], n_samples)

        # Put the higher-variance samples earlier in the ordering
        if left_var > right_var:
            new_result = left_result + right_result
        else:
            new_result = right_result + left_result

        new_n = left_n + right_n
        new_var = (left_n*left_var + right_n*right_var)/new_n

        return new_result, new_var, new_n

    return rec_cto(ch, n_s)


def order_rows(data_matrix, n_features=500):

    if n_features < data_matrix.shape[1]:
        vs = np.nanvar(data_matrix, axis=0)
        maxvar_idx = np.argsort(vs)[::-1][:n_features]
        data_matrix = data_matrix[:,maxvar_idx]

    print("ORDERING ROWS")
    row_clust = AgglomerativeClustering(compute_full_tree=True)
    print("\trunning hierarchical clustering")
    row_clust.fit(data_matrix)

    ordering, _, _ = cluster_to_ordering(row_clust.children_, data_matrix.shape[0], data_matrix)

    return ordering 


def order_cols(group_v, data_matrix):

    print("ORDERING COLUMNS")

    n_cols = data_matrix.shape[1]
    idx_v = np.arange(n_cols, dtype=int)
    
    for gp in np.unique(group_v):

        gp_idx = idx_v[group_v==gp]

        print("\trunning hierarchical clustering for group", gp)
        # Perform hierarchical clustering on the columns for this group
        col_clust = AgglomerativeClustering(compute_full_tree=True)
        col_clust.fit(data_matrix[:,gp_idx].transpose())

        gp_ordering, _, _ = cluster_to_ordering(col_clust.children_, len(gp_idx), data_matrix[:,gp_idx].transpose())
        idx_v[gp_idx] = gp_idx[gp_ordering]

    return idx_v.tolist() 
    


def load_hdf(input_hdf):

    with h5py.File(input_hdf, "r") as f_in:

        group_v = f_in["cancer_types"][:].astype(str)
        ctypes = sorted(list(np.unique(group_v)))
        ctype_to_int = {ct: idx for idx, ct in enumerate(ctypes)}
        ctype_to_int_v = np.vectorize(lambda x: ctype_to_int[x])
        group_v = ctype_to_int_v(group_v).astype(int) 
      
        data = f_in["data"][:,:] 
        #n_patients = sum([len(p) for p in patients_ls])
        #idx = f_in['index'][:]
        #
        #combined = np.empty((idx.shape[0], n_patients))
        #group_v = np.empty(n_patients, dtype=int)

        #leading = 0
        #lagging = 0
        #gp = 0
        #for pat, ct in zip(patients_ls, ctypes):
        #    leading += len(pat)
        #    combined[:, lagging:leading] = f_in[ct]['data'][:,:]
        #    group_v[lagging:leading] = gp
        #    lagging = leading
        #    gp += 1
 
    return group_v, ctypes, data 


def plot_heatmap(group_v, ctypes, arr, suptitle, out_png):

    group_l = [np.searchsorted(group_v, i) for i in range(group_v[-1]+1)]
    group_r = group_l[1:] + [len(group_v)]
    tick_locs = [0.5*(l + r) for l, r in zip(group_l, group_r)]

    fig = plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(nrows=2, ncols=2, 
                           height_ratios=[1,20], 
                           width_ratios=[20,1])
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    ax1 = fig.add_subplot(gs[0,0])
    ax1.matshow(np.reshape(group_v, (1,len(group_v))) % 2, cmap="binary", aspect="auto")
    ax1.set_xticks(tick_locs)
    ax1.set_xticklabels(ctypes, rotation=45)
    ax1.set_yticks([])

    ax2 = plt.subplot(gs[1,0])#, sharex=ax1)

    vmin = np.nanquantile(arr, 0.1)
    vmax = np.nanquantile(arr, 0.9)

    # symmetrize the bounds around 0
    symm_vmin = min(vmin, -vmax)
    symm_vmax = max(vmax, -vmin)

    palette = copy(plt.cm.bwr)
    palette.set_bad('gray')

    img = ax2.matshow(arr, cmap=palette, aspect="auto", norm=colors.Normalize(vmin=symm_vmin, vmax=symm_vmax))
    ax2.set_xticks([])
    ax2.set_yticks([])

    ax3 = plt.subplot(gs[1,1])

    fig.colorbar(img, cax=ax3)

    plt.suptitle(suptitle)
    plt.tight_layout()

    plt.savefig(out_png)

    return

if __name__=="__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="path to a specially formatted HDF file")
    parser.add_argument("out_png", help="path to the output PNG file")
    parser.add_argument("--save-col-ordering", type=str, default="")
    parser.add_argument("--save-row-ordering", type=str, default="")
    parser.add_argument("--use-col-ordering", type=str, default="")
    parser.add_argument("--use-row-ordering", type=str, default="")
    parser.add_argument("--log-transform", action="store_true")
    parser.add_argument("--standardize", action="store_true")

    args = parser.parse_args()

    group_v, ctypes, arr = load_hdf(args.data_file)

    if args.log_transform:
        arr = np.log( arr - np.nanmin(arr, axis=1, keepdims=True) + 1.0)
    
    if args.standardize:
        mus = np.nanmean( arr, axis=1, keepdims=True)
        sigmas = np.nanstd(arr, axis=1, keepdims=True)
        arr = (arr - mus) /sigmas


    # Use hierarchical clustering
    # to reorder the rows and columns
    if "" in (args.use_col_ordering, args.use_row_ordering):
        cluster_arr = arr.copy()
        cluster_arr[np.isnan(arr)] = 0.0
        
        if args.use_col_ordering == "": 
            col_ordering = order_cols(group_v, cluster_arr)
        else:
            with open(args.use_col_ordering, "r") as f:
                col_ordering = json.load(f)

        if args.use_row_ordering == "":
            row_ordering = order_rows(cluster_arr)
        else:
            with open(args.use_row_ordering, "r") as f:
                row_ordering = json.load(f)


    arr = arr[:, col_ordering]
    arr = arr[row_ordering, :]
    group_v = group_v[col_ordering]

    plot_heatmap(group_v, ctypes, arr, args.data_file, args.out_png)

    if args.save_col_ordering != "":
        with open(args.save_col_ordering, "w") as f:
            col_ordering = [int(idx) for idx in col_ordering]
            json.dump(col_ordering, f)
    if args.save_row_ordering != "":
        with open(args.save_row_ordering, "w") as f:
            row_ordering = [int(idx) for idx in row_ordering]
            json.dump(row_ordering, f)

