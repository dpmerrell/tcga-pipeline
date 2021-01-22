

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import h5py
import sys
import numpy as np


def load_hdf(input_hdf):

    with h5py.File(input_hdf, "r") as f_in:

        ctypes = [str(k) for k in f_in.keys() if k != "index"]
        patients_ls = [f_in[ct]['columns'][:] for ct in ctypes]
       
        n_patients = sum([len(p) for p in patients_ls])
        idx = f_in['index'][:]
        
        combined = np.empty((idx.size, n_patients)) 
        group_v = np.empty((1,n_patients), dtype=int)

        leading = 0
        lagging = 0
        gp = 0
        for pat, ct in zip(patients_ls, ctypes):
            leading += len(pat)
            combined[:, lagging:leading] = f_in[ct]['data'][:]
            group_v[0,lagging:leading] = gp
            lagging = leading
            gp += 1
 
    return group_v, ctypes, combined


def plot_heatmap(group_v, ctypes, arr, out_png):

    #group_v = group_v % 2

    group_l = [np.searchsorted(group_v[0,:], i) for i in range(group_v[0,-1]+1)]
    group_r = group_l[1:] + [len(group_v[0,:])]
    tick_locs = [0.5*(l + r) for l, r in zip(group_l, group_r)]

    fig = plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(nrows=2, ncols=2, 
                           height_ratios=[1,20], 
                           width_ratios=[20,1])
    fig.subplots_adjust(wspace=0.0, hspace=0.0)

    ax1 = fig.add_subplot(gs[0,0])
    ax1.matshow(group_v % 2, cmap="binary", aspect="auto")
    #ax1.xaxis.tick_top()
    ax1.set_xticks(tick_locs)
    ax1.set_xticklabels(ctypes, rotation=45)
    ax1.set_yticks([])

    ax2 = plt.subplot(gs[1,0])#, sharex=ax1)

    vmin = np.nanquantile(arr, 0.1)
    vmax = np.nanquantile(arr, 0.9)

    img = ax2.matshow(arr, cmap="bwr", aspect="auto", vmin=vmin, vmax=vmax)
    ax2.set_xticks([])
    ax2.set_yticks([])

    ax3 = plt.subplot(gs[1,1])

    fig.colorbar(img, cax=ax3)

    plt.suptitle(data_file)
    plt.tight_layout()

    plt.savefig(out_png)

    return

if __name__=="__main__":

    data_file = sys.argv[1]
    out_png = sys.argv[2]

    #M = 300
    #N = 8000
    #group_v = np.random.choice(range(10), (1,N))
    #group_v = np.sort(v, axis=1)
    #group_v = group_v % 2

    #arr = np.random.randn(M,N)
    #arr[50:100,3000:4000] = np.nan

    group_v, ctypes, arr = load_hdf(data_file)

    plot_heatmap(group_v, ctypes, arr, out_png)

