# grsn.py
# (c) David Merrell 2021-01
# 
# An implementation of Global Rank-invariant Set Normalization.
#
# Pelz et al., 2008. Global rank-invariant set normalization (GRSN) 
#     to reduce systematic distortions in microarray data. Bioinformatics.
#     https://dx.doi.org/10.1186%2F1471-2105-9-520
#
#
# We use a modifed version of the algorithm, described by Liu et al. (2014)
# 
# Liu et al., 2014. A Comprehensive Comparison of Normalization Methods for Loading
#     Control and Variance Stabilization of Reverse-Phase Protein Array Data.
#     Cancer informatics. https://dx.doi.org/10.4137%2FCIN.S13329
#


import numpy as np
import statsmodels.api as sm


def a_m_inv_transform(a,m):
    return a + 0.5*m, a - 0.5*m


def a_m_transform(x,y):
    return 0.5*(x+y), x-y


def correct_column(col, gris_idx, ref_vals):
    a, m = a_m_transform(col, ref_vals)

    a_smooth = np.linspace(np.nanmin(a), np.nanmax(a), 1001)
    m_smooth = sm.nonparametric.lowess(m[gris_idx], a[gris_idx], xvals=a_smooth)

    correction = np.interp(a, a_smooth, m_smooth)
    print(correction)
    m_corrected = m - correction 

    _, y_corrected = a_m_inv_transform(a, m_corrected)

    return y_corrected


def apply_corrections(arr, gris_idx, ref_vals):

    for j in range(arr.shape[1]):
        arr[:,j] = correct_column(arr[:,j], gris_idx, ref_vals)

    return arr


def compute_reference_values(arr):

    ref_values = np.empty(arr.shape[0])

    for i,_ in enumerate(ref_values):
        row = arr[i,:]
        ref_values[i] = row[np.logical_and(np.nanquantile(row, 0.25) <= row,
                                             row <= np.nanquantile(row,0.75))].mean()

    print(ref_values)
    return ref_values


def remove_greatest_variation(gris_idx, ranks, n_remove):

    variances = np.nanvar(ranks, axis=1)
    keep_idx = (-variances).argsort()[n_remove:]

    return np.sort(gris_idx[keep_idx])


def grsn(arr, set_size, iterations):

    # Transform to log space
    #arr = np.log2(arr)

    # Set up the Global Rank-Invariant Set
    gris_idx = np.arange(arr.shape[0], dtype=int)
    to_remove = arr.shape[0] - set_size
    chunk_size =  to_remove // iterations
    
    # Find the Global Rank-Invariant Set
    for _ in range(iterations - 1):
        ranks = arr[gris_idx,:].argsort(axis=0)
        gris_idx = remove_greatest_variation(gris_idx, ranks, chunk_size)
    # Remove the last chunk
    ranks = arr[gris_idx,:].argsort(axis=0)
    gris_idx = remove_greatest_variation(gris_idx, ranks, gris_idx.shape[0] - set_size)

    # Compute reference values
    ref_values = compute_reference_values(arr) 

    # Correct the data, using the reference values
    arr = apply_corrections(arr, gris_idx, ref_values)

    # Undo the log transform
    #return 2.0**arr
    print(arr)
    return arr


if __name__=="__main__":


    pass


