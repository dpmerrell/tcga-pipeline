
import pandas as pd
import sys


def assign_barcodes(data_file, barcoded_file):

    data_df = pd.read_csv(data_file, sep="\t")
    data_df.set_index(["gene"], inplace=True)
    
    barcode_df = pd.read_csv(barcoded_file, sep="\t", skiprows=[1])
    barcode_df.set_index(["Hybridization REF"], inplace=True)

    print("DATA DF:", data_df)
    print("BARCODE DF:", barcode_df)
    data_samples = [col[:15] for col in data_df.columns]
    barcode_samples = [col[:15] for col in barcode_df.columns]
    print("COLUMN SYMMETRIC DIFFERENCE:", set(data_samples).symmetric_difference(barcode_samples))
    print("DATA_COLS - BARCODE_COLS:", set(data_samples).difference(barcode_samples))
    print("BARCODE_COLS - DATA_COLS:", set(barcode_samples).difference(data_samples))

    sample_to_barcode = {col[:15]: col for col in barcode_df.columns}

    new_cols = [sample_to_barcode[pat] for pat in data_df.columns]

    data_df.columns = new_cols

    return data_df


if __name__=="__main__":


    data_file = sys.argv[1]
    barcoded_file = sys.argv[2]
    output_file = sys.argv[3]

    barcode_df = assign_barcodes(data_file, barcoded_file) 

    barcode_df.to_csv(output_file, sep="\t")
