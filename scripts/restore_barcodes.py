
import pandas as pd
import sys


def assign_barcodes(data_file, barcoded_file, backup_barcode_file):

    data_df = pd.read_csv(data_file, sep="\t")
    gene_col = data_df.columns[0]
    print("GENE COLUMN: ", gene_col)
    data_df.set_index(gene_col, inplace=True)

    barcodes = []
    with open(barcoded_file, "r") as f:
        barcodes = f.readline().strip().split("\t")[1:]
        
    barcode_samples = [bc[:15] for bc in barcodes]
    sample_to_barcode = {s: barcodes[i] for i, s in enumerate(barcode_samples)}

    if backup_barcode_file is not None:
        with open(backup_barcode_file, "r") as f_backup:
            backup_barcodes = f_backup.readline().strip().split("\t")[1:]
            backup_samples = [bc[:15] for bc in backup_barcodes]
            backup_sample_to_barcode = {s: backup_barcodes[i] for i, s in enumerate(backup_samples)}

            diff_samples = set(backup_samples).difference(barcode_samples)
            for ds in diff_samples:
                sample_to_barcode[ds] = backup_sample_to_barcode[ds]

    print("DATA COLUMNS:", data_df.columns)
    print("SAMPLE TO BARCODE:", sample_to_barcode)
    new_cols = [sample_to_barcode[pat] for pat in data_df.columns]
    data_df.columns = new_cols

    return data_df


if __name__=="__main__":


    data_file = sys.argv[1]
    barcoded_file = sys.argv[2]
    output_file = sys.argv[3]

    backup_barcode_file = None
    if len(sys.argv) > 4:
        backup_barcode_file = sys.argv[4]

    barcode_df = assign_barcodes(data_file, barcoded_file, backup_barcode_file) 

    barcode_df.to_csv(output_file, sep="\t")


