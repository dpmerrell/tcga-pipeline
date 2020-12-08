
import pandas as pd
import sqlite3
import argparse

def standardize_gene_id(gene_id):

    std_id = gene_id.split("|")[0]
    std_id = std_id.upper()

    return std_id


def standardize_patient_id(patient_id):
    std_id = patient_id.upper()
    std_id = std_id.replace("_","-")
    std_id = "-".join(std_id.split("-")[:3])
    return std_id


def is_valid_gene_id(gene_id):
    return gene_id.isalnum()


def read_mRNAseq_data(data_filename):

    gene_col = "gene"

    # standardize the gene IDs 
    df = pd.read_csv(data_filename, sep="\t")
    df.loc[:,gene_col] = df[gene_col].map(standardize_gene_id)

    # Discard the unidentified genes
    df = df.loc[df[gene_col].map(is_valid_gene_id),:]
    
    # index the df by gene
    df.set_index(gene_col, inplace=True)

    # standardize the patient IDs, then transpose the df
    df.columns = df.columns.map(standardize_patient_id)
    df = df.transpose()

    print(df)
 
    return df


def read_data(data_filename, data_type_str):

    if data_type_str=="mRNAseq_Preprocess":
        df = read_mRNAseq_data(data_filename)

    return df



def get_table_cols(conn, table_name):
    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info({})".format(table_name))
    res = cursor.fetchall()
    cols = [row[1] for row in res]
    cursor.close()
    return cols


def add_new_cols(conn, table_name, new_columns):

    c = conn.cursor()
    for col in new_columns:
        c.execute("ALTER TABLE {} ADD COLUMN {} NUMERIC".format(table_name, col))
        conn.commit()
    c.close()

    return


def row_to_tuple(df_row):

    return tuple(x if (not pd.isnull(x)) else None for x in df_row)


def add_data(df, db_path, data_type_str, log_path):

    # connect to DB
    conn = sqlite3.connect(db_path)

    # Add new columns to the table, if necessary
    existing_cols = get_table_cols(conn, data_type_str)
    new_cols = sorted(list(set(df.columns).difference(existing_cols)))
    if len(new_cols) > 0:
        add_new_cols(conn, data_type_str, new_cols)

    # Create a list of tuples to insert
    rows = [row_to_tuple(row) for row in df.iterrows()]

    # Build up the SQL command
    command = "INSERT INTO {} ("
    command += ",".join(df.columns.to_list())
    command += ") VALUES ("
    command += ",".join(["?"]*len(df.columns))
    command += ")"

    # Insert into the table
    c = conn.cursor()
    c.executemany(command, rows)
    c.close()

    conn.close()

    with open(log_path, "w") as f:
        f.write("success")
    return


if __name__=="__main__":

    parser = argparse.ArgumentParser("Take data from a text file and put it in the database")
    parser.add_argument("database", help="path to the sqlite database")
    parser.add_argument("data_file", help="path to the text file")
    parser.add_argument("cancer_type", help="string identifier for cancer type")
    parser.add_argument("data_type", help="string identifier for data type")
    parser.add_argument("log_file", help="path to an output log file")

    args = parser.parse_args()

    df = read_data(args.data_file, args.data_type) 

    add_data(df, args.database, args.data_type, args.log_file)



