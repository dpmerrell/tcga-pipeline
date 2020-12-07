
import sqlite3
import argparse


def initialize_db(db_filename, table_names):

    con = sqlite3.connect(db_filename)
    c = con.cursor()

    command = "CREATE TABLE IF NOT EXISTS Patients (PatientID TEXT PRIMARY KEY, CancerType TEXT)"
    c.execute(command)

    for tbn in table_names:
        command = "CREATE TABLE IF NOT EXISTS "+ tbn +" (PatientID TEXT PRIMARY KEY)"
        c.execute(command)
    con.commit()

    con.close()

    return

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Initialize a SQLite database for the TCGA data") 
    parser.add_argument('db_filename', help="Name of the output sqlite file")
    parser.add_argument('--table_names', help="the names of the tables to initialize", nargs="+")

    args = parser.parse_args()
    initialize_db(args.db_filename, args.table_names)
    print(args.table_names)
