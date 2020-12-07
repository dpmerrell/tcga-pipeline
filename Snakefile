
import os
import pathlib

configfile: "config.yaml"

PWD = str(pathlib.Path().absolute()) 

TIMESTAMP = config["timestamp"]
TIMESTAMP_ABBR = "".join(TIMESTAMP.split("_"))

DB_FILENAME = config["output_database"]

CANCER_TYPES = config["cancer_types"]
DATA_TYPE_DICT = config["data_types"]
DATA_TYPES = DATA_TYPE_DICT.keys()
DATA_EXCEPTIONS = config["exceptions"]
ZIPNAME_DECODER = {templates["zip_template"]: k for k, templates in DATA_TYPE_DICT.items()}
FNAME_DECODER = {templates["file_template"]: k for k, templates in DATA_TYPE_DICT.items()}
print(ZIPNAME_DECODER)

FIREHOSE_GET_URL = "http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip"

LOG_DIR = os.path.join(PWD, "logs")

def fill_template(template, ct):
    return template.format(cancer_type=ct,
		           timestamp=TIMESTAMP,
			   timestamp_abbr=TIMESTAMP_ABBR)


DATA_FILES = {ctype:{dtype: os.path.join(fill_template(DATA_TYPE_DICT[dtype]["zip_template"],ctype),
				         fill_template(DATA_TYPE_DICT[dtype]["file_template"],ctype)
	                        )\
		    for dtype in DATA_TYPES\
		    if (ctype not in DATA_EXCEPTIONS.keys() or dtype not in DATA_EXCEPTIONS[ctype])}\
             for ctype in CANCER_TYPES\
	         }
print(DATA_FILES)


rule all:
    input: 
        [os.path.join(LOG_DIR, ct+"."+dt+".log")\
            for ct in DATA_FILES.keys() for dt in DATA_FILES[ct].keys()]
        
#    input:
#        [DATA_FILES[ct][dt] for ct in DATA_FILES.keys() for dt in DATA_FILES[ct].keys()],
#	DB_FILENAME

rule add_data_to_db:
    input:
        db=DB_FILENAME,
	data_file=lambda w: DATA_FILES[w["cancer_type"]][w["data_type"]]
    output:
        os.path.join(LOG_DIR, "{cancer_type,[A-Z]+}.{data_type}.log")
    group: "add_data"
    shell:
        "python add_data_script.py {input.db} {input.data_file}"


rule initialize_db:
    output:
        DB_FILENAME
    shell:
        "python db_init_script.py {output} --table_names {DATA_TYPES}"


rule unzip_data:
    input:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}.tar.gz")
    output:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}",
		     "{output_file}.txt")
    group: "add_data"
    shell:
        "tar -xvzf {input} --directory {wildcards.run_type}__{TIMESTAMP}/{wildcards.cancer_type}/{TIMESTAMP_ABBR}/"


rule download_data:
    input:
        "firehose_get"
    output:
        os.path.join("{run_type}__{TIMESTAMP}",
	             "{cancer_type}", TIMESTAMP_ABBR,
                     "gdac.broadinstitute.org_{cancer_type}.{data_type}.Level_{data_level,[0-9]}.{TIMESTAMP_ABBR}00.0.0.tar.gz")
    shell:
        "./{input} -b -tasks {wildcards.data_type} {wildcards.run_type} {TIMESTAMP} {wildcards.cancer_type}"

# The CNA data is a little different from the others.
# We'll give it its own rule just for convenience.
rule download_cna:
    input:
        "firehose_get"
    output:
        os.path.join("analyses__"+TIMESTAMP,
	             "{cancer_type}", TIMESTAMP_ABBR,
                     "gdac.broadinstitute.org_{cancer_type}-TP.CopyNumber_Gistic2.Level_4.{TIMESTAMP_ABBR}00.0.0.tar.gz")
    shell:
        "./{input} -b -tasks CopyNumber_Gistic2 analyses {TIMESTAMP} {wildcards.cancer_type}"


rule unzip_firehose_get:
    input:
        "firehose_get_latest.zip"
    output:
        "firehose_get"
    shell:
        "unzip {input}"


rule download_firehose_get:
    output:
        "firehose_get_latest.zip"
    shell:
        "wget {FIREHOSE_GET_URL}"



