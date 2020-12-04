
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
FIREHOSE_DOWNLOAD_DIR = os.path.join(PWD,"stddata__{}".format(TIMESTAMP)) 

LOG_DIR = os.path.join(PWD, "logs")

DATA_FILES = {ctype:{dtype: os.path.join(FIREHOSE_DOWNLOAD_DIR,
	                         ctype, 
				 TIMESTAMP_ABBR,
                                 ".".join(["gdac.broadinstitute.org_"+ ctype,
					  DATA_TYPE_DICT[dtype]["zip_template"],
					  TIMESTAMP_ABBR+"00.0.0"]),
				 ctype+"."+DATA_TYPE_DICT[dtype]["file_template"]
	                        )\
		    for dtype in DATA_TYPES\
		    if (ctype not in DATA_EXCEPTIONS.keys() or dtype not in DATA_EXCEPTIONS[ctype])}\
             for ctype in CANCER_TYPES\
	         }
print(DATA_FILES)


rule all:
    input: 
        [os.path.join(LOG_DIR, ct+"."+DATA_TYPE_DICT[dt]["file_template"]+".log")\
			for ct in DATA_FILES.keys() for dt in DATA_FILES[ct].keys()]

    #input: [DATA_FILES[ct][dt] for ct in DATA_FILES.keys() for dt in DATA_FILES[ct].keys()]

rule add_data_to_db:
    input:
        db=DB_FILENAME,
	data_file=lambda w: os.path.join(FIREHOSE_DOWNLOAD_DIR,w["cancer_type"], TIMESTAMP_ABBR,
			    "gdac.broadinstitute.org_"+w["cancer_type"]+"."+\
		            DATA_TYPE_DICT[FNAME_DECODER[w["file_template"]]]["zip_template"]\
			        +"."+TIMESTAMP_ABBR+"00.0.0",
		            w["cancer_type"]+"."+w["file_template"])
    output:
        os.path.join(LOG_DIR, "{cancer_type,[A-Z]+}.{file_template}.log")
    group: "add_data"
    shell:
        "python add_data_script.py {input.db} {input.data_file}"


rule initialize_db:
    output:
        DB_FILENAME
    shell:
        "python db_init_script.py --output {output} --data-types {DATA_TYPES} --cancer-types {CANCER_TYPES}"


rule unzip_data:
    input:
        os.path.join(FIREHOSE_DOWNLOAD_DIR,"{cancer_type}", TIMESTAMP_ABBR,
		     "gdac.broadinstitute.org_{cancer_type}.{zip_template}.{TIMESTAMP_ABBR}00.0.0.tar.gz")
    params:
        zip_dir=os.path.join(FIREHOSE_DOWNLOAD_DIR,"{cancer_type}", TIMESTAMP_ABBR)
    output:
        os.path.join(FIREHOSE_DOWNLOAD_DIR,"{cancer_type}", TIMESTAMP_ABBR,
		     "gdac.broadinstitute.org_{cancer_type}.{zip_template}.{TIMESTAMP_ABBR}00.0.0",
		     "{cancer_type}.{file_template}")
    group: "add_data"
    shell:
        "tar -xvzf {input} --directory {params.zip_dir}"


rule download_data:
    input:
        "firehose_get"
    params:
        data_type=lambda w: ZIPNAME_DECODER[w["zip_template"]]
    output:
        os.path.join(FIREHOSE_DOWNLOAD_DIR,
	             "{cancer_type}", TIMESTAMP_ABBR,
                     "gdac.broadinstitute.org_{cancer_type}.{zip_template}.{TIMESTAMP_ABBR}00.0.0.tar.gz")
    shell:
        "./{input} -b -tasks {params.data_type} stddata {TIMESTAMP} {wildcards.cancer_type}"


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



