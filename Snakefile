
import os
import pathlib

configfile: "config.yaml"

PWD = str(pathlib.Path().absolute()) 

TIMESTAMP = config["timestamp"]
TIMESTAMP_ABBR = "".join(TIMESTAMP.split("_"))

CANCER_TYPES = config["cancer_types"]
DATA_TYPE_DICT = config["data_types"]
DATA_TYPES = list(DATA_TYPE_DICT.keys())
DATA_EXCEPTIONS = config["exceptions"]
#ZIPNAME_DECODER = {templates["zip_template"]: k for k, templates in DATA_TYPE_DICT.items()}
#FNAME_DECODER = {templates["file_template"]: k for k, templates in DATA_TYPE_DICT.items()}
#print(ZIPNAME_DECODER)

FIREHOSE_GET_URL = "http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip"


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
        [v for d in DATA_FILES.values() for k,v in d.items()] 


rule unzip_data:
    input:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}.tar.gz")
    output:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}",
		     "{output_file}.txt")
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
rule download_cna_data:
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



