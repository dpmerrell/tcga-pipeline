
import os
import pathlib

configfile: "config.yaml"

PWD = str(pathlib.Path().absolute()) 

TIMESTAMP = config["timestamp"]
TIMESTAMP_ABBR = "".join(TIMESTAMP.split("_"))

CANCER_TYPES = config["cancer_types"]
OMIC_DICT = config["data_types"]
OMIC_TYPES = list(OMIC_DICT.keys())
MISSING_FILES = config["missing_files"]
EXTRA_FILES = config["extra_files"]

CLINICAL_DICT = config['clinical_data_types']['Clinical_Pick_Tier1']

OMIC_HDF = config["omic_hdf"]
CLINICAL_HDF = config["clinical_hdf"]

def fill_template(template, ct):
    return template.format(cancer_type=ct,
		           timestamp=TIMESTAMP,
			   timestamp_abbr=TIMESTAMP_ABBR)


# Build a dictionary containing the names of all of the
# omic data CSV files we need to download
OMIC_CSVS = {ctype:{dtype: os.path.join(fill_template(OMIC_DICT[dtype]["zip_template"],ctype),
				        fill_template(OMIC_DICT[dtype]["file_template"],ctype)
	                        )\
		    for dtype in OMIC_TYPES\
		    if (ctype not in MISSING_FILES.keys() or dtype not in MISSING_FILES[ctype])}\
             for ctype in CANCER_TYPES\
	         }

# Add the files with nonstandard names
for ctype, d in EXTRA_FILES.items():
    for dtype, fname in d.items():
        OMIC_CSVS[ctype][dtype] = fill_template(fname, ctype)

# Build a list of all the clinical data
# CSV files we need to download
CLINICAL_CSVS = [os.path.join(fill_template(CLINICAL_DICT['zip_template'], ct),
	                       fill_template(CLINICAL_DICT['file_template'], ct)) for ct in CANCER_TYPES]

rule all:
    input:
        OMIC_HDF,
        CLINICAL_HDF


rule merge_all_omic_data:
    input:
        expand("temp_hdf/{cancer_type}.hdf", cancer_type=CANCER_TYPES)
    output:
        OMIC_HDF
    shell:
        "python merge_all_data.py --hdf-path-list {input} --group-name-list {CANCER_TYPES} --output {output}"


rule merge_clinical_data:
    input:
        CLINICAL_CSVS
    output:
        CLINICAL_HDF
    shell:
        "python merge_clinical_data.py {output} --csv-files {input} --cancer-types {CANCER_TYPES}"


def get_data_files(wc):
    srt_keys = sorted(list(OMIC_CSVS[wc['cancer_type']].keys()))
    return [OMIC_CSVS[wc['cancer_type']][k] for k in srt_keys]

def get_data_types(wc):
    return sorted(list(OMIC_CSVS[wc['cancer_type']].keys()))

rule merge_ctype_data:
    input:
        data_files=get_data_files
    output:
        temp("temp_hdf/{cancer_type}.hdf")
    params:
        data_types=get_data_types
    shell:
        "python merge_ctype_data.py {output} --csv-files {input} --data-types {params.data_types}"



rule unzip_data:
    input:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}."+TIMESTAMP_ABBR+"00.0.0.tar.gz")
    output:
        os.path.join("{run_type}__{TIMESTAMP}","{cancer_type}", TIMESTAMP_ABBR,
		     "{zip_dir}.{other_timestamp,[0-9]+}00.0.0",
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



