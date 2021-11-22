
import os
import pathlib

configfile: "config.yaml"

PWD = str(pathlib.Path().absolute()) 

SCRIPT_DIR = "scripts"

TIMESTAMP = config["timestamp"]
TIMESTAMP_ABBR = "".join(TIMESTAMP.split("_"))
FIREHOSE_GET_URL = config["firehose_get_url"]
CANCER_TYPES = config["cancer_types"]
OMIC_DICT = config["data_types"]
OMIC_TYPES = list(OMIC_DICT.keys())
MISSING_FILES = config["missing_files"]
EXTRA_FILES = config["extra_files"]
VIZ_OPTIONS = config["viz_options"]

WORKFLOWS = {k: v["workflow"] for k,v in OMIC_DICT.items() if v["workflow"] != "unused"}

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
        CLINICAL_HDF,
        expand("{ot}.png", ot=WORKFLOWS.keys())

def hdf_path(wc):
    return "temp/{}/{}.hdf".format(WORKFLOWS[wc["omic_type"]], wc["omic_type"])

def get_viz_options(wc):
    return VIZ_OPTIONS[wc["omic_type"]]

rule visualize_data:
    input:
        src="scripts/make_heatmap.py",
        hdf=hdf_path
    output:
        img="{omic_type}.png"
    params:
        opts=get_viz_options
    shell:
        "python {input.src} {input.hdf} {output.img} {params.opts}"

rule merge_all_omic_data:
    input:
        ["temp/{workflow}/{omic_type}.hdf".format(workflow=v, omic_type=k) for k,v in WORKFLOWS.items()]
    output:
        OMIC_HDF
    shell:
        "python scripts/merge_all_omic_data.py --hdf-path-list {input} --output {output}"



def get_cancer_types(wc):
    return sorted( [ctype for ctype in OMIC_CSVS.keys() if wc["omic_type"] in OMIC_CSVS[ctype].keys()]  )

def get_omic_csv_files(wc):
    ctypes = get_cancer_types(wc)
    return [OMIC_CSVS[ctype][wc["omic_type"]] for ctype in ctypes]


rule merge_default_omic_data:
    input:
        csv_files=get_omic_csv_files,
        script=os.path.join("scripts","merge_across_ctypes.py")
    params:
        ctypes=get_cancer_types
    output:
        "temp/default/{omic_type}.hdf"
    shell:
        "python {input.script} {wildcards.omic_type} {output} --csv-files {input.csv_files} --cancer-types {params.ctypes}"


MRNASEQ_BARCODE_CTYPES = get_cancer_types({"omic_type": "mRNAseq_Preprocess"})
rule merge_restore_barcode_data:
    input:
        restored=expand("temp/restore_barcode/{ctype}_{{omic_type}}.tsv", ctype=MRNASEQ_BARCODE_CTYPES),
        script=os.path.join("scripts", "merge_across_ctypes.py")
    params:
        ctypes=get_cancer_types
    output:
        "temp/restore_barcode/{omic_type}.hdf"
    shell:
        "python {input.script} {wildcards.omic_type} {output} --csv-files {input.restored} --cancer-types {params.ctypes}"



def get_mrnaseq_csv_file(wc):
    return OMIC_CSVS[wc["ctype"]]["mRNAseq_Preprocess"]

def get_barcoded_mrnaseq_file(wc):
    return OMIC_CSVS[wc["ctype"]]["Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"]

rule restore_mrnaseq_barcodes:
    input:
        script=os.path.join("scripts", "restore_barcodes.py"),
        mrnaseq_tsv=get_mrnaseq_csv_file,
        barcoded_tsv=get_barcoded_mrnaseq_file
    output:
        "temp/restore_barcode/{ctype}_mRNAseq_Preprocess.tsv"
    shell:
        "python {input.script} {input.mrnaseq_tsv} {input.barcoded_tsv} {output}"


def get_methylation_csv_file(wc):
    return OMIC_CSVS[wc["ctype"]]["Methylation_Preprocess"]

def get_barcoded_methylation_file(wc):
    return OMIC_CSVS[wc["ctype"]]["fingbob"]


MUT_CTYPES = get_cancer_types({"omic_type": "Mutation_Packager_Oncotated_Calls"})
rule merge_mutation_data:
    input:
        expand("temp/featurize_mutations/{ctype}.tsv", ctype=MUT_CTYPES),
        script=os.path.join("scripts","merge_across_ctypes.py")
    params:
        ctypes=get_cancer_types
    output:
        "temp/featurize_mutations/{omic_type}.hdf"
    shell:
        "python {input.script} {wildcards.omic_type} {output} --csv-files {input} --cancer-types {params.ctypes}"

def get_maf_csv_file(wc):
    return OMIC_CSVS[wc["ctype"]]["Mutation_Packager_Oncotated_Calls"]

def get_mutsig_csv_file(wc):
    return OMIC_CSVS[wc["ctype"]]["MutSigNozzleReport2CV"]

rule featurize_mutations:
    input:
        script=os.path.join("scripts","featurize_mutations.py"),
        annotations_txt=get_maf_csv_file,
        mutsig=get_mutsig_csv_file
    output:
        "temp/featurize_mutations/{ctype}.tsv"
    shell:
        "python {input.script} {input.annotations_txt} {input.mutsig} {output}"


rule merge_clinical_data:
    input:
        CLINICAL_CSVS
    output:
        CLINICAL_HDF
    shell:
        "python scripts/merge_clinical_data.py {output} --csv-files {input} --cancer-types {CANCER_TYPES}"



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


rule download_analyses:
    input:
        "firehose_get"
    output:
        os.path.join("analyses__"+TIMESTAMP,
	             "{cancer_type}", TIMESTAMP_ABBR,
                     "gdac.broadinstitute.org_{cancer_type}-{tumor_type,[A-Z]+}.{analysis_id}.Level_4.{TIMESTAMP_ABBR}00.0.0.tar.gz")
    shell:
        "./{input} -b -tasks {wildcards.analysis_id} analyses {TIMESTAMP} {wildcards.cancer_type}"


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



