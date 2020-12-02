


rule unzip_data:
    input:
    output:
    shell:

rule download_data:
    input:
        "firehose_get"
    output:
        
    shell:
        "{input} -tasks preprocess clinical segmented_scna_*_hg19 rppa_annotatewithgene stddata latest"

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
        "wget http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip"
