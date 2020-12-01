# tcga-database
Downloads TCGA data and puts it into a nice relational database.

## The Cancer Genome Atlas ([TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga))

* TCGA includes multi-omic and clinical data, collected from 11,000+ patients divided into 38 cancer types.
* It is a really compelling dataset for supervised or unsupervised machine learning tasks.
* You can access TCGA data via [FireBrowse](http://firebrowse.org/) or the [`firehose_get`](https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844333139/Download) command line tool.
* Unfortunately, these tools (FireBrowse/Firehose) are inconvenient for the uninitiated. 
    - Unless you have a lot of bioinformatics background knowledge, it can be very difficult to understand which data you should actually download.
    - Once you've downloaded the data, it exists as a bunch of zipped CSV files.

This repository contains a complete workflow for (i) downloading TCGA data and (ii) putting it into a nice SQLite database with a sensible structure.

## Setup and execution

## Database structure

## Some bioinformatics details

## License
