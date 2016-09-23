# MSEA Plotting scripts

These scripts are designed to function as a pipeline, starting with an annotated vcf and ending with plots. They are based on: [MSEA: detection and quantification of mutation hotspots through mutation set enrichment analysis][article] and greatly extend the source code found on [github][mseagit].

Follow prepare_data.R for the full workflow, but here it is in short:

  - Download gene reference for the genome version of your data
  - Create domain data for exonic plots
  - Create transcription factor binding site data for promoter plots
  - Process vcf to create data formatted for input into MSEA.clust
  - Run MSEA.clust.R
  - Use the output to generate your desired plots with make.msea.plot.R

# Requirements

These scripts have been run on Python 2.7.12 and R 3.3.0 under 64-bit CentOS Linux 7 Core.

The only odd R package requirement is that `ggplot2` must be version 1.0.1 or older.

The scripts use relative paths, please maintain the directory structure.

[article]: <http://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0489-9>
[mseagit]: <https://github.com/bsml320/MSEA>
