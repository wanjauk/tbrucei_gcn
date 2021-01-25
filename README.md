# *Trypanosoma brucei* gene co-expression network
This repository contains data analysis procedure and results for *T. brucei* Gene Co-expression Network study.

Preprint: https://doi.org/10.21203/rs.3.rs-78868/v1  
Research publication: https://doi.org/10.1186/s13071-021-04597-6

The analysis workflow is documented on [00-data_analysis_workflow.pdf](https://github.com/wanjauk/tbrucei_gcn/blob/master/scripts/analysis/00-data_analysis_workflow.pdf) in `scripts/analysis` directory.

## Project directory structure

```
   
|-- README.md: This readme file
|
|-- environment.yml: contains libraries and packages used in this analysis.
|
|-- data
|    |-- intermediate: data generated during the intermediate steps of the
|    |   analysis. This is output data of one step that is used as input of
|    |   another step.
|    |
|    |-- raw: holds the raw data used in the analysis including metadata files,
|    |   and links to download the data from European Nucleotide Archive (ENA).
|
|-- results
|    |-- figures: png files generated from scripts/analysis and scripts/figures
|    |
|    |-- tables: tables generated from scripts/analysis
|    |    |-- includes csv, xlsx, txt, pdf files.
|    
|-- scripts
     |-- analysis: scripts to do all analyses
     |    |-- includes an Rmd file and its pdf output that documents the analysis workflow
     |
     |-- exploration: Rmd documents and scripts with exploratory analyses.
     |
     |-- figures: code used to produce some of the figures in results/figures
     |
     |-- utils: custom functions used in the analysis
```
