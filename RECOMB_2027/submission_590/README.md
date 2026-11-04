# Where Biostatistics meets Bioinformatics: Transparent and Trusted Uncertainty Analysis of Biomedical HT Study Data with Fuzzifier*

The scripts and Makefiles were written for data processing, fuzzification and the generation of the specific FuzzSys reports.

## Input Data

The data originate from an unpublished work. Therefore, currently there are no count matrices available. Sixty-four genes were sequenced in a targeted manner in three data sets, namely:

1. Samples collected from neutrophils of 155 adults, which were categorised into five groups according to BMI.

2. Samples collected from neutrophils of 107 kids, which were categorised into four groups according to the weight condition of the kids and their mothers.

3. Samples collected from T cells of 109 kids, which were categorised into four groups according to the weight condition of the kids and their mothers.

## Dependencies

| Python packages | R packages |
| :---: | :---: |
| matplotlib | DESeq2 |
| numpy | edgeR |
| pandas | pheatmap |
| scipy | RColorBrewer |
| seaborn | readxl |
|  | tidyverse |

## Running the pipeline

1. Adjust parameter settings in `Makefile`. Specify programming language and path to the `FuzzSys` directory.

2. Use `make prepare` to extract counts and metadata from the Excel tables, perform log2-transformation on the raw counts.

3. Use `make run_test` to perform all seven hypothesis tests and generate reports, respectively.


