# Feature Selection Pipeline

The pipeline uses miRNA-seq and RNA-seq data downloaded from The Cancer Genome Atlas Project (TCGA). Our main purpose is to demonstrate the different methods of applying fuzzification on expression data. It not only models various value distributions based on fitted Gaussian distribution and quantitative classifiation, but also provides inputs for down-stream analyses such as data-based feature selection.


# Input Data

Nine cancer types containing at least 30 pairs of primary tumour and solid normal tissue samples collected from the same patient under identical sample submitter ID were selected, namely:

| abbreviation | full name |
| :---: | :--- |
| BRCA | breast invasive carcinoma |
| HNSC | head and neck squamous cell carcinoma |
| KIRC | kidney renal clear cell carcinoma |
| KIRP | kidney renal papillary cell carcinoma |
| LIHC | liver hepatocellular carcinoma |
| LUSC | lung squamous cell carcinoma |
| PRAD | prostate adenocarcinoma |
| STAD | stomach adenocarcinoma |
| THCA | thyroid carcinoma |


# Dependencies

| Python package | R package |
| :---: | :---: |
| flowsets | dplyr |
|matplotlib| GenomicRanges |
| numpy | ggplot2 |
| pandas | ggVennDiagram |
| polars | Matrix |
| scipy | TCGAbiolinks |
| seaborn |  |


# Running the Pipeline

1. Adjust parameter settings in `Makefile`. Specify both programming languages, path to the `FuzzSys` directory and base directory.

2. Call `make download` to download miRNA-seq and RNA-seq data from TCGA.

3. Call `make DESeq2` for differential expression analysis as well as fuzzification of DESeq2 log2 fold change and significance.

4. Call `make raw_log2FC` for the pipeline starting from raw counts. TO investigate the effect of DESeq2 normalization, call `make norm_log2FC` afterwards.

5. Use `paper_visualization.py` for generating figures from the paper. The index `FIGURE_ID_IN_PAPER` corresponds to the figure label from paper, i.e. `F1` for figure 1 in main text and `A2a` for subfigure a of figure 2 in appendix.

```
python paper_visualization.py --dir BASE_DIRECTORY --index FIGURE_ID_IN_PAPER
```

6. Change the comments on top of the `Makefile` to run the pipeline for the other data set.

