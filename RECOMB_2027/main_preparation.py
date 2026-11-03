import os
import json
import argparse
import numpy as np
import pandas as pd

### python main_preparation.py --mtx rawValueMatrix --metadata metadata --config config --output outputDirectory


def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--mtx", type = str, required = True, help = "Raw value matrix (TSV)")
    parser.add_argument ("--metadata", type = str, required = True, help = "Metadata containing clustering and comparison columns (TSV)")
    parser.add_argument ("--config", type = str, required = True, help = "Config file for fold change calculation (JSON)")
    parser.add_argument ("--output", type = str, required = True, help = "Output directory for splitted matrices")
    args = parser.parse_args ()

    mtx = pd.read_csv (args.mtx, index_col = 0, sep = "\t")
    metadata = pd.read_csv (args.metadata, index_col = None, sep = "\t")
    if metadata.columns[0] == "Unnamed: 0":
        metadata = metadata.rename (columns = {"Unnamed: 0": "index"})
    with open (args.config) as f:
        config = json.load (f)

    indexCol = config.get ("metadata_index_column", "index"); clusterCol = config.get ("metadata_cluster_column", "cluster")
    numeratorCol = config["numeratorCol"]; denominatorCol = config["denominatorCol"]
    pseudoCount = config.get ("pseudo_count", 0); pseudoCount = 0 if not isinstance (pseudoCount, (float, int)) else pseudoCount
    mtx += pseudoCount
    if (mtx < 0).any (axis = None):
        raise ValueError ("Negative values in matrix!")

    with np.errstate (divide = "ignore", invalid = "ignore"):
        numerator = np.log2 (pd.DataFrame (mtx[metadata[numeratorCol]].values, index = mtx.index, columns = metadata[indexCol]))
        denominator = np.log2 (pd.DataFrame (mtx[metadata[denominatorCol]].values, index = mtx.index, columns = metadata[indexCol]))
    rawLog2FC = numerator - denominator

    if not os.path.exists (args.output):
        os.makedirs (args.output, exist_ok = True)
    rawLog2FC.to_csv (os.path.join (args.output, "paired_log2FC.tsv"), sep = "\t")

    if not os.path.exists (os.path.join (args.output, "numerator_submatrices")):
        os.makedirs (os.path.join (args.output, "numerator_submatrices"), exist_ok = True)
    if not os.path.exists (os.path.join (args.output, "denominator_submatrices")):
        os.makedirs (os.path.join (args.output, "denominator_submatrices"), exist_ok = True)
    for cluster in set (metadata[clusterCol]):
        tmp = numerator[list (set (metadata.loc[metadata[clusterCol] == cluster, indexCol]))]
        tmp.to_csv (os.path.join (args.output, "numerator_submatrices", f"submatrix_{cluster}.tsv"), sep = "\t")
        tmp = denominator[list (set (metadata.loc[metadata[clusterCol] == cluster, indexCol]))]
        tmp.to_csv (os.path.join (args.output, "denominator_submatrices", f"submatrix_{cluster}.tsv"), sep = "\t")



if __name__ == "__main__":
    main ()


