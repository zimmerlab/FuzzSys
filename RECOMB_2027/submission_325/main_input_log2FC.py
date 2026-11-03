import os
import json
import argparse
import pandas as pd
from functools import reduce

# python main_input_log2FC.py --log2FC log2FC_directory --metadata metadata --config config --output outputFile


def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--log2FC", type = str, required = True, help = "Directory for log2FC fuzzy values")
    parser.add_argument ("--metadata", type = str, required = True, help = "Metadata containing clustering information (TSV)")
    parser.add_argument ("--config", type = str, required = True, help = "Config file containing detailed parameters (JSON)")
    parser.add_argument ("--output", type = str, required = True, help = "Output file name for FlowSets input (TSV)")
    args = parser.parse_args ()

    with open (args.config) as f:
        config = json.load (f)
    idxSep = config.get ("index_separator", "__"); colSep = config.get ("series_separator", "*"); log2FC_FS = config["log2FC_fuzzy_sets"]
    indexCol = config.get ("metadata_index_column", "index"); clusterCol = config.get ("metadata_cluster_column", "cluster")

    metadata = pd.read_csv (args.metadata, index_col = None, sep = "\t")
    if metadata.columns[0] == "Unnamed: 0":
        metadata = metadata.rename (columns = {"Unnamed: 0": "index"})
    clustering = pd.Series (metadata[clusterCol].values, index = metadata[indexCol].values); allClusters = sorted (set (clustering))
    
    foldChangeDF = list ()
    for cluster in allClusters:
        sampleList = clustering[clustering == cluster].index
        log2FC = [pd.read_csv (os.path.join (args.log2FC, f"fuzzyValues_{sample}.tsv"), index_col = 0, sep = "\t") for sample in sampleList]
        log2FC = (reduce (lambda x, y: x.add (y, fill_value = 0), log2FC) / len (log2FC))[log2FC_FS]
        foldChangeDF.append (log2FC.rename (index = {idx: f"{idx}{idxSep}{cluster}" for idx in log2FC.index},
                                            columns = {col: f"{col}{colSep}log2FC" for col in log2FC.columns}))
    foldChangeDF = pd.concat (foldChangeDF, axis = 0).reset_index (names = "feature")

    outputDir = os.path.abspath (os.path.dirname (args.output))
    if not os.path.exists (outputDir):
        os.makedirs (outputDir, exist_ok = True)
    foldChangeDF.to_csv (args.output, index = None, sep = "\t")



if __name__ == "__main__":
    main ()


