import os
import argparse
import numpy as np
import pandas as pd

# python main_select_DESeq2.py --data rawValueDirectory --fuzzy fuzzyValueDirectory --clusterCol metadataClusterColumn --output outputDirectory


def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--data", type = str, required = True, help = "Directory for DESeq2 results (TSV)")
    parser.add_argument ("--fuzzy", type = str, required = True, help = "Directory for fuzzy values of DESeq2 results")
    parser.add_argument ("--clusterCol", type = str, required = False, default = "cluster", help = "Column name for clusters in metadata")
    parser.add_argument ("--output", type = str, required = True, help = "Output directory for cluster-specific and marker features")
    args = parser.parse_args ()

    allClusters = sorted (set (pd.read_csv (os.path.join (args.data, "metadata.tsv"), index_col = None, sep = "\t")[args.clusterCol]))
    avgLog2FC = pd.read_csv (os.path.join (args.data, "DESeq2_log2FC.tsv"), index_col = 0, sep = "\t")
    dir = os.path.join (args.fuzzy, "log2FC", "fuzzy_values"); FS = ["--", "-", "o", "+", "++"]
    log2FC_FV = np.array ([pd.read_csv (os.path.join (dir, f"fuzzyValues_{feature}.tsv"), index_col = 0, sep = "\t").loc[allClusters, FS].to_numpy ()
                           for feature in avgLog2FC.index])
    padj = pd.read_csv (os.path.join (args.data, "DESeq2_padj.tsv"), index_col = 0, sep = "\t")
    dir = os.path.join (args.fuzzy, "padj", "fuzzy_values"); FS = ["o", "*", "**", "***", "****"]
    padj_FV = np.array ([pd.read_csv (os.path.join (dir, f"fuzzyValues_{feature}.tsv"), index_col = 0, sep = "\t").loc[allClusters, FS].to_numpy ()
                         for feature in padj.index])
    
    os.makedirs (args.output, exist_ok = True)
    
    specific = pd.concat ([pd.concat ([pd.DataFrame ({"feature": avgLog2FC.index[(log2FC_FV[:, idx, 0] == 1) & (padj_FV[:, idx, -1] == 1)],
                                                      "cluster": allClusters[idx], "regulation": "DOWN"}),
                                       pd.DataFrame ({"feature": padj.index[(log2FC_FV[:, idx, -1] == 1) & (padj_FV[:, idx, -1] == 1)],
                                                      "cluster": allClusters[idx], "regulation": "UP"})],
                                      axis = 0, ignore_index = True)
                           for idx in range (len (allClusters))],
                          axis = 0, ignore_index = True)
    specific = specific.sort_values (["cluster", "feature", "regulation"]).reset_index (drop = True)
    specific.to_csv (os.path.join (args.output, "DESeq2_cluster-specific_features.tsv"), index = None, sep = "\t")
    markers = list ()
    for feature in set (specific["feature"]):
        featureInfo = specific.loc[specific["feature"] == feature]
        if featureInfo.shape[0] == 1:
            markers.append (featureInfo)
        else:
            num = featureInfo.value_counts ("regulation"); least = num.index[-1]
            if (num[least] == 1):
                markers.append (featureInfo.loc[featureInfo["regulation"] == least])
    markers = pd.concat (markers, axis = 0).sort_values (["cluster", "regulation", "regulation"]).reset_index (drop = True)
    markers.to_csv (os.path.join (args.output, "DESeq2_marker_features.tsv"), index = None, sep = "\t")



if __name__ == "__main__":
    main ()


