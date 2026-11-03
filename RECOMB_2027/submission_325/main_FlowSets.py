import os
import json
import argparse
import pandas as pd
import polars as pl
from flowsets import *

# python main_FlowSets.py --input input_directory --config config --output outputDirectory


def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--input", type = str, required = True, help = "Directory for FlowSets inputs")
    parser.add_argument ("--config", type = str, required = True, help = "Config file containing detailed parameters (JSON)")
    parser.add_argument ("--output", type = str, required = True, help = "Output directory for FlowSets plots and cluster-specific features")
    args = parser.parse_args ()

    class fuzzy:
        def __init__ (self, fuzzy_variables):
            self.terms = {fuzzy_var: list () for fuzzy_var in fuzzy_variables}

    with open (args.config) as f:
        config = json.load (f)
    numName = config.get ("numerator_label", "numerator"); denName = config.get ("denominator_label", "denominator")
    idxSep = config.get ("index_separator", "__"); colSep = config.get ("series_separator", "*")
    exp_FS = config["expression_fuzzy_sets"]; log2FC_FS = config["log2FC_fuzzy_sets"]
    cutoff = config["flow_membership_cutoff"]; comp_cutoff = max (0.05, (1 - cutoff) / 10)

    if not os.path.exists (args.output):
        os.makedirs (args.output, exist_ok = True)

    valuesDF = pd.concat ([pd.read_csv (os.path.join (args.input, "concept_is_fitted.tsv"), index_col = 0, sep = "\t"),
                           pd.read_csv (os.path.join (args.input, "concept_fit_score.tsv"), index_col = 0, sep = "\t"),
                           pd.read_csv (os.path.join (args.input, "numerator_input.tsv"), index_col = 0, sep = "\t"),
                           pd.read_csv (os.path.join (args.input, "denominator_input.tsv"), index_col = 0, sep = "\t"),
                           pd.read_csv (os.path.join (args.input, "log2FC.tsv"), index_col = 0, sep = "\t")],
                          axis = 1).reset_index ()
    series2name = tuple ([("is_fitted", "is_fitted"), ("fit", "fit"), (numName, numName), (denName, denName), ("log2FC", "log2FC")])
    mfFuzzy = {"is_fitted": fuzzy (["FALSE", "TRUE"]), "fit": fuzzy (["BAD", "MEDIUM", "GOOD"]),
               numName: fuzzy (exp_FS), denName: fuzzy (exp_FS), "log2FC": fuzzy (log2FC_FS)}
    fa = FlowAnalysis (pl.from_pandas (valuesDF), "feature", series2name, mfFuzzy, colSep)
    fa.plot_flows (seriesFontsize = 12, classFontsize = 10, figsize = (14, 8), outfile = os.path.join (args.output, "FlowSets"))
    
    strictDown = fa.flow_finder (["?", "?", "<", "?"], minLevels = ["TRUE", "MEDIUM", "NO_EXP", "LOW", "-INF"],
                                 maxLevels = ["TRUE", "GOOD", "high", "HIGH", "-"], verbose = False)
    downFlows = fa.flow_finder (["?", "?", "<=", "?"], minLevels = ["TRUE", "MEDIUM", "LOW", "LOW", "-INF"],
                                maxLevels = ["TRUE", "GOOD", "HIGH", "HIGH", "-"], verbose = False)
    downFlows |= fa.flow_finder (["?", "?", "<", "?"], minLevels = ["TRUE", "MEDIUM", "NO_EXP", "LOW", "-INF"],
                                 maxLevels = ["TRUE", "GOOD", "NOISE", "HIGH", "-"], verbose = False)
    strictUp = fa.flow_finder (["?", "?", ">", "?"], minLevels = ["TRUE", "MEDIUM", "LOW", "NO_EXP", "+"],
                               maxLevels = ["TRUE", "GOOD", "HIGH", "high", "+INF"], verbose = False)
    upFlows = fa.flow_finder (["?", "?", ">=", "?"], minLevels = ["TRUE", "MEDIUM", "LOW", "LOW", "+"],
                              maxLevels = ["TRUE", "GOOD", "HIGH", "HIGH", "+INF"], verbose = False)
    upFlows |= fa.flow_finder (["?", "?", ">", "?"], minLevels = ["TRUE", "MEDIUM", "LOW", "NO_EXP", "+"],
                               maxLevels = ["TRUE", "GOOD", "HIGH", "NOISE", "+INF"], verbose = False)
    
    highlight = {f: "black" for f in strictUp}; highlight.update ({f: "steelblue" for f in downFlows})
    fa.plot_flows (use_edges = strictUp | downFlows, specialColors = highlight, seriesFontsize = 12, classFontsize = 10, title = "DOWN REGULATION",
                   figsize = (14, 8), outfile = os.path.join (args.output, "FlowSets_DOWN"))
    highlight = {f: "black" for f in strictDown}; highlight.update ({f: "firebrick" for f in upFlows})
    fa.plot_flows (use_edges = strictDown | upFlows, specialColors = highlight, seriesFontsize = 12, classFontsize = 10, title = "UP REGULATION",
                   figsize = (14, 8), outfile = os.path.join (args.output, "FlowSets_UP"))
    
    memberships = pd.DataFrame ({"STRICT_DOWN": fa.calc_coarse_flow_memberships (use_edges = strictDown).to_pandas ().set_index ("feature")["membership"],
                                 "STRICT_UP": fa.calc_coarse_flow_memberships (use_edges = strictUp).to_pandas ().set_index ("feature")["membership"],
                                 "DOWN": fa.calc_coarse_flow_memberships (use_edges = downFlows).to_pandas ().set_index ("feature")["membership"],
                                 "UP": fa.calc_coarse_flow_memberships (use_edges = upFlows).to_pandas ().set_index ("feature")["membership"]})
    memberships = memberships.round (5).reset_index (); memberships[["feature", "cluster"]] = memberships["feature"].str.split (idxSep, expand = True)
    memberships = memberships[["feature", "cluster", "STRICT_DOWN", "STRICT_UP", "DOWN", "UP"]]
    memberships.to_csv (os.path.join (args.output, "FlowSets_memberships.tsv"), index = None, sep = "\t")
    memberships["regulation"] = memberships[["DOWN", "UP"]].idxmax (axis = 1); memberships["FV"] = memberships[["DOWN", "UP"]].max (axis = 1)
    comp = "STRICT_" + memberships[["DOWN", "UP"]].idxmin (axis = 1)
    memberships["complement"] = [memberships.loc[idx, comp[idx]] for idx in memberships.index]

    specific = memberships.loc[(memberships["FV"] > cutoff) & (memberships["complement"] < comp_cutoff), ["feature", "cluster", "regulation", "FV"]]
    specific = specific.sort_values (["cluster", "feature", "regulation"]).reset_index (drop = True)
    specific.to_csv (os.path.join (args.output, "data_cluster-specific_features.tsv"), index = None, sep = "\t")

    markers = list ()
    for feature in set (specific["feature"]):
        featureInfo = specific.loc[specific["feature"] == feature].sort_values ("FV")
        if featureInfo.shape[0] == 1:
            markers.append (featureInfo)
        else:
            num = featureInfo.value_counts ("regulation"); least = num.index[-1]; worst = featureInfo["regulation"].iloc[0]
            if (num[least] == 1) and (least != worst):
                markers.append (featureInfo.loc[featureInfo["regulation"] == least])
    markers = pd.concat (markers, axis = 0).sort_values (["cluster", "regulation", "regulation"]).reset_index (drop = True)
    markers.to_csv (os.path.join (args.output, "data_marker_features.tsv"), index = None, sep = "\t")



if __name__ == "__main__":
    main ()


