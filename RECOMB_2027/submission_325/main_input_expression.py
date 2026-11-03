import os
import json
import argparse
import numpy as np
import pandas as pd

# python main_input_expression.py --expression expression_directory --metadata metadata --config config --output outputDirectory


def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--expression", type = str, required = True, help = "Directory for numerator and denominator fuzzy values")
    parser.add_argument ("--metadata", type = str, required = True, help = "Metadata containing clustering information (TSV)")
    parser.add_argument ("--config", type = str, required = True, help = "Config file containing detailed parameters (JSON)")
    parser.add_argument ("--output", type = str, required = True, help = "Output directory for FlowSets inputs (TSV)")
    args = parser.parse_args ()

    with open (args.config) as f:
        config = json.load (f)
    numName = config.get ("numerator_label", "numerator"); denName = config.get ("denominator_label", "denominator")
    idxSep = config.get ("index_separator", "__"); colSep = config.get ("series_separator", "*")
    exp_FS = config["expression_fuzzy_sets"]

    metadata = pd.read_csv (args.metadata, index_col = None, sep = "\t")
    if metadata.columns[0] == "Unnamed: 0":
        metadata = metadata.rename (columns = {"Unnamed: 0": "index"})
    allClusters = sorted (set (metadata[config.get ("metadata_cluster_column", "cluster")]))
    
    numDir = os.path.join (args.expression, "numerator"); denDir = os.path.join (args.expression, "denominator")
    numerator = dict (); denominator = dict (); fitScore = list (); fitStatus = list ()
    for cluster in allClusters:
        featureList = sorted (set ([file[12:-4] for file in os.listdir (os.path.join (numDir, f"submatrix_{cluster}", "fuzzy_values"))
                                    if file.endswith (".tsv")]) &
                              set ([file[12:-4] for file in os.listdir (os.path.join (denDir, f"submatrix_{cluster}", "fuzzy_values"))
                                    if file.endswith (".tsv")]))
        for feature in featureList:
            numFV = pd.read_csv (os.path.join (numDir, f"submatrix_{cluster}", "fuzzy_values", f"fuzzyValues_{feature}.tsv"),
                                 index_col = 0, sep = "\t").mean (axis = 0)
            numFV = numFV[exp_FS].rename (index = {idx: f"{idx}{colSep}{numName}" for idx in numFV.index})
            denFV = pd.read_csv (os.path.join (denDir, f"submatrix_{cluster}", "fuzzy_values", f"fuzzyValues_{feature}.tsv"),
                                 index_col = 0, sep = "\t").mean (axis = 0)
            denFV = denFV[exp_FS].rename (index = {idx: f"{idx}{colSep}{denName}" for idx in denFV.index})
            numerator[f"{feature}{idxSep}{cluster}"] = numFV; denominator[f"{feature}{idxSep}{cluster}"] = denFV
        summary = pd.read_csv (os.path.join (numDir, f"submatrix_{cluster}", "evaluations", "summary.tsv"), index_col = 0, sep = "\t").reset_index (names = "feature")
        tmp = pd.read_csv (os.path.join (denDir, f"submatrix_{cluster}", "evaluations", "summary.tsv"), index_col = 0, sep = "\t").reset_index (names = "feature")
        summary = summary.merge (tmp, on = ["feature", "individual_concept"], how = "outer", suffixes = ("_numerator", "_denominator"))
        with np.errstate (divide = "ignore", invalid = "ignore"):
            summary["deviation"] = np.log2 (summary["deviation_numerator"] / summary["deviation_denominator"]).round (3)
        summary = summary[["feature", "deviation", "individual_concept"]]
        summary["feature"] = summary["feature"] + f"{idxSep}{cluster}"
        mu = summary["deviation"].replace ([-np.inf, np.inf], np.nan).dropna ().median ()
        std1 = summary.loc[summary["deviation"] < mu, "deviation"].replace ([-np.inf, np.inf], np.nan).dropna ().std ()
        std2 = summary.loc[summary["deviation"] > mu, "deviation"].replace ([-np.inf, np.inf], np.nan).dropna ().std ()
        sigma = np.sqrt (std1 ** 2 + std2 ** 2)
        score = pd.DataFrame ({"feature": summary["feature"], f"MEDIUM{colSep}fit": np.exp (-(summary["deviation"] - mu) ** 2 / (2 * sigma ** 2))})
        score[f"BAD{colSep}fit"] = (summary["deviation"] < mu).astype (int) * (1 - score[f"MEDIUM{colSep}fit"])
        score[f"GOOD{colSep}fit"] = (summary["deviation"] > mu).astype (int) * (1 - score[f"MEDIUM{colSep}fit"])
        fitScore.append (score[["feature", f"BAD{colSep}fit", f"MEDIUM{colSep}fit", f"GOOD{colSep}fit"]].set_index ("feature"))
        isFitted = pd.DataFrame ({"feature": summary["feature"], f"FALSE{colSep}is_fitted": 1 - summary["individual_concept"].astype (int),
                                  f"TRUE{colSep}is_fitted": summary["individual_concept"].astype (int)})
        fitStatus.append (isFitted.set_index ("feature"))
    numerator = pd.DataFrame.from_dict (numerator, orient = "index").reset_index (names = "feature")
    denominator = pd.DataFrame.from_dict (denominator, orient = "index").reset_index (names = "feature")
    fitScore = pd.concat (fitScore, axis = 0).reset_index (); fitStatus = pd.concat (fitStatus, axis = 0).reset_index ()

    if not os.path.exists (args.output):
        os.makedirs (args.output, exist_ok = True)
    numerator.to_csv (os.path.join (args.output, "numerator_input.tsv"), index = None, sep = "\t")
    denominator.to_csv (os.path.join (args.output, "denominator_input.tsv"), index = None, sep = "\t")
    fitScore.to_csv (os.path.join (args.output, "concept_fit_score.tsv"), index = None, sep = "\t")
    fitStatus.to_csv (os.path.join (args.output, "concept_is_fitted.tsv"), index = None, sep = "\t")



if __name__ == "__main__":
    main ()


