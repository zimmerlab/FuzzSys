import os
import re
import json
import argparse
import itertools
import subprocess
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from flowsets import *
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

# python paper_visualization.py --dir baseDirectory --index plotNumber

os.makedirs ("visualizations", exist_ok = True)


def parseIndex (argIndex):
    if re.search (r"^[AF][0-9]+$", argIndex):
        prefix = argIndex[0]; idx = int (argIndex[1:]); suffix = ""
    elif re.search (r"^[AF][0-9]+[a-z]$", argIndex):
        prefix = argIndex[0]; idx = int (argIndex[1:-1]); suffix = argIndex[-1]
    else:
        raise ValueError ("Figure index: F1 for figure 1 in main text; A1 for figure 1 in appendix.")
    return prefix, idx, suffix



def generateFA (inputDir, padjDir):
    class fuzzy:
        def __init__ (self, fuzzy_variables):
            self.terms = {fuzzy_var: list () for fuzzy_var in fuzzy_variables}
    exp_FS = ["NO_EXP", "NOISE", "LOW", "low", "MEDIUM", "high", "HIGH"]
    log2FC_FS = ["NA", "-INF", "--", "-", "o", "+", "++", "+INF"]
    padj_FS = ["o", "*", "**", "***", "****"]
    series2name = tuple ([("is_fitted", "is_fitted"), ("fit", "fit"), ("tumour", "tumour"), ("normal", "normal"),
                          ("log2FC", "log2FC"), ("significance", "significance")])
    mfFuzzy = {"is_fitted": fuzzy (["FALSE", "TRUE"]), "fit": fuzzy (["BAD", "MEDIUM", "GOOD"]),
               "tumour": fuzzy (exp_FS), "normal": fuzzy (exp_FS), "log2FC": fuzzy (log2FC_FS),
               "significance": fuzzy (padj_FS)}
    df = pd.concat ([pd.read_csv (os.path.join (inputDir, "concept_is_fitted.tsv"), index_col = 0, sep = "\t"),
                     pd.read_csv (os.path.join (inputDir, "concept_fit_score.tsv"), index_col = 0, sep = "\t"),
                     pd.read_csv (os.path.join (inputDir, "numerator_input.tsv"), index_col = 0, sep = "\t"),
                     pd.read_csv (os.path.join (inputDir, "denominator_input.tsv"), index_col = 0, sep = "\t"),
                     pd.read_csv (os.path.join (inputDir, "log2FC.tsv"), index_col = 0, sep = "\t")],
                    axis = 1).reset_index ()
    df = df.rename (columns = {col: col.replace ("*", "+") for col in df.columns[1:]}); padjDF = list ()
    for file in os.listdir (padjDir):
        feature = file[12:-4]
        try:
            FV = pd.read_csv (os.path.join (padjDir, f"fuzzyValues_{feature}.tsv"), index_col = 0, sep = "\t")[padj_FS]
        except FileNotFoundError:
            continue
        FV = FV.rename (index = {idx: f"{feature}__{idx}" for idx in FV.index},
                        columns = {col: f"{col}+significance" for col in padj_FS})
        padjDF.append (FV.reset_index (names = "feature"))
    padjDF = pd.concat (padjDF, axis = 0).set_index ("feature").loc[df["feature"]]
    df = pd.concat ([df.set_index ("feature"), padjDF], axis = 1).reset_index ()
    FA = FlowAnalysis (pl.from_pandas (df), "feature", series2name, mfFuzzy, "+")
    flows = {"strictDOWN": FA.flow_finder (["?", "?", "<", "?", "?"],
                                           minLevels = ["TRUE", "MEDIUM", "NO_EXP", "LOW", "-INF", "*"],
                                           maxLevels = ["TRUE", "GOOD", "high", "HIGH", "-", "****"],
                                           verbose = False),
             "down": FA.flow_finder (["?", "?", "<=", "?", "?"],
                                     minLevels = ["TRUE", "MEDIUM", "LOW", "LOW", "-INF", "*"],
                                     maxLevels = ["TRUE", "GOOD", "HIGH", "HIGH", "-", "****"],
                                     verbose = False),
             "strictUP": FA.flow_finder (["?", "?", ">", "?", "?"],
                                         minLevels = ["TRUE", "MEDIUM", "LOW", "NO_EXP", "+", "*"],
                                         maxLevels = ["TRUE", "GOOD", "HIGH", "high", "+INF", "****"],
                                         verbose = False),
             "up": FA.flow_finder (["?", "?", ">=", "?", "?"],
                                   minLevels = ["TRUE", "MEDIUM", "LOW", "LOW", "+", "*"],
                                   maxLevels = ["TRUE", "GOOD", "HIGH", "HIGH", "+INF", "****"],
                                   verbose = False)}
    flows["down"] |= FA.flow_finder (["?", "?", "<", "?", "?"], minLevels = ["TRUE", "MEDIUM", "NO_EXP", "LOW", "-INF", "*"],
                                     maxLevels = ["TRUE", "GOOD", "NOISE", "HIGH", "-", "****"], verbose = False)
    flows["up"] |= FA.flow_finder (["?", "?", ">", "?", "?"], minLevels = ["TRUE", "MEDIUM", "LOW", "NO_EXP", "+", "*"],
                                   maxLevels = ["TRUE", "GOOD", "HIGH", "NOISE", "+INF", "****"], verbose = False)
    return FA, flows



def plotVolcano (dataDir, fuzzyDir, outputPath = ""):
    log2FC = pd.read_csv (os.path.join (dataDir, "DESeq2_log2FC.tsv"), index_col = 0, sep = "\t").clip (lower = -10, upper = 10)
    with open (os.path.join (fuzzyDir, "concepts_log2FC", "concepts_detailed.json")) as f:
        log2FC_concept = json.load (f)["DEFAULT"]
    padj = pd.read_csv (os.path.join (dataDir, "DESeq2_padj.tsv"), index_col = 0, sep = "\t")
    with open (os.path.join (fuzzyDir, "concepts_padj", "concepts_detailed.json")) as f:
        padj_concept = json.load (f)["DEFAULT"]
    xRange = [np.floor (log2FC.mask (~np.isfinite (log2FC)).min (axis = None, skipna = True)),
              np.ceil (log2FC.mask (~np.isfinite (log2FC)).max (axis = None, skipna = True))]
    yRange = [np.floor (padj.mask (~np.isfinite (padj)).min (axis = None, skipna = True)) - 0.05,
              np.ceil (padj.mask (~np.isfinite (padj)).max (axis = None, skipna = True)) + 1]
    concept = [log2FC_concept[key][0] for key in ["--", "-", "o", "+", "++"]]
    mu, sigma = concept[2]; log2FC_cutoff = [xRange[0], (concept[0][2] + concept[0][3]) / 2]
    xValues = np.linspace (concept[1][2], concept[3][1], 1000)
    values = pd.DataFrame ({"x": xValues, "left": ((xValues - concept[1][3]) / (concept[1][2] - concept[1][3])).clip (min = 0),
                            "middle": np.exp (-(xValues - mu) ** 2 / (2 * sigma ** 2)),
                            "right": ((xValues - concept[3][0]) / (concept[3][1] - concept[3][0])).clip (min = 0)})
    values["diff1"] = values["left"] - values["middle"]; values["diff2"] = values["middle"] - values["right"]
    log2FC_cutoff += [values.loc[values["diff1"].abs ().sort_values ().index].iloc[:2, 0].mean (),
                      values.loc[values["diff2"].abs ().sort_values ().index].iloc[:2, 0].mean ()]
    log2FC_cutoff += [(concept[4][0] + concept[4][1]) / 2, xRange[1]]
    concept = np.array ([padj_concept[key][0] for key in ["o", "*", "**", "***", "****"]])
    padj_cutoff = [yRange[0]] + concept[:-1, -2:].mean (axis = 1).tolist () + [yRange[1]]
    pltData = log2FC.reset_index ().melt (id_vars = "index", value_name = "log2FC")
    pltData = pltData.merge (padj.reset_index ().melt (id_vars = "index", value_name = "padj"), on = ["index", "variable"], how = "inner")
    pltData["class"] = ""; coords = dict ()
    for (i, j) in itertools.product (range (5), range (5)):
        mask = (pltData["log2FC"] > log2FC_cutoff[i]) & (pltData["log2FC"] < log2FC_cutoff[i + 1]) & \
            (pltData["padj"] > padj_cutoff[j]) & (pltData["padj"] < padj_cutoff[j + 1])
        pltData.loc[mask, "class"] = f"{i}_{j}"
        coords[f"{i}_{j}"] = [(log2FC_cutoff[i] + log2FC_cutoff[i + 1]) / 2, (padj_cutoff[j] + padj_cutoff[j + 1]) / 2]
    labels = pltData.value_counts ("class"); labels = {key: f"{(labels.get (key, 0) / pltData.shape[0]):.2%}\n({labels.get (key, 0)})"
                                                       for key in coords.keys ()}
    pltData = log2FC.reset_index ().melt (id_vars = "index", value_name = "log2FC")
    pltData = pltData.merge (padj.reset_index ().melt (id_vars = "index", value_name = "padj"), on = ["index", "variable"], how = "inner")
    fig, ax = plt.subplots (figsize = (10, 6))
    sns.scatterplot (pltData, x = "log2FC", y = "padj", color = "silver", legend = None, ax = ax)
    ax.set_xlim (xRange); ax.set_ylim (yRange)
    for val in log2FC_cutoff[1:-1]:
        ax.axvline (val, color = "black", linestyle = "dashed")
    for val in padj_cutoff[1:-1]:
        ax.axhline (val, color = "black", linestyle = "dashed")
    for key in coords.keys ():
        ax.text (*coords[key], labels[key], fontdict = {"size": 6, "ha": "center", "weight": "bold"})
    ax.set_xlabel ("DESeq2 average log2 fold change", size = 10); ax.set_ylabel ("-log10 (DESeq2 corrected p-value)", size = 10)
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def _getLines (concept, xRange):
    lines = list (); curves = list (); names = list (); handles = list ()
    const = {"-infinity": xRange[0], "-inf": xRange[0],
             "+infinity": xRange[1], "+inf": xRange[1], "infinity": xRange[1], "inf": xRange[1],
             "nan": np.nan, "na": np.nan, "zero": 0}
    cutoffs = [concept.get ("MIN-NOISE", xRange[0]), concept.get ("MAX-NOISE", xRange[1])]
    cutoffs[0] = const.get (cutoffs[0].lower (), xRange[0]) if isinstance (cutoffs[0], str) else cutoffs[0]
    cutoffs[1] = const.get (cutoffs[1].lower (), xRange[1]) if isinstance (cutoffs[1], str) else cutoffs[1]
    if concept.get ("number_fuzzy_sets", 0) == 0 or cutoffs[0] >= cutoffs[1]:
        return lines, curves, names, handles
    for key in concept:
        if key in ["number_fuzzy_sets", "label_values", "MIN-NOISE", "MAX-NOISE"]:
            continue
        params, typeFS, color, _ = concept[key]
        names.append (key); handles.append (Line2D ([0], [0], color = color, linewidth = 2))
        if typeFS == "Gaussian":
            if params[1] > 0:
                xValues = np.linspace (*cutoffs, 1000)
                yValues = np.exp (-(xValues - params[0]) ** 2 / (2 * params[1] ** 2))
                curves.append ([xValues, yValues, color])
            else:
                continue
        elif typeFS == "trapezoidal":
            if cutoffs[1] <= params[0] and params[0] != params[1]:
                continue
            elif cutoffs[1] > params[0] and cutoffs[1] < params[1]:
                y_cutoffs = [(cutoffs[0] - params[0]) / (params[1] - params[0]),
                             (cutoffs[1] - params[0]) / (params[1] - params[0])]
                lines += [(max (cutoffs[0], params[0]), cutoffs[1]),
                          (max (y_cutoffs[0], 0), y_cutoffs[1]), color]
            elif cutoffs[1] >= params[1] and cutoffs[1] <= params[2]:
                if cutoffs[0] < params[1]:
                    y_cutoffs = [0 if params[0] == params[1] else (cutoffs[0] - params[0]) / (params[1] - params[0]), 1]
                    lines += [(max (cutoffs[0], params[0]), params[1]),
                              (max (y_cutoffs[0], 0), 1), color,
                              (params[1], cutoffs[1]), (1, 1), color]
                else:
                    lines += [(cutoffs[0], cutoffs[1]), (1, 1), color]
            else:
                if cutoffs[0] < params[1]:
                    y_cutoffs = [0 if params[0] == params[1] else (cutoffs[0] - params[0]) / (params[1] - params[0]),
                                 0 if params[2] == params[3] else (cutoffs[1] - params[3]) / (params[2] - params[3])]
                    lines += [(max (cutoffs[0], params[0]), params[1]), (max (y_cutoffs[0], 0), 1), color,
                              (params[1], params[2]), (1, 1), color,
                              (params[2], min (cutoffs[1], params[3])), (1, max (y_cutoffs[1], 0)), color]
                elif cutoffs[0] >= params[1] and cutoffs[0] <= params[2]:
                    y_cutoffs = [1, 0 if params[2] == params[3] else (cutoffs[1] - params[3]) / (params[2] - params[3])]
                    lines += [(cutoffs[0], params[2]), (1, 1), color,
                              (params[2], min (cutoffs[1], params[3])),
                              (1, max (y_cutoffs[1], 0)), color]
                else:
                    if params[2] == params[3]:
                        return
                    y_cutoffs = [(cutoffs[0] - params[3]) / (params[2] - params[3]),
                                 (cutoffs[1] - params[3]) / (params[2] - params[3])]
                    lines += [(cutoffs[0], min (cutoffs[1], params[3])),
                              (y_cutoffs[0], max (y_cutoffs[1], 0)), color]
        else:
            raise ValueError
    return lines, curves, names, handles



def plotConcept (mtx, concept, axis_label, outputPath = ""):
    values = mtx.mask (~np.isfinite (mtx)).melt ()["value"].dropna ()
    xRange = [np.floor (values.min ()), np.ceil (values.max ())]
    lines, curves, names, handles = _getLines (concept, xRange)
    fig, ax = plt.subplots (figsize = (8, 4))
    ax.hist (values, bins = 100, color = "silver"); ax.set_xlim (xRange)
    ax.set_xlabel (axis_label, size = 10); ax.set_ylabel ("number of values", size = 10)
    ax2 = ax.twinx (); ax2.set_ylim ((0, 1.05)); ax2.set_ylabel ("fuzzy value", size = 10)
    ax2.plot (*lines, linewidth = 2)
    for c in curves:
        ax2.plot (c[0], c[1], color = c[2], linewidth = 2)
    ax2.legend (handles, names, facecolor = "white")
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotScatter (baseDir, selection, allClusters,  outputPath = ""):
    base = plt.cm.viridis; colors = base (np.linspace (0, 1, 256)); colors[0] = np.array ([1, 1, 1, 1])
    cmap = mcolors.LinearSegmentedColormap.from_list ("vlag", [(0, "blue"), (0.25, "lightblue"), (0.5, "snow"), (0.75, "pink"), (1, "red")])
    fct = pd.Series ([0, -2, -2, -1, 0, 1, 2, 2], index = ["NA", "-INF", "--", "-", "o", "+", "++", "+INF"]) / 2
    features = set (selection.sort_values ("FV").groupby ("cluster").tail (5).reset_index ()["feature"])
    dummyData = pd.DataFrame ({"feature": "_dummy", "cluster": allClusters, "color": -1, "size": 0})
    dummyData.iloc[1, 2] = 1; dummyData.iloc[1, 3] = 1; dummyData["label"] = "NULL"
    pltData = pd.read_csv (os.path.join (baseDir, "FlowSets_input", "raw_counts", "log2FC.tsv"), index_col = None, sep = "\t")
    pltData = pltData.rename (columns = {col: col.replace ("*log2FC", "") for col in pltData.columns[1:]})
    pltData = (pltData.set_index ("feature") * fct).sum (axis = 1).reset_index (name = "color")
    pltData[["feature", "cluster"]] = pltData["feature"].str.split ("__", expand = True)
    pltData = pltData[["feature", "cluster", "color"]]
    size = pd.DataFrame (columns = allClusters, dtype = float)
    for cluster in allClusters:
        summary = pd.read_csv (os.path.join (baseDir, "FV_raw_expression", "numerator", f"submatrix_{cluster}", "evaluations", "summary.tsv"),
                               index_col = 0, sep = "\t").reset_index (names = "feature")
        tmp = pd.read_csv (os.path.join (baseDir, "FV_raw_expression", "denominator", f"submatrix_{cluster}", "evaluations", "summary.tsv"),
                           index_col = 0, sep = "\t").reset_index (names = "feature")
        summary = summary.merge (tmp, on = ["feature", "individual_concept"], how = "outer", suffixes = ("_numerator", "_denominator"))
        with np.errstate (divide = "ignore", invalid = "ignore"):
            size[cluster] = pd.Series (np.log2 (summary["deviation_numerator"] / summary["deviation_denominator"]).round (3).values,
                                       index = summary["feature"].values)
    size = size.reset_index (names = "feature").melt (id_vars = "feature", var_name = "cluster", value_name = "size")
    pltData = pltData.merge (size, on = ["feature", "cluster"], how = "inner")
    pltData = pltData.merge (selection[["feature", "cluster"]].rename (columns = {"cluster": "label"}), on = "feature", how = "left")\
        .replace (np.nan, "NULL")
    pltData["label"] = "[" + pltData["label"] + "] " + pltData["feature"]
    pltData = pd.concat ([pltData.loc[pltData["feature"].isin (features)].sort_values (["cluster", "label"]), dummyData],
                         axis = 0, ignore_index = True)
    fig, ax = plt.subplots (figsize = (7.5, 12.5))
    sns.scatterplot (pltData, x = "cluster", y = "label", hue = "color", palette = cmap,
                     size = "size", sizes = (0, 150), size_norm = (-2, 5), ax = ax)
    ax.set_ylim ((-0.5, len (features) - 0.5)); ax.set_xlabel (""); ax.set_ylabel (""); ax.get_legend ().remove (); ax.invert_yaxis ()
    sm = plt.cm.ScalarMappable (cmap = cmap, norm = plt.Normalize (-1, 1)); sm.set_array (list ())
    ax.figure.colorbar (sm, ax = ax, location = "bottom", orientation = "horizontal", pad = 0.05, label = "level of regulation")
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotFoldChange (baseDir, fileLog2F, clustering, outputPath = ""):
    dataSpecific = pd.read_csv (os.path.join (baseDir, "raw_counts", "data_cluster-specific_features.tsv"),
                                index_col = None, sep = "\t")
    testSpecific = pd.read_csv (os.path.join (baseDir, "DESeq2_cluster-specific_features.tsv"),
                                index_col = None, sep = "\t")
    testSpecific["FV"] = 1
    common = dataSpecific.merge (testSpecific, on = ["feature", "cluster", "regulation"], how = "outer",
                                 suffixes = ("_data", "_DESeq2"))
    log2FC = pd.read_csv (fileLog2F, index_col = 0, sep = "\t")
    average = log2FC.mask (~np.isfinite (log2FC))
    average = pd.DataFrame ({cluster: average[clustering[clustering == cluster].index].median (axis = 1, skipna = True)
                             for cluster in sorted (set (clustering))})
    average = average.reset_index (names = "feature").melt (id_vars = "feature", var_name = "cluster",
                                                            value_name = "average")
    common = common.merge (average, on = ["feature", "cluster"], how = "left")
    fig, ax = plt.subplots (figsize = (8, 4))
    ax.hist (common.dropna ()["average"], bins = 50, label = "by both data and DESeq2")
    ax.hist (common.loc[np.isnan (common["FV_DESeq2"]), "average"], bins = 50, color = "palegreen",
             alpha = 0.8, label = "only by data")
    ax.hist (common.loc[np.isnan (common["FV_data"]), "average"], bins = 50, color = "silver",
             alpha = 0.6, label = "only by DESeq2")
    ax.set_xlabel ("paired log2 fold change", size = 10); ax.set_ylabel ("number of values", size = 10)
    ax.legend (loc = "upper right"); fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotRidgeline (feature, baseDir, outputPath = ""):
    clustering = pd.read_csv (os.path.join (baseDir, "data", "metadata.tsv"), index_col = None, sep = "\t")
    clustering = pd.Series (clustering["context"].values, index = clustering["sample"].values)
    allClusters = sorted (set (clustering))
    allMethods = ["raw log2FC", "norm. log2FC", "DESeq2"]
    candidates = pd.read_csv (os.path.join (baseDir, "results", "DESeq2_cluster-specific_features.tsv"), index_col = None, sep = "\t")
    candidates["method"] = 2
    tmp = pd.read_csv (os.path.join (baseDir, "results", "raw_counts", "data_cluster-specific_features.tsv"), index_col = None, sep = "\t")
    tmp = tmp.drop ("FV", axis = 1); tmp["method"] = 0; candidates = pd.concat ([candidates, tmp], axis = 0, ignore_index = True)
    tmp = pd.read_csv (os.path.join (baseDir, "results", "norm_counts", "data_cluster-specific_features.tsv"), index_col = None, sep = "\t")
    tmp = tmp.drop ("FV", axis = 1); tmp["method"] = 1; candidates = pd.concat ([candidates, tmp], axis = 0, ignore_index = True)
    candidates = candidates.sort_values ("feature").set_index ("feature"); featureData = candidates.loc[feature]
    clusters = [featureData["cluster"]] if isinstance (featureData, pd.Series) else list (set (featureData["cluster"]))
    rawFC = pd.read_csv (os.path.join (baseDir, "data", "raw_counts", "paired_log2FC.tsv"), index_col = 0, sep = "\t")
    normFC = pd.read_csv (os.path.join (baseDir, "data", "norm_counts", "paired_log2FC.tsv"), index_col = 0, sep = "\t")
    rawFC = rawFC.loc[feature].replace ([-np.inf, np.inf], np.nan); normFC = normFC.loc[feature].replace ([-np.inf, np.inf], np.nan)
    log2FC = pd.read_csv (os.path.join (baseDir, "data", "DESeq2_log2FC.tsv"), index_col = 0, sep = "\t").loc[feature].replace ([-np.inf, np.inf], np.nan)
    padj = 10 ** -pd.read_csv (os.path.join (baseDir, "data", "DESeq2_padj.tsv"), index_col = 0, sep = "\t").loc[feature]
    rawTumor = {C: pd.read_csv (os.path.join (baseDir, "data", "raw_counts", "numerator_submatrices", f"submatrix_{C}.tsv"),
                                index_col = 0, sep = "\t").loc[feature]
                for C in allClusters}
    rawNormal = {C: pd.read_csv (os.path.join (baseDir, "data", "raw_counts", "denominator_submatrices", f"submatrix_{C}.tsv"),
                                 index_col = 0, sep = "\t").loc[feature]
                 for C in allClusters}
    normTumor = {C: pd.read_csv (os.path.join (baseDir, "data", "norm_counts", "numerator_submatrices", f"submatrix_{C}.tsv"),
                                 index_col = 0, sep = "\t").loc[feature]
                 for C in allClusters}
    normNormal = {C: pd.read_csv (os.path.join (baseDir, "data", "norm_counts", "denominator_submatrices", f"submatrix_{C}.tsv"),
                                  index_col = 0, sep = "\t").loc[feature]
                  for C in allClusters}
    rawParams = dict (); rawMin = list (); rawMax = list (); rawMinCutoffs = dict (); rawMaxCutoffs = dict ()
    normParams = dict (); normMin = list (); normMax = list (); normMinCutoffs = dict (); normMaxCutoffs = dict ()
    for cluster in sorted (set (clustering)):
        with open (os.path.join (baseDir, "FV_raw_expression", "concepts", f"submatrix_{cluster}", "concepts_detailed.json")) as f:
            concepts = json.load (f)
        default = concepts.get ("DEFAULT", dict ()).get ("MEDIUM", [list (), "", "", 0])
        rawParams[cluster] =  concepts.get (feature, dict ()).get ("MEDIUM", default)[0]
        rawMin.append (min ([rawTumor[cluster].min (skipna = True), rawNormal[cluster].min (skipna = True)]))
        rawMax.append (max ([rawTumor[cluster].max (skipna = True), rawNormal[cluster].max (skipna = True)]))
        tmp = concepts.get (feature, dict ()).get ("MIN-NOISE", -np.inf)
        rawMinCutoffs[cluster] = -np.inf if isinstance (tmp, str) else tmp
        tmp = concepts.get (feature, dict ()).get ("MAX-NOISE", -np.inf)
        rawMaxCutoffs[cluster] = np.inf if isinstance (tmp, str) else tmp
        with open (os.path.join (baseDir, "FV_norm_expression", "concepts", f"submatrix_{cluster}", "concepts_detailed.json")) as f:
            concepts = json.load (f)
        default = concepts.get ("DEFAULT", dict ()).get ("MEDIUM", [list (), "", "", 0])
        normParams[cluster] = concepts.get (feature, dict ()).get ("MEDIUM", default)[0]
        normMin.append (min ([normTumor[cluster].min (skipna = True), normNormal[cluster].min (skipna = True)]))
        normMax.append (max ([normTumor[cluster].max (skipna = True), normNormal[cluster].max (skipna = True)]))
        tmp = concepts.get (feature, dict ()).get ("MIN-NOISE", -np.inf)
        normMinCutoffs[cluster] = -np.inf if isinstance (tmp, str) else tmp
        tmp = concepts.get (feature, dict ()).get ("MAX-NOISE", -np.inf)
        normMaxCutoffs[cluster] = np.inf if isinstance (tmp, str) else tmp
    rawRangeExp = [np.floor (min (rawMin)), np.ceil (max (rawMax))]; normRangeExp = [np.floor (min (normMin)), np.ceil (max (normMax))]
    rawRangeFC = [np.floor (rawFC.dropna ().min ()), np.ceil (rawFC.dropna ().max ())]
    normRangeFC = [np.floor (normFC.dropna ().min ()), np.ceil (normFC.dropna ().max ())]
    rawDummy = np.linspace (*rawRangeExp, 1000); normDummy = np.linspace (*normRangeExp, 1000)
    fig, axs = plt.subplots (len (allClusters), 4, sharex = False, sharey = False, layout = "constrained", figsize = (12, 1.5 * len (allClusters)))
    for idx in range (len (allClusters)):
        cluster = allClusters[idx]; sampleList = clustering[clustering == cluster].index
        if isinstance (featureData, pd.Series):
            methods = [featureData["method"]] if featureData["cluster"] == allClusters[idx] else list ()
        else:
            methods = featureData.loc[featureData["cluster"] == allClusters[idx], "method"].tolist ()
        if not np.isnan (rawTumor[cluster]).all ():
            axs[idx, 0].hist (rawTumor[cluster].dropna (), bins = 20, color = "firebrick", alpha = 0.6)
        if not np.isnan (rawNormal[cluster]).all ():
            axs[idx, 0].hist (rawNormal[cluster].dropna (), bins = 20, color = "steelblue", alpha = 0.6)
        axs[idx, 0].set_xlim (rawRangeExp)
        if rawMinCutoffs[cluster] > rawRangeExp[0]:            axs[idx, 0].axvline (rawMinCutoffs[cluster], c = "darkgray", linestyle = "dashed")
        if rawMaxCutoffs[cluster] < rawRangeExp[1]:
            axs[idx, 0].axvline (rawMaxCutoffs[cluster], c = "darkgray", linestyle = "dashed")
        ax = axs[idx, 0].twinx (); ax.set_ylim ((0, 1.05)); mu, sigma = rawParams[cluster]
        if np.isnan (mu) or np.isnan (sigma):
            ax.plot (list (), list (), color = "black", label = r"$\mu$: " + f"{mu:.2f}\n" + r"$\sigma$: " + f"{sigma:.2f}")
        else:
            ax.plot (rawDummy, np.exp (-(rawDummy - mu) ** 2 / (2 * sigma ** 2)), color = "black",
                     label = r"$\mu$: " + f"{mu:.2f}\n" + r"$\sigma$: " + f"{sigma:.2f}")
        ax.legend (loc = "upper left", title = cluster, fontsize = 8, title_fontsize = 9, handlelength = 1)
            
        if not np.isnan (rawFC[sampleList]).all ():
            axs[idx, 1].hist (rawFC[sampleList].dropna (), bins = 20, color = "silver")
        axs[idx, 1].set_xlim (rawRangeFC)
        axs[idx, 1].axvline (log2FC[cluster], color = "black", linestyle = "dashed")
        for i in range (len (allMethods)):
            if i in methods:
                axs[idx, 1].scatter (list (), list (), s = 20, c = "black", edgecolors = "black", marker = "o", label = allMethods[i])
            else:
                axs[idx, 1].scatter (list (), list (), s = 20, c = "white", edgecolors = "black", marker = "o", label = allMethods[i])
        axs[idx, 1].legend (loc = "upper left", fontsize = 8)
        if not np.isnan (normTumor[cluster]).all ():
            axs[idx, 2].hist (normTumor[cluster].dropna (), bins = 20, color = "firebrick", alpha = 0.6)
        if not np.isnan (normNormal[cluster]).all ():
            axs[idx, 2].hist (normNormal[cluster].dropna (), bins = 20, color = "steelblue", alpha = 0.6)
        axs[idx, 2].set_xlim (normRangeExp)
        if normMinCutoffs[cluster] > normRangeExp[0]:
            axs[idx, 2].axvline (normMinCutoffs[cluster], c = "darkgray", linestyle = "dashed")
        if normMaxCutoffs[cluster] < normRangeExp[1]:
            axs[idx, 2].axvline (normMaxCutoffs[cluster], c = "darkgray", linestyle = "dashed")
        ax = axs[idx, 2].twinx (); ax.set_ylim ((0, 1.05)); mu, sigma = normParams[cluster]
        if np.isnan (mu) or np.isnan (sigma):
            ax.plot (list (), list (), color = "black", label = r"$\mu$: " + f"{mu:.2f}\n" + r"$\sigma$: " + f"{sigma:.2f}")
        else:
            ax.plot (normDummy, np.exp (-(normDummy - mu) ** 2 / (2 * sigma ** 2)), color = "black",
                     label = r"$\mu$: " + f"{mu:.2f}\n" + r"$\sigma$: " + f"{sigma:.2f}")
        ax.legend (loc = "upper left", title = cluster, fontsize = 8, title_fontsize = 9, handlelength = 1)
        if not np.isnan (normFC[sampleList]).all ():
            axs[idx, 3].hist (normFC[sampleList].dropna (), bins = 20, color = "silver")
        axs[idx, 3].set_xlim (normRangeFC)
        axs[idx, 3].axvline (log2FC[cluster], color = "black", linestyle = "dashed")
        axs[idx, 3].text (0.05, 0.72,
                          f"DE log2FC: {log2FC[cluster]:.2f}\nadj. p-val.: {padj[cluster]:.2f}",
                          transform = axs[idx, 3].transAxes, fontsize = 8,
                          bbox = {"facecolor": "white", "edgecolor": "silver", "alpha": 0.8})
        if cluster in clusters:
            axs[idx, 0].set_facecolor ("lavender"); axs[idx, 1].set_facecolor ("lavender")
            axs[idx, 2].set_facecolor ("lavender"); axs[idx, 3].set_facecolor ("lavender")
    axs[idx, 0].set_xlabel ("log2 (raw count)\nprimary tumour: red\nsolid normal tissue: blue", size = 10)
    axs[idx, 1].set_xlabel ("log2 fold change of raw values\ndashed line from DESeq2", size = 10)
    axs[idx, 2].set_xlabel ("log2 (norm. count)\nprimary tumour: red\nsolid normal tissue: blue", size = 10)
    axs[idx, 3].set_xlabel ("log2 fold change of norm. values\ndashed line from DESeq2", size = 10)
    fig.supylabel ("number of sample pairs", size = 10); fig.suptitle (feature, size = 12.5)
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def violinDepth (baseDir, allClusters, clustering, minLevel = -np.inf, maxLevel = np.inf, outputPath = ""):
    tumor = pd.concat ([pd.read_csv (os.path.join (baseDir, "numerator_submatrices", f"submatrix_{cluster}.tsv"),
                                     index_col = 0, sep = "\t")
                        for cluster in allClusters], axis = 1)
    normal = pd.concat ([pd.read_csv (os.path.join (baseDir, "denominator_submatrices", f"submatrix_{cluster}.tsv"),
                                      index_col = 0, sep = "\t")
                         for cluster in allClusters], axis = 1)
    pltData = pd.concat ([pd.DataFrame ({"depth": (2 ** tumor.mask ((tumor <= minLevel) | (tumor >= maxLevel), -np.inf))\
                                         .astype (int).sum (axis = 0),
                                         "type": "primary tumour"}).reset_index (),
                          pd.DataFrame ({"depth": (2 ** normal.mask ((normal <= minLevel) | (normal >= maxLevel), -np.inf))\
                                         .astype (int).sum (axis = 0),
                                         "type": "solid normal tissue"}).reset_index ()],
                         axis = 0, ignore_index = True)
    pltData["cluster"] = clustering[pltData["index"]].values
    fig, ax = plt.subplots (figsize = (6, 4))
    sns.violinplot (pltData, x = "depth", y = "cluster", order = sorted (set (clustering)), hue = "type",
                    palette = {"primary tumour": "firebrick", "solid normal tissue": "steelblue"},
                    inner = "quart", log_scale = 10, legend = None, ax = ax)
    ax.set_xlabel ("sequencing depth [primary tumour: red; solid normal tissue: blue]", size = 10); ax.set_ylabel ("")
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def histExpression (baseDir, allClusters, axis_label, minLevel = -np.inf, maxLevel = np.inf, outputPath = ""):
    tumor = pd.concat ([pd.read_csv (os.path.join (baseDir, "numerator_submatrices", f"submatrix_{cluster}.tsv"),
                                     index_col = 0, sep = "\t")
                        for cluster in allClusters], axis = 1)
    normal = pd.concat ([pd.read_csv (os.path.join (baseDir, "denominator_submatrices", f"submatrix_{cluster}.tsv"),
                                      index_col = 0, sep = "\t")
                         for cluster in allClusters], axis = 1)
    fig, ax = plt.subplots (figsize = (6, 3))
    ax.hist (tumor.mask (~np.isfinite (tumor)).melt ()["value"].dropna (), bins = 100, color = "firebrick")
    ax.hist (normal.mask (~np.isfinite (normal)).melt ()["value"].dropna (), bins = 100,
             alpha = 0.6, color = "steelblue")
    if np.isfinite (minLevel):
        ax.axvline (minLevel, color = "black", linestyle = "dashed")
    if np.isfinite (maxLevel):
        ax.axvline (maxLevel, color = "black", linestyle = "dashed")
    ax.set_xlabel (f"{axis_label} [primary tumour: red; solid normal tissue: blue]", size = 10)
    ax.set_ylabel ("number of expression values", size = 10); fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotCardinalityDE (dir, allClusters, allSets, axis_label, outputPath = ""):
    allFV = np.array ([pd.read_csv (os.path.join (dir, file), index_col = 0, sep = "\t").loc[allClusters].to_numpy ()
                       for file in os.listdir (dir) if file.startswith ("fuzzyValues_") and file.endswith (".tsv")])
    mainFS = allFV.argmax (axis = 2)
    pctMainFS = pd.DataFrame.from_dict ({allSets[idx]: (mainFS == idx).mean (axis = 0) for idx in range (len (allSets))},
                                        orient = "index", columns = allClusters)
    pctMainFS = pctMainFS.reset_index ().melt (id_vars = "index")
    fig, ax = plt.subplots (figsize = (8, 4))
    sns.lineplot (pctMainFS, x = "index", y = "value", hue = "variable", hue_order = allClusters,
                  style = "variable", markers = True, dashes = False, ax = ax)
    ax.set_xlabel (""); ax.set_ylabel (f"percent of {axis_label}", size = 10)
    ax.legend (loc = "upper left", title = None, ncol = 3); fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotCardinalityFC (dir, clustering, allSets, axis_label, outputPath = ""):
    percent = dict ()
    for cluster in sorted (set (clustering)):
        counts = pd.Series (0, index = allSets); n = 0
        for sample in clustering[clustering == cluster].index:
            file = f"fuzzyValues_{sample}.tsv"
            num = pd.read_csv (os.path.join (dir, file), index_col = 0, sep = "\t"); dim = num.shape[0]
            num = num.idxmax (axis = 1).value_counts ()
            counts += pd.Series ({FS: num.get (FS, 0) for FS in allSets}); n += 1
        percent[cluster] = counts / (dim * n)
    percent = pd.DataFrame (percent).T.reset_index ().melt (id_vars = "index")
    fig, ax = plt.subplots (figsize = (8, 4))
    sns.lineplot (percent, x = "variable", y = "value", hue = "index", hue_order = sorted (set (clustering)),
                  style = "index", markers = True, dashes = False, ax = ax)
    ax.set_xlabel (""); ax.set_ylabel (f"percent of {axis_label}\nwith highest memberships", size = 10)
    ax.legend (loc = "upper left", title = None, ncol = 3); fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotCardinality (dir, clustering, allSets, axis_label, outputPath = ""):
    percent = dict ()
    for cluster in sorted (set (clustering)):
        subdir = os.path.join (dir, f"submatrix_{cluster}", "fuzzy_values")
        counts = pd.Series (0, index = allSets); n = 0
        for file in os.listdir (subdir):
            if (not file.startswith ("fuzzyValues_")) or (not file.endswith (".tsv")):
                continue
            num = pd.read_csv (os.path.join (subdir, file), index_col = 0, sep = "\t").idxmax (axis = 1).value_counts ()
            counts += pd.Series ({FS: num.get (FS, 0) for FS in allSets}); n += 1
        percent[cluster] = counts / ((clustering == cluster).sum () * n)
    percent = pd.DataFrame (percent).T.reset_index ().melt (id_vars = "index")
    fig, ax = plt.subplots (figsize = (8, 4))
    sns.lineplot (percent, x = "variable", y = "value", hue = "index", hue_order = sorted (set (clustering)),
                  style = "index", markers = True, dashes = False, ax = ax)
    ax.set_xlabel (""); ax.set_ylabel (f"percent of {axis_label}\nwith highest memberships", size = 10)
    ax.legend (loc = "upper left", title = None, ncol = 3); fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotFeatureMu (dir, allClusters, axis_label, outputPath = ""):
    fit = pd.DataFrame (columns = allClusters, dtype = float); default = list ()
    for cluster in allClusters:
        with open (os.path.join (dir, f"submatrix_{cluster}", "concepts_detailed.json")) as f:
            concepts = json.load (f)
        default.append (concepts["DEFAULT"]["MEDIUM"][0][0])
        fit[cluster] = pd.Series ([concepts[x].get ("MEDIUM", [[np.nan]])[0][0] for x in concepts.keys ()],
                                  index = concepts.keys ()).drop ("DEFAULT")
    fig, ax = plt.subplots (figsize = (8, 4))
    sns.violinplot (fit.melt (), x = "value", y = "variable", color = "silver", inner = "quart", legend = None, ax = ax)
    for idx in range (len (default)):
        ax.plot (default[idx], idx, c = "red", marker = "x")
        ax.annotate (default[idx], (default[idx], idx - 0.25))
    ax.set_xlabel (r"fitted $\mu$ for " + axis_label + r" (cross: $\mu$ in default fuzzy concept)", size = 10); ax.set_ylabel ("")
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotSampleMu (file, clustering, axis_label, outputPath = ""):
    with open (file) as f:
        concepts = json.load (f)
    params = pd.DataFrame ({"mu": [concepts[x].get ("o", [[np.nan]])[0][0] for x in clustering.index],
                            "cluster": clustering}, index = clustering.index).sort_values ("cluster")
    fig, ax = plt.subplots (figsize = (8, 4))
    ax.axvline (0, color = "red", linestyle = "dashed", alpha = 0.8)
    sns.violinplot (params, x = "mu", y = "cluster", color = "silver", inner = "quart", ax = ax)
    ax.set_xlabel (r"fitted $\mu$ for " + axis_label + " (red line: zero)", size = 10); ax.set_ylabel ("")
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def plotFitFlowFV (baseDir, regulation, color, outputPath = ""):
    rawFV = pd.read_csv (os.path.join (baseDir, "raw_counts", "FlowSets_memberships.tsv"), index_col = None, sep = "\t")
    rawFV["regulation"] = rawFV[["DOWN", "UP"]].idxmax (axis = 1); rawFV["FV"] = rawFV[["DOWN", "UP"]].max (axis = 1)
    normFV = pd.read_csv (os.path.join (baseDir, "norm_counts", "FlowSets_memberships.tsv"), index_col = None, sep = "\t")
    normFV["regulation"] = normFV[["DOWN", "UP"]].idxmax (axis = 1); normFV["FV"] = normFV[["DOWN", "UP"]].max (axis = 1)
    memberships = rawFV[["feature", "cluster", "regulation", "FV"]].merge (normFV[["feature", "cluster", "regulation", "FV"]],
                                                                           on = ["feature", "cluster"], how = "outer",
                                                                           suffixes = ("_raw", "_norm")).replace (np.nan, 0)
    FV = memberships.loc[(memberships["regulation_raw"] == regulation) & (memberships["regulation_norm"] == regulation)]
    xValues = np.linspace (0, 1, 1000); allClusters = sorted (set (memberships["cluster"]))
    fig, axs = plt.subplots (3, 3, sharex = False, sharey = False, figsize = (9, 9))
    for i in range (3):
        for j in range (3):
            cluster = allClusters[3 * i + j]; tmp = FV.loc[FV["cluster"] == cluster]
            res = stats.linregress (x = tmp["FV_raw"], y = tmp["FV_norm"])
            sns.scatterplot (tmp, x = "FV_raw", y = "FV_norm", color = color, legend = None, ax = axs[i, j])
            axs[i, j].plot (xValues, res.intercept + res.slope * xValues, color = "black", linestyle = "dashed",
                            label = f"y={res.intercept:.2f}+{res.slope:.2f}x\n" + r"$r^2$: " + f"{(res.rvalue ** 2):.2f}")
            axs[i, j].set_xlabel (""); axs[i, j].set_ylabel (""); axs[i, j].set_title (cluster, size = 12)
            axs[i, j].set_xlim ((0, 1)); axs[i, j].set_ylim ((0, 1)); axs[i, j].legend (loc = "upper left")
    fig.supxlabel ("flow membership from raw log2 fold change", size = 10)
    fig.supylabel ("flow membership from norm. log2 fold change", size = 10)
    fig.tight_layout ()
    if outputPath == "":
        plt.show ()
    else:
        plt.savefig (outputPath); plt.close ()



def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--dir", required = True, type = str, help = "Base directory")
    parser.add_argument ("--index", required = True, type = str, help = "Figure index in paper")
    args = parser.parse_args ()
    
    os.chdir (args.dir); prefix, figureIndex, suffix = parseIndex (args.index)

    miRNA_cluster = pd.read_csv ("./miRNA/data/metadata.tsv", index_col = None, sep = "\t")
    miRNA_cluster = pd.Series (miRNA_cluster["context"].values, index = miRNA_cluster["sample"].values)
    RNA_cluster = pd.read_csv ("./RNA/data/metadata.tsv", index_col = None, sep = "\t")
    RNA_cluster = pd.Series (RNA_cluster["context"].values, index = RNA_cluster["sample"].values)
    allClusters = sorted (set (miRNA_cluster) & set (RNA_cluster))
    exp_FS = ["NO_EXP", "NOISE", "LOW", "low", "MEDIUM", "high", "HIGH"]
    log2FC_FS = ["NA", "-INF", "--", "-", "o", "+", "++", "+INF"]
    padj_FS = ["o", "*", "**", "***", "****"]

    if (prefix == "F" and figureIndex == 3) or figureIndex == 5 or (prefix == "A" and figureIndex == 11):
        regFA, flows = generateFA ("./miRNA/FlowSets_input/raw_counts/", "./miRNA/FV_DESeq2/padj/fuzzy_values/")
        tarFA, _ = generateFA ("./RNA/FlowSets_input/raw_counts/", "./RNA/FV_DESeq2/padj/fuzzy_values/")

    if prefix == "F":
        if figureIndex == 1:
            plotVolcano ("./RNA/data/", "./RNA/FV_DESeq2/", outputPath = "./visualizations/F2_volcano.png")
            log2FC = pd.read_csv ("./RNA/data/DESeq2_log2FC.tsv", index_col = 0, sep = "\t").clip (lower = -10, upper = 10)
            with open ("./RNA/FV_DESeq2/concepts_log2FC/concepts_detailed.json") as f:
                log2FC_concept = json.load (f)["DEFAULT"]
            padj = pd.read_csv ("./RNA/data/DESeq2_padj.tsv", index_col = 0, sep = "\t")
            with open ("./RNA/FV_DESeq2/concepts_padj/concepts_detailed.json") as f:
                padj_concept = json.load (f)["DEFAULT"]
            plotConcept (log2FC, log2FC_concept, "DESeq2 average log2 fold change", outputPath = "./visualizations/F2_DESeq2_log2FC.png")
            plotConcept (padj, padj_concept, "-log10 (DESeq2 corrected p-value)", outputPath = "./visualizations/F2_DESeq2_padj.png")
        if figureIndex == 2:
            if suffix == "a" or suffix == "":
                regFA.plot_flows (classFontsize = 9, seriesFontsize = 10, figsize = (10, 6),
                                  outfile = "./visualizations/F3a_FlowSets_miRNA")
            if suffix == "b" or suffix == "":
                tarFA.plot_flows (classFontsize = 9, seriesFontsize = 10, figsize = (10, 6),
                                  outfile = "./visualizations/F3b_FlowSets_RNA")
        if figureIndex == 3:
            if suffix == "a" or suffix == "":
                markers = pd.read_csv ("./miRNA/results/raw_counts/data_marker_features.tsv", index_col = None, sep = "\t")
                plotScatter ("./miRNA/", markers, allClusters,  outputPath = "./visualizations/F3a_miRNA_markers.png")
                markers = pd.read_csv ("./RNA/results/raw_counts/data_marker_features.tsv", index_col = None, sep = "\t")
                plotScatter ("./RNA/", markers, allClusters, outputPath = "./visualizations/F3b_RNA_markers.png")
        if figureIndex == 4:
            if suffix == "a" or suffix == "":
                plotFoldChange ("./miRNA/results/", "./miRNA/data/raw_counts/paired_log2FC.tsv", miRNA_cluster,
                                outputPath = "./visualizations/F4a_miRNA_raw_comparison.png")
            if suffix == "b" or suffix == "":
                plotFoldChange ("./RNA/results/", "./RNA/data/raw_counts/paired_log2FC.tsv", RNA_cluster,
                                outputPath = "./visualizations/F4b_RNA_raw_comparison.png")
        if figureIndex == 5:
            if suffix == "a" or suffix == "":
                plotRidgeline ("hsa-mir-221", "./miRNA/", outputPath = "./visualizations/F5a_ridgeline_hsa-mir-221.png")
            if suffix == "b" or suffix == "":
                plotRidgeline ("KIT", "./RNA/", outputPath = "./visualizations/F5b_ridgeline_KIT.png")
    elif prefix == "A":
        if figureIndex == 1:
            if suffix == "a" or suffix == "":
                violinDepth ("./miRNA/data/raw_counts/", allClusters, miRNA_cluster, outputPath = "./visualizations/A1a_miRNA_raw_depth.png")
            if suffix == "b" or suffix == "":
                histExpression ("./miRNA/data/raw_counts/", allClusters, "log2 (raw count)", minLevel = 3, maxLevel = np.inf,
                                outputPath = "./visualizations/A1b_miRNA_raw_expression.png")
            if suffix == "c" or suffix == "":
                violinDepth ("./miRNA/data/norm_counts/", allClusters, miRNA_cluster, outputPath = "./visualizations/A1c_miRNA_norm_depth.png")
            if suffix == "d" or suffix == "":
                histExpression ("./miRNA/data/norm_counts/", allClusters, "log2 (norm. count)", minLevel = 3, maxLevel = np.inf,
                                outputPath = "./visualizations/A1d_miRNA_norm_expression.png")
            if suffix == "e" or suffix == "":
                violinDepth ("./RNA/data/raw_counts/", allClusters, RNA_cluster, outputPath = "./A1e_RNA_raw_depth.png")
            if suffix == "f" or suffix == "":
                histExpression ("./RNA/data/raw_counts/", allClusters, "log2 (raw count)", minLevel = 3, maxLevel = np.inf,
                                outputPath = "./visualizations/A1f_RNA_raw_expression.png")
            if suffix == "g" or suffix == "":
                violinDepth ("./RNA/data/norm_counts/", allClusters, RNA_cluster, outputPath = "./visualizations/A1g_RNA_norm_depth.png")
            if suffix == "h" or suffix == "":
                histExpression ("./RNA/data/norm_counts/", allClusters, "log2 (norm. count)", minLevel = 3, maxLevel = np.inf,
                                outputPath = "./visualizations/A1h_RNA_norm_expression.png")
        if figureIndex == 2:
            if suffix == "a" or suffix == "":
                plotCardinalityDE ("./miRNA/FV_DESeq2/padj/fuzzy_values/", allClusters, ["NULL", "NA"] + padj_FS,
                                   "-log10 (DESeq2 corrected p-value)", outputPath = "./visualizations/A2a_miRNA_DESeq2_padj_cardinality.png")
            if suffix == "b" or suffix == "":
                plotCardinalityDE ("./RNA/FV_DESeq2/padj/fuzzy_values/", allClusters, ["NULL", "NA"] + padj_FS,
                                   "-log10 (DESeq2 corrected p-values)", outputPath = "./visualizations/A2b_RNA_DESeq2_padj_cardinality.png")
        if figureIndex == 3:
            if suffix == "a" or suffix == "":
                plotCardinalityFC ("./miRNA/FV_raw_log2FC/final_fuzzy_values/", miRNA_cluster, log2FC_FS, "raw log2FC",
                                   outputPath = "./visualizations/A3a_miRNA_raw_log2FC_cardinality.png")
            if suffix == "b" or suffix == "":
                plotCardinalityFC ("./RNA/FV_raw_log2FC/final_fuzzy_values/", RNA_cluster, log2FC_FS, "raw log2FC",
                                   outputPath = "./visualizations/A3b_RNA_raw_log2FC_cardinality.png")
            if suffix == "c" or suffix == "":
                plotCardinality ("./miRNA/FV_raw_expression/numerator/", miRNA_cluster, exp_FS, "raw counts",
                                 outputPath = "./visualizations/A3c_miRNA_raw_tumor_cardinality.png")
            if suffix == "d" or suffix == "":
                plotCardinality ("./RNA/FV_raw_expression/numerator/", RNA_cluster, exp_FS, "raw counts",
                                 outputPath = "./visualizations/A3d_RNA_raw_tumor_cardinality.png")
            if suffix == "e" or suffix == "":
                plotCardinality ("./miRNA/FV_raw_expression/denominator/", miRNA_cluster, exp_FS, "raw counts",
                                 outputPath = "./visualizations/A3e_miRNA_raw_normal_cardinality.png")
            if suffix == "f" or suffix == "":
                plotCardinality ("./RNA/FV_raw_expression/denominator/", RNA_cluster, exp_FS, "raw counts",
                                 outputPath = "./visualizations/A3f_RNA_raw_normal_cardinality.png")
        if figureIndex == 4:
            if suffix == "a" or suffix == "":
                plotFeatureMu ("./miRNA/FV_raw_expression/concepts/", allClusters, "miRNA raw expression",
                               outputPath = "./visualizations/A4a_miRNA_raw_fitted.png")
            if suffix == "b" or suffix == "":
                plotFeatureMu ("./RNA/FV_raw_expression/concepts/", allClusters, "RNA raw expression",
                               outputPath = "./visualizations/A4b_RNA_raw_fitted.png")
            if suffix == "c" or suffix == "":
                plotSampleMu ("./miRNA/FV_raw_log2FC/concepts_detailed.json", miRNA_cluster, "miRNA raw log2 fold change",
                              outputPath = "./visualizations/A4c_miRNA_raw_log2FC_fitted.png")
            if suffix == "d" or suffix == "":
                plotSampleMu ("./RNA/FV_raw_log2FC/concepts_detailed.json", RNA_cluster, "RNA raw log2 fold change",
                              outputPath = "./visualizations/A4d_RNA_raw_log2FC_fitted.png")
        if figureIndex == 5:
            suffixList = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
            for i in range (len (allClusters)):
                if suffix == suffixList[i] or suffix == "":
                    subprocess.call (["Rscript", "venn_visualization.R",
                                      "./miRNA/results/raw_counts/data_cluster-specific_features.tsv",
                                      "./miRNA/results/norm_counts/data_cluster-specific_features.tsv",
                                      "./miRNA/results/DESeq2_cluster-specific_features.tsv",
                                      allClusters[i], "DOWN",
                                      f"./visualizations/A5{suffixList[i]}_{allClusters[i]}_DOWN"])
        if figureIndex == 6:
            suffixList = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
            for i in range (len (allClusters)):
                if suffix == suffixList[i] or suffix == "":
                    subprocess.call (["Rscript", "venn_visualization.R",
                                      "./miRNA/results/raw_counts/data_cluster-specific_features.tsv",
                                      "./miRNA/results/norm_counts/data_cluster-specific_features.tsv",
                                      "./miRNA/results/DESeq2_cluster-specific_features.tsv",
                                      allClusters[i], "UP",
                                      f"./visualizations/A6{suffixList[i]}_{allClusters[i]}_UP"])
        if figureIndex == 7:
            suffixList = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
            for i in range (len (allClusters)):
                if suffix == suffixList[i] or suffix == "":
                    subprocess.call (["Rscript", "venn_visualization.R",
                                      "./RNA/results/raw_counts/data_cluster-specific_features.tsv",
                                      "./RNA/results/norm_counts/data_cluster-specific_features.tsv",
                                      "./RNA/results/DESeq2_cluster-specific_features.tsv",
                                      allClusters[i], "DOWN",
                                      f"./visualizations/A7{suffixList[i]}_{allClusters[i]}_DOWN"])
        if figureIndex == 8:
            suffixList = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
            for i in range (len (allClusters)):
                if suffix == suffixList[i] or suffix == "":
                    subprocess.call (["Rscript", "venn_visualization.R",
                                      "./RNA/results/raw_counts/data_cluster-specific_features.tsv",
                                      "./RNA/results/norm_counts/data_cluster-specific_features.tsv",
                                      "./RNA/results/DESeq2_cluster-specific_features.tsv",
                                      allClusters[i], "UP",
                                      f"./visualizations/A8{suffixList[i]}_{allClusters[i]}_UP"])
        if figureIndex == 9:
            if suffix == "a" or suffix == "":
                plotFitFlowFV ("./miRNA/results/", "DOWN", "steelblue", outputPath = "./visualizations/A9a_miRNA_DOWN_comparison.png")
            if suffix == "b" or suffix == "":
                plotFitFlowFV ("./miRNA/results/", "UP", "firebrick", outputPath = "./visualizations/A9b_miRNA_UP_comparison.png")
            if suffix == "c" or suffix == "":
                plotFitFlowFV ("./RNA/results/", "DOWN", "steelblue", outputPath = "./visualizations/A9c_RNA_DOWN_comparison.png")
            if suffix == "d" or suffix == "":
                plotFitFlowFV ("./RNA/results/", "UP", "firebrick", outputPath = "./visualizations/A9d_RNA_UP_comparison.png")
        if figureIndex == 10:
            if suffix == "a" or suffix == "":
                regFA.plot_flows (genes = ["hsa-mir-221__THCA"], specialColors = {f: "firebrick" for f in flows["up"]},
                                  classFontsize = 9, seriesFontsize = 10, figsize = (10, 6),
                                  outfile = "./visualizations/A10a_hsa-mir-221")
            if suffix == "b" or suffix == "":
                tarFA.plot_flows (genes = ["KIT__THCA"], specialColors = {f: "steelblue" for f in flows["down"]},
                                  classFontsize = 9, seriesFontsize = 10, figsize = (10, 6),
                                  outfile = "./visualizations/A10b_KIT")
    else:
        raise ValueError



if __name__ == "__main__":
    main ()


