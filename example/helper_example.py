import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plotVolcano (log2FC, padj, log2FC_concept, padj_concept):
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
        mask = (pltData["log2FC"] > log2FC_cutoff[i]) & (pltData["log2FC"] < log2FC_cutoff[i + 1]) & (pltData["padj"] > padj_cutoff[j]) & (pltData["padj"] < padj_cutoff[j + 1])
        pltData.loc[mask, "class"] = f"{i}_{j}"; coords[f"{i}_{j}"] = [(log2FC_cutoff[i] + log2FC_cutoff[i + 1]) / 2, (padj_cutoff[j] + padj_cutoff[j + 1]) / 2]
    labels = pltData.value_counts ("class"); labels = {key: f"{(labels.get (key, 0) / pltData.shape[0]):.2%}\n({labels.get (key, 0)})" for key in coords.keys ()}
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
    ax.set_title (f"{log2FC.shape[0]} features in {log2FC.shape[1]} comparisons", size = 12)
    fig.tight_layout (); plt.show ()



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



def plotConcept (mtx, concept, axis_label):
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
    fig.tight_layout (); plt.show ()
    

