import os
import argparse
import numpy as np
import pandas as pd

# python main_fuzzyRule.py --paired pairedDirectory --expression expressionDirectory --rules fuzzyRules \
#                          --fuzzySets fuzzySetsForFC --output outputDirectory


def readFV (dir, nameDict = dict ()):
    allFV = list ()
    featureList = nameDict["feature"] if "feature" in nameDict else sorted ([file[12:-4] for file in os.listdir (dir)])
    for feature in featureList:
        FV = pd.read_csv (os.path.join (dir, f"fuzzyValues_{feature}.tsv"), index_col = 0, sep = "\t")
        if "sample" in nameDict and "FS" in nameDict:
            allFV.append (FV.loc[nameDict["sample"], nameDict["FS"]].to_numpy ())
        else:
            allFV.append (FV.to_numpy ())
    allFV = np.array (allFV); nameDict = {"feature": featureList, "sample": list (FV.index), "FS": list (FV.columns)}
    return allFV, nameDict



def main ():
    parser = argparse.ArgumentParser ()
    parser.add_argument ("--paired", type = str, required = True, help = "Directory of fuzzy values for paired log2FC")
    parser.add_argument ("--expression", type = str, required = True, help = "Directory of fuzzy values for expression values")
    parser.add_argument ("--rules", type = str, required = True, help = "Table of fuzzy rules (CSV)")
    parser.add_argument ("--fuzzySets", type = str, required = True, help = "String of ordered list for log2FC fuzzy variables")
    parser.add_argument ("--output", type = str, required = True, help = "Ouput directory for fuzzy fold changes")
    args = parser.parse_args ()

    if not os.path.exists (args.output):
        os.makedirs (args.output, exist_ok = True)

    rules = pd.read_csv (args.rules, index_col = 0, sep = ",").replace (np.nan, "NA").replace (np.inf, "+INF").replace (-np.inf, "-INF")
    allSets = args.fuzzySets.split (","); print (allSets); print (rules)
    mask = {FS: rules.to_numpy () == FS for FS in allSets}
    for dir in set (os.listdir (os.path.join (args.expression, "numerator"))) & set (os.listdir (os.path.join (args.expression, "denominator"))):
        numDir = os.path.join (args.expression, "numerator", dir, "fuzzy_values")
        denDir = os.path.join (args.expression, "denominator", dir, "fuzzy_values")
        nameDict = {"feature": [file[12:-4] for file in sorted (set (os.listdir (numDir)) & set (os.listdir (denDir))) if file.endswith (".tsv")],
                    "FS": list (rules.index)}
        numeratorFV, nameDict = readFV (numDir, nameDict = nameDict)
        denominatorFV, _ = readFV (denDir, nameDict = {"feature": nameDict["feature"], "sample": nameDict["sample"], "FS": list (rules.columns)})
        prod = np.einsum ("ijk, ijl -> ijkl", numeratorFV, denominatorFV)
        fuzzyFC = np.array ([[[prod[i, j, :, :][mask[FS]].sum () for FS in allSets] for j in range (prod.shape[1])] for i in range (prod.shape[0])])
        for idx in range (len (nameDict["sample"])):
            sample = nameDict["sample"][idx]
            origFV = pd.read_csv (os.path.join (args.paired, f"fuzzyValues_{sample}.tsv"), index_col = 0, sep = "\t")
            newFV = pd.DataFrame (fuzzyFC[:, idx, :], index = nameDict["feature"], columns = allSets)[origFV.columns].round (3)
            idxMask = (newFV.sum (axis = 1) != 0); origFV.loc[idxMask] = newFV.loc[idxMask]
            origFV.to_csv (os.path.join (args.output, f"fuzzyValues_{sample}.tsv"), sep = "\t")



if __name__ == "__main__":
    main ()


