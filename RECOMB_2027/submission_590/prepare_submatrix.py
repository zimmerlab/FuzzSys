import os
import argparse
import pandas as pd

# python prepare_submatrix.py --data dataDirectory --category categoryNameList --output outputFile


parser = argparse.ArgumentParser ()
parser.add_argument ("--data", type = str, required = True, help = "Directory with concatenated raw value matrix and metadata")
parser.add_argument ("--category", type = str, required = True, help = "List of comma-separated category names")
parser.add_argument ("--output", type = str, required = True, help = "Output file path")
args = parser.parse_args ()

matrix = pd.read_csv (os.path.join (args.data, "raw_expression.tsv"), index_col = 0, sep = "\t")
metadata = pd.read_csv (os.path.join (args.data, "metadata.tsv"), index_col = None, sep = "\t")

category = [x.strip () for x in args.category.split (",")]
submatrix = matrix[metadata.loc[metadata["category"].isin (category), "patient"]]

os.makedirs (os.path.dirname (args.output), exist_ok = True)
submatrix.to_csv (args.output, sep = "\t")


