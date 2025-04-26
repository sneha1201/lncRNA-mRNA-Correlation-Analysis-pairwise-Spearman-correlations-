import pandas as pd
import numpy as np
import sys
from scipy.stats import spearmanr
from multiprocessing import Pool, cpu_count

# ------------- Arguments from user ----------------
if len(sys.argv) < 5:
    print("Usage: python parallel_correlation.py <lncRNA_file.csv> <mRNA_file.csv> <DEG_lncRNA_file.txt> <DEG_mRNA_file.txt> <num_cores>")
    sys.exit(1)

lncRNA_file = sys.argv[1]  # lncRNA raw counts file
mRNA_file = sys.argv[2]    # mRNA raw counts file
lncRNA_deg_file = sys.argv[3]  # txt file with lncRNA DEG IDs
mRNA_deg_file = sys.argv[4]    # txt file with mRNA DEG IDs
num_cores = int(sys.argv[5])  # Number of cores for parallel processing

# ------------- Load Data ----------------
lnc_df = pd.read_csv(lncRNA_file, index_col=0)
mRNA_df = pd.read_csv(mRNA_file, index_col=0)

# ------------- Load DEG IDs ----------------
with open(lncRNA_deg_file, 'r') as f:
    lnc_deg_ids = [line.strip() for line in f.readlines()]

with open(mRNA_deg_file, 'r') as f:
    mRNA_deg_ids = [line.strip() for line in f.readlines()]

# ------------- Preprocess ----------------
# Transpose so rows are samples, columns are genes
lnc_df = lnc_df.T
mRNA_df = mRNA_df.T

# Keep only common samples
common_samples = lnc_df.index.intersection(mRNA_df.index)
lnc_df = lnc_df.loc[common_samples]
mRNA_df = mRNA_df.loc[common_samples]

# Filter data based on DEG IDs
lnc_df = lnc_df.loc[:, lnc_df.columns.isin(lnc_deg_ids)]
mRNA_df = mRNA_df.loc[:, mRNA_df.columns.isin(mRNA_deg_ids)]

# Remove constant genes (zero variance)
lnc_df = lnc_df.loc[:, lnc_df.std() != 0]
mRNA_df = mRNA_df.loc[:, mRNA_df.std() != 0]

# ------------- Define Correlation Function ----------------
def compute_correlation(lnc_mRNA_pair):
    lnc_col, mRNA_col = lnc_mRNA_pair
    corr, pval = spearmanr(lnc_df[lnc_col], mRNA_df[mRNA_col])
    return (lnc_col, mRNA_col, corr, pval)

# ------------- Create All Pairs ----------------
all_pairs = [(lnc, mrna) for lnc in lnc_df.columns for mrna in mRNA_df.columns]

# ------------- Parallel Computation ----------------
print(f"Using {num_cores} cores...")
with Pool(processes=num_cores) as pool:
    results = pool.map(compute_correlation, all_pairs)

# ------------- Format and Filter ----------------
results_df = pd.DataFrame(results, columns=["lncRNA_ID", "mRNA_ID", "correlation", "p_value"])
filtered = results_df[(results_df["correlation"].abs() >= 0.9) & (results_df["p_value"] < 0.05)]

# ------------- Save Output ----------------
filtered.to_csv("lncRNA_mRNA_correlation_filtered.csv", index=False)
results_df.to_csv("lncRNA_mRNA_correlation_all.csv", index=False)
print("Filtered correlations saved to 'lncRNA_mRNA_correlation_filtered.csv'")

