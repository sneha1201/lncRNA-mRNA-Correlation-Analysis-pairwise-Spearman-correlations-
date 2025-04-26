# lncRNA-mRNA-Correlation-Analysis-pairwise-Spearman-correlations-

# 🧬 lncRNA–mRNA Correlation Analysis (Parallelized)

This script computes **pairwise Spearman correlations** between differentially expressed lncRNAs and mRNAs using **parallel processing**. It identifies significantly correlated lncRNA–mRNA pairs—both **positively and negatively correlated**—based on user-defined statistical thresholds.

---

## 🚀 Features

- Fast correlation computation using **multiprocessing**
  
- Filters lncRNA–mRNA pairs based on:
  
  - Correlation magnitude (|r| ≥ 0.9 → strong positive or negative)
    
  - Statistical significance (p-value < 0.05)
- Automatically excludes constant-expression genes (zero variance)
  
- Outputs both full and filtered correlation results

---

## 📁 Input Files

| File                     | Description                                |
|--------------------------|--------------------------------------------|
| `lncRNA_file.csv`        | Raw count matrix of lncRNAs (genes × samples) |
| `mRNA_file.csv`          | Raw count matrix of mRNAs (genes × samples)  |
| `DEG_lncRNA_file.txt`    | List of differentially expressed lncRNA gene IDs (one per line) |
| `DEG_mRNA_file.txt`      | List of differentially expressed mRNA gene IDs (one per line)   |

> ⚠️ Count files should have **gene IDs as rows** and **sample names as columns**.  
> The script transposes them internally to have samples as rows and genes as columns.

---

## 📦 Requirements

- Python 3.6+
- `pandas`
- `numpy`
- `scipy`

Install dependencies with:

```bash
pip install pandas numpy scipy

🖥️ Usage

python parallel_correlation.py <lncRNA_file.csv> <mRNA_file.csv> <DEG_lncRNA_file.txt> <DEG_mRNA_file.txt> <num_cores>

🔹 Example

python parallel_correlation.py lnc_counts.csv mrna_counts.csv lnc_deg.txt mrna_deg.txt 8

This runs the script using 8 CPU cores.
🧪 Output Files

    lncRNA_mRNA_correlation_all.csv
    → Contains all lncRNA–mRNA correlation results

    lncRNA_mRNA_correlation_filtered.csv
    → Contains only the filtered significant pairs where:

        |Spearman correlation| ≥ 0.9 (i.e., r ≥ 0.9 or r ≤ -0.9)

        p-value < 0.05

📌 Notes

    The script only compares lncRNA–mRNA pairs where both genes are differentially expressed.

    Genes with zero variance across all samples are excluded.

    Designed to handle large datasets efficiently using parallelization.
