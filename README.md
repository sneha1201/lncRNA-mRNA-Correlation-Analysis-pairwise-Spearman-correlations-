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
