# COPD EWAS Analysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

# -----------------------------
# 1. Load Data
# -----------------------------


# Load GEO methylation site data
raw_data = pd.read_csv("data/GSE326391_sites.csv")

print("Original dataset shape:")
print(raw_data.shape)

print("Few initial columns:")
print(raw_data.columns[:20])

print("Few initial rows:")
print(raw_data.head())


# -----------------------------
# 2. Prepare Methylation Matrix
# -----------------------------

metadata_col = ["ID", "Chromosome", "Start", "End", "Strand"]
sample_col = [col for col in raw_data.columns if col not in metadata_col]

print("Number of samples:")
print(len(sample_col))

print("Number of methylation sites:")
print(raw_data.shape[0])

# Create simple phenotype labels (temporary)
phenotype_data = pd.DataFrame({
    "sample_id": sample_col,
    "COPD_status": [1 if i < len(sample_col)//2 else 0 for i in range(len(sample_col))]
})

print("Phenotype data:")
print(phenotype_data.head())


# -----------------------------
# 3. Transpose Data
# -----------------------------

# Use CpG IDs as feature names
methylation_matrix = raw_data[sample_col].copy()
methylation_matrix.index = raw_data["ID"]

# Transpose so rows = samples, columns = CpG sites
methylation_matrix = methylation_matrix.T

# Add sample IDs as column
methylation_matrix["sample_id"] = methylation_matrix.index

# Merge with phenotype
analysis_data = phenotype_data.merge(methylation_matrix, on="sample_id")

print("Analysis dataset shape:")
print(analysis_data.shape)

print("First few rows of analysis data:")
print(analysis_data.head())

# -----------------------------
# 4. Association Testing
# -----------------------------

results = []

# IMPORTANT: limit features for speed
feature_columns = analysis_data.columns[2:1000]  # initial 1000 CpG sites

for feature in feature_columns:
    try:
        y = analysis_data["COPD_status"]
        X = analysis_data[[feature]]

        X = sm.add_constant(X)

        model = sm.Logit(y, X)
        result = model.fit(disp=False)

        beta = result.params[feature]
        p_value = result.pvalues[feature]

        results.append([feature, beta, p_value])

    except:
        continue

# Convert to DataFrame
results_df = pd.DataFrame(results, columns=["Feature", "Beta", "P_Value"])

# Sort by significance
results_df = results_df.sort_values("P_Value")

print("Top results:")
print(results_df.head())


# -----------------------------
# 5. Manhattan Plot
# -----------------------------

results_df.to_csv("methylation_results.csv", index=False)

results_df["minus_log10_p"] = -np.log10(results_df["P_Value"])
results_df["CpG_index"] = range(1, len(results_df) + 1)

plt.figure(figsize=(12, 6))
plt.scatter(results_df["CpG_index"], results_df["minus_log10_p"], s=10)

plt.axhline(-np.log10(0.05), linestyle="--")
plt.xlabel("CpG Site Index")
plt.ylabel("-log10(P-value)")

plt.title("Manhattan Plot for COPD Methylation Association")
plt.tight_layout()
plt.savefig("manhattan_plot.png")
plt.close()

print("Manhattan plot has been saved as manhattan_plot.png")


# -----------------------------
# 6. QQ Plot
# -----------------------------

# Remove the missing p-values then sort 
observed_pv = results_df["P_Value"].dropna().sort_values()

# The expected p-values
expected_pv = np.arange(1, len(observed_pv) + 1) / (len(observed_pv) + 1)

# Convert into -log10 scale
obs_log = -np.log10(observed_pv)
exp_log = -np.log10(expected_pv)

# QQ plot made
plt.figure(figsize=(6, 6))
plt.scatter(exp_log, obs_log, s=10)

# diagonal line
max_val = max(exp_log.max(), obs_log.max())
plt.plot([0, max_val], [0, max_val])

plt.xlabel("Expected -log10(P-value)")
plt.ylabel("Observed -log10(P-value)")
plt.title("COPD Methylation Analysis - QQ Plot")

plt.tight_layout()
plt.savefig("qq_plot.png")
plt.close()

print("QQ plot has been saved as qq_plot.png")


# -----------------------------
# 7. Export of top CpG Sites
# -----------------------------

# Choosing the top 20 CpG sites regarding smallest p-values
top_sites = results_df.head(20)
top_sites.to_csv("top_cpg_sites.csv", index=False)

# Will save CpG sites that contain suggestive significance
sig_sites = results_df[results_df["P_Value"] < 0.05]
sig_sites.to_csv("significant_cpg_sites.csv", index=False)

print("Top CpG sites has been saved as top_cpg_sites.csv")
print("CpG sites that contain suggestive significance has been saved as significant_cpg_sites.csv")
print("Amount of significant CpG sites:", len(sig_sites))