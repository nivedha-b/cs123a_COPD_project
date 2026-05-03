# COPD GWAS Analysis

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
res_df = pd.DataFrame(results, columns=["Feature", "Beta", "P_Value"])

# Sort by significance
results_df = res_df.sort_values("P_Value")

print("Top results:")
print(results_df.head())



# # -----------------------------
# # OLD - Basic Quality Control
# # -----------------------------

# # Remove rows with missing phenotype information
# required_columns = ["COPD_status", "age", "sex", "BMI", "smoking_status", "PC1", "PC2", "PC3"]
# data = data.dropna(subset=required_columns)

# # Identify SNP columns
# snp_columns = []

# for column_name in genotype_data.columns:
#     if column_name != "sample_id":
#         snp_columns.append(column_name)

# clean_snp_columns = []

# for snp in snp_columns:
#     call_rate = data[snp].notna().mean()

#     allele_frequency = data[snp].mean() / 2
#     minor_allele_frequency = min(allele_frequency, 1 - allele_frequency)

#     if call_rate >= 0.95 and minor_allele_frequency >= 0.01:
#         clean_snp_columns.append(snp)

# print("Number of SNPs before QC:")
# print(len(snp_columns))

# print("Number of SNPs after QC:")
# print(len(clean_snp_columns))

# qc_report = pd.DataFrame({
#     "total_samples_after_cleaning": [data.shape[0]],
#     "total_snps_before_qc": [len(snp_columns)],
#     "total_snps_after_qc": [len(clean_snp_columns)]
# })

# qc_report.to_csv("qc_report.csv", index=False)

# # -----------------------------
# # OLD - Run GWAS Logistic Regression
# # -----------------------------

# results = []

# covariates = ["age", "sex", "BMI", "smoking_status", "PC1", "PC2", "PC3"]

# for snp in clean_snp_columns:
#     analysis_columns = ["COPD_status", snp] + covariates
#     snp_data = data[analysis_columns].dropna()

#     y = snp_data["COPD_status"]
#     x = snp_data[[snp] + covariates]

#     x = sm.add_constant(x)

#     try:
#         model = sm.Logit(y, x)
#         result = model.fit(disp=False)

#         beta = result.params[snp]
#         p_value = result.pvalues[snp]
#         odds_ratio = np.exp(beta)

#         results.append([snp, beta, odds_ratio, p_value])

#     except:
#         print("Could not run model for:")
#         print(snp)

# gwas_results = pd.DataFrame(
#     results,
#     columns=["SNP", "Beta", "Odds_Ratio", "P_Value"]
# )

# gwas_results = gwas_results.sort_values("P_Value")

# gwas_results.to_csv("gwas_results.csv", index=False)

# print("Top GWAS results:")
# print(gwas_results.head(10))

# # -----------------------------
# # OLD - Manhattan Plot
# # -----------------------------

# gwas_results["minus_log10_p"] = -np.log10(gwas_results["P_Value"])
# gwas_results["SNP_index"] = range(1, len(gwas_results) + 1)

# plt.figure(figsize=(12, 6))
# plt.scatter(gwas_results["SNP_index"], gwas_results["minus_log10_p"], s=10)

# plt.axhline(-np.log10(5e-8), linestyle="--")
# plt.xlabel("SNP Index")
# plt.ylabel("-log10(P-value)")
# plt.title("Manhattan Plot for COPD GWAS")
# plt.tight_layout()
# plt.savefig("manhattan_plot.png")
# plt.close()

# # -----------------------------
# # OLD - QQ Plot
# # -----------------------------

# observed_p_values = gwas_results["P_Value"].dropna().sort_values()
# expected_p_values = np.arange(1, len(observed_p_values) + 1) / (len(observed_p_values) + 1)

# observed_values = -np.log10(observed_p_values)
# expected_values = -np.log10(expected_p_values)

# plt.figure(figsize=(6, 6))
# plt.scatter(expected_values, observed_values, s=10)

# maximum_value = max(expected_values.max(), observed_values.max())
# plt.plot([0, maximum_value], [0, maximum_value])

# plt.xlabel("Expected -log10(P-value)")
# plt.ylabel("Observed -log10(P-value)")
# plt.title("QQ Plot for COPD GWAS")
# plt.tight_layout()
# plt.savefig("qq_plot.png")
# plt.close()

# # -----------------------------
# # OLD - Top SNP Annotation
# # -----------------------------

# known_copd_genes = {
#     "FAM13A": "Rho GTPase signaling",
#     "HHIP": "Hedgehog signaling",
#     "IREB2": "Iron regulation",
#     "AGER": "Inflammation and lung disease",
#     "MMP1": "Matrix remodeling and emphysema",
#     "MMP12": "Matrix remodeling and emphysema"
# }

# top_snps = gwas_results.head(20).copy()

# # Placeholder annotation column
# # In a real project, this would come from VEP, FUMA, or GTEx.
# top_snps["Possible_Gene"] = "Needs annotation"
# top_snps["Possible_Pathway"] = "Needs pathway lookup"

# top_snps.to_csv("top_snps_for_annotation.csv", index=False)

# print("Pipeline complete.")
# print("Files created:")
# print("qc_report.csv")
# print("gwas_results.csv")
# print("manhattan_plot.png")
# print("qq_plot.png")
# print("top_snps_for_annotation.csv")
