# cs123a_COPD_project
This project focuses on conducting a genomic as well as epigenomic analysis of the susceptibility of Chronic Obstructive Pulmonary Disease (COPD) based on publicly available datasets.

This project mainly implements an Epigenome - Wide Association Study (EWAS) pipeline via Gmethylation datasets in order to find and observe methylation patterns from smoker associated samples associated with COPD.

The second portion of this project will cover the COPD genetic literature as well as GWAS statistical summaries to further verify the biological perspective of these observations.

The dataset (csv) used for EWAS can be accessed and downloaded at the following link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326391
The dataset csv should be added to a folder titled "data".

To run the pipeline and obtain the outputs, simply run the file "copd_ewas_analysis.py" on an IDE, such as VSCode.

## EWAS Pipeline
This pipeline contains:

- Preprocessing and loading of Data
- Preparation of the Methylation matrix 
- Association Testing
- Manhattan Plot 
- QQ Plot
- Export of Top CpG Sites


The pipeline outputs generate:
- 'manhattan_plot.png'
- 'qq_plot.png'
- 'top_cpg_sites.csv'
- 'significant_cpg_sites.csv'

---

## Biological Interpretation

The EWAS analysis identified multiple CpG methylation sites linked to COPD status in smoker associated samples.

These changes in methylation may reflect biological mechanisms involved in:

- Lung tissue damage associated with smoking exposure
- Chronic inflammation
- Immune response dysregulation
- Oxidative stress
- Remodeling of airway

Alongside the primary focuses on methylation analysis, COPD literature as well as published
GWAS summary statistics were also observed to fortify biological interpretation of observations.

Known COPD-associated genes discussed in prior GWAS studies include:

- HHIP
- FAM13A
- IREB2
- AGER
- CHRNA3/5

---

## Limitations

Individual-level GWAS genotype datasets such as COPDGene
and ECLIPSE are controlled-access through dbGaP which were not accessible
within the project timeline.

Due to this limitation, this project has instead implemented an EWAS-based methylation analysis based on publicly available GEO datasets in addition to incorporating published GWAS summary statistics for fortifying observations.

