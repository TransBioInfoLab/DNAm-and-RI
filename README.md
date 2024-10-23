# **DNA Methylation Signature of a Lifestyle-based Resilience Index for Cognitive Health**

Wei Zhang, David Lukacsovich, Juan I. Young, Lissette Gomez, Michael A. Schmidt, Eden R. Martin, Brian W. Kunkle, Xi Chen, Deirdre M. O’Shea, James E. Galvin, Lily Wang

## Discription

Cognitive resilience (CR) contributes to the variability in risk for developing and progressing in Alzheimer’s disease (AD) among individuals. Beyond genetics, recent studies highlight the critical role of lifestyle factors in enhancing CR and delaying cognitive decline. DNA methylation (DNAm), an epigenetic mechanism influenced by both genetic and environmental factors, including CR-related lifestyle factors, offers a promising pathway for understanding the biology of CR. We studied DNAm changes associated with the Resilience Index (RI), a composite measure of lifestyle factors, using blood samples from the Healthy Brain Initiative (HBI) cohort. After corrections for multiple comparisons, our analysis identified 19 CpGs and 24 differentially methylated regions significantly associated with the RI, adjusting for covariates age, sex, *APOE* *ε4*, and immune cell composition. The RI-associated methylation changes are significantly enriched in pathways related to lipid metabolism, synaptic plasticity, and neuroinflammation, and highlight the connection between cardiovascular health and cognitive function. Furthermore, we developed a Methylation-based Resilience Score (MRS) that successfully predicted future cognitive decline, independent of age, sex, *APOE* *ε4*, years of education, baseline diagnosis, and baseline MMSE score, in an external dataset from the Alzheimer’s Disease Neuroimaging Initiative (ADNI). Our findings are particularly relevant for a better understanding of epigenetic architecture underlying cognitive resilience. Importantly, the significant association between baseline MRS and future cognitive decline demonstrated that DNAm could be a predictive marker for linking lifestyle factors to AD, laying the foundation for future studies on personalized AD prevention.

### 1. Preprocessing of DNA methylation data

The DNA methylation samples from HBI study were measured using the Illumina Infinium MethylationEPIC v2.0 Beadchips at the Center for Genomic Technology (CGT), John P. Hussman Institute for Human Genomics (HIHG). Quality control of both probes and samples were performed. 

| File and folder                                              | Description                        |
| ------------------------------------------------------------ | ---------------------------------- |
| [code/01_preprocessing/01_read_methyl_data.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/01_preprocessing/01_read_methyl_data.Rmd) | Reading DNA methylation data       |
| [code/01_preprocessing/02_preprocess.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/01_preprocessing/02_preprocess.Rmd) | Preprocessing DNA methylation data |

The preprocessing of DNAm samples from the ADNI study can be found in the repository ([Link](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia)) of our previous study[^1].

### 2.  **Statistical analyses to identify DNA methylation significantly associated with Resilience Index** 

For each CpG, we fitted a robust linear model (RLM) with DNAm M-values as the outcome, the RI as the primary independent variable, and relevant covariates including age, sex, diagnosis, APOE ε4 allele count, and the first two principal components of immune cell type proportions (which accounted for 90.7% of the variance in all estimated cell-type proportions). For region-based meta-analysis, we used the comb-p method[^2].

| File and folder                                              | Description                        |
| ------------------------------------------------------------ | ---------------------------------- |
| [code/02_resilience/resilience_model.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/02_resilience/resilience_model.Rmd) | Individual CpGs analysis using RLM |
| [code/03_combp/annotation.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/03_combp/annotation.Rmd) | Annotation of combp DMR results    |

To further validate the robustness of our results, we conducted sensitivity analyses by fitting the robust linear models separately to CN and MCI subjects. 

| File and folder                                              | Description                                |
| ------------------------------------------------------------ | ------------------------------------------ |
| [code/02_resilience/resilience_cn.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/02_resilience/resilience_cn.Rmd) | Individual CpGs analysis with CN subjects  |
| [code/02_resilience/resilience_mci.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/02_resilience/resilience_mci.Rmd) | Individual CpGs analysis with MCI subjects |

### 3. Pathway analysis

To identify biological pathways enriched with significant DNA methylation differences, we used the methylRRA function in the methylGSA R package[^3].

| File and folder                                              | Description      |
| ------------------------------------------------------------ | ---------------- |
| [code/04_pathway_analysis/pathway.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/04_pathway_analysis/pathway.Rmd) | Pathway analysis |

### 4. Integrative analyses with gene expression, genetic variants, and brain-to-blood correlations

To evaluate the effect of DNA methylation on the expression of nearby genes, we overlapped our dementia-associated CpGs, including both significant individual CpGs and those located within DMRs, with eQTm analysis results in Supplementary Tables 2 and 3 of Yao et al (2021)[^4].

To assess the correlation of dementia-associated CpGs and DMRs methylation levels in blood and brain samples, we used the London cohort, which consisted of 69 samples with matched PFC and blood samples.

For correlation and overlap with genetic susceptibility loci, we searched for mQTLs in the blood using the GoDMC database[^5] (http://mqtldb.godmc.org.uk/downloads).

| File and folder                                              | Description                                       |
| ------------------------------------------------------------ | ------------------------------------------------- |
| [code/05_check_overlap/check_overlap.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/05_check_overlap/check_overlap.Rmd) | Check overlapping CpGs with eQTm analysis results |
| [code/06_brain_blood_corr/brain_blood_corr.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/06_brain_blood_corr/brain_blood_corr.Rmd) | Brain-to-blood correlations using London cohort   |
| [code/07_mqtl/coloc_code_mQTL_gwasBellenguez.R](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/07_mqtl/coloc_code_mQTL_gwasBellenguez.R) | mQTLs analysis                                    |

### 5. Out-of-sample validation on the ADNI data using the Methylation-based Resilience Scores (MRS)

To assess disease progression, we performed out-of-sample validation on the ADNI data using the MRS by summing the methylation M-values for the CpGs with non-zero weights and adding the intercept obtained from the final model. We then conducted Cox proportional hazards regression analyses on the ADNI dataset to evaluate the association between MRS and disease progression.

| File and folder                                              | Description  |
| ------------------------------------------------------------ | ------------ |
| [code/08_mrs/mrs.Rmd](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/05_check_overlap/check_overlap.Rmd) | MRS analysis |

## For reproducible research

To perform the analysis, begin by installing the packages found in `00_utility/session_info.R` ([Link to the script](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/00_utility/session_info.R)). Then, load the auxiliary functions in the rest of the R scripts from the `00_utility` folder ([Link to the folder](https://github.com/TransBioInfoLab/DNAm-and-RI/blob/main/code/00_utility)). Follow the sequence provided in the Description to conduct the analysis.

## Acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

## Reference 

[^1]: Zhang, W. *et al.* Blood DNA Methylation Signature for Incident Dementia: Evidence from Longitudinal Cohorts. *Alzheimer's & Dementia*, In Review (2024).
[^2]: Pedersen, B.S., Schwartz, D.A., Yang, I.V. & Kechris, K.J. Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values. *Bioinformatics* **28**, 2986-8 (2012).
[^3]: Ren, X. & Kuan, P.F. methylGSA: a Bioconductor package and Shiny app for DNA methylation data length bias adjustment in gene set testing. *Bioinformatics* **35**, 1958-1959 (2019)
[^4]: Yao, C. *et al.* Epigenome-wide association study of whole blood gene expression in Framingham Heart Study participants provides molecular insight into the potential role of CHRNA5 in cigarette smoking-related lung diseases. *Clin Epigenetics* **13**, 60 (2021).
[^5]: Min, J.L. *et al.* Genomic and phenotypic insights from an atlas of genetic effects on DNA methylation. *Nat Genet* **53**, 1311-1321 (2021).
