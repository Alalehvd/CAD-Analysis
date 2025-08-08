# CAD-Analysis
 Machine learning–driven exploration of genetic and omics data to identify biomarkers and develop predictive tools for Canine Atopic Dermatitis.


---

## 📊 Data
The primary dataset originates from:
- **Forman et al., Dryad PLINK dataset** on CAD
- Breed-specific case–control SNP data for >20 dog breeds
- Reference genome: **CanFam3.1**

> **Note:** Raw genotype data is not included in this repository due to licensing. Please follow the original publication’s data access guidelines.

---

## 🧬 Methodology
1. **Data Cleaning & QC**
   - PLINK-based filtering (MAF, missingness, Hardy–Weinberg equilibrium)
2. **Statistical Analysis**
   - Case–control association testing per breed
   - Bonferroni/FDR correction
3. **Annotation**
   - Mapping SNPs to genes (Ensembl, UCSC)
4. **Pathway Analysis**
   - Functional enrichment (KEGG, GO)
5. **Machine Learning**
   - Feature selection from SNP data
   - Model training (e.g., XGBoost, Logistic Regression)
   - Evaluation (AUC, accuracy, feature importance)

---

## 🛠️ Tools & Dependencies
- **PLINK 1.9/2.0**
- **Python 3.10+** (pandas, scikit-learn, xgboost, matplotlib, seaborn)
- **R 4.0+** (tidyverse, clusterProfiler, WGCNA, qqman)

---

## 📜 Citation
If you use this code or methodology, please cite:
- Forman OP, et al. (Year). *Title of the CAD study*. Journal Name. DOI.

---

## 🤝 Contributing
Contributions are welcome! Please open an issue or submit a pull request.

---

## 📧 Contact
Dr. Alaleh V.D.  
Veterinary Medicine & AI Research  
✉️ [your.email@example.com]  
🌐 [LinkedIn Profile URL]
