# CAD-Analysis
Machine learning–driven exploration of genetic and omics data to identify biomarkers and develop predictive tools for **Canine Atopic Dermatitis (CAD)**.

---

## 📊 Data
The primary dataset originates from:

- **Forman et al.**, Dryad PLINK dataset on CAD  
- Breed-specific case–control SNP data for **20+ dog breeds**  
- Reference genome: **CanFam3.1**  

> **Note:** Raw genotype data is **not included** in this repository due to licensing restrictions. Please refer to the original publication’s data access guidelines.

---

## 🧬 Methodology

1. **Data Cleaning & QC**  
   - PLINK-based filtering (MAF, missingness, Hardy–Weinberg equilibrium)

2. **Statistical Analysis**  
   - Case–control association testing per breed  
   - Multiple testing correction (Bonferroni/FDR)

3. **Annotation**  
   - SNP-to-gene mapping (Ensembl, UCSC resources)

4. **Pathway Analysis**  
   - Functional enrichment: **GO** and **KEGG** (dog-specific, `org.Cf.eg.db`, `cfa`)

5. **Machine Learning**  
   - Feature selection from SNP data  
   - Model training (e.g., **XGBoost**, Logistic Regression)  
   - Performance evaluation (AUC, accuracy, feature importance)

---

## 🛠️ Tools & Dependencies

**Core Tools**
- PLINK **1.9 / 2.0**
- **Python 3.10+**: `pandas`, `scikit-learn`, `xgboost`, `matplotlib`, `seaborn`
- **R 4.0+**: `tidyverse`, `clusterProfiler`, `WGCNA`, `qqman`

**Reference Data**
- CanFam3.1 reference genome

---

## 📜 Citation
If you use this repository, please cite:

- Forman OP, *et al.* (Year). *Title of the CAD study*. Journal Name. DOI.

---

## 🤝 Contributing
Contributions are welcome!  
- Open an issue for questions or suggestions  
- Submit a pull request for fixes or enhancements

---

## 📧 Contact
**Dr. Alaleh V.D.**  
Veterinary Medicine & AI Research  
✉️ [your.email@example.com]  
🌐 [LinkedIn Profile URL]