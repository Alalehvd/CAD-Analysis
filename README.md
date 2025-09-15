# CAD-Analysis
Transcriptomic and functional enrichment analysis of **Canine Atopic Dermatitis (CAD)** using gene expression data from [NCBI GEO: GSE39278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39278).  
The goal is to identify differentially expressed genes (DEGs), perform functional enrichment (GO, KEGG), and explore biological pathways involved in CAD pathogenesis.

---

## 📊 Data
**Source:**  
- **GSE39278** — Gene expression profiles of lesional and non-lesional skin from dogs with CAD and healthy controls.
- Platform: **Affymetrix Canine Genome Array**  
- Species: *Canis lupus familiaris* (dog)  
- Reference genome: **CanFam3.1**

> **Note:** Raw CEL files are not included in this repository. Please download directly from GEO using the link above.

---

## 🧬 Methodology

1. **Data Acquisition & Preprocessing**  
   - Download raw CEL files from GEO  
   - Preprocessing using `affy` (RMA normalization, background correction)  
   - Probe-to-gene annotation using `org.Cf.eg.db`

2. **Differential Expression Analysis**  
   - Conducted with **limma**  
   - Pairwise contrasts:  
     - Lesional CAD vs Healthy skin  
     - Non-lesional CAD vs Healthy skin  
     - Lesional CAD vs Non-lesional CAD  
   - Significance thresholds: `adj.P.Val < 0.05`, `|log2FC| > 1`
   - Significance thresholds: `adj.P.Val < 0.1`, `|log2FC| > 0.58`

3. **Functional Enrichment Analysis**  
   - **GO** (Biological Process, Cellular Component, Molecular Function) — *dog-specific (`OrgDb = org.Cf.eg.db`)*  
   - **KEGG** — *dog-specific (`organism = "cfa"`)*
   - Analysis performed for all DEGs and significant DEGs separately

4. **Visualization**  
   - Volcano plots  
   - Heatmaps of top DEGs  
   - KEGG/GO barplots & dotplots (`enrichplot`)  
   - PCA plots of expression profiles

---

## 🛠️ Tools & Dependencies

**R Packages**
- `affy`
- `limma`
- `clusterProfiler`
- `org.Cf.eg.db`
- `AnnotationDbi`
- `enrichplot`
- `ggplot2`
- `pheatmap`
- `EnhancedVolcano`

**Other**
- R ≥ 4.0.0  
- Platform annotation data: Affymetrix Canine Genome Array

---

## 🤝 Contributing
Contributions are welcome!  
- Please open an **issue** for bug reports or questions  
- Submit a **pull request** for new analyses or improvements

---

## 📧 Contact
**Dr. Alaleh Vazifehdoost**  
Veterinary Medicine & AI Research  
✉️ alalehvd@gmail.com  
🌐 https://www.linkedin.com/in/alaleh-vazifehdoost/
