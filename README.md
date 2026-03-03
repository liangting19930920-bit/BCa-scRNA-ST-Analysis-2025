# Scripts for data analysis and plotting in "Integrin α5β1-Mediated Multicellular Crosstalk in the Tumor Microenvironment Drives Bladder Cancer Progression and Reveals Targetable Vulnerabilities"

Scripts for statistics, bioinformatics analysis, and plotting figures in "Integrin α5β1-Mediated Multicellular Crosstalk in the Tumor Microenvironment Drives Bladder Cancer Progression and Reveals Targetable Vulnerabilities", published in iMetaOmics (2026).

## Repository Structure

All analytical R scripts used in this study are provided below:

* `01_ScBLCA_integration_pipeline.R`: scRNA-seq integration and global annotation.
* `02_Epithelial_Downstream_Analysis.R`: Epithelial sub-clustering, stemness scoring, and pseudotime trajectory.
* `03_Fibroblast_Analysis.R`: CAF sub-clustering and functional enrichment.
* `04_Myeloid_Downstream_Analysis.R`: Myeloid cells trajectory and SCENIC TF networks.
* `05_Endothelial_Downstream_Analysis.R`: Endothelial cells metabolic reprogramming and PROGENy.
* `06_Spatial_Transcriptomics_Pipeline.R`: ST data integration, scoring, and CellChat spatial networks.
* `07_Spatial_CellTrek_Mapping.R`: Single-cell to spatial mapping and distance (kdist) analysis.
* `08_Spatial_MIA_Integration.R`: Multimodal Intersection Analysis (MIA) for functional regions.
* `09_Virtual_Knockout_ITGA5_ITGB1.R`: scGRN construction and virtual knockout simulation.

*(Note: The processed single-cell and spatial transcriptomics data generated in this study have been deposited in the OMIX database of the China National Center for Bioinformation (CNCB) under the accession number OMIX015404, which is publicly accessible at [https://ngdc.cncb.ac.cn/omix/preview/6lcfgVxv](https://ngdc.cncb.ac.cn/omix/preview/6lcfgVxv). Output results and figures are provided in the Supplementary Materials of the paper.)

## Environment

All analyses were performed in `R` (version >= 4.1.0). Key packages include `Seurat` (v4.4.0), `harmony`, `Monocle2`, `CellChat`, `CellTrek`, and `scTenifoldKnk`.

## Citation

If you use these scripts or data in your research, please cite the paper below:

> Ting Liang#, Wuwu Xu#, Hu Fang#, Lu Fu# ,......., Guangzhi Li*, Song Wu*, Tao Tao*. Integrin α5β1-Mediated Multicellular Crosstalk in the Tumor Microenvironment Drives Bladder Cancer Progression and Reveals Targetable Vulnerabilities. 2026. iMetaOmics. 
