# FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study

**Short Title:** FAP in malignant transformation in NF1  
**Keywords:** Neurofibromatosis Type I, FAP, Fibroblast Activation Protein, Molecular Imaging, PET, MPNST

---

## Overview

This repository contains the analysis scripts and supporting code used in the study titled:

**"FAP Expression as a Marker of Malignancy Enabling In-Vivo Imaging in NF1-Associated Peripheral Nerve Tumors: A Multimodal and Translational Study"**

This work investigates Fibroblast Activation Protein (FAP) as a diagnostic biomarker and theranostic target distinguishing malignant peripheral nerve sheath tumors (MPNSTs) from benign neurofibromas in Neurofibromatosis Type I (NF1) patients. Using bulk, spatial, and single-cell transcriptomic datasets combined with immunohistochemistry and clinical PET/CT imaging data, we demonstrate the elevated expression of FAP in MPNSTs, supporting its potential as a clinical imaging target.

The manuscsript is currently under review.

---
## Repository Structure and Scripts

This repository includes five R scripts that replicate the key transcriptomic analyses presented in the manuscript:

| Script Name                             | Description                                                                                      |
|---------------------------------------|------------------------------------------------------------------------------------------------|
| `01_Bulk_RNA_EBioMed2023.R`           | Analysis of bulk RNA expression data from the EBioMedicine 2023 dataset (GSE241224), comparing FAP expression between neurofibromas and MPNSTs. |
| `02_Bulk_RNA_MolOncol2015.R`          | Analysis of bulk RNA expression data from Molecular Oncology 2015 dataset (GSE66743) for validation of FAP expression patterns. |
| `03_Spatial_Transcritpomics_MPNST_NF_NeuroOncology_2025.R` | Processing and analysis of spatial transcriptomics data (10x Genomics Visium) from NF1-associated peripheral nerve tumors, focusing on spatial FAP expression and co-localization with tumor and stromal markers. |
| `04_scRNA_MPNST_SciAdv2022.R`         | Single-cell RNA-sequencing analysis of NF1-associated MPNST samples (GSE179033), including clustering, marker identification, and differentiation trajectory analysis. |
| `05_TCGA_FAP_Sarcoma_PanCancerAtlas.R` | Analysis of FAP gene expression across sarcoma subtypes in the TCGA PanCancer Atlas dataset. |



---

## References

For more information on the datasets that were used, please refer to the following manuscripts:
- Comprehensive and Integrated Genomic Characterization of Adult Soft Tissue Sarcomas. Cell. 2017 Nov 2;171(4):950-965.e28. doi: 10.1016/j.cell.2017.10.014
Suppiah S, Mansouri S, Mamatjan Y, et al. Multiplatform molecular profiling uncovers two subgroups of malignant peripheral nerve sheath tumors with distinct therapeutic vulnerabilities. Nat Commun. 2023;14(1):2696. Published 2023 May 10. doi:10.1038/s41467-023-38432-6
- Bremer J, Franco P, Menstell JA, Tey S, Zajt KK, Popzhelyazkova K, Nolte K, Schlegel J, Pedro MT, Osterloh A, Delev D, Hohenhaus M, Scholz C, Schnell O, Beck J, Weis J, Heiland DH. Spatially resolved transcriptomics of benign and malignant peripheral nerve sheath tumors. Neuro Oncol. 2025 Jan 23:noaf016. doi: 10.1093/neuonc/noaf016
- Høland M, Berg KCG, Eilertsen IA, Bjerkehagen B, Kolberg M, Boye K, Lingjærde OC, Guren TK, Mandahl N, van den Berg E, Palmerini E, Smeland S, Picci P, Mertens F, Sveen A, Lothe RA. Transcriptomic subtyping of malignant peripheral nerve sheath tumours highlights immune signatures, genomic profiles, patient survival and therapeutic targets. EBioMedicine. 2023 Nov;97:104829. doi: 10.1016/j.ebiom.2023.104829
- Kolberg M, Høland M, Lind GE, Ågesen TH, Skotheim RI, Hall KS, Mandahl N, Smeland S, Mertens F, Davidson B, Lothe RA. Protein expression of BIRC5, TK1, and TOP2A in malignant peripheral nerve sheath tumours--A prognostic test after surgical resection. Mol Oncol. 2015 Jun;9(6):1129-39. doi: 10.1016/j.molonc.2015.02.005
- Wu LMN, Zhang F, Rao R, Adam M, Pollard K, Szabo S, Liu X, Belcher KA, Luo Z, Ogurek S, Reilly C, Zhou X, Zhang L, Rubin J, Chang LS, Xin M, Yu J, Suva M, Pratilas CA, Potter S, Lu QR. Single-cell multiomics identifies clinically relevant mesenchymal stem-like cells and key regulators for MPNST malignancy. Sci Adv. 2022 Nov 4;8(44):eabo5442. doi: 10.1126/sciadv.abo5442

---

## Contact

For questions or collaborations, please contact:  
**Dr. med. Nic G. Reitsam**  
**Pathology, Faculty of Medicine, University of Augsburg, Augsburg, Germany**  
**Email:** nic.reitsam@uk-augsburg.de
