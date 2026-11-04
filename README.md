# Methylation Consistency Across Individuals

A framework for identifying tissue-specific and conserved DNA methylation patterns in cattle and sheep.

## Overview

This repository contains the analysis code for the study:

*Methylation consistency across individuals reveals tissue-specific and conserved methylation patterns in cattle and sheep.*

The study introduces the **Methylation Consistency Index (MCI)** to quantify the consistency of DNA methylation states across individuals and investigate tissue-specific and cross-species methylation patterns.

Conventional DNA methylation analyses primarily focus on mean methylation levels and differences between biological groups. Methylation consistency across individuals represents a complementary dimension of methylome organization, providing additional insights into tissue-specific epigenetic patterns and their conservation across species.

This repository provides a reproducible analytical framework for DNA methylation site identification, adaptive filtering, methylation-level distribution analysis, MCI quantification, tissue-specific analysis, and cross-species conservation analysis.

## Repository Structure

The repository is organized into seven major analysis modules:

| Directory | Description |
|-----------|-------------|
| `01Ref_Pipeline_Prepare/` | Reference preparation and pipeline setup |
| `02Methy_identification/` | DNA methylation site identification |
| `03Adaptive_filtering/` | Adaptive filtering of reliable methylation sites |
| `04ML_Distribution/` | Analysis of methylation-level distributions |
| `05MCI_index/` | Methylation Consistency Index (MCI) calculation and analysis |
| `06MCI_Tissues/` | Tissue-specific methylation consistency analysis |
| `07Species/` | Cross-species comparison and conservation analysis between cattle and sheep |

## Analytical Workflow

The analytical framework consists of the following steps:

1. **Reference and pipeline preparation**  
   Prepare reference resources and inputs required for DNA methylation analysis.

2. **DNA methylation site identification**  
   Identify methylated CpG sites from whole-genome bisulfite sequencing (WGBS) data.

3. **Adaptive filtering**  
   Apply adaptive filtering strategies to identify reliable methylation sites.

4. **Methylation-level distribution analysis**  
   Characterize methylation-level distributions across samples, tissues, and species.

5. **MCI quantification**  
   Quantify the consistency of methylation states across individuals using the Methylation Consistency Index (MCI).

6. **Tissue-specific analysis**  
   Identify highly consistent methylation states associated with tissue identity.

7. **Cross-species conservation analysis**  
   Investigate conserved methylation patterns between cattle and sheep.

## Data

This study analyzed whole-genome bisulfite sequencing (WGBS) methylomes from cattle and sheep.

Sample information, data sources, and accession numbers are provided in the manuscript and its supplementary materials.


## Software and Requirements

The analyses were performed using computational tools and statistical packages described in the manuscript and the corresponding scripts.

Please refer to the individual analysis scripts for specific software requirements, parameters, and execution details.

## Usage

The analysis modules are organized according to the analytical workflow.

Users can navigate to the corresponding module directories to access the scripts for each analysis step.

For reproducibility, please refer to the software versions, parameters, and data resources specified in the manuscript and supplementary materials.

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions regarding the code or methodology, please contact:

**Shang-Qian Xie**  
Email: your_email@missouri.edu
