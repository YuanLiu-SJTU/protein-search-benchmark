# Fold classification

## Overview

**Fold_classification** is a benchmark for evaluating protein fold-level similarity search methods on the **CATH-S20** dataset. This repository contains code and scripts used to perform **all-versus-all protein structure search** and compute the performance of different search paradigms.

The benchmark focuses on comparing **structure alignment-based**, **structure representation-based**, and **sequence-based** search approaches, highlighting how search paradigm affects retrieval performance.

---

## Dataset

- **CATH-S20**: a representative subset of the CATH database, containing proteins classified at the fold level.  
- The dataset is used to perform **all-versus-all searches**, ensuring a comprehensive evaluation of remote homology detection capabilities.

---

## Evaluation Metric

The primary metric used is:

**Sensitivity up to the first false positive**

- For each query protein, the search results are ranked according to the method-specific score.  
- Sensitivity is computed as the fraction of true positive hits retrieved **before the first false positive** occurs.  
- This metric captures the ability of a method to prioritize structurally similar proteins at the top of the ranking list, which is critical for fold-level search tasks.

---

## Benchmark Procedure

1. **All-versus-all search**:  
   Each protein in CATH-S20 is used as a query against the rest of the dataset.  

2. **Methods evaluated**:  
   - **Structure alignment-based**: GTalign, TM-align, Dali, Foldseek  
   - **Representation-based**: GraSR, TMvec, PLMSearch, DHR, FoldExplorer  
   - **Sequence alignment-based**: BLAST, jackhmmer, Diamond, MMseqs

3. **Performance calculation**:  
   For each query, compute the **Sensitivity up to the first false positive**. Aggregate results across all queries to compare methods.

4. **Analysis**:  
   - Examine sensitivity trends for each method.  
   - Identify which search paradigms are more effcient on fold classification tasks.

## Usage

Run the sensitivity computation for a single search result TSV file:

```bash
python compute_sensitivity.py --search_result search_result_file.tsv
