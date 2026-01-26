# Functional Evaluation of Protein Search Results

This folder contains scripts for evaluating the **functional consistency** of protein search results, complementing structural or sequence similarity analyses. Functional consistency is assessed using Gene Ontology Molecular Function (GO-MF) annotations.

---

## Background

Structural or sequence similarity does not necessarily imply functional relatedness. To assess the biological relevance of protein search results, we evaluated functional consistency on a **SwissProt dataset filtered to <20% sequence identity**, reducing redundancy and emphasizing remote functional relationships.

---

## Scripts

### `compute_function_BMA.py`

This script computes **GO semantic similarity scores** between query proteins and their retrieved Top-k hits using the **Best Match Average (BMA)** metric.  

Features:

- Accepts a single **search result TSV file** with query-target pairs.  
- Performs **random padding** to ensure each query has top-k targets (default `topk=10`).  
- Computes **BMA GO-MF similarity** using `go-basic.obo` and `GO.json` annotation files.  
- Generates a **boxplot of similarity scores** for all query-target pairs.  

#### Usage

```bash
python compute_function_BMA.py --search_result search_result_file.tsv --output_dir figures
```

## Requirements

- **Python version:** 3.8 or higher  
- **Python packages:**  
  - `numpy`  
  - `pandas`  
  - `matplotlib`  
  - `seaborn`  
  - `tqdm`  
  - `goatools`  

- **Data files:**  
  - `GO.json` – Gene Ontology annotations for all proteins  
  - `go-basic.obo` – Gene Ontology DAG file
