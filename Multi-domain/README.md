### lDDT Evaluation Script

This Folder includes a modified version of the `lddt` protein distance score script originally developed by DeepMind Technologies and adapted by Martin Steinegger for alignment comparison.

- **Original source**: [Foldseek-analysis GitHub](https://github.com/steineggerlab/foldseek-analysis)
- **License**: Apache License 2.0  
  (See [LICENSE](https://www.apache.org/licenses/LICENSE-2.0) for details)
- **Modifications**:  
  - Adjusted to support batch evaluation from a json alignment file.  
  - Uses residue indices directly for alignments instead of parsing CIGAR strings.  
  - Computes **average lDDT** over aligned residues instead of per-residue scores.

**Usage Example**:

```bash
python compute_lddt.py \
    --alignment alignment_results.tsv \
    --query /path/to/query_pdbs \
    --target /path/to/target_pdbs
