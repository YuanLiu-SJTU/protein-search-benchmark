# Impact of Predicted Structure Quality on Search Performance

Recent advances in protein structure prediction, exemplified by **AlphaFold** and related methods, have enabled large-scale proteome-wide structural databases. Consequently, querying against predicted structures has become increasingly important. Understanding how **prediction uncertainty** influences search performance is critical for proper interpretation and method design.

## Benchmark Setup

- **Dataset**: Proteome-wide predicted structures stratified by local confidence (**pLDDT** scores).
- **Evaluation**: Functional relevance assessed via **GO semantic similarity** of top-10 search hits.
- **pLDDT stratification**:  
  1. **High confidence**: pLDDT ≥ 90  
  2. **Moderate confidence**: 70 ≤ pLDDT < 90  
  3. **Low confidence**: pLDDT < 70  

## Notes

- Evaluation uses the **same GO semantic similarity pipeline** as the `Functional_consistency` and `IDP` folders.  
