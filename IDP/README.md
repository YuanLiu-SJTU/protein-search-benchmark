# Search Performance on Intrinsically Disordered Proteins (IDPs)

Intrinsically disordered proteins (IDPs) and regions (IDRs) play critical roles in cellular regulation, signaling, and molecular recognition, and constitute a substantial fraction of proteomes. Unlike well-folded proteins, IDPs lack stable tertiary structures under physiological conditions and often exert their functions through transient, context-dependent interactions. As a result, identifying functionally related IDPs remains an important and challenging task for protein search methods.

## Benchmark Dataset

- **Query set**: 534 experimentally validated IDP sequences from **DisProt**, with pairwise sequence identity <20%.
- **Target set**: SwissProt human dataset containing 6955 proteins, also filtered to <20% sequence identity.
- **Functional evaluation**: GO semantic similarity computed for the **top 10 query-target pairs**.

This dataset ensures low redundancy and emphasizes remote functional relationships among IDPs.

## Notes

- Evaluation code is **identical to the functional consistency analysis** used for SwissProt structured proteins (see `functional_consistency` folder).
- Only the dataset differs; no additional code is required.
- The benchmark demonstrates the need for methods and evaluation frameworks **tailored to disordered proteins**.
