## version 0.2.2

---
### Updates.
Manuscript accepted in [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac179/6553658?login=true).


### Bugs fixed.
-prepro_na: An error is displayed if any features contains more than 20%
of missing data.

## version 0.2.1

---

### Bugs fixed.
-cov.adj: in this update cov.adj allows to have missing values in covariates.
Also, when used outside of $moss$, the omic block to adjust can also contain
missing values.

### New function.
-moss_venn.

## version 0.2.0

---


### Extra features included in this version.

- Option to tune degrees of sparsity and tSNE maps in parallel.
- Option to adjust omic blocks for covariates effects.
- Option to select degrees of sparsity by 'liberal' or 'conservative' methods.
- Mean-Imputation of missing data.
- Summaries of features selected by omic block.
- Features signatures by cluster of subjects.
- Bi-cluster representations using 'ComplexHeatmap'.

### New functions.
- moss_select.
- moss_signatures.
- moss_heatmap.
- cov_adj.
