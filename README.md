# codondiffR

An R package for the calculation and visualisation of codon usage metrics in user-supplied protein-coding nucleotide sequences, and the comparison with those of reference taxa. The package is designed primarily for the analysis of viral open reading frames (ORFs), with the intention of identifying potential host taxa based on overall codon usage similarity, but it can be used for the comparison of any protein-coding input sequence(s).

Pre-defined codon usage statistics for reference taxa come from the latest release of the Codon Usage Table Database made by [Athey et al. (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28865429), and the mean codon frequency usage difference (MCUFD) metric is calculated as described in [Stedman et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23308027).

## Getting Started

To install codondiffR, run the following commands in R:
```{r}
# install.packages("devtools")
devtools::install_github("adamd3/codondiffR")
```
