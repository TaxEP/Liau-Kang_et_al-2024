## Disparity and sampling analysis for _Mimosa_ pollen grains morphology

**1. morphospace.R script**

This script outputs the figure 8 of the manuscript (morphospace):
- As described in the main text, the original dataset (pollen_data-mimosa.csv) included 139 taxa in total and 11 variables. However, we excluded four variables (number of pores, pore diameter, and thickness of sexine an exine) that lacked more than 60% of representation to avoid missing data biases.
- Then, we discretized continuous characters using the Freedman-Diaconis (F-D) rule for histograms (function _discretizeDF_) from the _arules_ package. This procedure was necessary due to restrictions for joint use of continuous and categorical data in currently available software.
- Afterwards, a distance matrix using the Maximum Observable Distance (MORD) with the _Claddis_ package was calculated, and then used as the input in a Principal Coordinate Analysis (PCoA). We resolved intraspecific polymorphism by selecting the state that contributed the smallest pairwise distance and, to avoid negative-eigenvalues, we applied a Cailliez correction to the PCoA.

**2. sampling.R script**

This script supply the current phylogenetic coverage of pollen morphology (Vasconcelos et al. 2018 _Mimosa_ phylogeny) and the percentages of pollen data for present for the main _Mimosa_ clades. Such informations were presented in the "Notes on pollen disparity" section of the Discussion.
