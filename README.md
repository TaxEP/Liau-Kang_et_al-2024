## Disparity and sampling analysis for _Mimosa_ pollen grains morphology

**1. analysis.R script**

**1.1. Section "Sampling Analysis"**
This section provides an overview of the current phylogenetic coverage of pollen morphology (Vasconcelos et al. 2018 Mimosa phylogeny) and the percentages of pollen data present for the main Mimosa clades. This information was previously presented in the "Notes on pollen disparity" section of the Discussion.

**1.2. Section "Tree Plot"**
This section outputs the annotated phylogeny in Figure 1 of the manuscript. It maps the current pollen sampling in the Mimosa genus, highlighting phylogeny regions that lack pollen descriptions.

**1.3. Section "Disparity analysis"**
This section outputs graphs used to generate figure 8 of the manuscript:
- As described in the main text, the original dataset (**pollen_data-mimosa.csv**) included 139 taxa in total and 13 variables. However, we excluded two variables for lack of variation (grain cohesion and aperture type - this are not even present in **pollen_data-mimosa.csv**), and four other (number of pores, pore diameter, and thickness of sexine an exine) that lacked more than 60% of representation to avoid missing data biases.
- Then, to infer the morphospace, we first discretized continuous characters using the Freedman-Diaconis (F-D) rule for histograms (function _discretizeDF_) from the _arules_ package. This procedure was necessary due to restrictions for joint use of continuous and categorical data in currently available software.
- Afterwards, a distance matrix using the Maximum Observable Distance (MORD) with the _Claddis_ package was calculated, and then used as the input in a Principal Coordinate Analysis (PCoA). We resolved intraspecific polymorphism by selecting the state that contributed the smallest pairwise distance and, to avoid negative-eigenvalues, we applied a Cailliez correction to the PCoA.
