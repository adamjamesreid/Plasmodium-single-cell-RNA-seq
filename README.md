# Plasmodium-single-cell-RNA-seq

This repository contains files related to the publication Single cell RNA-seq reveals hidden transcriptional variation in malaria parasites by Reid et al. 

This paper describes the development and application of a Smart-seq2-based single-cell RNA-seq approach to the study of Plasmodium parasite transcriptomes

The files here include read counts, meta data and R code for our key datasets describing the late asexual cycle and mature gametocytes of the human malaria parasite Plasmodium falciparum and the rodent malaria model parasite Plasmodium berghei

Plasmodium berghei mixed blood stages (trophozoites, schizonts, male and female gametocytes)

- PbM_analysis.R is the R code to filter and normalise the data
- PbM_counts.txt is the raw count data
- PbM_meta.txt is the meta data describing the samples/cells
- berg.desc contains the product descriptions for the genes

Plasmodium falciparum late asexual stages (trophozoites and schizonts)

- PfAsex_analysis.R is the R code to filter and normalise the data
- PfAsex_counts.txt is the raw count data
- PfAsex_meta.txt is the meta data describing the samples/cells
- fal.desc contains the product descriptions for the genes

Plasmodium falciparum mature gametocytes

- PfGam_analysis.R is the R code to filter and normalise the data
- PfGam_counts.txt is the raw count data
- PfGam_meta.txt is the meta data describing the samples/cells
- fal.desc contains the product descriptions for the genes
