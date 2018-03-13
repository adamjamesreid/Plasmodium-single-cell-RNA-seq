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

Other code

- Fourier_code.tar.gz contains the code used to identify the peak of expression for cycling genes

  To unpack use:
    tar -xvf Fourier_code.tar.gz
    
  The resulting README.rtf file explains how to compile the code and what the input/outputs are for the find_periodicity  script
  
- stage_rankmethod.pl was used to identify the most likely cell cycle stage for each single cell based on bulk RNA-seq data

  It takes as input:
  
  1. A matrix of normalised gene expression values for your cells with cells in columns and genes in rows (first column header should be 'id')
  2. A matrix of normalised gene expression values for the reference data in the same format as (1) (see example file Derisi_3d7_smoothed.txt.redPct)
  3. An expression cutoff value for the query data e.g. only genes with an expression level above this value will be included for a particular cell. The default is 10 and is suggested for linear FPKM values. For lscran values we recommend 3.
  4. An optional file with a mapping between the sample name in the reference data matrix and another useful name. This was used because in the DeRisi data, some samples with different names relate to the same timepoint (see example file Derisi_3d7_tp_mapping.txt)
  
  The output is:
  
  Column 1: Sample name
  Column 2: Best prediction
  Column 3: Number of genes used in the comparison with bulk
  Column 4: Spearmans r
  Column 5: Standard deviation
  Column 6: String of Spearman's r values for each bulk reference sample
