# TE.pipeline
## Introduction
TE.pipeline is used for analysis of TE loci across samples.

<img width="1733" height="1196" alt="image" src="https://github.com/user-attachments/assets/14892636-035a-462c-8374-3f0a5ba59dd8" />
Due to the imperfectness of the existing TE quantitative methods, there are often situations where it is impossible to achieve site-level quantification or the calculations are overly complex. Therefore, based on the traditional alignment methods and considering the positional relationship of genes and TE in the genome, we have developed a new quantitative process and subsequent analysis procedure for TE. The specific analysis process consists of several steps. RNA-seq reads were re-aligned to the mouse reference genome using STAR to generate intermediate BAM files. TE annotation file was obtained from the Dfam database (https://www.dfam.org/home). TE expression matrix were then quantified using FeatureCounts. 

## Getting Started
