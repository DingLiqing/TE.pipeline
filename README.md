# TE.pipeline
## Overview
TE.pipeline is used for analysis of TE loci across samples.

<img width="1733" height="1196" alt="image" src="https://github.com/user-attachments/assets/14892636-035a-462c-8374-3f0a5ba59dd8" />
Due to the imperfectness of the existing TE quantitative methods, there are often situations where it is impossible to achieve site-level quantification or the calculations are overly complex. Therefore, based on the traditional alignment methods and considering the positional relationship of genes and TE in the genome, we have developed a new quantitative process and subsequent analysis procedure for TE. The specific analysis process consists of several steps. RNA-seq reads were re-aligned to the mouse reference genome using STAR to generate intermediate BAM files. TE annotation file was obtained from the Dfam database (https://www.dfam.org/home). TE expression matrix were then quantified using FeatureCounts. 

## Requirements
1. Linux: (tested on Ubuntu 16.04.7)
2. R with tidyr, dply, and Matrix installed (only required for QC)
3. STAR v2.7.11 (does not work with earlier versions, might work with newer versions but not tested)
4. Subread v2.0.3 (does not work with earlier versions)
We are planning on creating a docker and/or singularity container for this pipeline as well, but have not done so yet.

## Contributors
TE.pipeline was developed by Liqing Ding. Please contact Liqing Ding (liqing_ding@sina.com) for any questions or suggestions.
