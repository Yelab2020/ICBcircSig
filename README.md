# ICBcircSig
ICBcircSig is to calculate candidate circRNA and construct circRNA score to predict ICB response. We provided all main codes of calling circRNA
of four well known method, identify differential expressed circRNAs and filtering candidate circRNA then build the circRNA score. The main process was described below.

## 1. Identification of circRNAs
Four tools, including CIRI2, find_circ, CircExplorer2, and CircRNA_finder were applied to identify circRNA with default settings. After FastQC for assessment
of the data quality, reads that passed thresholds were aligned to reference genome (GRCh38) using hisat275 with the default setting to obtain mapped and unmapped
reads in bam files. Unmapped reads were retrieved by samtools from bam files, and unmapped reads in fastq format were done by bedtools bamtofastq.
We employed each program to identify circRNAs with default parameters and annotated with gencode_v28. 

## 2. Identification differential expressed circRNAs
To identify differentially expressed circRNAs between responders and non-responders samples in pre-treatment (PRE) and early during treatment (EDT), respectively, 
a linear mixed-effects model(LME) which allows for nested random effects (each individual sample) and considers for  potential confounding factors was utilized and 
executed by the lme program in the nlme R package.

## 3. filtering candidate circRNA and build the circRNA score for ICB
We utilized a machine learning-based algorithm as previously described to construct ICBcircSig. Briefly, (i) we performed univariate survival analysis to identify 
prognosis relevant circRNAs by assessing the association of progression-free survival (PFS) and the expression of circRNAs; (ii) Based on LASSO Cox regression model, 
cv.glmnet function in R package glmnet was used. We first set seed 123, and deviance to measure loss to use for cross-validation and 5 folds, to develop the LASSO 
Cox regression model. Then we filtered circRNAs with lambda coefficient > 0 to retain the optimal combination from circRNAs in (i). The final signature, named “ICBcirSig”,
include significant circRNAs (p < 0.05) by multivariate cox analysis of circRNAs in (ii). (iii) The ICBcircSig score of each sample was built through the following 
equations based on the expression value and multivariate Cox regression coefficient (1.001 ∗ circ-TMTC3 + 1.048 ∗ circ-FAM117B).

