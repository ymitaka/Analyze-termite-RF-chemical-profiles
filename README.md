# How to use

Introduction

   These text and R command files are used for automatically comparing lipid or peptide profiles, which were analyzed by using LC-MS/MS, between king and queen foods or among midgut contents of kings, queens, soldiers, and workers in the Japanese subterranean termite <i>Reticulitermes speratus</i>. 
   The .R files included in this repository not only outputs a Mean-EPI (Excess proportion index) plot comparing the intensity of each ion in royal foods between KF and QF, but also a heatmap comparing the intensity of ions detected in midgut contents between castes, which are used in our research article (Eisuke et al. Under review).

Preparation in advance

   As a preliminary preparation, a group of .raw format files containing the analysis results of each sample saved by Xcalibur (the LC-MS/MS control software, Thermo Fisher) and a group of .txt files containing a list of precursor ions that could be detected up to product ion in each sample that were output separately are required. The following files must be prepared. The former .raw files are registered in DDBJ MetaboBank (BioProject: PRJDB14286, ).The latter .txt files are contained in this repository. 

