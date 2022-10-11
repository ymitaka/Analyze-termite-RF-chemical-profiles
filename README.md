<h2>Read me</h2>

<h3>Introduction</h3>

&ensp;&ensp; These text and R command files are used for automatically comparing lipid or peptide profiles, which were analyzed by using LC-MS/MS, between king food (KF) and queen food (QF) or among midgut contents of kings, queens, soldiers, and workers in the Japanese subterranean termite <i>Reticulitermes speratus</i>. <br>
&ensp;&ensp; The .R files included in this repository not only outputs a Mean-EPI (Excess proportion index) plot comparing the intensity of each ion in royal foods between KF and QF, but also a heatmap comparing the intensity of ions detected in midgut contents between castes, which are used in our research article (Eisuke et al. Under review). 
<br><br>

<h3>Preparation in advance</h3>

&ensp;&ensp; As a preliminary preparation, a group of .raw format files containing the analysis results of each sample saved by Xcalibur (the LC-MS/MS control software, Thermo Fisher) is required. These .raw files are registered in DDBJ MetaboBank (BioProject: PRJDB14286, ), and therefore, please download these files before starting to use the codes in this repository. 
<br><br>

<h3>How to use</h3>

<h4>STEP 1. Merging the lists of candidates for product-ion-detected (MS2-detected) precursor ions</h4>
&ensp;&ensp; First, download all the .R files and folders in this repository, open a Terminal.app on your mac (on Windows, use another suitable app) and set the directory to the "MS2-detected-lipid-lists" folder. Subsequently, run the shell scripts written in "UNIX-shell-script-(Lipids)" file to output new text files named "RFtoPK_lip.txt", "RFtoSQ_lip.txt", and "MidgutAll_lip.txt". "RFtoPK_lip.txt" and "RFtoSQ_lip.txt" contain the list of MS2-detected precursor ions common among termite colonies in KF and QF, respectively. "MidgutAll_lip.txt" contains the list of all MS2-detected precursor ions in midgut contents of all castes of all colonies. <br>
&ensp;&ensp; Perform the same operation for peptide ion candidates.

<h4>STEP 2. Generate plots and heatmaps</h4>
&ensp;&ensp; 
