<h2>Read me</h2>

<h3>Introduction</h3>

&ensp;&ensp; These text and R code files are used for automatically comparing lipid or peptide profiles, which were analyzed by using LC-MS/MS, between king food (KF) and queen food (QF) or among midgut contents of kings, queens, soldiers, and workers in the Japanese subterranean termite <i>Reticulitermes speratus</i>. <br>
&ensp;&ensp; The .R files included in this repository not only outputs a Mean-EPI (Excess proportion index) plot comparing the intensity of each ion in royal foods between KF and QF, but also a heatmap comparing the intensity of ions detected in midgut contents between castes, which are used in our research article <b>(Eisuke et al. 2023 PNAS Nexus DOI: <a href="https://doi.org/10.1093/pnasnexus/pgad222">10.1093/pnasnexus/pgad222</a>)</b>. 
<br><br>

<h3>Preparation in advance</h3>

&ensp;&ensp; As a preliminary preparation, a group of .txt files beggining with the name "Lipid spectra list MS1" (lipid data) or "Specta list MS1" (peptide data), which contains the analysis results of each sample exported from mzMine2 software, is also required. These .txt files are registered in DDBJ MetaboBank (BioProject: <b>PRJDB14286</b>, ). So, please download not only all the .R files and folders (beggining with the name "MS2-detected") in this repository but also the .txt files containing LC-MS/MS analysis results from DDBJ server, and then save these files/folders in the same single folder. <br>
&ensp;&ensp; In this automatic analyses, you will need to use <b>UNIX shell scripts</b> and <b>R software</b> (<a href="https://www.r-project.org/">The R project for Statistical Computing</a>). To run the shell scripts, you will need to install GNU grip (ggrep) by means of Homebrew on Terminal.app (for Mac/Linux users) or on other suitable software (for Windows users). Also, to run the R codes, you will need to install four additional packages "dplyr", "reshape2", "ggplot2", and "gplots" on your computer. 
<br><br>

<h3>How to use</h3>

First, as an example, a plot and heat map will be created using lipid data.<br>
<ol>
<li>Open a Terminal.app on your mac (on Windows, use another suitable app) and set the directory to the <b>"MS2-detected-lipid-lists" folder</b>. Subsequently, run the shell scripts written in <b>"UNIX-shell-script-(Lipids)"</b> file to output new text files named <b>"RFtoPK_lip.txt", "RFtoSQ_lip.txt", and "MidgutAll_lip.txt"</b>. (RFtoPK_lip.txt and RFtoSQ_lip.txt contain the lists of product-ion-detected (MS2-detected) precursor ions common among termite colonies in KF and QF, respectively. MidgutAll_lip.txt contains the list of all MS2-detected precursor ions in midgut contents of all castes of all colonies.)<br></li>

<li>Move (or copy and paste) these 3 files (RFtoPK_lip.txt, RFtoSQ_lip.txt, and MidgutAll_lip.txt) to a folder one level above. This folder should contain not only these three text files but also the .R files, and a group of text files downloaded from DDBJ, which contains lipid analysis results.<br></li>

<li>Open "<b>RF MC Lipid profile comparison (for Github).R</b>" file and run the code in order from top to bottom. The execution of this code will save in a folder a group of csv files containing data sets in the intermediate stages of processing, as well as csv files containing data sets used directly for plotting and drawing heatmaps.<br></li>

<li>A plot and a heatmap will be generated and displayed in a new window on the R software. Save them in PDF (or JPEG, PNG, and so on) format.</li>
</ol>

Perform the same operation for peptide ion candidates. <br>
