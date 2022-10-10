### Royal food analysis from LC-MS/MS results

## Peptides


### Comparisons of peptide or lipid profile between king food (KF) and queen food (QF) 

## 1. Extract m/z and the strongest peak intensity values of each precursor ion
## (If CSV files have already been exported, you can skip the following commands.)

p_kf6 <- read.table("Spectra List MS1 - RF6 To PK.txt", header=T)
p_qf6 <- read.table("Spectra List MS1 - RF6 To SQ.txt", header=T)
p_kf7 <- read.table("Spectra List MS1 - RF7 To PK.txt", header=T)
p_qf7 <- read.table("Spectra List MS1 - RF7 To SQ.txt", header=T)
colnames(p_kf6) <- c("mz", "Intensity")
colnames(p_qf6) <- c("mz", "Intensity")
colnames(p_kf7) <- c("mz", "Intensity")
colnames(p_qf7) <- c("mz", "Intensity")
p_mc_k6 <- read.table("Specta list MS1 - Midgut6 PK.txt", header=T)
p_mc_q6 <- read.table("Specta list MS1 - Midgut6 SQ.txt", header=T)
p_mc_s6 <- read.table("Specta list MS1 - Midgut6 Sol.txt", header=T)
p_mc_w6 <- read.table("Specta list MS1 - Midgut6 W.txt", header=T)
p_mc_k7 <- read.table("Specta list MS1 - Midgut7 PK.txt", header=T)
p_mc_q7 <- read.table("Specta list MS1 - Midgut7 SQ.txt", header=T)
p_mc_s7 <- read.table("Specta list MS1 - Midgut7 Sol.txt", header=T)
p_mc_w7 <- read.table("Specta list MS1 - Midgut7 W.txt", header=T)
colnames(p_mc_k6) <- c("mz", "Intensity")
colnames(p_mc_q6) <- c("mz", "Intensity")
colnames(p_mc_s6) <- c("mz", "Intensity")
colnames(p_mc_w6) <- c("mz", "Intensity")
colnames(p_mc_k7) <- c("mz", "Intensity")
colnames(p_mc_q7) <- c("mz", "Intensity")
colnames(p_mc_s7) <- c("mz", "Intensity")
colnames(p_mc_w7) <- c("mz", "Intensity")

# Ignore the information of retention time
p_kf6 <- p_kf6[!p_kf6$mz =="Name:", ]
p_qf6 <- p_qf6[!p_qf6$mz =="Name:", ]
p_kf7 <- p_kf7[!p_kf7$mz =="Name:", ]
p_qf7 <- p_qf7[!p_qf7$mz =="Name:", ]
p_mc_k6 <- p_mc_k6[!p_mc_k6$mz =="Name:", ]
p_mc_q6 <- p_mc_q6[!p_mc_q6$mz =="Name:", ]
p_mc_s6 <- p_mc_s6[!p_mc_s6$mz =="Name:", ]
p_mc_w6 <- p_mc_w6[!p_mc_w6$mz =="Name:", ]
p_mc_k7 <- p_mc_k7[!p_mc_k7$mz =="Name:", ]
p_mc_q7 <- p_mc_q7[!p_mc_q7$mz =="Name:", ]
p_mc_s7 <- p_mc_s7[!p_mc_s7$mz =="Name:", ]
p_mc_w7 <- p_mc_w7[!p_mc_w7$mz =="Name:", ]

# m/z values and intensities are rounded off to 2 and 4 decimal places, respectively
p_kf6[, 1] <- round(as.numeric(p_kf6[, 1]), digits=4)
p_kf6[, 2] <- round(as.numeric(p_kf6[, 2]), digits=4)
p_qf6[, 1] <- round(as.numeric(p_qf6[, 1]), digits=4)
p_qf6[, 2] <- round(as.numeric(p_qf6[, 2]), digits=4)
p_kf7[, 1] <- round(as.numeric(p_kf7[, 1]), digits=4)
p_kf7[, 2] <- round(as.numeric(p_kf7[, 2]), digits=4)
p_qf7[, 1] <- round(as.numeric(p_qf7[, 1]), digits=4)
p_qf7[, 2] <- round(as.numeric(p_qf7[, 2]), digits=4)
p_mc_k6[, 1] <- round(as.numeric(p_mc_k6[, 1]), digits=4)
p_mc_k6[, 2] <- round(as.numeric(p_mc_k6[, 2]), digits=4)
p_mc_q6[, 1] <- round(as.numeric(p_mc_q6[, 1]), digits=4)
p_mc_q6[, 2] <- round(as.numeric(p_mc_q6[, 2]), digits=4)
p_mc_s6[, 1] <- round(as.numeric(p_mc_s6[, 1]), digits=4)
p_mc_s6[, 2] <- round(as.numeric(p_mc_s6[, 2]), digits=4)
p_mc_w6[, 1] <- round(as.numeric(p_mc_w6[, 1]), digits=4)
p_mc_w6[, 2] <- round(as.numeric(p_mc_w6[, 2]), digits=4)
p_mc_k7[, 1] <- round(as.numeric(p_mc_k7[, 1]), digits=4)
p_mc_k7[, 2] <- round(as.numeric(p_mc_k7[, 2]), digits=4)
p_mc_q7[, 1] <- round(as.numeric(p_mc_q7[, 1]), digits=4)
p_mc_q7[, 2] <- round(as.numeric(p_mc_q7[, 2]), digits=4)
p_mc_s7[, 1] <- round(as.numeric(p_mc_s7[, 1]), digits=4)
p_mc_s7[, 2] <- round(as.numeric(p_mc_s7[, 2]), digits=4)
p_mc_w7[, 1] <- round(as.numeric(p_mc_w7[, 1]), digits=4)
p_mc_w7[, 2] <- round(as.numeric(p_mc_w7[, 2]), digits=4)

# Sort m/z values in an ascending order and intensity in a descending order
library(dplyr)
p_kf6 <- arrange(p_kf6, mz, desc(Intensity))
p_qf6 <- arrange(p_qf6, mz, desc(Intensity))
p_kf7 <- arrange(p_kf7, mz, desc(Intensity))
p_qf7 <- arrange(p_qf7, mz, desc(Intensity))
p_mc_k6 <- arrange(p_mc_k6, mz, desc(Intensity))
p_mc_q6 <- arrange(p_mc_q6, mz, desc(Intensity))
p_mc_s6 <- arrange(p_mc_s6, mz, desc(Intensity))
p_mc_w6 <- arrange(p_mc_w6, mz, desc(Intensity))
p_mc_k7 <- arrange(p_mc_k7, mz, desc(Intensity))
p_mc_q7 <- arrange(p_mc_q7, mz, desc(Intensity))
p_mc_s7 <- arrange(p_mc_s7, mz, desc(Intensity))
p_mc_w7 <- arrange(p_mc_w7, mz, desc(Intensity))

# Remove duplicated m/z values
p_kf6n <- distinct(p_kf6, mz, .keep_all=T)
p_qf6n <- distinct(p_qf6, mz, .keep_all=T)
p_kf7n <- distinct(p_kf7, mz, .keep_all=T)
p_qf7n <- distinct(p_qf7, mz, .keep_all=T)
p_mc_k6n <- distinct(p_mc_k6, mz, .keep_all=T)
p_mc_q6n <- distinct(p_mc_q6, mz, .keep_all=T)
p_mc_s6n <- distinct(p_mc_s6, mz, .keep_all=T)
p_mc_w6n <- distinct(p_mc_w6, mz, .keep_all=T)
p_mc_k7n <- distinct(p_mc_k7, mz, .keep_all=T)
p_mc_q7n <- distinct(p_mc_q7, mz, .keep_all=T)
p_mc_s7n <- distinct(p_mc_s7, mz, .keep_all=T)
p_mc_w7n <- distinct(p_mc_w7, mz, .keep_all=T)

# Cut off the rows that intensity is less than 1000
p_kf6n <- subset(p_kf6n, Intensity > 1000)
p_qf6n <- subset(p_qf6n, Intensity > 1000)
p_kf7n <- subset(p_kf7n, Intensity > 1000)
p_qf7n <- subset(p_qf7n, Intensity > 1000)
p_mc_k6n <- subset(p_mc_k6n, Intensity > 1000)
p_mc_q6n <- subset(p_mc_q6n, Intensity > 1000)
p_mc_s6n <- subset(p_mc_s6n, Intensity > 1000)
p_mc_w6n <- subset(p_mc_w6n, Intensity > 1000)
p_mc_k7n <- subset(p_mc_k7n, Intensity > 1000)
p_mc_q7n <- subset(p_mc_q7n, Intensity > 1000)
p_mc_s7n <- subset(p_mc_s7n, Intensity > 1000)
p_mc_w7n <- subset(p_mc_w7n, Intensity > 1000)

# Export csv files
write.csv(p_kf6n, "Peptides MS1 mz and peak 3 - KF colony6.csv", row.names=F)
write.csv(p_qf6n, "Peptides MS1 mz and peak 3 - QF colony6.csv", row.names=F)
write.csv(p_kf7n, "Peptides MS1 mz and peak 3 - KF colony7.csv", row.names=F)
write.csv(p_qf7n, "Peptides MS1 mz and peak 3 - QF colony7.csv", row.names=F)
write.csv(p_mc_k6n, "Peptides MS1 mz and peak 3 - MC K colony6.csv", row.names=F)
write.csv(p_mc_q6n, "Peptides MS1 mz and peak 3 - MC Q colony6.csv", row.names=F)
write.csv(p_mc_s6n, "Peptides MS1 mz and peak 3 - MC S colony6.csv", row.names=F)
write.csv(p_mc_w6n, "Peptides MS1 mz and peak 3 - MC W colony6.csv", row.names=F)
write.csv(p_mc_k7n, "Peptides MS1 mz and peak 3 - MC K colony7.csv", row.names=F)
write.csv(p_mc_q7n, "Peptides MS1 mz and peak 3 - MC Q colony7.csv", row.names=F)
write.csv(p_mc_s7n, "Peptides MS1 mz and peak 3 - MC S colony7.csv", row.names=F)
write.csv(p_mc_w7n, "Peptides MS1 mz and peak 3 - MC W colony7.csv", row.names=F)

#########################################################################


## 2. Merge the above CSV files containing m/z and intensity values of each precursor ion
## (If CSV files have already been exported, you can skip the following commands.)

# Read CSV files
p_kf6 <- read.csv("Peptides MS1 mz and peak 3 - KF colony6.csv")
p_qf6 <- read.csv("Peptides MS1 mz and peak 3 - QF colony6.csv")
p_kf7 <- read.csv("Peptides MS1 mz and peak 3 - KF colony7.csv")
p_qf7 <- read.csv("Peptides MS1 mz and peak 3 - QF colony7.csv")
p_mc_k6 <- read.csv("Peptides MS1 mz and peak 3 - MC K colony6.csv")
p_mc_q6 <- read.csv("Peptides MS1 mz and peak 3 - MC Q colony6.csv")
p_mc_s6 <- read.csv("Peptides MS1 mz and peak 3 - MC S colony6.csv")
p_mc_w6 <- read.csv("Peptides MS1 mz and peak 3 - MC W colony6.csv")
p_mc_k7 <- read.csv("Peptides MS1 mz and peak 3 - MC K colony7.csv")
p_mc_q7 <- read.csv("Peptides MS1 mz and peak 3 - MC Q colony7.csv")
p_mc_s7 <- read.csv("Peptides MS1 mz and peak 3 - MC S colony7.csv")
p_mc_w7 <- read.csv("Peptides MS1 mz and peak 3 - MC W colony7.csv")

# Add the information about food type and colony
p_kf6n <- cbind(Food=c("KF"), Colony=c("MT552"), p_kf6)
p_qf6n <- cbind(Food=c("QF"), Colony=c("MT552"), p_qf6)
p_kf7n <- cbind(Food=c("KF"), Colony=c("MT548"), p_kf7)
p_qf7n <- cbind(Food=c("QF"), Colony=c("MT548"), p_qf7)
p_mc_k6n <- cbind(Food=c("MCK"), Colony=c("MT552"), p_mc_k6)
p_mc_q6n <- cbind(Food=c("MCQ"), Colony=c("MT552"), p_mc_q6)
p_mc_s6n <- cbind(Food=c("MCS"), Colony=c("MT552"), p_mc_s6)
p_mc_w6n <- cbind(Food=c("MCW"), Colony=c("MT552"), p_mc_w6)
p_mc_k7n <- cbind(Food=c("MCK"), Colony=c("MT548"), p_mc_k7)
p_mc_q7n <- cbind(Food=c("MCQ"), Colony=c("MT548"), p_mc_q7)
p_mc_s7n <- cbind(Food=c("MCS"), Colony=c("MT548"), p_mc_s7)
p_mc_w7n <- cbind(Food=c("MCW"), Colony=c("MT548"), p_mc_w7)

# Merge all data
d <- rbind(p_kf6n, p_qf6n, p_kf7n, p_qf7n, p_mc_k6n, p_mc_q6n, p_mc_s6n, p_mc_w6n, p_mc_k7n, p_mc_q7n, p_mc_s7n, p_mc_w7n)

library(dplyr)
library(reshape2)
d2n <- dcast(d, mz ~ Food + Colony, sum)

write.csv(d2n, "Peptides MS1 raw mz and peak list v2.csv", row.names=F)	#This file will be used not only for comparing chemical profile between king food and queen food but also for comparing midgut content's chemical profile among castes.

############################################################################


## 3. Find MS2-interpretable precursor ions
## (If the dataset (CSV file) for generating a plot have already been exported, you can skip this section.)

d2n <- read.csv("Peptides MS1 raw mz and peak list v2.csv")
d2n2 <- d2n[, c(1:3, 12, 13)]	#Extract KF & QF data
library(dplyr)
library(reshape2)

# Make the data frame of MS2-detected precursor ions 
ms2_kf <- read.table("RFtoPK_pep.txt", header=F)	# KF common
ms2_qf <- read.table("RFtoSQ_pep.txt", header=F)	# QF common

ms2_kf2 <- cbind(Food=c("KF"), ms2_kf)
ms2_qf2 <- cbind(Food=c("QF"), ms2_qf)
ms2 <- rbind(ms2_kf2, ms2_qf2)
colnames(ms2) <- c("Food", "mz")

# Extract MS2-detected precursor ions
d3 <- mutate(d2n2, KF_mean=(KF_MT548 + KF_MT552)/2, QF_mean=(QF_MT548+QF_MT552)/2, Diff=QF_mean - KF_mean, Average=(KF_mean + QF_mean)/2, logA=log10(Average), EPI=Diff/(KF_mean + QF_mean),  Bias= ifelse(KF_mean/QF_mean >= 2, "KF-biased", ifelse(QF_mean/KF_mean >= 2, "QF-biased", "Common to KF & QF")))
d3 <-cbind(d3, MS2=0)

for(i in 1:nrow(ms2)){
	row.no <- grep(ms2$mz[i], d3$mz)
	d3$MS2[row.no] <- 1
	}
d4 <- d3[d3$MS2==1, ]

d5n <- mutate(d4, mz.01=round(mz, digits=2))

# Sort m/z values in an ascending order and intensity in a descending order
d5n <- arrange(d5n, mz.01, desc(Average))
# Remove duplicated m/z values
d5 <- distinct(d5n, mz.01, .keep_all=T)
d5 <- d5[d5$Average > 0, ]

# At this time, the numer of ions is 172.

write.csv(d5, "Peptides MS2 raw mz and peak list v4.csv", row.names=F)

############################################################################


## 4. Generage the KF vs QF plot 

# Mean-EPI plot
d5 <- read.csv("Peptides MS2 raw mz and peak list v4.csv")
d5$Bias <- factor(d5$Bias, levels=c("KF-biased", "Common to KF & QF", "QF-biased"))
 library(ggplot2)
 p <- ggplot(d5, aes(x=EPI, y=logA, colour=Bias))　
 p <- p + geom_point(size=2) + theme_bw() + xlim(-1.1, 1.1) + 
 			geom_vline(xintercept=0.33, linetype="dashed") + 
 			geom_vline(xintercept=-0.33, linetype="dashed") + 
  			labs(x="EPI", y="log10(Average)", title="Peptides") + 
 			scale_colour_manual(values=c("deepskyblue", "grey60", "orange"))
 dev.new(width=7, height=4)
 plot(p)

############################################################################
############################################################################


### Heatmap for comparing mean peak intensity among midgut contents

## 1. Find MS2-interpretable precursor ions in midgut contents of all castes
## (If you finished generating a list of precursor ions in midgut contents of all castes, you can skip this section.)

d2n <- read.csv("Peptides MS1 raw mz and peak list v2.csv")
d2n <- d2n[, c(1, 4:11)]
library(dplyr)
library(reshape2)

# Extract MCK (or MCQ, MCS, MCW) ions detected in both colonies
d3 <- d2n[(d2n$MCK_MT548 != 0 & d2n$MCK_MT552 != 0) | (d2n$MCQ_MT548 != 0 & d2n$MCQ_MT552 !=0) | (d2n$MCS_MT548 != 0 & d2n$MCS_MT552 != 0) | (d2n$MCW_MT548 != 0 & d2n$MCW_MT552 != 0), ]

# m/z values is rounded down to the second decimal place
d3 <- mutate(d3, MS2=0, mz1=trunc(mz*10)/10, Duplicate=0)	

# Sort m/z values in an ascending order
d3 <- arrange(d3, mz1)
# Remove duplicated m/z values
d3 <- distinct(d3, mz1, .keep_all=T)

# Make the data frame of MS2-detected precursor ions 
ms2_mc <- read.table("MidgutAll_pep.txt", header=F)	# MS2-detected ions in midgut contents of all castes 
colnames(ms2_mc) <- c("mz")		# 28644 ions

# Extract MS2-detected precursor ions
for(i in 1:nrow(ms2_mc)){
	row.no <- grep(ms2_mc$mz[i], d3$mz)
	d3$MS2[row.no] <- 1
	}
d3 <- d3[d3$MS2==1, ]

# Calculate mean intensities
d3 <- mutate(d3, MCK_mean = (MCK_MT548 + MCK_MT552)/2, 
		MCQ_mean = (MCQ_MT548 + MCQ_MT552)/2, 
		MCS_mean = (MCS_MT548 + MCS_MT552)/2, 
		MCW_mean = (MCW_MT548 + MCW_MT552)/2)

write.csv(d3, "Peptides MS2 raw mz and peak list (MC) v2.csv", row.names=F)

############################################################################


## 2. Generate heatmap

d3 <- read.csv("Peptides MS2 raw mz and peak list (MC) v2.csv")
ion.values <- d3[, 15]
d4 <- as.matrix(d3[, 13:16])		# Exclude KF & QF data
samples <- c("King", "Queen", "Soldier", "Worker")
dimnames(d4) <- list(ion.values, samples)
d5 <- t(scale(t(d4), center=F))		# Scaling (SD=1, Mean is not zero)
library(gplots)
mycol <- colorpanel(128, low="white", high="red2")	# Manual color setting 
dev.new(width=5, height=10)
par(oma=c(3,1,2,1), mar=c(0,0,0,0))
heatmap.2(d5, Colv=F, Rowv=T, dendrogram="row", col=mycol, key=T, keysize=0.01, key.par=list(cex=0.5), density.info="none", trace="none", main="Peptides", lwid=c(3, 7), lhei=c(1.1, 11), cexRow=0.1, cexCol=0.8)

