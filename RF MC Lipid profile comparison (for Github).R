### Royal food analysis from LC-MS/MS results

## Lipids


### Comparisons of peptide or lipid profile between king food (KF) and queen food (QF) 

## 1. Extract m/z and the strongest peak intensity values of each precursor ion
## (If CSV files have already been exported, you can skip this section)

l_kf6 <- read.table("Lipid spectra list MS1 - RF6 To PK.txt", header=T)
l_qf6 <- read.table("Lipid spectra list MS1 - RF6 To SQ.txt", header=T)
l_kf7 <- read.table("Lipid spectra list MS1 - RF7 To PK.txt", header=T)
l_qf7 <- read.table("Lipid spectra list MS1 - RF7 To SQ.txt", header=T)
colnames(l_kf6) <- c("mz", "Intensity")
colnames(l_qf6) <- c("mz", "Intensity")
colnames(l_kf7) <- c("mz", "Intensity")
colnames(l_qf7) <- c("mz", "Intensity")
l_mc_k6 <- read.table("Lipid spectra list MS1 - Midgut6 PK.txt", header=T)
l_mc_q6 <- read.table("Lipid spectra list MS1 - Midgut6 SQ.txt", header=T)
l_mc_s6 <- read.table("Lipid spectra list MS1 - Midgut6 Sol.txt", header=T)
l_mc_w6 <- read.table("Lipid spectra list MS1 - Midgut6 W.txt", header=T)
l_mc_k7 <- read.table("Lipid spectra list MS1 - Midgut7 PK.txt", header=T)
l_mc_q7 <- read.table("Lipid spectra list MS1 - Midgut7 SQ.txt", header=T)
l_mc_s7 <- read.table("Lipid spectra list MS1 - Midgut7 Sol.txt", header=T)
l_mc_w7 <- read.table("Lipid spectra list MS1 - Midgut7 W.txt", header=T)
colnames(l_mc_k6) <- c("mz", "Intensity")
colnames(l_mc_q6) <- c("mz", "Intensity")
colnames(l_mc_s6) <- c("mz", "Intensity")
colnames(l_mc_w6) <- c("mz", "Intensity")
colnames(l_mc_k7) <- c("mz", "Intensity")
colnames(l_mc_q7) <- c("mz", "Intensity")
colnames(l_mc_s7) <- c("mz", "Intensity")
colnames(l_mc_w7) <- c("mz", "Intensity")

# Ignore the information of retention time
l_kf6 <- l_kf6[!l_kf6$mz =="Name:", ]
l_qf6 <- l_qf6[!l_qf6$mz =="Name:", ]
l_kf7 <- l_kf7[!l_kf7$mz =="Name:", ]
l_qf7 <- l_qf7[!l_qf7$mz =="Name:", ]
l_mc_k6 <- l_mc_k6[!l_mc_k6$mz =="Name:", ]
l_mc_q6 <- l_mc_q6[!l_mc_q6$mz =="Name:", ]
l_mc_s6 <- l_mc_s6[!l_mc_s6$mz =="Name:", ]
l_mc_w6 <- l_mc_w6[!l_mc_w6$mz =="Name:", ]
l_mc_k7 <- l_mc_k7[!l_mc_k7$mz =="Name:", ]
l_mc_q7 <- l_mc_q7[!l_mc_q7$mz =="Name:", ]
l_mc_s7 <- l_mc_s7[!l_mc_s7$mz =="Name:", ]
l_mc_w7 <- l_mc_w7[!l_mc_w7$mz =="Name:", ]

# m/z values and intensities are rounded off to 4 and 4 decimal places, respectively
l_kf6[, 1] <- round(as.numeric(l_kf6[, 1]), digits=4)
l_kf6[, 2] <- round(as.numeric(l_kf6[, 2]), digits=4)
l_qf6[, 1] <- round(as.numeric(l_qf6[, 1]), digits=4)
l_qf6[, 2] <- round(as.numeric(l_qf6[, 2]), digits=4)
l_kf7[, 1] <- round(as.numeric(l_kf7[, 1]), digits=4)
l_kf7[, 2] <- round(as.numeric(l_kf7[, 2]), digits=4)
l_qf7[, 1] <- round(as.numeric(l_qf7[, 1]), digits=4)
l_qf7[, 2] <- round(as.numeric(l_qf7[, 2]), digits=4)
l_mc_k6[, 1] <- round(as.numeric(l_mc_k6[, 1]), digits=4)
l_mc_k6[, 2] <- round(as.numeric(l_mc_k6[, 2]), digits=4)
l_mc_q6[, 1] <- round(as.numeric(l_mc_q6[, 1]), digits=4)
l_mc_q6[, 2] <- round(as.numeric(l_mc_q6[, 2]), digits=4)
l_mc_s6[, 1] <- round(as.numeric(l_mc_s6[, 1]), digits=4)
l_mc_s6[, 2] <- round(as.numeric(l_mc_s6[, 2]), digits=4)
l_mc_w6[, 1] <- round(as.numeric(l_mc_w6[, 1]), digits=4)
l_mc_w6[, 2] <- round(as.numeric(l_mc_w6[, 2]), digits=4)
l_mc_k7[, 1] <- round(as.numeric(l_mc_k7[, 1]), digits=4)
l_mc_k7[, 2] <- round(as.numeric(l_mc_k7[, 2]), digits=4)
l_mc_q7[, 1] <- round(as.numeric(l_mc_q7[, 1]), digits=4)
l_mc_q7[, 2] <- round(as.numeric(l_mc_q7[, 2]), digits=4)
l_mc_s7[, 1] <- round(as.numeric(l_mc_s7[, 1]), digits=4)
l_mc_s7[, 2] <- round(as.numeric(l_mc_s7[, 2]), digits=4)
l_mc_w7[, 1] <- round(as.numeric(l_mc_w7[, 1]), digits=4)
l_mc_w7[, 2] <- round(as.numeric(l_mc_w7[, 2]), digits=4)

# Sort m/z values in an ascending order and intensity in a descending order
library(dplyr)
l_kf6 <- arrange(l_kf6, mz, desc(Intensity))
l_qf6 <- arrange(l_qf6, mz, desc(Intensity))
l_kf7 <- arrange(l_kf7, mz, desc(Intensity))
l_qf7 <- arrange(l_qf7, mz, desc(Intensity))
l_mc_k6 <- arrange(l_mc_k6, mz, desc(Intensity))
l_mc_q6 <- arrange(l_mc_q6, mz, desc(Intensity))
l_mc_s6 <- arrange(l_mc_s6, mz, desc(Intensity))
l_mc_w6 <- arrange(l_mc_w6, mz, desc(Intensity))
l_mc_k7 <- arrange(l_mc_k7, mz, desc(Intensity))
l_mc_q7 <- arrange(l_mc_q7, mz, desc(Intensity))
l_mc_s7 <- arrange(l_mc_s7, mz, desc(Intensity))
l_mc_w7 <- arrange(l_mc_w7, mz, desc(Intensity))

# Remove duplicated m/z values
l_kf6n <- distinct(l_kf6, mz, .keep_all=T)
l_qf6n <- distinct(l_qf6, mz, .keep_all=T)
l_kf7n <- distinct(l_kf7, mz, .keep_all=T)
l_qf7n <- distinct(l_qf7, mz, .keep_all=T)
l_mc_k6n <- distinct(l_mc_k6, mz, .keep_all=T)
l_mc_q6n <- distinct(l_mc_q6, mz, .keep_all=T)
l_mc_s6n <- distinct(l_mc_s6, mz, .keep_all=T)
l_mc_w6n <- distinct(l_mc_w6, mz, .keep_all=T)
l_mc_k7n <- distinct(l_mc_k7, mz, .keep_all=T)
l_mc_q7n <- distinct(l_mc_q7, mz, .keep_all=T)
l_mc_s7n <- distinct(l_mc_s7, mz, .keep_all=T)
l_mc_w7n <- distinct(l_mc_w7, mz, .keep_all=T)

# Cut off the rows that intensity is less than 1000
l_kf6n <- subset(l_kf6n, Intensity > 1000)
l_qf6n <- subset(l_qf6n, Intensity > 1000)
l_kf7n <- subset(l_kf7n, Intensity > 1000)
l_qf7n <- subset(l_qf7n, Intensity > 1000)
l_mc_k6n <- subset(l_mc_k6n, Intensity > 1000)
l_mc_q6n <- subset(l_mc_q6n, Intensity > 1000)
l_mc_s6n <- subset(l_mc_s6n, Intensity > 1000)
l_mc_w6n <- subset(l_mc_w6n, Intensity > 1000)
l_mc_k7n <- subset(l_mc_k7n, Intensity > 1000)
l_mc_q7n <- subset(l_mc_q7n, Intensity > 1000)
l_mc_s7n <- subset(l_mc_s7n, Intensity > 1000)
l_mc_w7n <- subset(l_mc_w7n, Intensity > 1000)

# Export csv files
write.csv(l_kf6n, "Lipids MS1 mz and peak 4 - KF colony6.csv", row.names=F)
write.csv(l_qf6n, "Lipids MS1 mz and peak 4 - QF colony6.csv", row.names=F)
write.csv(l_kf7n, "Lipids MS1 mz and peak 4 - KF colony7.csv", row.names=F)
write.csv(l_qf7n, "Lipids MS1 mz and peak 4 - QF colony7.csv", row.names=F)
write.csv(l_mc_k6n, "Lipids MS1 mz and peak 4 - MC K colony6.csv", row.names=F)
write.csv(l_mc_q6n, "Lipids MS1 mz and peak 4 - MC Q colony6.csv", row.names=F)
write.csv(l_mc_s6n, "Lipids MS1 mz and peak 4 - MC S colony6.csv", row.names=F)
write.csv(l_mc_w6n, "Lipids MS1 mz and peak 4 - MC W colony6.csv", row.names=F)
write.csv(l_mc_k7n, "Lipids MS1 mz and peak 4 - MC K colony7.csv", row.names=F)
write.csv(l_mc_q7n, "Lipids MS1 mz and peak 4 - MC Q colony7.csv", row.names=F)
write.csv(l_mc_s7n, "Lipids MS1 mz and peak 4 - MC S colony7.csv", row.names=F)
write.csv(l_mc_w7n, "Lipids MS1 mz and peak 4 - MC W colony7.csv", row.names=F)

##################################################################################


## 2. Merge the above CSV files containing m/z and intensity values of each precursor ion
## (If the merged CSV files have already been exported, you can skip this section)

# Read CSV files
l_kf6 <- read.csv("Lipids MS1 mz and peak 4 - KF colony6.csv")
l_qf6 <- read.csv("Lipids MS1 mz and peak 4 - QF colony6.csv")
l_kf7 <- read.csv("Lipids MS1 mz and peak 4 - KF colony7.csv")
l_qf7 <- read.csv("Lipids MS1 mz and peak 4 - QF colony7.csv")
l_mc_k6 <- read.csv("Lipids MS1 mz and peak 4 - MC K colony6.csv")
l_mc_q6 <- read.csv("Lipids MS1 mz and peak 4 - MC Q colony6.csv")
l_mc_s6 <- read.csv("Lipids MS1 mz and peak 4 - MC S colony6.csv")
l_mc_w6 <- read.csv("Lipids MS1 mz and peak 4 - MC W colony6.csv")
l_mc_k7 <- read.csv("Lipids MS1 mz and peak 4 - MC K colony7.csv")
l_mc_q7 <- read.csv("Lipids MS1 mz and peak 4 - MC Q colony7.csv")
l_mc_s7 <- read.csv("Lipids MS1 mz and peak 4 - MC S colony7.csv")
l_mc_w7 <- read.csv("Lipids MS1 mz and peak 4 - MC W colony7.csv")


# Add the information about food type and colony
l_kf6n <- cbind(Food=c("KF"), Colony=c("MT552"), l_kf6)
l_qf6n <- cbind(Food=c("QF"), Colony=c("MT552"), l_qf6)
l_kf7n <- cbind(Food=c("KF"), Colony=c("MT548"), l_kf7)
l_qf7n <- cbind(Food=c("QF"), Colony=c("MT548"), l_qf7)
l_mc_k6n <- cbind(Food=c("MCK"), Colony=c("MT552"), l_mc_k6)
l_mc_q6n <- cbind(Food=c("MCQ"), Colony=c("MT552"), l_mc_q6)
l_mc_s6n <- cbind(Food=c("MCS"), Colony=c("MT552"), l_mc_s6)
l_mc_w6n <- cbind(Food=c("MCW"), Colony=c("MT552"), l_mc_w6)
l_mc_k7n <- cbind(Food=c("MCK"), Colony=c("MT548"), l_mc_k7)
l_mc_q7n <- cbind(Food=c("MCQ"), Colony=c("MT548"), l_mc_q7)
l_mc_s7n <- cbind(Food=c("MCS"), Colony=c("MT548"), l_mc_s7)
l_mc_w7n <- cbind(Food=c("MCW"), Colony=c("MT548"), l_mc_w7)

# Merge all data
d <- rbind(l_kf6n, l_qf6n, l_kf7n, l_qf7n, l_mc_k6n, l_mc_q6n, l_mc_s6n, l_mc_w6n, l_mc_k7n, l_mc_q7n, l_mc_s7n, l_mc_w7n)

library(dplyr)
library(reshape2)
d2n <- dcast(d, mz ~ Food + Colony, sum)

write.csv(d2n2, "Lipids MS1 raw mz and peak list v3.csv", row.names=F)	#This file will be used not only for comparing chemical profile between king food and queen food but also for comparing midgut content's chemical profile among castes.

##################################################################################


## 3. Find MS2-interpretable precursor ions in royal food
## (If the dataset (CSV file) for generating a plot have already been exported, you can skip this section.)

d2n <- read.csv("Lipids MS1 raw mz and peak list v3.csv")
d2n2 <- d2n[, c(1:3, 12, 13)]  #Extract only KF & QF data
library(dplyr)
library(reshape2)

# Make the data frame of MS2-detected precursor ions 
ms2_kf <- read.table("RFtoPK_lip.txt", header=F)	# KF common
ms2_qf <- read.table("RFtoSQ_lip.txt", header=F)	# QF common

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

# Find the same compounds having different adduct ions (H+, NH4+, Na+)
# Atomic weight: H+: 1.0073, NH4+: 18.0338, Na+: 22.9893
d5n <- mutate(d4, mz.01=round(mz, digits=2), Duplicate=0)
for(i in 1:length(d5n)){
	for(j in 1:length(d5n)){
		# Search H+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 1.0073, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
		# Search NH4+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 18.0338, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
		# Search Na+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 22.9893, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
	}
}

# Sort m/z values in an ascending order and intensity in a descending order
d5n <- arrange(d5n, mz.01, desc(Average))
# Remove duplicated m/z values
d5 <- distinct(d5n, mz.01, .keep_all=T)

# Estimate noise ions based on the percentage of absolute value of mass defect to integer part of m/z value
d6 <- mutate(d5, Contaminant=ifelse(100*(mz - trunc(mz)) / trunc(mz) < 0.062 | 100*(mz - trunc(mz)) / trunc(mz) > 0.1, 1, 0))
d6 <- d6[d6$Contaminant==0, ]
d6 <- d6[d6$Average > 0, ]

write.csv(d6, "Lipids MS2 mz and peak list v3.csv", row.names=F)

##################################################################################


## 4. Generate the KF vs QF plot 

# Mean-EPI plot
d6 <- read.csv("Lipids MS2 mz and peak list v3.csv")
d6$Bias <- factor(d6$Bias, levels=c("KF-biased", "Common to KF & QF", "QF-biased"))
 library(ggplot2)
 p <- ggplot(d6, aes(x=EPI, y=logA, colour=Bias))
 p <- p + geom_point(size=2) + theme_bw() + xlim(-1.1, 1.1) + 
 			geom_vline(xintercept=0.33, linetype="dashed") + 
 			geom_vline(xintercept=-0.33, linetype="dashed") + 
  			labs(x="EPI", y="log10(Average)", title="Lipids") + 
 			scale_colour_manual(values=c("deepskyblue", "grey60", "orange"))
 dev.new(width=7, height=4)
 plot(p)

##################################################################################
##################################################################################


### Heatmap for comparing mean peak intensity among midgut contents

## 1. Find MS2-interpretable precursor ions in midgut contents
## (If you finished generating a list of precursor ions in midgut contents of all castes, you can skip this section.)

d2n <- read.csv("Lipids MS1 raw mz and peak list v3.csv")
library(dplyr)
library(reshape2)

# Extract MCK (or MCQ, MCS, MCW) ions detected in both colonies
d3 <- d2n[(d2n$MCK_MT548 != 0 & d2n$MCK_MT552 != 0) | (d2n$MCQ_MT548 != 0 & d2n$MCQ_MT552 !=0) | (d2n$MCS_MT548 != 0 & d2n$MCS_MT552 != 0) | (d2n$MCW_MT548 != 0 & d2n$MCW_MT552 != 0), ]
d3 <- d3[, c(1, 4:11)]   #Remove KF & QF intensity information


# Make the data frame of MS2-detected precursor ions 
ms2_mc <- read.table("MidgutAll_lip.txt", header=F)	# All MS2-detected ions in midgut contents of all castes
colnames(ms2_mc) <- c("mz")	

# Extract MS2-detected precursor ions
d3 <-cbind(d3, MS2=0)

for(i in 1:nrow(ms2_mc)){
	row.no <- grep(ms2_mc$mz[i], d3$mz)
	d3$MS2[row.no] <- 1
	}
d4 <- d3[d3$MS2==1, ]

# Find the same compounds having different adduct ions (H+, NH4+, Na+)
# Atomic weight: H+: 1.0073, NH4+: 18.0338, Na+: 22.9893
d5n <- mutate(d4, mz.01=round(mz, digits=2), Duplicate=0)
for(i in 1:length(d5n)){
	for(j in 1:length(d5n)){
		# Search H+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 1.0073, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
		# Search NH4+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 18.0338, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
		# Search Na+ added ions
		ifelse(all.equal(d5n$mz.01[j], d5n$mz.01[i] + 22.9893, tolerance=0.01), d5n$Duplicate[j] == 1, 0)
	}
}

# Sort m/z values in an ascending order
d5n <- arrange(d5n, mz.01)
# Remove duplicated m/z values
d5 <- distinct(d5n, mz.01, .keep_all=T)

# Estimate noise ions based on the percentage of absolute value of mass defect to integer part of m/z value
d6 <- mutate(d5, Contaminant=ifelse(100*(mz - trunc(mz)) / trunc(mz) < 0.062 | 100*(mz - trunc(mz)) / trunc(mz) > 0.1, 1, 0))
d6 <- d6[d6$Contaminant==0, ]

write.csv(d6, "Lipids MS2 mz and peak list (MC) v3.csv", row.names=F)

##################################################################################


## 2.  Generate heatmap

d6 <- read.csv("Lipids MS2 mz and peak list (MC) v3.csv")
library(dplyr)
library(reshape2)
d6n <- mutate(d6, MCK_mean = (MCK_MT548 + MCK_MT552)/2, 
		MCQ_mean = (MCQ_MT548 + MCQ_MT552)/2, 
		MCS_mean = (MCS_MT548 + MCS_MT552)/2, 
		MCW_mean = (MCW_MT548 + MCW_MT552)/2)

ion.values <- d6n[, 1]
d6n2 <- as.matrix(d6n[, 14:17])
samples <- c("King", "Queen", "Soldier", "Worker")
dimnames(d6n2) <- list(ion.values, samples)
d7 <- t(scale(t(d6n2), center=F))		# Scaling (SD=1, Mean is not zero)

library(gplots)
mycol <- colorpanel(128, low="white", high="red2")	# Manual color setting 
dev.new(width=5, height=10)
res.heat <- heatmap.2(d7, Colv=F, Rowv=T, dendrogram="row", col=mycol, key=T, keysize=0.01, key.par=list(cex=0.5), density.info="none", trace="none", main="Lipids", lwid=c(3, 7), lhei=c(1, 11), cexRow=0.1, cexCol=0.8)


# Find caste-biased ions in midgut contents
row.den <- as.hclust(res.heat$rowDendrogram)
km <- cutree(row.den, h=2.0)
d9 <- as.data.frame(t(res.heat$carpet))
d10 <- d9[order(rownames(d9)), ]
d10 <- cbind(d6, d10, clust.no=km)
rownames(d10) <- c()
# 3: no difference among castes, 4: King-biaed, 1: Soldier-biased, 2: Queen-biased, 5: Worker-biased
d10 <- cbind(d10, Bias.in.MC="")
for(i in 1:nrow(d10)){
	if(d10$clust.no[i] == 4){d10$Bias.in.MC[i] <- "King-biased"} else {
		if(d10$clust.no[i] == 1){d10$Bias.in.MC[i] <- "Soldier-biased"} else {
			if(d10$clust.no[i] == 2){d10$Bias.in.MC[i] <- "Queen-biased"} else{
				if(d10$clust.no[i] == 5){d10$Bias.in.MC[i] <- "Worker-biased"} else {}
			}
		}
	}
}

write.csv(d10, "Lipids MS2 heatmap clustering result v5.csv", row.names=F)

