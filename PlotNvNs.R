#Plots
setwd('~/Desktop/NvsN/');
library(ggplot2)
lines10 <- data.frame(y = 1:5000)
lines10$x = 10/lines10$y
lines8 <- data.frame(y = 1:5000)
lines8$x = 8/lines8$y
#SNV
All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo <- read.delim("AllSample_MergedResults_SNPs.txt")
data.exonic <- All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo[All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo$Exonic == "exonic-nonsynonymous", ]
# Raw filter for somaticness - Ratio >= 5
data.exonic <- data.exonic[which(as.numeric(as.character(data.exonic[,'Ratio'])) >= 5),];
allVarPlot <-ggplot(data.exonic, aes(x=as.numeric(VF_T), y=as.numeric(DP_T),color=VarType))+ geom_point(size = 3) + geom_line(data=lines10, aes(x = x, y = y),color="blue") + xlim(0, 1) + geom_vline(xintercept = 0.05,color="green") + geom_hline(yintercept = 20,color="red")
allVarPlot <- allVarPlot + facet_grid(Panel~Hotspot) + theme_bw(base_size = 20, base_family = "Helvetica") +xlab("VF") + ylab("DP") + ggtitle("Normal vs Normal SNP Analysis")
ggsave(file="SNP_AD10.png",width=20, height=10,dpi=400)
allVarPlot <-ggplot(data.exonic, aes(x=as.numeric(VF_T), y=as.numeric(DP_T),color=VarType))+ geom_point(size = 3) + geom_line(data=lines8, aes(x = x, y = y),color="blue") + xlim(0, 1) + geom_vline(xintercept = 0.02,color="green") + geom_hline(yintercept = 20,color="red")
allVarPlot <- allVarPlot + facet_grid(Panel~Hotspot) + theme_bw(base_size = 20, base_family = "Helvetica") +xlab("VF") + ylab("DP") + ggtitle("Normal vs Normal SNP Analysis")
ggsave(file="SNP_AD8.png",width=20, height=10,dpi=400)

#indel
indel.data <- read.delim("AllSample_MergedResults_Indels.txt")
# Raw filter for somaticness - Ratio >= 5
indel.data <- indel.data[which(as.numeric(as.character(indel.data[,'Ratio'])) >= 5),];
# Raw filter for artifacts - Occurance in other normals <= 20%
indel.data <- indel.data[which(as.numeric(as.character(indel.data[,'OccuranceInNormals'])) <= 4),];
allVarPlot <-ggplot(indel.data, aes(x=as.numeric(VF_T), y=as.numeric(DP_T),color=VarType))+ geom_point(size = 3) + geom_line(data=lines10, aes(x = x, y = y),color="blue") + xlim(0, 1) + geom_vline(xintercept = 0.05,color="green") + geom_hline(yintercept = 20,color="red")
allVarPlot <- allVarPlot + facet_grid(Panel~Hotspot) + theme_bw(base_size = 20, base_family = "Helvetica") +xlab("VF") + ylab("DP") + ggtitle("Normal vs Normal Indel Analysis")
ggsave(file="Indel_AD10.png",width=20, height=10,dpi=400)
allVarPlot <-ggplot(indel.data, aes(x=as.numeric(VF_T), y=as.numeric(DP_T),color=VarType))+ geom_point(size = 3) + geom_line(data=lines8, aes(x = x, y = y),color="blue") + xlim(0, 1) + geom_vline(xintercept = 0.02,color="green") + geom_hline(yintercept = 20,color="red")
allVarPlot <- allVarPlot + facet_grid(Panel~Hotspot) + theme_bw(base_size = 20, base_family = "Helvetica") +xlab("VF") + ylab("DP") + ggtitle("Normal vs Normal Indel Analysis")
ggsave(file="Indel_AD8.png",width=20, height=10,dpi=400)




