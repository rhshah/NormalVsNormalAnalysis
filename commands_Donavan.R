setwd('~/Desktop/NvsN/');
library('ggplot2');
All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo <- read.delim("AllSample_MergedResults_SNPs.txt")
data.exonic <- All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo[All_SNP_WithAlleleDepth_annotated_exonicOnly_withHotspotInfo$Exonic == "exonic-nonsynonymous", ]
# Raw filter for somaticness - Ratio >= 5
data.exonic <- data.exonic[which(as.numeric(as.character(data.exonic[,'Ratio'])) >= 5),];

make.ad.curve <- function(dat,ad){
	vf.max <- floor(max(as.numeric(dat[,'VF']))*10+0.5)/10;
	x <- seq(0.02,vf.max,0.02);
	y <- ad/x;
	return(cbind('x'=x,'y'=y));
}	

makeScatterPlot <- function(dat,my.title){
	if(my.title != 'Hotspot'){
		vf.max <- floor(max(as.numeric(dat[,'VF_T']))*10+0.5)/10;
		dp.max <- max(as.numeric(dat[,'DP_T']));
	}else{
		vf.max=0.25;
		dp.max = 1000;
	}
	plot(x=dat[,'VF_T'],y=dat[,'DP_T'],pch=19,col='grey20',xlab='VF',xlim=c(0,vf.max),ylim=c(0,dp.max),ylab='DP',main="",cex=0.6);
}

dp_snv.hs <- 20;
ad_snv.hs <- 8;
vf_snv.hs <- 0.02;
dp_snv.reg <- 20;
ad_snv.reg <- 10;
vf_snv.reg <- 0.05;

data.exonic.hs <- data.exonic[which(data.exonic[,'Hotspot'] == 'Hotspot'),];
data.exonic.reg <- data.exonic[which(data.exonic[,'Hotspot'] != 'Hotspot'),];

idx.hs <- which(as.numeric(as.character(data.exonic.hs[,'DP_T'])) >= dp_snv.hs &
                as.numeric(as.character(data.exonic.hs[,'AD_T'])) >= ad_snv.hs &
                as.numeric(as.character(data.exonic.hs[,'VF_T'])) >= vf_snv.hs);
length(idx.hs);
                
idx.reg <- which(as.numeric(as.character(data.exonic.reg[,'DP_T'])) >= dp_snv.reg &
                as.numeric(as.character(data.exonic.reg[,'AD_T'])) >= ad_snv.reg &
                as.numeric(as.character(data.exonic.reg[,'VF_T'])) >= vf_snv.reg);
length(idx.reg);
                
pdf('snv_hotspot.pdf',height=6,width=6);
par(mfrow=c(1,1));
makeScatterPlot(data.exonic[which(data.exonic[,'Hotspot'] == 'Hotspot'),],'Hotspot');
abline(v=0.02,lty=2,col='red');
#abline(h=20,lty=2,col='red');
h <- cbind('x'=seq(0.01,0.25,0.01),'y'=8/seq(0.01,0.25,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='red',lty=2);
dev.off();

pdf('snv_nonhotspot.pdf',height=6,width=6);
par(mfrow=c(1,1));
makeScatterPlot(data.exonic[which(data.exonic[,'Hotspot'] != 'Hotspot'),],'NotHotspot');
abline(v=0.05,lty=2,col='red');
#abline(h=20,lty=2,col='red');
h <- cbind('x'=seq(0.01,0.25,0.01),'y'=10/seq(0.01,0.25,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='red',lty=2);
dev.off();

indel.data <- read.delim("AllSample_MergedResults_Indels.txt")
# Raw filter for somaticness - Ratio >= 5
indel.data <- indel.data[which(as.numeric(as.character(indel.data[,'Ratio'])) >= 5),];
# Raw filter for artifacts - Occurance in other normals <= 20%
indel.data <- indel.data[which(as.numeric(as.character(indel.data[,'OccuranceInNormals'])) <= 4),];

dp_indel.hs <- 20;
ad_indel.hs <- 8;
vf_indel.hs <- 0.02;
dp_indel.reg <- 20;
ad_indel.reg <- 10;
vf_indel.reg <- 0.05;

indel.data.hs <- indel.data[which(indel.data[,'Hotspot'] == 'Hotspot'),];
indel.data.reg <- indel.data[which(indel.data[,'Hotspot'] != 'Hotspot'),];

idx.hs <- which(as.numeric(as.character(indel.data.hs[,'DP_T'])) >= dp_indel.hs &
                as.numeric(as.character(indel.data.hs[,'AD_T'])) >= ad_indel.hs &
                as.numeric(as.character(indel.data.hs[,'VF_T'])) >= vf_indel.hs);
length(idx.hs);
                
idx.reg <- which(as.numeric(as.character(indel.data.reg[,'DP_T'])) >= dp_indel.reg &
                as.numeric(as.character(indel.data.reg[,'AD_T'])) >= ad_indel.reg &
                as.numeric(as.character(indel.data.reg[,'VF_T'])) >= vf_indel.reg);
length(idx.reg);

# pdata.hs - combine data.exonic.hs [SNVs] and indel.data.hs [Indels]
pdata.hs <- rbind(cbind(data.exonic.hs[,1:8],'Type'=rep('SNV',nrow(data.exonic.hs))),
                  cbind(indel.data.hs[,1:8],'Type'=rep('Indel',nrow(indel.data.hs))));

pdata.reg <- rbind(cbind(data.exonic.reg[,1:8],'Type'=rep('SNV',nrow(data.exonic.reg))),
                  cbind(indel.data.reg[,1:8],'Type'=rep('Indel',nrow(indel.data.reg))));

pdf('FPs_NormalvsNormal.pdf',height=8,width=8);
par(mfrow=c(1,1));
vf.max=0.30;
dp.max = 2000;
idx <- which(pdata.hs[,'Type'] == 'SNV');
plot(x=pdata.hs[idx,'VF_T'],y=pdata.hs[idx,'DP_T'],pch=19,col='red',xlim=c(0,vf.max),ylim=c(0,dp.max),xlab='',ylab='',main="",cex=0.6, cex.lab=1.4,cex.axis=1.3);
par(new=T);
idx <- which(pdata.hs[,'Type'] == 'Indel');
plot(x=pdata.hs[idx,'VF_T'],y=pdata.hs[idx,'DP_T'],pch=19,col='blue',xlim=c(0,vf.max),ylim=c(0,dp.max),xlab='Variant Frequency',ylab='Unique sequence coverage (X)',main="FPs in first-tier locations",cex=0.6,cex.lab=1.4,cex.axis=1.3);
abline(v=vf_snv.hs,lty=2,col='black');
abline(h=dp_snv.hs,lty=2,col='black');
h <- cbind('x'=seq(0.01,0.25,0.01),'y'=ad_snv.hs/seq(0.01,0.25,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='black',lty=2);
legend(x='topright',legend=c('SNV','Indel'),col=c('red','blue'),pch=19,bty='n');

vf.max=0.30;
dp.max = 2000;
idx <- which(pdata.reg[,'Type'] == 'SNV');
plot(x=pdata.reg[idx,'VF_T'],y=pdata.reg[idx,'DP_T'],pch=19,col='red',xlim=c(0,vf.max),ylim=c(0,dp.max),xlab='',ylab='',main="",cex=0.5, cex.lab=1.4,cex.axis=1.3);
par(new=T);
idx <- which(pdata.reg[,'Type'] == 'Indel');
plot(x=pdata.reg[idx,'VF_T'],y=pdata.reg[idx,'DP_T'],pch=19,col='blue',xlim=c(0,vf.max),ylim=c(0,dp.max),xlab='Variant Frequency',ylab='Unique sequence coverage (X)',main="FPs in second-tier locations",cex=0.5,cex.lab=1.4,cex.axis=1.3);
abline(v=vf_snv.reg,lty=2,col='black');
abline(h=dp_snv.reg,lty=2,col='black');
h <- cbind('x'=seq(0.01,0.25,0.01),'y'=ad_snv.reg/seq(0.01,0.25,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='black',lty=2);
legend(x='topright',legend=c('SNV','Indel'),col=c('red','blue'),pch=19,bty='n');
dev.off();

idx.reg <- which(as.numeric(as.character(pdata.reg[,'DP_T'])) >= dp_snv.reg &
                as.numeric(as.character(pdata.reg[,'AD_T'])) >= ad_snv.reg &
                as.numeric(as.character(pdata.reg[,'VF_T'])) >= vf_snv.reg);
length(idx.reg);


# Previous plots:


indel.data.20prcntOccurrence <- indel.data[indel.data$OccuranceInNormals <= 4, ]
indel.data.20prcntOccurrence <- indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Ratio >= 5, ]

pdf('indel_hotspot.pdf',height=6,width=6);
par(mfrow=c(1,1));
makeScatterPlot(indel.data.20prcntOccurrence[which(indel.data.20prcntOccurrence[,'Hotspot'] == 'Hotspot'),],'Hotspot');
abline(v=0.02,lty=2,col='red');
#abline(h=20,lty=2,col='red');
h <- cbind('x'=seq(0.01,0.25,0.01),'y'=8/seq(0.01,0.25,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='red',lty=2);
dev.off();

pdf('indel_nonhotspot.pdf',height=6,width=6);
par(mfrow=c(1,1));
makeScatterPlot(indel.data.20prcntOccurrence[which(indel.data.20prcntOccurrence[,'Hotspot'] != 'Hotspot'),],'NotHotspot');
abline(v=0.05,lty=2,col='red');
#abline(h=20,lty=2,col='red');
h <- cbind('x'=seq(0.01,0.4,0.01),'y'=10/seq(0.01,0.4,0.01));
lines(x=as.numeric(h[,'x']),y=as.numeric(h[,'y']),col='red',lty=2);
dev.off();


ggplot(data.exonic[data.exonic$Hotspot == "Hotspot", ], aes(x = DP_T, y = VF_T)) + geom_line(data=lines10, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 4000), breaks = c(seq(0, 4000, by = 500))) + scale_y_continuous(limits = c(0, 0.20), breaks = c(seq(0, 0.20, by = 0.01)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("SNV_OnlyExonic_Hotspots_AD10_VF5.pdf", width = 14, height = 14)
ggplot(data.exonic[data.exonic$Hotspot == "NotHotspot", ], aes(x = DP_T, y = VF_T)) + geom_line(data=lines10, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 4000), breaks = c(seq(0, 4000, by = 500))) + scale_y_continuous(limits = c(0, 0.20), breaks = c(seq(0, 0.20, by = 0.01)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("SNV_OnlyExonic_NON-Hotspots_AD10_VF5.pdf", width = 14, height = 14)
ggplot(data.exonic[data.exonic$Hotspot == "Hotspot", ], aes(x = DP_T, y = VF_T)) + geom_line(data=lines8, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 4000), breaks = c(seq(0, 4000, by = 500))) + scale_y_continuous(limits = c(0, 0.20), breaks = c(seq(0, 0.20, by = 0.01)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("SNV_OnlyExonic_Hotspots_AD8_VF2.pdf", width = 14, height = 14)
ggplot(data.exonic[data.exonic$Hotspot == "NotHotspot", ], aes(x = DP_T, y = VF_T)) + geom_line(data=lines8, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 4000), breaks = c(seq(0, 4000, by = 500))) + scale_y_continuous(limits = c(0, 0.20), breaks = c(seq(0, 0.20, by = 0.01)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("SNV_OnlyExonic_NON-Hotspots_AD8_VF2.pdf", width = 14, height = 14)

table(data.exonic$Hotspot)
dim(data.exonic)
dim(subset(data.exonic, AD_T>= 10 & VF_T >= 0.05))
dim(subset(data.exonic, AD_T>= 8 & VF_T >= 0.02))

indel.data <- read.delim("AllSample_MergedResults_Indels.txt")
indel.data.20prcntOccurrence <- indel.data[indel.data$OccuranceInNormals <= 4, ]
indel.data.20prcntOccurrence <- indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Ratio >= 5, ]



ggplot(indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Hotspot == "Hotspot", ], aes(x = DP_T, y = VF_T, color = Exonic)) + geom_line(data=lines10, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 6200), breaks = c(seq(0, 6200, by = 500))) + scale_y_continuous(limits = c(0, 0.3), breaks = c(seq(0, 0.3, by = 0.02)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("INDEL_OnlyExonic_Hotspots_AD10_VF5.pdf", width = 14, height = 14)
ggplot(indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Hotspot == "NotHotspot", ], aes(x = DP_T, y = VF_T, color = Exonic)) + geom_line(data=lines10, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 6200), breaks = c(seq(0, 6200, by = 500))) + scale_y_continuous(limits = c(0, 0.3), breaks = c(seq(0, 0.3, by = 0.02)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("INDEL_OnlyExonic_NON-Hotspots_AD10_VF5.pdf", width = 14, height = 14)

ggplot(indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Hotspot == "Hotspot", ], aes(x = DP_T, y = VF_T, color = Exonic)) + geom_line(data=lines8, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 6200), breaks = c(seq(0, 6200, by = 500))) + scale_y_continuous(limits = c(0, 0.3), breaks = c(seq(0, 0.3, by = 0.02)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("INDEL_OnlyExonic_Hotspots_AD8_VF2.pdf", width = 14, height = 14)
ggplot(indel.data.20prcntOccurrence[indel.data.20prcntOccurrence$Hotspot == "NotHotspot", ], aes(x = DP_T, y = VF_T, color = Exonic)) + geom_line(data=lines8, aes(x = DP, y = VF, group = Filter), color = "red") + ggtitle("") + scale_x_continuous(limits = c(0, 6200), breaks = c(seq(0, 6200, by = 500))) + scale_y_continuous(limits = c(0, 0.3), breaks = c(seq(0, 0.3, by = 0.02)))+ geom_point(alpha = 0.7,size = 4, color = "steelblue") + my_theme(xAngle=0, hJust=0.5)
ggsave("INDEL_OnlyExonic_NON-Hotspots_AD8_VF2.pdf", width = 14, height = 14)

table(indel.data.20prcntOccurrence$Hotspot)
dim(indel.data.20prcntOccurrence)
dim(subset(indel.data.20prcntOccurrence, AD_T >= 10 & VF_T >= 0.05))
dim(subset(indel.data.20prcntOccurrence, AD_T >= 8 & VF_T >= 0.02))