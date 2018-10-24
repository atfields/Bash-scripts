library('stringr')
library('data.table')

#import data
idepth<-read.table(file="out.idepth", header=TRUE)
imiss<-read.table(file="out.imiss", header=TRUE)
ldepth<-read.table(file="out.ldepth.mean", header=TRUE)
lmiss<-read.table(file="out.lmiss", header=TRUE)
frq<-read.table(file="out.frq", header=TRUE, fill=TRUE, row.names=NULL)
het<-read.table(file="out.het", header=TRUE)
colnames(het)<-c("INDV", "O.HOM.", "E.HOM.", "N_SITES", "Fis")
lqual<-read.table(file="out.lqual", header=TRUE)
rel<-read.table(file="out.relatedness", header=TRUE)
pop<-read.table(file="popmap", header=FALSE)
colnames(pop)<-c("INDV", "POP")
isam<-read.table(file="out.samples", header=TRUE, fill=TRUE)
hardy<-read.table(file="out.hwe.alter", header=TRUE)

#combine data
indv<-merge(idepth, imiss, by="INDV")
indv<-merge(indv, het, by="INDV")
indv<-merge(indv, pop, by="INDV")
indv<-merge(indv, isam, by="INDV")
#temp<-str_split_fixed(indv[,1], "\\.", 3)
#colnames(temp)<-c("Sample", "Index", "Lib")
#indv<-cbind(indv,temp)
loci<-merge(ldepth, lmiss, by.x=c("CHROM", "POS"), by.y=c("CHR", "POS") )
loci<-merge(loci, lqual, by=c("CHROM", "POS") )
loci<-merge(loci, frq, by=c("CHROM", "POS") )
loci<-merge(loci, hardy, by=c("CHROM", "POS"))
#Data Manipulation
loci$P_ADJ<-p.adjust(loci$P_HET_EXCESS, method="BH")
top<-data.frame()
tmp2<-matrix(ncol=2, nrow=dim(loci)[1])
i=0
j=0

print("fitting the line")
while(xor(xor(length(which(tmp2[,2]>0)) < 20,
      length(which(tmp2[,2]<0))/length(tmp2[,2]) < 0.99),
      sum(is.na(tmp2)) > 0)) {
  top<-loci[which(loci$MEAN_DEPTH>max(loci$MEAN_DEPTH, na.rm=T)*(0.95-(j*0.05)) & loci$QUAL>max(loci$QUAL, na.rm=T)*(0.95-(i*0.01))),]
    while(dim(top)[1]==0) {
#    print(c(i,j))
    i<-i+1
    top<-loci[which(loci$MEAN_DEPTH>max(loci$MEAN_DEPTH, na.rm=T)*(0.95-(j*0.05)) & loci$QUAL>max(loci$QUAL, na.rm=T)*(0.95-(i*0.01))),]
    top
    i<-i+1
    if(i>80){j<-j+1; i=0}
  }
#  print(c(i,j))
  top<-top[order(top$QUAL),]
  for(k in 1:dim(top)[1]) {
#    print(paste("top trial", k, sep=" "))
    tmp<-matrix(c(0,0,top$MEAN_DEPTH[k],top$QUAL[k]),byrow=T, nrow=2)
    test<-lm(tmp[,2]~tmp[,1])
    tmp2<-cbind(loci$MEAN_DEPTH*test$coefficients[2]+test$coefficients[1], loci$QUAL-loci$MEAN_DEPTH*test$coefficients[2]+test$coefficients[1])
    if((length(which(tmp2[,2]>0)) > 20) &
       (length(which(tmp2[,2]<0))/length(tmp2[,2]) > 0.99)) break
  }
  i<-i+1
  if(i>99){j<-j+1; i=0}
  if(j>20) break
}

loci<-cbind(loci, loci$MEAN_DEPTH*test$coefficients[2]+test$coefficients[1])
loci<-cbind(loci, loci$QUAL-loci[,17])
colnames(loci)[colnames(loci)=="loci$MEAN_DEPTH * test$coefficients[2] + test$coefficients[1]"] <- "DEPTH_EST"
colnames(loci)[colnames(loci)=="loci$QUAL - loci[, 17]"] <- "RESID"
out.20<-loci[which(loci$RESID>max(loci$RESID*0.2)),1:2]
out.10<-loci[which(loci$RESID>max(loci$RESID*0.1)),1:2]

rel.fil<-rel[which(rel[,1]!=rel[,2]),]

#Output bad loci
write.table(out.20, file="SNP_Qual_outlier.20.txt", quote=F, col.names=F, row.names=F)
write.table(out.10, file="SNP_Qual_outlier.10.txt", quote=F, col.names=F, row.names=F)

#Output tables
write.csv(indv, file="indv.csv", row.names=F)
write.csv(loci, file="loci.csv", row.names=F)

#plotting
##barplots
tiff(file ="barplots.Pop.tif", width = 1500, height = 720, res =200)
par(mar=c(6.5,4.1,4.1,2.1))
par(mfrow=c(1,3))
t<-table(indv$POP)
boxplot(Fis ~ POP, data=indv, ylab="Fis", names=paste(levels(indv$POP), "(n=",t, ")"), las = 2)
boxplot(F_MISS ~ POP, data=indv, ylab="% Missing per individual", names=paste(levels(indv$POP), "(n=",t, ")"), las = 2)
boxplot(MEAN_DEPTH ~ POP, data=indv, ylab="Mean Depth", names=paste(levels(indv$POP), "(n=",t, ")"), las = 2)
dev.off()

tiff(file ="barplots.Lib.tif", width = 1500, height = 720, res =200)
par(mar=c(6.5,4.1,4.1,2.1))
par(mfrow=c(1,3))
t<-table(indv$Lib)
boxplot(Fis ~ Lib, data=indv, ylab="Fis", names=paste(levels(indv$Lib), "(n=",t, ")"), las = 2)
boxplot(F_MISS ~ Lib, data=indv, ylab="% Missing per individual", names=paste(levels(indv$Lib), "(n=",t, ")"), las = 2)
boxplot(MEAN_DEPTH ~ Lib, data=indv, ylab="Mean Depth", names=paste(levels(indv$Lib), "(n=",t, ")"), las = 2)
dev.off()

tiff(file ="barplots.Index.tif", width = 1500, height = 720, res =200)
par(mar=c(6.5,4.1,4.1,2.1))
par(mfrow=c(1,3))
t<-table(indv$Index)
boxplot(Fis ~ Index, data=indv, ylab="Fis", names=paste(levels(indv$Index), "(n=",t, ")"), las = 2)
boxplot(F_MISS ~ Index, data=indv, ylab="% Missing per individual", names=paste(levels(indv$Index), "(n=",t, ")"), las = 2)
boxplot(MEAN_DEPTH ~ Index, data=indv, ylab="Mean Depth", names=paste(levels(indv$Index), "(n=",t, ")"), las = 2)
dev.off()

tiff(file="relatedness.tif", res=200, height = 720, width = 720)
hist(as.numeric(rel.fil[,3]), ylab="Relatedness", xlab="Frequency", main="")
dev.off()

##Individuals Histograms
tiff(file ="indv.hist.tif", width = 1500, height = 720, res =200)
par(mfrow=c(1,3))
hist(indv$F_MISS, xlim=c(min(floor(indv$F_MISS*10)/10, na.rm=T),max(ceiling(10*indv$F_MISS)/10, na.rm=T)), xlab="% Missing data per individual", ylab="Counts", main = NULL, breaks=nclass.FD(indv$F_MISS)*2)
hist(indv$Fis, xlim=c(min(floor(indv$Fis*10)/10, na.rm=T),max(ceiling(10*indv$Fis)/10, na.rm=T)), xlab="Fis", ylab="Counts", main = NULL, breaks=nclass.FD(indv$Fis)*2)
hist(indv$MEAN_DEPTH, xlim=c(min(floor(indv$MEAN_DEPTH*10)/10, na.rm=T),max(ceiling(10*indv$MEAN_DEPTH)/10, na.rm=T)), xlab="Mean read depth per individual", ylab="Counts", main = NULL, breaks=nclass.FD(indv$MEAN_DEPTH)*2)
dev.off()

##Loci Histograms
tiff(file ="loci.hist.tif", width = 1500, height = 1500, res =200)
par(mfrow=c(2,2))
hist(loci$F_MISS, xlim=c(min(floor(loci$F_MISS*10)/10, na.rm=T),max(ceiling(10*loci$F_MISS)/10, na.rm=T)), xlab="% Missing data per locus", ylab="Counts", main = NULL, breaks=nclass.Sturges(loci$F_MISS)*2)
hist(loci$MEAN_DEPTH, xlim=c(min(floor(loci$MEAN_DEPTH*10)/10, na.rm=T),max(ceiling(10*loci$MEAN_DEPTH)/10, na.rm=T)), xlab="Mean read depth per locus", ylab="Counts", main = NULL, breaks=nclass.Sturges(loci$MEAN_DEPTH)*2)
a<-table(loci$CHROM)
hist(a, xlab="SNPs per Contig", ylab="Counts", main = NULL, breaks=nclass.Sturges(a)*2)
a<-as.data.table(table(loci$CHROM))
b<-aggregate(loci$MEAN_DEPTH,by=list(loci$CHROM), mean)
c<-merge(a,b,by.x ="V1", by.y= "Group.1")
names(c)<-c("CHROM", "COUNT", "MEAN_DEPTH")
plot(MEAN_DEPTH ~ COUNT, data =c, xlab ="SNPs per locus", ylab = "Mean Depth" )
dev.off()

##Loci Heterozygosity
tiff(file ="loci.He.tif", width = 1500, height = 1500, res =200)
par(mfrow=c(2,2))
hist(loci$PER_HET, xlim=c(0,1), ylab="Counts", xlab="Percent Heterozygosity", main=NULL, breaks=50)
plot(PER_HET ~ MEAN_DEPTH, data=loci, ylab="Percent Heterozygosity", xlab="Mean depth per locus", ylim=c(0,1))
plot(PER_HET ~ F_MISS, data=loci, ylab="Percent Heterozygosity", xlab="Percent Missing Data", ylim=c(0,1))
plot(PER_HET ~ P_HET_EXCESS, data=loci, ylab="Percent Heterozygosity", xlab="p-value of excess", ylim=c(0,1))
dev.off()

HWE.loci<-loci[which(loci$P_ADJ<0.01),1:2]
write.table(HWE.loci, "Het_Excess.01.txt", quote=F, row.names=F, col.names=F)
HWE.loci<-loci[which(loci$P_ADJ<0.05),1:2]
write.table(HWE.loci, "Het_Excess.05.txt", quote=F, row.names=F, col.names=F)

#Scatterplots
library("RColorBrewer")
getPalette=colorRampPalette(brewer.pal(8,"Dark2"))
palette(getPalette(max(length(levels(indv$POP)),length(levels(indv$Lib)))))

tiff(file ="scatter.POP.tif", width = 1000, height = 2160, res =200)
par(mfrow=c(3,2))

plot(F_MISS~MEAN_DEPTH,xlim =c(floor(min(indv$MEAN_DEPTH, na.rm=T)/100)*100,ceiling(max(indv$MEAN_DEPTH, na.rm=T)/100)*100), data=indv, ylab="% missing data", xlab="Mean depth per individual", bty="n", pch =20, col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5), cex=0.5,pch=16)

plot(F_MISS~MEAN_DEPTH, data=loci, ylab="% missing data", xlab="Mean depth per locus", bty="n", pch =20)

plot(F_MISS~Fis, data=indv,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylab="% missing data", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=1:length(indv$POP), cex=0.5,pch=16)

plot(MEAN_DEPTH~Fis,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylim =c(0,round(max(indv$MEAN_DEPTH/100))*100), data=indv, ylab="mean depth per individual", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=1:length(indv$POP), cex=0.5,pch=16)

plot(QUAL~MEAN_DEPTH, data=loci, ylab="SNP quality", xlab="Mean depth per locus", pch =20, bty="n")
points(loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="blue")
points(loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="red")
abline(test, col="red")
abline(a=max(loci$RESID)*0.2, b=test$coefficients[2], col="blue")
abline(a=max(loci$RESID)*0.1, b=test$coefficients[2], col="blue", lty=3)

dev.off()

tiff(file ="scatter.LIB.tif", width = 1000, height = 2160, res =200)
par(mfrow=c(3,2))
plot(F_MISS~MEAN_DEPTH, xlim =c(floor(min(indv$MEAN_DEPTH, na.rm=T)/100)*100,ceiling(max(indv$MEAN_DEPTH, na.rm=T)/100)*100), data=indv, ylab="% missing data", xlab="Mean depth per individual", bty="n", pch =20, col=adjustcolor(getPalett:
HWE.loci<-loci[which(loci$P_ADJ<0.01),1:2]
write.table(HWE.loci, "Het_Excess.01.txt", quote=F, row.names=F, col.names=F)
HWE.loci<-loci[which(loci$P_ADJ<0.05),1:2]
write.table(HWE.loci, "Het_Excess.05.txt", quote=F, row.names=F, col.names=F)

#Scatterplots
library("RColorBrewer")
getPalette=colorRampPalette(brewer.pal(8,"Dark2"))
palette(getPalette(max(length(levels(indv$POP)),length(levels(indv$Lib)))))

tiff(file ="scatter.POP.tif", width = 1000, height = 2160, res =200)
par(mfrow=c(3,2))

plot(F_MISS~MEAN_DEPTH,xlim =c(floor(min(indv$MEAN_DEPTH, na.rm=T)/100)*100,ceiling(max(indv$MEAN_DEPTH, na.rm=T)/100)*100), data=indv, ylab="% missing data", xlab="Mean depth per individual", bty="n", pch =20, col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5), cex=0.5,pch=16)

plot(F_MISS~MEAN_DEPTH, data=loci, ylab="% missing data", xlab="Mean depth per locus", bty="n", pch =20)

plot(F_MISS~Fis, data=indv,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylab="% missing data", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=1:length(indv$POP), cex=0.5,pch=16)

plot(MEAN_DEPTH~Fis,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylim =c(0,round(max(indv$MEAN_DEPTH/100))*100), data=indv, ylab="mean depth per individual", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$POP))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$POP),col=1:length(indv$POP), cex=0.5,pch=16)

plot(QUAL~MEAN_DEPTH, data=loci, ylab="SNP quality", xlab="Mean depth per locus", pch =20, bty="n")
points(loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="blue")
points(loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="red")
abline(test, col="red")
abline(a=max(loci$RESID)*0.2, b=test$coefficients[2], col="blue")
abline(a=max(loci$RESID)*0.1, b=test$coefficients[2], col="blue", lty=3)

dev.off()

tiff(file ="scatter.LIB.tif", width = 1000, height = 2160, res =200)
par(mfrow=c(3,2))
plot(F_MISS~MEAN_DEPTH, xlim =c(floor(min(indv$MEAN_DEPTH, na.rm=T)/100)*100,ceiling(max(indv$MEAN_DEPTH, na.rm=T)/100)*100), data=indv, ylab="% missing data", xlab="Mean depth per individual", bty="n", pch =20, col=adjustcolor(getPalette(length(levels(indv$Lib))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$Lib),col=1:length(indv$Lib), cex=0.5,pch=16)

plot(F_MISS~MEAN_DEPTH, data=loci, ylab="% missing data", xlab="Mean depth per locus", bty="n", pch =20)

plot(F_MISS~Fis, data=indv,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylab="% missing data", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$Lib))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$Lib),col=1:length(indv$Lib), cex=0.5,pch=16)

plot(MEAN_DEPTH~Fis,xlim =c(floor(min(indv$Fis, na.rm=T)),ceiling(max(indv$Fis, na.rm=T))), ylim =c(0,round(max(indv$MEAN_DEPTH/100))*100), data=indv, ylab="mean depth per individual", xlab="Fis per individual", pch =20, bty="n", col=adjustcolor(getPalette(length(levels(indv$Lib))),alpha=0.5))
legend("topright", inset=0.025, legend=levels(indv$Lib),col=1:length(indv$Lib), cex=0.5, pch=16)

plot(QUAL~MEAN_DEPTH, data=loci, ylab="SNP quality", xlab="Mean depth per locus", pch =20, bty="n")
points(loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.1+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="blue")
points(loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="MEAN_DEPTH")],loci[which(loci$RESID>(max(loci$RESID)*0.2+test$coefficients[2]*loci$MEAN_DEPTH)),which(names(loci)=="QUAL")], pch=16, col="red")
abline(test, col="red")
abline(a=max(loci$RESID)*0.2, b=test$coefficients[2], col="blue")
abline(a=max(loci$RESID)*0.1, b=test$coefficients[2], col="blue", lty=3)

dev.off()
