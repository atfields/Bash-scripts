library('stringr')
library('tidyr')

print("Loading data")
opt<-read.table(file="opt.012", head=F)
opt.names<-read.table(file="opt.012.indv", head=F)
opt.loci<-read.table(file="opt.012.pos", head=F)
samples<-read.table(file="out.samples",head=T, fill=T)

#dim(opt)
#dim(opt.names)
#dim(opt.loci)
opt<-cbind(opt.names,opt)
rownames(opt)<-opt[,1]
opt<-opt[,3:dim(opt)[2]]
opt<-as.data.frame(cbind(paste(opt.loci[,1],opt.loci[,2],sep="_"),t(opt)))
rownames(opt)<-opt[,1]
opt<-as.data.frame(t(opt[,2:dim(opt)[2]]))
opt[opt==-1]<-NA
#opt[1:7,1:7]

tmp<-samples[grep("cat",samples[,1], invert=T),]
tmp2<-table(tmp[,2])
dup.names<-names(tmp2[tmp2==2])
out<-vector()
geno.error<-data.frame(row.names = c("1","2","3","H1","H2","Het"))
geno.tmp<-list()
print("Comparing duplicate samples")
ptm <- proc.time()
for (i in dup.names){
        print(paste("Processing ", i))
        test<-opt[grep(i,rownames(opt), value=T),]
        for (j in 1:dim(test)[2]){
                if(is.na(test[1,j])| is.na(test[2,j])){out[j]="NA"} else
                if(test[1,j]==0 & test[2,j]==0){out[j]="H1"} else
                if(test[1,j]==1 & test[2,j]==1){out[j]="Het"} else
                if(test[1,j]==2 & test[2,j]==2){out[j]="H2"} else
                if(test[1,j]==1 & test[2,j]!=1){out[j]="1"} else
                if(test[1,j]!=1 & test[2,j]==1){out[j]="2"}
                else {out[j]="3"}
}
geno.tmp[[which(dup.names==i)]]<-out
tmp<-t(t(table(out[out!="NA"])/length(out[out!="NA"])))
geno.error<-cbind(geno.error,tmp[match(rownames(geno.error),rownames(tmp))])
}
colnames(geno.error)<-dup.names
print("Time to process duplicates (mins)")
print(round((proc.time() - ptm)/60,digits=3))

geno.error_tab<-as.data.frame(t(geno.error))*100
geno.error_tab$Good<-round(apply(data.matrix(geno.error_tab[,4:6]),1,function(x) sum(x,na.rm=T)),digits=3)
geno.error_tab$Bad<-round(apply(data.matrix(geno.error_tab[,1:3]),1,function(x) sum(x,na.rm=T)), digits=3)

print("Percent of Genotyping")
print(geno.error_tab)

geno.table<-matrix(unlist(geno.tmp), ncol=dim(opt)[2],byrow=T)
geno.table<-as.data.frame(cbind(dup.names,geno.table))
colnames(geno.table)<-c("Sample",colnames(opt))

geno.td<-geno.table %>% gather(POS,CODE,2:dim(geno.table)[2])
tmp<-table(geno.td[geno.td$CODE==3,1],geno.td[geno.td$CODE==3,2])
print("Number of positions with different homozygous calls in a sample")
print(table(colSums(tmp)))

tmp<-table(geno.td[geno.td$CODE==1|geno.td$CODE==2,1],geno.td[geno.td$CODE==1|geno.td$CODE==2,2])
print("Number of positions with a homozygote and a heterozygote")
print(table(colSums(tmp)))

write.table(geno.error_tab, file="geno.error.txt", row.names=T, quote=F)
write.table(geno.table, file="geno.table.txt", row.names=F, quote=F)
