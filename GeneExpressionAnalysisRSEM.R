rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Catalyst Project. Gene Expression Analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: July, 2016
############################################################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")


#Work Directory
setwd("/Users/Pinedasans/Catalyst/Data/GeneExpression")

#################################
### The RSEM normalized data ###
################################
df.gene.expr <- read.table("/Users/Pinedasans/Catalyst/Data/GeneExpression/gtex_RSEM_Hugo_norm_count",header=T)

##Preparing the data to work in R
##The columns are stored as factors. We need to convert to numeric by columns
df.gene.expr.2 <- df.gene.expr[2:nrow(df.gene.expr),2:ncol(df.gene.expr)]
xx <- apply (df.gene.expr.2,2,as.character)
gene.expr <- apply(xx,2,as.numeric)
colnames.gene<-apply(df.gene.expr[1,2:ncol(df.gene.expr)],2,as.character)
rownames(gene.expr) <- df.gene.expr[2:nrow(df.gene.expr),1]
colnames(gene.expr) <- colnames.gene ##This is the gene expression matrix where columns are samples and rows are genes

##To annotate the samples
annotate.sample <- read.table("TissueAnnotation.txt",sep="\t",header=T)
#Convert to kidney==1, the others==0
annotate.sample$kidney <- ifelse(annotate.sample$SMTS!="Kidney","0","1")
annotate.sample$vessels <- ifelse(annotate.sample$SMTS!="Blood Vessel","0","1")
##To find the annotated individuals in the annotation file
id.sample <- match(colnames(gene.expr),annotate.sample$SAMPID)
id.na <- colnames(gene.expr)[which(is.na(id.sample)==TRUE)]
id.sample.na <- lapply(1:length(id.na), function(x) grep(id.na[x],annotate.sample$SAMPID))       

id.na <- which(is.na(id.sample)==TRUE)
for (i in 1:length(id.na)){
  if (length(id.sample.na[[i]])>1) {
    id.sample[id.na[i]] <- id.sample.na[[i]][1]
  } else {
    id.sample[id.na[i]] <- id.sample.na[[i]]
  }
} 
table(annotate.sample[id.sample,3]) ##Kidney 28, other 7824
my.annotate.sample <- annotate.sample[id.sample,]


###To annotate the genes (only for protein_coding genes) #CH38 version #19797 #To find this I have used the script GeneExprAnnotation.sh
annotate.gene <- read.table("coding_genes.Ch38.txt",sep="\t")
colnames(annotate.gene) <- c("chr" , "start" , "end" , "gene_id" , "gene_type", "gene_name")
id.gene <- match(rownames(gene.expr),annotate.gene$gene_name) ### 19,725 protein-coding genes

my.annotate.gene <- annotate.gene[na.omit(id.gene),] #To consider only the protein-coding genes
id.gene <- match(my.annotate.gene$gene_name,rownames(gene.expr))
gene.expr.coding <- gene.expr[id.gene,] 

save(gene.expr.coding,my.annotate.sample,my.annotate.gene,file="GeneExprRSEM.Rdata")


###############
##### Differential expression analysis using three approaches
#################
load("GeneExprRSEM.Rdata")

##The genes that are non-zero for all the samples
gene.expr.nonzero <- gene.expr.coding[which(rowSums(gene.expr.coding)!=0),] ## 19,649 genes that have non zero values

boxplot(t(gene.expr.coding)~my.annotate.sample$SMTS,col=rainbow(32))




#The genes that are non-zero in kidney samples
gene.expr.kidney <- gene.expr.nonzero[which(rowSums(gene.expr.nonzero[,which(my.annotate.sample[,3]==1)])!=0),] #19,183 genes that are some expression in the kidney
#The genes that are non-zero in blood vessels
gene.expr.vessel <- gene.expr.nonzero[which(rowSums(gene.expr.nonzero[,which(my.annotate.sample[,4]==1)])!=0),] #19,595 genes

#############
### 1.kidney differential expressed genes kidney vs. all
###############
##Apply a non-parametric test (Wilcoxon rank test) to find differences considering kidney vs. all
p.value <- sapply(1:nrow(gene.expr.kidney), function(x) wilcox.test(gene.expr.kidney[x,]~ my.annotate.sample[ ,3])$p.value)
##To find the mean amd sd by gene within each gorup
mean.kidney <- apply(gene.expr.kidney[,which(my.annotate.sample$kidney==1)],1,mean)
sd.kidney <- apply(gene.expr.kidney[,which(my.annotate.sample$kidney==1)],1,sd)
mean.others <- apply(gene.expr.kidney[,which(my.annotate.sample$kidney==0)],1,mean)
sd.others <- apply(gene.expr.kidney[,which(my.annotate.sample$kidney==0)],1,sd)

write.table(cbind(p.value,mean.kidney,mean.others,sd.kidney,sd.others),file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEKidney.txt")

##DE genes
results.DE.genes<-read.table(file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEKidney.txt")
results.DE.genes$p.value.BY <- p.adjust(results.DE.genes$p.value,method="BY") #7,405
results.DE.genes.order <- results.DE.genes[order(results.DE.genes$p.value.BY),]
results.order.DE.up <- results.DE.genes.order[which(as.numeric(results.DE.genes.order[,2])-as.numeric(results.DE.genes.order[,3])>"1"),] 
results.order.DE.sign.up<-results.order.DE.up[which(results.order.DE.up$p.value.BY<0.01),]##1,817 genes upregulated and DE
results.order.DE.sign.up[which(results.order.DE.sign.up[,2]-results.order.DE.sign.up[,3]>1),] #1817

write.table(results.order.DE.sign.up,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.sign.upregulated.kidney.txt",row.names = T,col.names = T)


#############
### 1. Blood vessels differential expressed genes blood vessels vs. all
###############
##Apply a non-parametric test (Wilcoxon rank test) to find differences considering kidney vs. all
p.value <- sapply(1:nrow(gene.expr.vessel), function(x) wilcox.test(gene.expr.vessel[x,]~ my.annotate.sample[ ,4])$p.value)
##To find the mean amd sd by gene within each gorup
mean.vessel <- apply(gene.expr.vessel[,which(my.annotate.sample$vessels==1)],1,mean)
sd.vessel <- apply(gene.expr.vessel[,which(my.annotate.sample$vessels==1)],1,sd)
mean.others <- apply(gene.expr.vessel[,which(my.annotate.sample$vessels==0)],1,mean)
sd.others <- apply(gene.expr.vessel[,which(my.annotate.sample$vessels==0)],1,sd)

write.table(cbind(p.value,mean.vessel,mean.others,sd.vessel,sd.others),file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEvessels.txt")

##DE genes
results.DE.genes<-read.table(file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEvessels.txt")
results.DE.genes$p.value.BY <- p.adjust(results.DE.genes$p.value,method="BY") 
results.DE.genes.order <- results.DE.genes[order(results.DE.genes$p.value.BY),]

results.order.DE.up <- results.DE.genes.order[which(as.numeric(results.DE.genes.order[,2])-as.numeric(results.DE.genes.order[,3])>"1"),] 
results.order.DE.sign.up<-results.order.DE.up[which(results.order.DE.up$p.value.BY<0.01),]##1,689
results.order.DE.sign.up[which(results.order.DE.sign.up[,2]-results.order.DE.sign.up[,3]>1),] ##1,688

write.table(results.order.DE.sign.up,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.sign.upregulated.vessels.txt",row.names = T,col.names = T)



#############
### 3. To find kidney differential expressed genes kidney vs. each tissue
###############
####

my.annotate.sample$SAMPID<-as.character(my.annotate.sample$SAMPID)
splitpop <- strsplit(my.annotate.sample$SAMPID,"_")
my.annotate.sample$SAMPID<-unlist(lapply(splitpop, "[", 1))

samples<-seq(1,32,1)
samples<-samples[-c(1,4,10,13,15)]
gene.nonzero<-NULL
for (i in samples){
  print(i)
  my.annotate.sample.filter<-my.annotate.sample[which(as.numeric(my.annotate.sample$SMTS)==i | my.annotate.sample$SMTS=="Kidney"),]
  id.filter <- match(my.annotate.sample.filter$SAMPID,colnames(gene.expr.kidney))
  gene.expr.specific<-gene.expr.kidney[,na.omit(id.filter)]
  assign(paste("p.value.",i,sep=""),sapply(1:nrow(gene.expr.specific), function(x) wilcox.test(gene.expr.specific[x,]~ my.annotate.sample.filter[ ,3])$p.value))
  assign(paste("mean.others",i,sep=""),sapply(1:nrow(gene.expr.specific), function(x) mean(gene.expr.specific[x,my.annotate.sample.filter[,3]==0])))
}

mean.others.total <- cbind(mean.others2,mean.others3,mean.others5,mean.others6,mean.others7,mean.others8,mean.others9,mean.others11,mean.others12
                           ,mean.others14,mean.others16,mean.others17,mean.others18,mean.others19,mean.others20,mean.others21,mean.others22,mean.others23
                           ,mean.others24,mean.others25,mean.others26,mean.others27,mean.others28,mean.others29,mean.others30,mean.others31,mean.others32)
rownames(mean.others.total)<-rownames(gene.expr.kidney)
colnames(mean.others.total)<-names(table(my.annotate.sample$SMTS)[samples])
write.table(mean.others.total,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/mean.others.kidney.specific.txt")

for (i in samples){
  assign(paste("p.value.BY.",i,sep=""),p.adjust(get(paste("p.value.",i,sep=""))))
} 

p.value.total<-cbind(p.value.BY.2,p.value.BY.3,p.value.BY.5,p.value.BY.6,p.value.BY.7,p.value.BY.8,p.value.BY.9,p.value.BY.11,p.value.BY.12
                     ,p.value.BY.14,p.value.BY.16,p.value.BY.17,p.value.BY.18,p.value.BY.19,p.value.BY.20,p.value.BY.21,p.value.BY.22,p.value.BY.23
                     ,p.value.BY.24,p.value.BY.25,p.value.BY.26,p.value.BY.27,p.value.BY.28,p.value.BY.29,p.value.BY.30,p.value.BY.31,p.value.BY.32)

rownames(p.value.total)<-rownames(gene.expr.kidney)
colnames(p.value.total)<-names(table(my.annotate.sample$SMTS)[samples])
write.table(p.value.total,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/p.value.DE.kidney.specific.txt")


######
p.value.total<-read.table("/Users/Pinedasans/Catalyst/Results/GeneExpr/p.value.DE.kidney.specific.txt")

###DE genes kidney vs. Each
tissue.sign.gene<-apply(p.value.total,1,function(x)length(which(x<0.01)))
table(tissue.sign.gene)
genes.kidney.27<-names(tissue.sign.gene[which(tissue.sign.gene==27)])
genes.kidney.26<-names(tissue.sign.gene[which(tissue.sign.gene==26)])
genes.kidney.each<-c(genes.kidney.26,genes.kidney.27)
write.table(genes.kidney.each,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/genes.kidney.each.txt")

###DE genes kidney vs. each up-regulated
mean.others.up <- -(mean.others.total - mean.kidney)
tissue.up<-apply(mean.others.up,1,function(x)length(which(x>0)))
genes.up.kidney.27<-names(tissue.up[which(tissue.up==27)])
id.merge<-match(genes.kidney.each,genes.up.kidney.27)
genes.up.kidney.each[na.omit(id.merge)] ##101 genes up-regulates DE for each other tissue
write.table(genes.up.kidney.each[na.omit(id.merge)] ,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/genes.kidney.each.txt")


#############
### 3. To find blood vessels differential expressed genes kidney vs. each tissue
###############

my.annotate.sample$SAMPID<-as.character(my.annotate.sample$SAMPID)
splitpop <- strsplit(my.annotate.sample$SAMPID,"_")
my.annotate.sample$SAMPID<-unlist(lapply(splitpop, "[", 1))

samples<-seq(1,32,1)
samples<-samples[-c(1,4,6,10,13)]
for (i in samples){
  print(i)
  my.annotate.sample.filter<-my.annotate.sample[which(as.numeric(my.annotate.sample$SMTS)==i | my.annotate.sample$SMTS=="Blood Vessel"),]
  id.filter <- match(my.annotate.sample.filter$SAMPID,colnames(gene.expr.vessel))
  gene.expr.specific<-gene.expr.vessel[,na.omit(id.filter)]
  #assign(paste("p.value.",i,sep=""),sapply(1:nrow(gene.expr.specific), function(x) wilcox.test(gene.expr.specific[x,]~ my.annotate.sample.filter[ ,4])$p.value))
  assign(paste("mean.others",i,sep=""),sapply(1:nrow(gene.expr.specific), function(x) mean(gene.expr.specific[x,my.annotate.sample.filter[,4]==0])))
}

mean.others.total <- cbind(mean.others2,mean.others3,mean.others5,mean.others7,mean.others8,mean.others9,mean.others11,mean.others12
                           ,mean.others14,mean.others15,mean.others16,mean.others17,mean.others18,mean.others19,mean.others20,mean.others21,mean.others22,mean.others23
                           ,mean.others24,mean.others25,mean.others26,mean.others27,mean.others28,mean.others29,mean.others30,mean.others31,mean.others32)

rownames(mean.others.total)<-rownames(gene.expr.vessel)
colnames(mean.others.total)<-names(table(my.annotate.sample$SMTS)[samples])
write.table(mean.others.total,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/mean.others.vessels.specific.txt")

for (i in samples){
  assign(paste("p.value.BY.",i,sep=""),p.adjust(get(paste("p.value.",i,sep=""))))
} 

paste("p.value.BY.",samples,sep="")
p.value.total<-cbind(p.value.BY.2,p.value.BY.3,p.value.BY.5,p.value.BY.7,p.value.BY.8,p.value.BY.9,p.value.BY.11,p.value.BY.12
                    ,p.value.BY.14,p.value.BY.15,p.value.BY.16,p.value.BY.17,p.value.BY.18,p.value.BY.19,p.value.BY.20,p.value.BY.21,p.value.BY.22,p.value.BY.23
                    ,p.value.BY.24,p.value.BY.25,p.value.BY.26,p.value.BY.27,p.value.BY.28,p.value.BY.29,p.value.BY.30,p.value.BY.31,p.value.BY.32)

rownames(p.value.total)<-rownames(gene.expr.vessel)
colnames(p.value.total)<-names(table(my.annotate.sample$SMTS)[samples])
write.table(p.value.total,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/p.value.DE.vessel.specific.txt")

######
p.value.total<-read.table("/Users/Pinedasans/Catalyst/Results/GeneExpr/p.value.DE.vessel.specific.txt")

tissue.sign.gene<-apply(p.value.total,1,function(x)length(which(x<0.01)))
table(tissue.sign.gene)
genes.vessels.27<-names(tissue.sign.gene[which(tissue.sign.gene==27)])
genes.vessels.26<-names(tissue.sign.gene[which(tissue.sign.gene==26)])
genes.vessels.each<-c(genes.vessels.26,genes.vessels.27)
write.table(genes.vessels.each,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/genes.vessels.each.txt")

###DE genes kidney vs. each up-regulated
mean.others.up <- -(mean.others.total - mean.vessel)
tissue.up<-apply(mean.others.up,1,function(x)length(which(x>0)))
genes.up.vessels.27<-names(tissue.up[which(tissue.up==27)])
id.merge<-match(genes.vessels.each,genes.up.vessels.27)
genes.up.vessels.27[na.omit(id.merge)] ##101 genes up-regulates DE for each other tissue
write.table(genes.up.vessels.27[na.omit(id.merge)] ,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/genes.vessels.each.txt")




#############
### 2. Highly expressed in kidney mean(gene)>mean(all_genes)+1.SD(all_genes)
###############
results.DEKidney<-read.table(file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEKidney.txt")
gene.expr.kidney.samples<-gene.expr.kidney[,which(my.annotate.sample$kidney==1)]
mean(gene.expr.kidney.samples) #7.383752
sd(gene.expr.kidney.samples) #4.095384
high.expr<-NULL
for (i in 1:nrow(gene.expr.kidney.samples)){
  high.expr[i]<-mean(gene.expr.kidney.samples[i,])>7.4+(1*4.1)
}
gene.kidney.high.expr<-results.DEKidney[which(high.expr==TRUE),] # 1,326 genes
gene.kidney.high.expr[order(gene.kidney.high.expr$mean.kidney),]
write.table(gene.kidney.high.expr,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/gene.kidney.high.expr.txt")



#############
### 2. Highly expressed in blood vessels mean(gene)>mean(all_genes)+1.SD(all_genes)
###############
results.DEVessel<-read.table(file="/Users/Pinedasans/Catalyst/Results/GeneExpr/results.DEvessels.txt")
gene.expr.vessel.samples<-gene.expr.vessel[,which(my.annotate.sample$vessels==1)]
mean(gene.expr.vessel.samples) # 7.119112
sd(gene.expr.vessel.samples) #4.285212
high.expr<-NULL
for (i in 1:nrow(gene.expr.vessel.samples)){
  high.expr[i]<-mean(gene.expr.vessel.samples[i,])>7.1+(1*4.3)
}
gene.vessel.high.expr<-results.DEVessel[which(high.expr==TRUE),] # 2,482 genes
gene.vessel.high.expr[order(gene.vessel.high.expr$mean.vessel),]
write.table(gene.vessel.high.expr,file="/Users/Pinedasans/Catalyst/Results/GeneExpr/gene.vessel.high.expr.txt")



