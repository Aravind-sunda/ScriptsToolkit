library(dplyr)
library(edgeR)
library(gplots)
require(DESeq)
library(ggplot2)

Genewisecounts<- read.delim("Genecode_genes_pdxrc1.txt",sep = "\t")
dim(Genewisecounts)

Genewisecounts<-Genewisecounts %>% distinct(Genes, .keep_all = TRUE)
dim(Genewisecounts)

rownames(Genewisecounts) <- Genewisecounts[,1]
head(Genewisecounts)

substring(colnames(Genewisecounts),1,6)
colnames(Genewisecounts)<-substring(colnames(Genewisecounts),1,6)
head(Genewisecounts)

Genewisecounts[,1] <- NULL
head(Genewisecounts)

colData <- read.csv("pheno_data.csv",row.names = 1)
ct_group <- colData$condition



y<- DGEList(Genewisecounts,group=ct_group, genes = Genewisecounts[,1,drop=F])
options(digits = 3)
y$samples
head(y$counts,3)
keep<-rowSums(cpm(y)>0.5)>=2
table(keep)
y<-y[keep,,keep.lib.sizes=F]
y<- calcNormFactors(y)


y1 <- estimateCommonDisp(y, verbose=T)
y1 <- estimateTagwiseDisp(y1)
names(y1)
plotBCV(y1)
design<- model.matrix( ~0+y$samples$group)
colnames(design)
colnames(design)<-levels(y$samples$group)
design
y2 <- estimateGLMCommonDisp(y, design )
y2 <- estimateGLMTrendedDisp(y2, design, method = "power")
y2 <- estimateGLMTagwiseDisp(y2, design)
plotBCV(y2)

et12 <- exactTest(y1,pair = c(1,2))
topTags(et12, n=10)
hist(et12$table[,"PValue"],breaks = 50, main="Histogram- p-values for exact Test", xlab="p-values")
res_exact <- topTags(et12, n=nrow(et12$table))
head(res_exact$table)
head(res_exact$comparison)
dim(res_exact$table)
write.csv(res_exact$table,"DE_exacttest_Genecode_6.csv")
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12 <- rownames(y1)[as.logical(de1)]
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
fit <- glmFit(y2, design)
lrt12 <- glmLRT(fit, contrast=c(-1,1))
topTags(lrt12, n=10)
hist(lrt12$table[,"PValue"],breaks = 50, main="Histogram- p-values for GLM Test", xlab="p-values")
res_glm <- topTags(lrt12, n=nrow(lrt12$table))
dim(res_glm)
write.csv(res_glm$table,"DE_GLM_test_Genecode_6.csv")

de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
summary(de2)
design
de2tags12 <- rownames(y2)[as.logical(de2)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")

volcano_exact <- cbind(res_exact$table$logFC, -log10(res_exact$table$FDR))
colnames(volcano_exact) <- c("logFC","negLogPval")
head(volcano_exact)
plot(volcano_exact, pch=19)
#add ggplot

#ggplot(volcano_exact, aes(x=volcano_exact$logFC, y=volcano_exact$negLogPval))+geom_point()

logy<- cpm(y, log=T, prior.count = 1)
head(logy)
#sely<- logy[rownames(res_exact$table)[res_exact$table$FDR<0.05 & abs(res_exact$table$logFC)>5],]
sely<- logy[rownames(res_exact$table)[res_exact$table$PValue<0.05 & abs(res_exact$table$logFC)>1],]

log_sely <- t(scale(t(sely)))
head(log_sely)
dim(log_sely)


col.pan <- colorpanel(100,"navy","darkgrey","yellow")
heatmap.2(log_sely, col= col.pan,Rowv = T, scale= "none", trace = "none", dendrogram = "column", labRow = F, cexCol = 1, srtCol = 45)
pca <- prcomp(t(logy))

plot(pca$x[,1],pca$x[,2],pch=19,xlab="PC1",ylab="PC2",ylim = c(-150,150))
text(pca$x[,1],pca$x[,2], labels = colnames(logy))
#pos 4= left, 2=right 1 =down etc
text(pca$x[,1],pca$x[,2], labels = colnames(logy),cex=0.6, pos=1)
summary(pca)
dev.off()



sely<- logy[rownames(res_exact$table)[res_exact$table$PValue<0.05 & abs(res_exact$table$logFC)>1.5],]
dim(sely)
log_sely <- t(scale(t(sely)))
head(log_sely)
col.pan <- colorpanel(100,"navy","darkgrey","yellow")
heatmap.2(log_sely, col= col.pan,Rowv = T, scale= "none", trace = "none", dendrogram = "column", labRow = F, cexCol = 1, srtCol = 45)



glm_y <- logy[rownames(res_glm$table)[res_glm$table$FDR<0.05 & abs(res_glm$table$logFC)>5],]
log_glm <- t(scale(t(glm_y)))
col.pan <- colorpanel(100,"darkgreen","lightblue","darkred")
heatmap.2(log_glm, col= col.pan,Rowv = T, scale= "none", trace = "none", dendrogram = "column", labRow = F, cexCol = 1, srtCol = 45)


results <- as.data.frame(res_glm$table)
head(results)
results <- results[,-1]
dim(results)

library(ggplot2)
results$threshold = as.factor(abs(results$logFC)>2 & results$FDR < 0.05)
ggplot(results, aes(x=logFC, y= -log10(FDR), colour=threshold))+geom_point()+theme_bw()+labs(title = "Volcano Plot based on GLMtest\nCre vs tumors")
results$name <-row.names(results)

ggplot(results, aes(x=logFC, y= -log10(FDR), colour=threshold))+geom_point()+theme_bw()+labs(title = "Volcano Plot based on GLMtest\n Normal vs MpBC")+geom_text(data=subset(results,abs(results$logFC) > 10 & results$FDR <0.001),aes(x=logFC, y= -log10(FDR),label=name),check_overlap = F,nudge_x = 1, hjust=0.5,size=2)

