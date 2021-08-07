# data pre-process & different analysis
library(DESeq2)
library(stats)
library(pheatmap)
datanew <- read.table(file.choose(), header=TRUE)
data1 <- data[,-16]
rownames(data1) <- data1$X
data1 <- data1[,-1]

sampleNames <- c("C55","C56","C5","S511","S513","S514","S517",
                 "S519","S521","S523","S524","S527","S528",
                 "S52","S54","S57","S59")
colnames(data1) <- sampleNames
write.csv(data1,file = "none_s531.csv",quote = TRUE)
######
rownames(data) <- data[,1]
data <- data[,-1]
countData <- data
#
data <- data[rowSums(data)>2,]
head(data)
sampleNames
database <- data.frame(name=sampleNames, condition=c(rep("Non.COVID19",3),rep("COVID19",14)))
rownames(database) <- sampleNames

## 
dds <- DESeqDataSetFromMatrix(data, colData=database, design= ~ condition)

## 
dds <- DESeq(dds)
#####
countdata <-counts(dds,normalized = TRUE)

normalized_counts_mad <- apply(countdata, 1, mad)

# 
write.table(countdata, file="DESeq2.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
# 
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

write.table(rlogMat, file="DESeq2.normalized.rlog.xls", quote=F, sep="\t", row.names=T, col.names=T)
##
# 
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
# 
heatmap(pearson_cor,col = allcolour)
###
annotation_col = data.frame(PatienType = factor(rep(c("Non.COVID19","COVID19"),c(3,14))))
rownames(annotation_col) = rownames(pearson_cor)
ann_colors <- list(PatienType = c(Non.COVID19="#00CCCC",COVID19="#FF9999"))
pheatmap(pearson_cor,scale = "none",
         cluster_cols = TRUE,cluster_rows = TRUE,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         cellwidth = 20, cellheight = 20,angle_col = "45",fontsize = 15,fontsize_col = 12,fontsize_row = 12)


write.table(pearson_cor, file="pearson_correlative.xls", quote=F, sep="\t", row.names=T, col.names=T)
write.csv(pearson_cor,file = "SFIG1B.csv")
#
p1 <- plotPCA(rld)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  annotate('text', label ='C55', x =-22, y =3.0,size=4, colour ='#00CCCC') +
  annotate('text', label ='C56', x =-28, y =3,size=4, colour ='#00CCCC') +
  annotate('text', label ='C5', x =-29, y =-3,size=4, colour ='#00CCCC') +
  annotate('text', label ='S511', x =0, y =-18,size=4, colour ='red') +
  annotate('text', label ='S513', x =-4, y =-12,size=4, colour ='red') +
  annotate('text', label ='S514', x =4.9, y =4.17,size=4, colour ='red') +
  annotate('text', label ='S517', x =13.6, y =-10.7,size=4, colour ='red') +
  annotate('text', label ='S519', x =7.2, y =14.67,size=4, colour ='red') +
  annotate('text', label ='S521', x =3, y =-1.9,size=4, colour ='red') +
  annotate('text', label ='S523', x =8.4, y =-14.55,size=4, colour ='red') +
  annotate('text', label ='S524', x =5.2, y =10.3,size=4, colour ='red') +
  annotate('text', label ='S527', x =4, y =19.8,size=4, colour ='red') +
  annotate('text', label ='S528', x =9, y =-7,size=4, colour ='red') +
  annotate('text', label ='S52', x =2.5, y =-10.5,size=4, colour ='red') +
  annotate('text', label ='S54', x =15, y =2.7,size=4, colour ='red') +
  annotate('text', label ='S57', x =-2, y =8.5,size=4, colour ='red') +
  annotate('text', label ='S59', x =0, y =15.6,size=4, colour ='red') + 
  theme(legend.text=element_text(size=12))+theme(legend.title = element_text(size = 15))
p1
write.csv(p1$data,file = "FIG1A.csv")
###
# 
sampleA = "control"
sampleB = "treat"
contrastv <- c("condition",sampleA,sampleB)
res <- results(dds,contrast = contrastv)
res

# 
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]

if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)

# 
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)
# 
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
head(res)
# 
res <- cbind(ID=rownames(res), as.data.frame(res))
res$baseMean <- rowMeans(cbind(baseA, baseB))
# 
res$padj[is.na(res$padj)] <- 1
# 
res <- res[order(res$pvalue),]
head(res)
write.csv(res, "diff_express.csv")

library(ggplot2)
library(ggrepel)
####################################################
volcano <- res
# 
volcano$significant <-NULL
volcano$significant <-ifelse(volcano$log2FoldChange >= 1,"up","unchanged")
volcano$significant[volcano$log2FoldChange <= -1] = "down"
volcano$label[volcano$padj < 10e-15 & volcano$significant != "unchanged"] <- rownames(volcano)[volcano$padj < 10e-15 & volcano$significant != "unchanged"]
##############

###########################################################
colnames(gene.df) <- c("ID","SYMBOL")
volcano.new <- merge(volcano,gene.df,by = "ID",all = TRUE)
#
volcano.new<- volcano.new[!duplicated(volcano.new$SYMBOL),]
volcano.new$label <-NULL
volcano.new$significant <-NULL
volcano.new$Significant <-ifelse(volcano.new$log2FoldChange >= 1&volcano.new$padj<=0.05,"Non-Covid 19","Unchanged")
volcano.new$Significant[volcano.new$log2FoldChange <= -1&volcano.new$padj<=0.05] = "Covid 19"
volcano.new$label <- ""
volcano.new$label[volcano.new$padj < 5e-10 & volcano.new$Significant != "Unchanged"] <- volcano.new$SYMBOL[volcano.new$padj < 5e-10 & volcano.new$significant != "Unchanged"]

volcano.new[8864,11] <- ""
##############ggplot2
volcano.new_plot <- ggplot(volcano.new,aes(log2FoldChange,-1*log10(padj)))+ 
  geom_point(aes(color = padj),size = 1) +
  geom_jitter(aes(color = padj),alpha=0.4)+
  # 
  geom_text_repel(aes(label = label,color = padj),direction = "both",nudge_x = -5,nudge_y = 1.2,
                  ylim = c(0,7),xlim = c(-15,10)) +  
  # 
  labs(x = expression(log[2](FoldChange)), y = expression(-log[10](FDR))) +
  scale_color_gradientn(colours = heat.colors(12))+
  # log10(0.5)~1.3
  geom_hline(yintercept = 1.0, linetype = 4) +#
  # log2(2)=1,log2(0.5)=-1#
  geom_vline(xintercept = c(-1, 1), linetype = 4) + #  
  
  ylim(0, 7) +xlim(-10,10) + #
  theme_bw() + #
  annotate('text', label =' Non-Covid 19', x =8.5, y =7,size=6, colour ='orange') +
  annotate('text', label =' Covid 19', x =-9, y =7,size=6, colour ='red') +
  ###
  theme(axis.title.x = element_text(size = 15,  color = "black", 
                                    face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+
  theme(axis.title.y = element_text(size = 15,  color = "black", 
                                    face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+
  theme(legend.text = element_text(colour = 'black', angle = 0, size = 10.5, 
                                   hjust = 0, vjust = 0.5, face = 'plain'))+
  theme(legend.title=element_blank())+
  theme(legend.position = "none")+
  theme(legend.background=element_rect())+
  theme(legend.text=element_text(size=13))+theme(legend.title = element_text(size = 17))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) #
volcano.new_plot

write.csv(volcano.new, "volcan.new_plot.csv")
pdf("volcano.new.pdf",width=10,height=10)
p<-volcano.new_plot
print(p)
dev.off()
###########################################

##########
volcano.new$label <- NULL
diffs <- na.omit(volcano.new)
diffs <- diffs[order(diffs$log2FoldChange),]
head(diffs)
diffs <- diffs[!duplicated(diffs$SYMBOL),]
rownames(diffs) <- diffs$SYMBOL

DEG_up <-diffs[diffs$significant == "up" & diffs$padj < 0.05,]
DEG_down <-diffs[diffs$significant == "down" & diffs$padj < 0.05,]
############
library(AnnotationHub)	
library(org.Hs.eg.db)   
library(clusterProfiler) 
library(topGO) 
library(pathview)  
library(Rgraphviz)  
library(DOSE)
library(GO.db)
library(GSEABase)
#
DEG1 <- DEG_up$ID #
DEG.gene_symbolup = as.character(DEG1)#
#
DEG_up_entrize_id <- mapIds(x= org.Hs.eg.db,
                            keys = DEG.gene_symbolup,
                            keytype = "ENSEMBL",
                            column = "ENTREZID")

DEG2 <- DEG_down$ID #
DEG.gene_symboldown = as.character(DEG2)#
#
DEG_down_entrize_id <- mapIds(x= org.Hs.eg.db,
                              keys = DEG.gene_symboldown,
                              keytype = "ENSEMBL",
                              column = "ENTREZID")
#
DEG_up_enter_id <- na.omit(DEG_up_entrize_id)
DEG_down_enter_id <- na.omit(DEG_down_entrize_id)
write.csv(gene.df, "allgenesymbol.csv")

#####3##enrichment--analysis#####################
DEG_up_enter_id <- as.vector(DEG_up_enter_id)

#up
up.erich.go.BP = enrichGO(gene = DEG1,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          ont = "BP",  #
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = T)
write.csv(up.erich.go.BP@result,"up_go_bp.csv")

#########################
#down
down.erich.go.BP = enrichGO(gene = DEG2,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENSEMBL",
                            ont = "BP",  # 
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.05,
                            readable = T)
write.csv(down.erich.go.BP@result,"down_go_bp.csv")

#########
genes_diff <- as.data.frame(rlogMat)
genes_diff$ID <- rownames(genes_diff)
genes_diff <- merge(genes_diff,diffs,by = "ID")
genes_diff <- genes_diff[genes_diff$significant != "unchanged",]
genes_diff <- genes_diff[genes_diff$padj <= 0.05,]
genes_diff <- genes_diff[,-1]
rownames(genes_diff) = genes_diff$SYMBOL
genes_diff <- genes_diff[order(genes_diff$log2FoldChange),]
genes_diff <- genes_diff[,1:17]

write.csv(genes_diff,file = "genes_diffs_heatmap.csv")
######################################################
#######################################
genelist <-volcano.new[,c(1,5)]
genelist <- genelist[order(genelist$log2FoldChange),]
gseGO(geneList, ont = "BP", OrgDb, keyType = "ENTREZID", 
      exponent = 1, nPerm = 1000, minGSSize = 10, maxGSSize = 500, 
      pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = TRUE, 
      seed = FALSE, by = "fgsea")

##########################################
library(pheatmap)
library(showtext)
library(ggplot2)
library(gtable)
library(grid)
library(ggplot2)
library(ggpubr)
library(scatterplot3d)
library(showtext)
library(pheatmap)
showtext_auto(enable = TRUE)
show_point_shapes()
########################################pca
merge <- read.csv(file.choose())
rownames(merge) <- merge$ID
merge <- merge[,-1]

###################################
merge2 <- merge[,1:12]
all_pca <- t(merge2)
all_pca <- scale(all_pca)
all_pca <- prcomp(all_pca)
label_merge <- rep(c('liver','heart','spleen','kidney'),c(6,2,2,2))
label_gender <-c("female","female","male","male","male","male",
                 "female","male",
                 "female","male",
                 "female","male")

all_pcs <- data.frame(all_pca$x,Comments = label_merge,gender = label_gender)

###############################
merge3<- merge2[,c(-3,-4)]
all_pca <- t(merge3)
all_pca <- scale(all_pca)
all_pca <- prcomp(all_pca)
label_merge <- rep(c('liver','heart','spleen','kidney'),c(4,2,2,2))
label_gender <-c("female","female","male","male",
                 "female","male",
                 "female","male",
                 "female","male")

all_pcs <- data.frame(all_pca$x,Comments = label_merge,gender = label_gender)

#3)

percentage<-round(all_pca$sdev / sum(all_pca$sdev) *100,2)
percentage<-paste(colnames(all_pcs),"(", paste(as.character(percentage),"%",")", sep=""))

mytheme <- theme(plot.title=element_text(face="plain",family = " Arial-Black",
                                         size="24", color="black"),
                 
                 axis.title=element_text(face="plain",family = " Arial-Black",
                                         size=24, color="black"),
                 
                 axis.text=element_text(face="plain", size=18,family = " Arial-Black",
                                        color="black"),
                 
                 panel.background=element_rect(fill="white",size = 1, color="black"),
                 
                 panel.grid.major.y=element_line(color="white", linetype=1),
                 
                 panel.grid.minor.y=element_line(color="white",linetype=1),
                 #legend.title = element_text(face="plain",family = " Arial-Black",
                 #size=24, color="black"),
                 legend.title = element_blank(),
                 
                 panel.grid.minor.x=element_line(color="white", linetype=1,),
                 
                 panel.grid.major.x=element_line(color="white", linetype=1),
                 
                 axis.text.x = element_text(angle = 0,face = "plain",family = "Arial-Black",
                                            size = 18,color = "black"),
                 axis.text.y = element_text(face = "plain",family = "Arial-Black",
                                            size = 18,color = "black" ),
                 axis.ticks.x  = element_blank(),
                 legend.key.size = unit(.3,"inches"),
                 legend.position="top")+
  theme(plot.margin=unit(rep(3,4),'lines'))+
  theme(panel.background = element_rect(fill = "white", colour = "black", 
                                        size = 2)) +
  theme(legend.text=element_text(size=24,family = " Arial-Black",face = "plain"))


p1 <- ggplot(all_pcs,aes(x=PC1,y=PC2,color=Comments,shape = gender))+
  #stat_ellipse(level = 0.95, show.legend = F) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  geom_point() +
  theme_bw() +
  #labs(title="all PCA Clustering",
  #subtitle=" PC1 and PC2 principal components ",
  #caption="Source: DSP_Lung_all")+
  scale_colour_manual(values = c("liver"='red', 
                                 "kidney"='#9900FF',
                                 #'lung' = "#00ba38",
                                 'heart'= 'orange',
                                 'spleen'='blue'))+#
  scale_shape_manual(values = c(17,16)) + #
  geom_point(size = 8)+
  # scale_y_continuous(expand = c(0,0), limits=c(-5, 5))+
  # scale_x_continuous(expand = c(0,0), limits=c(-20, 10))+
  annotate('text', label ='226f', x =-102, y =30,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='306s', x =-24, y =30,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='306z', x =-28, y =65,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='312g', x =100, y =90,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='226', x =-28, y =-50,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='229', x =30, y =-70,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='226', x =-65, y =-7,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='229', x =100, y =-60,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='306s', x =13, y =-30,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  annotate('text', label ='229', x =10, y =25,size=5, colour ='#669900',family = "Arial-Black",fontface = "bold") +
  
  theme(panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="#333333"))+mytheme
p1
write.csv(p1$data,file = "FIG1D2.csv")
pdf(file = "multi-pca.pdf",width = 10,height = 10)
p1
dev.off()

################################barplot#################
data <- as.data.frame(barplot_nc)
barplot_data <- as.data.frame(barplot_data)

mytheme <- theme(plot.title=element_text(face="plain",family = " Arial-Black",
                                         size="24", color="black"),
                 
                 axis.title=element_text(face="plain",family = " Arial-Black",
                                         size=24, color="black"),
                 
                 axis.text=element_text(face="plain", size=18,family = " Arial-Black",
                                        color="black"),
                 
                 panel.background=element_rect(fill="white",size = 1, color="black"),
                 
                 panel.grid.major.y=element_line(color="grey", linetype=1),
                 
                 panel.grid.minor.y=element_line(color="grey",linetype=1),
                 #legend.title = element_text(face="plain",family = " Arial-Black",
                 #size=24, color="black"),
                 legend.title = element_blank(),
                 
                 panel.grid.minor.x=element_line(color="grey", linetype=1,),
                 
                 panel.grid.major.x=element_line(color="grey", linetype=1),
                 
                 axis.text.x = element_text(angle = 0,face = "plain",family = "Arial-Black",
                                            size = 18,color = "black"),
                 axis.text.y = element_text(face = "plain",family = "Arial-Black",
                                            size = 18,color = "black" ),
                 axis.ticks.x  = element_blank(),
                 legend.key.size = unit(.3,"inches"),
                 legend.position="right")+
  theme(plot.margin=unit(rep(3,4),'lines'))+
  theme(legend.text=element_text(size=24,family = " Arial-Black",face = "plain"))


##
data.order <- c( "PRF1","IFNGR1","GZMK","IL5","IL6","IL6R","IL7R","IL23R",
                 "CCR1","CCR5","CCR7","CXCL2","CXCL3","CXCL16","CCL5","CCL20","CCL3L1")
data.order <- factor(data.order,levels = data.order)
####
p1 <- ggplot(data = barplot_data,mapping = aes(x =Gene,y=barplot_data$Log2fc,label = organ))+
  geom_bar(stat = "identity",
           width = 0.9,aes(fill=organ),
           position = position_dodge(0.8),color = 'black')+
  theme_bw()+##set background
  scale_y_continuous(expand = c(0,0), limits=c(-10, 10))+
  labs(title="", x=expression('Gene Name'), y="log2FoldChange")+
  scale_fill_manual( values = c('Liver'='red',
                                'Kidney'='#9900FF',
                                'Lung'="#00ba38",
                                'Heart'= 'orange',
                                'Spleen'='blue')) +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(data.order))) +
  theme(panel.grid.minor =element_blank())+#
  labs(title= "Diffgene Co-expression")  +
  theme(axis.ticks = element_blank())+    ## 
  mytheme

p1
write.csv(barplot_data,file = "FIG1E.csv")
pdf(file = "DIFFGENE-CO-EXPRESS.pdf",width = 8,height = 8)
p1
dev.off()
#######################################################
library(ggplot2)
library(RColorBrewer)
library(showtext)
library(ggplot2)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
library(dplyr)
library(Cairo)
library(ggplot2)
library(gtable)
library(grid)
showtext_auto(enable = TRUE)
##############
mytheme <- theme(plot.title=element_text(face="plain",family = " Arial-Black",
                                         size="24", color="black"),
                 
                 axis.title=element_text(face="plain",family = " Arial-Black",
                                         size=24, color="black"),
                 
                 axis.text=element_text(face="plain", size=18,family = " Arial-Black",
                                        color="black"),
                 
                 panel.background=element_rect(fill="white",size = 1, color="black"),
                 
                 panel.grid.major.y=element_line(color="grey", linetype=1),
                 
                 panel.grid.minor.y=element_line(color="grey",linetype=1),
                 legend.title = element_text(face="plain",family = " Arial-Black",
                                             size=24, color="black"),
                 #legend.title = element_blank(),
                 
                 panel.grid.minor.x=element_line(color="grey", linetype=1,),
                 
                 panel.grid.major.x=element_line(color="grey", linetype=1),
                 
                 axis.text.x = element_text(angle = 0,face = "plain",family = "Arial-Black",
                                            size = 18,color = "black"),
                 axis.text.y = element_text(face = "plain",family = "Arial-Black",
                                            size = 18,color = "black" ),
                 axis.ticks.x  = element_blank(),
                 legend.key.size = unit(.3,"inches"),
                 legend.position="right")+
  theme(plot.margin=unit(rep(3,4),'lines'))+
  theme(legend.text=element_text(size=24,family = " Arial-Black",face = "plain"))

figure1f <- as.data.frame(goplot_nc)

order2 <- c('liver','heart','kidney','lung','spleen')
order2 <- factor(order2,levels = order2)
######################buble plot################################################33
p1 <-ggplot(data = figure1f,mapping = aes(x = Description, y = organ,size = GeneRatio)) +
  geom_point(aes(colour =log),show.legend = T) +
  coord_flip()+
  #scale_colour_distiller(palette = "Set2")+####set color Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  #scale_color_gradientn(colours = heat.colors(10))+
  scale_color_gradient2(low = "yellow",
                        mid = "orange",
                        high = "red")+
  #scale_size_continuous(range=c(0,25))+
  #scale_size_continuous(breaks = c(5,10,15,20,25))+
  #scale_y_continuous(breaks=seq(0,0.3,5))+
  scale_x_discrete(limits = rev(levels(order1))) +
  scale_y_discrete(limits = rev(levels(order2))) +#
  theme_bw()+##set background
  labs(title="Gene Ontology enrichment", x=" ", y=" ")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+#
  mytheme
p1  
pdf(file = "Figure1F.pdf",width = 12,height = 8)
p1
dev.off()
################################barplot################
order1 <- c(  'positive regulation of lymphocyte activation',
              'positive regulation of T cell activation', 
              'interleukin-6 production',
              'response to virus',
              'T cell activation involved in immune response',
              'positive regulation of cell killing',
              'positive regulation of interferon-gamma production',
              'positive regulation of adaptive immune response',          
              'chemokine production',
              'positive regulation of natural killer cell chemotaxis')
order1 <- factor(order1,levels = order1)


p2 <- ggplot(data = figure1f,mapping = aes(x = Description, y = log,label = organ))+
  geom_bar(stat = "identity",
           width = 0.9,aes(fill=organ),
           position = position_dodge(0.8),color = 'black')+
  theme_bw()+##set background
  scale_y_continuous(expand = c(0,0), limits=c(0, 25))+
  labs(title="", x=expression('Pathway'), y="-log10pvalue")+
  scale_fill_manual( values = c('liver'='red',
                                'kidney'='#9900FF',
                                'lung'="#00ba38",
                                'heart'= 'orange',
                                'spleen'='blue')) +
  scale_color_gradient2(low = "yellow",
                        mid = "orange",
                        high = "red")+
  coord_flip()+
  scale_x_discrete(limits = rev(levels(order1))) +
  theme(panel.grid.minor =element_blank())+
  labs(title= "Multi-pathways")  +
  theme(axis.ticks = element_blank())+    
  mytheme

p2
write.csv(figure1f,file = "FIG1F.csv")
pdf(file = "GO-ANALYSIS.pdf",width = 12,height = 8)
p2
dev.off()
########################################


####pre-process
table <- read.table(file = file.choose(),header = TRUE)
table <- as.data.frame(table[,c(1,7:14)])
colnames(table) <- c("Geneid","226_liver","226_heart","226_spleen","229_heart","229_kidney","229_spleen","306L_liver","306S_kidney")
########################
library(edgeR)
library(DESeq2)
library(stats)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(scatterplot3d)
library(showtext)
library(pheatmap)
library(AnnotationHub)	
library(org.Hs.eg.db)   
library(clusterProfiler) 
library(topGO) 
library(pathview)  
library(Rgraphviz) 
library(DOSE)
library(GO.db)
library(GSEABase)
###################
rownames(table) <- table$Geneid
table <- table[,-1]
data <- table

data <- data[rowSums(data)>2,]

###################
data_heart <- data[,c(2,4)]
data_heart <- data_heart[rowSums(data_heart)>2,]

data_spleen <- data[,c(3,6)]
data_spleen<- data_spleen[rowSums(data_spleen)>2,]

data_kidney <- data[,c(8,5)]
data_kidney<- data_kidney[rowSums(data_kidney)>2,]

group <-c("female","male")

y_heart <- DGEList(counts = data_heart,group = group)
y_spleen <- DGEList(counts = data_spleen,group = group)
y_kidney <- DGEList(counts = data_kidney,group = group)
#######normalized
y_heart <- calcNormFactors(y_heart)
y_spleen <- calcNormFactors(y_spleen)
y_kidney <- calcNormFactors(y_kidney)
###diff gene express analysis
human_bcv <- 0.4
heart_bcv <- y_heart
spleen_bcv <- y_spleen
kidney_bcv <- y_kidney
et_heart <- exactTest(heart_bcv,dispersion = human_bcv^2)
et_spleen<- exactTest(spleen_bcv,dispersion = human_bcv^2)
et_kidney <- exactTest(kidney_bcv,dispersion = human_bcv^2)

###diff gene number
gene_heart <- decideTestsDGE(et_heart, p.value = 0.05, lfc = 0)
summary(gene_heart)
gene_heart <- as.data.frame(et_heart$table)
gene_heart <- data.frame(gene_heart,rownames(gene_heart))
gene_heart <- gene_heart[order(gene_heart$logFC,decreasing = T),]
colnames(gene_heart) <- c("logFC","logCPM","PValue","ID")
#####
gene_spleen <- decideTestsDGE(et_spleen, p.value = 0.05, lfc = 0)
summary(gene_spleen)
gene_spleen <- as.data.frame(et_spleen$table)
gene_spleen <- data.frame(gene_spleen,rownames(gene_spleen))
gene_spleen <- gene_spleen[order(gene_spleen$logFC,decreasing = T),]
colnames(gene_spleen) <- c("logFC","logCPM","PValue","ID")
####
gene_kidney <- decideTestsDGE(et_kidney, p.value = 0.05, lfc = 0)
summary(gene_kidney)
gene_kidney <- as.data.frame(et_kidney$table)
gene_kidney <- data.frame(gene_kidney,rownames(gene_kidney))
gene_kidney <- gene_kidney[order(gene_kidney$logFC,decreasing = T),]
colnames(gene_kidney) <- c("logFC","logCPM","PValue","ID")

###############################################
#####
DEG.all <- rownames(gene_heart)
DEG.all = as.character(DEG.all)#
#
DEG.all<- mapIds(x= org.Hs.eg.db,
                 keys = DEG.all,
                 keytype = "ENSEMBL",
                 column = "SYMBOL")
#
DEG.all.symbol <- na.omit(DEG.all)
DEG.all.symbol <- as.data.frame(DEG.all.symbol)
DEG.all.symbol <- data.frame(DEG.all.symbol,rownames(DEG.all.symbol))
colnames(DEG.all.symbol) <- c("SYMBOL","ID")
#################################################################################3
merge_heart <- merge(gene_heart,DEG.all.symbol,by = "ID")
write.csv(merge_heart,file = "merge_heart.csv")
merge_spleen <- merge(gene_spleen,DEG.all.symbol,by = "ID")
write.csv(merge_spleen,file = "merge_spleen.csv")
merge_kidney <- merge(gene_kidney,DEG.all.symbol,by = "ID")
write.csv(merge_kidney,file = "merge_kidney.csv")
######################################################################################
heart.go <- rownames(gene_heart)[gene_heart$logFC > 0 & gene_heart$PValue < 0.05] #
spleen.go <- rownames(gene_spleen)[gene_spleen$logFC > 0 & gene_spleen$PValue < 0.05] #
kidney.go <- rownames(gene_kidney)[gene_kidney$logFC > 0 & gene_kidney$PValue < 0.05] #
liver.go <- res$ID[res$log2FoldChange > 0 & res$pvalue < 0.05] #
lung.go <- res$ID[res$log2FoldChange > 0 & res$pvalue < 0.05] #
#####
liver.go = enrichGO(gene =liver.go,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont = "BP",  # 
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = T)
dotplot(liver.go,showCategory = 10,title = "EnrichmentGO")
write.csv(liver.go@result,"liver_go.csv")
####################################################################################################################################
##########################
liver <- read.csv(file = file.choose())
liver <- as.data.frame(liver)
colnames(liver) <- c("ID","219X","229D","306S","306Z","308Y","312G")
data_liver <- data[,c(1,7)]
data_liver <- data.frame(data_liver,rownames(data_liver))
colnames(data_liver) <- c("226F","306L","ID")
merge_liver <- merge(data_liver,liver,by = "ID")

merge_liver <- merge_liver[,c(1,2,6,3,4,5,7,9,8)]
merge_liver <- merge_liver[,-9]
##
##
data_liver <- merge_liver
rownames(data_liver)<- data_liver$ID
data_liver <- data_liver[,-1]
data_liver <- data_liver[,-3]
data_liver <- data_liver[,c(-3,-4)]
sampleNames <- colnames(data_liver)
# 
data_liver <- data_liver[rowSums(data_liver)>2,]
head(data_liver)
sampleNames
database <- data.frame(name=sampleNames, condition=c(rep("female",2),rep("male",2)))
rownames(database) <- sampleNames
###################################################to pre-handle lung
lung <- read.csv(file = file.choose())
lung <- lung[,c(1,5,6,7,8,9,10,11,12,13,16,18)]
colnames(lung) <- c("ID","219","229D","304W","306L","306S","306X","306Z","307C","308C","226","218")
lung <- lung[,c(1,4,6,7,9,10,11,12,2,3,5,8)]
rownames(lung)<- lung$ID
lung <- as.data.frame(lung[,-1])
# 
data_lung <- lung[rowSums(lung)>2,]
head(data_lung)
sampleNames <- colnames(lung)
database <- data.frame(name=sampleNames, condition=c(rep("female",7),rep("male",4)))
rownames(database) <- sampleNames


## 
dds <- DESeqDataSetFromMatrix(data_liver, colData=database, design= ~ condition)

## 
dds <- DESeq(dds)
#####
countdata <-counts(dds,normalized = TRUE)

normalized_counts_mad <- apply(countdata, 1, mad)

# 
write.table(countdata, file="liver.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
# 
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

##样品层级聚类分析，判断样品的相似性和组间组内差异
# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
# 层级聚类
#增加Time,CellType分组信息
heatmap(pearson_cor)
###
annotation_col = data.frame(PatienType = factor(rep(c("female","male"),c(2,2))))
rownames(annotation_col) = rownames(pearson_cor)
ann_colors <- list(PatienType = c(female="#00CCCC",male="#FF9999"))
pheatmap(pearson_cor,scale = "none",
         cluster_cols = TRUE,cluster_rows = TRUE,
         annotation_col = annotation_col,annotation_colors = ann_colors,
         cellwidth = 20, cellheight = 20,angle_col = "45",fontsize = 15,fontsize_col = 12,fontsize_row = 12)


write.table(pearson_cor, file="pearson_correlative.xls", quote=F, sep="\t", row.names=T, col.names=T)
#
p1 <- plotPCA(rld)+mytheme
p1
###
sampleA = "female"
sampleB = "male"
contrastv <- c("condition",sampleA,sampleB)
res <- results(dds,contrast = contrastv)
res

# 
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleA]

if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)

# 
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)
# 
res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
head(res)
# 
res <- cbind(ID=rownames(res), as.data.frame(res))
res$baseMean <- rowMeans(cbind(baseA, baseB))
# 
res$padj[is.na(res$padj)] <- 1
# 
res <- res[order(res$pvalue),]
head(res)
merge_liver2 <- merge(res,DEG.all.symbol,by = "ID")
write.csv(merge_liver2, "merge_liver2.csv")
#####################
data_liver2 <- data.frame(data_liver,rownames(data_liver))
colnames(data_liver2)<- c("226F_liver","306S_liver","219x_liver","229D_liver","306Z_liver","312G_liver","ID")

data_heart2 <- data.frame(data_heart,rownames(data_heart))
colnames(data_heart2)<- c("226_heart" ,"229_heart","ID")

data_spleen2 <- data.frame(data_spleen,rownames(data_spleen))
colnames(data_spleen2)<- c("226_spleen" ,"229_spleen","ID")

data_kidney2 <- data.frame(data_kidney,rownames(data_kidney))
colnames(data_kidney2)<- c("306S_kidney" ,"229_kidney","ID")

data_lung2 <- data.frame(data_lung,rownames(data_lung))
colnames(data_lung2)<- c("304W_lung", "306S_lung" ,"306X_lung", "307C_lung" ,
                         "308C_lung" ,"226_lung" , "218_lung" , "219_lung" , "229D_lung" ,"306L_lung", "306Z_lung","ID")

merge1 <- merge(data_liver2,data_heart2,by = "ID")
merge2 <- merge(merge1,data_spleen2,by = "ID")
merge3 <- merge(merge2,data_kidney2,by = "ID")
merge4 <- merge(merge3,data_lung2,by = "ID")
write.csv(merge4,file = "merge.csv")























