## salmon results_v2
setwd('~/Projects/JQ1/')

## deseq2
library(tximport)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggvenn)

raw_count <- read.table("raw_count.txt",sep="\t",header=T)
id_gene <- raw_count[,c(1,2)]

row.names(raw_count) <- raw_count[,1]
raw_count <- raw_count[,-c(1,2)]
sample_data <- data.frame(row.names = colnames(raw_count),
                          condition = factor(c("control","control","control",
                                        "JQ1","JQ1","JQ1",
                                        "control","control","control",
                                        "JQ1","JQ1","JQ1"),
                                        levels = c("control","JQ1")),
                          genotype = factor(c(rep("WT",6),rep("NF2null",6)),
                                            levels=c("WT","NF2null")))
sample_data$name <- paste(sample_data$condition,
                          sample_data$genotype,
                          c("1","2","3",
                            "1","2","3",
                            "1","2","3",
                            "1","2","3"),sep="_")

# nf2null only
nf2null_only <- raw_count[,c(7:12)]
nf2null_sample <- data.frame(row.names = colnames(nf2null_only),
                             condition = factor(c("control","control","control",
                                                  "JQ1","JQ1","JQ1"),
                                                levels = c("control","JQ1")))
# WT only
wt_only <- raw_count[,c(1:6)]
wt_sample <- data.frame(row.names = colnames(wt_only),
                             condition = factor(c("control","control","control",
                                                  "JQ1","JQ1","JQ1"),
                                                levels = c("control","JQ1")))

#================================================
# 1. How does JQ1 inhibition impact NF2-null SCs?
#================================================
# model with nf2null cell only (6 samples)
dds_nf2null <- DESeqDataSetFromMatrix(countData = nf2null_only,
                              colData = nf2null_sample,
                              design = ~ condition)
dds_nf2null <- DESeq(dds_nf2null)

resultsNames(dds_nf2null)
res_1 <- results(dds_nf2null, contrast = c("condition","JQ1","control"))

# lfc shrinking
lfcs_NF2null <- lfcShrink(dds_nf2null,coef="condition_JQ1_vs_control", res=res_1)

#================================================
# 1.2 How does JQ1 inhibition impact wt SCs?
#================================================
# model with nf2null cell only (6 samples)
dds_wt <- DESeqDataSetFromMatrix(countData = wt_only,
                                      colData = wt_sample,
                                      design = ~ condition)
dds_wt <- DESeq(dds_wt)

resultsNames(dds_wt)
res_wt <- results(dds_wt, contrast = c("condition","JQ1","control"))

# lfc shrinking
lfcs_wt <- lfcShrink(dds_wt,coef="condition_JQ1_vs_control", res=res_wt)

#================================================
# 2. Explore differential response of hSCs +/- NF2 to JQ1.
#================================================
# model with all 12 samples
# the condition (JQ1 treatment) effect for genotype II (NF2null)
# the main effect+ the interaction term
# (the extra condition effect in genotype II compared to genotype I)
dds <- DESeqDataSetFromMatrix(countData = raw_count,
                              colData = sample_data,
                              design = ~ genotype+condition+genotype:condition)
dds <- DESeq(dds)

# the condition (JQ1 treatment) effect for genotype I (WT) (the main effect)
# res_WT <- results_v2(dds, contrast=c("condition","JQ1","control"))
# lfcs_WT<-lfcShrink(dds, coef="condition_JQ1_vs_control",res = res_WT)

# the interaction term, answering: is the condition effect "different" across genotypes?
res_interaction <- results(dds,name = "genotypeNF2null.conditionJQ1")
lfcs_interaction<-lfcShrink(dds, coef = "genotypeNF2null.conditionJQ1",res=res_interaction)

#================================================
# 3. Gene expression changes in hSCs +/- NF2.
#================================================
# the genotype (NF2null) effect for condition I (control) (the main effect)
res_NF2null_WT_control <- results(dds,contrast=c("genotype", "NF2null","WT"))
lfcs_NF2null_WT_control <- lfcShrink(dds,coef = "genotype_NF2null_vs_WT",res=res_NF2null_WT_control)

# the genotype (NF2null) effect for condition II (JQ1)
# the main effect+ the interaction term
# (the extra genotype effect in condition II compared to condition I)
res_NF2null_WT_JQ1 <- results(dds,list( c("genotype_NF2null_vs_WT","genotypeNF2null.conditionJQ1") ))

dds_temp <- DESeqDataSetFromMatrix(countData = raw_count,
                                   colData = sample_data,
                                   design = ~ condition+genotype:condition)
dds_temp <- DESeq(dds_temp)
results_v2Names(dds_temp)
res_test <- results(dds_temp, name=c("conditionJQ1.genotypeNF2null"))

# lfc shrinking
lfcs_NF2null_WT_JQ1 <- lfcShrink(dds_temp,coef="conditionJQ1.genotypeNF2null",res=res_test)

#================================================
# pca & correlation
#================================================
library(corpcor)

vsd_nf2null <- counts(dds_nf2null, normalized=TRUE)
vsd_nf2null <- assay(rlog(dds_nf2null, blind=FALSE))
vsd <- assay(rlog(dds, blind=FALSE))
svd_tbl <- wt.scale(t(vsd),center=TRUE,scale=TRUE)
svd_run <- fast.svd(svd_tbl)
ds <- round(svd_run$d^2/sum(svd_run$d^2)*100,2)

# PCA
pdf("results_v2/PCA.pdf", width = 8,height = 6)
par(mfrow=c(1,2),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
plot(svd_run$u[,1],svd_run$u[,2],col=as.integer(as.factor(sample_data[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC2: ',ds[2],'%'),main='')
text(svd_run$u[,1],svd_run$u[,2],sample_data[,2],cex=0.6,pos=1)
legend('topright',levels(as.factor(sample_data[,1])),col=as.integer(as.factor(levels(as.factor(sample_data[,1]))))+1,pch=20,bty='n',cex=0.6)
plot(svd_run$u[,1],svd_run$u[,3],col=as.integer(as.factor(sample_data[,1]))+1,pch=20,
     xlab=paste0('PC1: ',ds[1],'%'),ylab=paste0('PC3: ',ds[3],'%'),main='')
text(svd_run$u[,1],svd_run$u[,3],sample_data[,2],cex=0.6,pos=1)
dev.off()

# correlation
pdf("results_v2/correlation.pdf", width = 8,height = 8)
par(mfrow=c(4,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:4){
    for(j in 1:3){
        x <- (i-1)*3+j
        y <- (i-1)*3+j+1
        if(y > i*3) y <- y-3
        cat(x,y,'\n')
        #plot(vsd[,x],vsd[,y],pch=20,col=rgb(0,0,0,alpha=0.1),xlab=sample_data$name[x],ylab=sample_data$name[y])
        plot(log2(assay(dds)[,x]+1),log2(assay(dds)[,y]+1),pch=20,col=rgb(0,0,0,alpha=0.1),
             xlab=sample_data$name[x],
             ylab=sample_data$name[y])
        abline(a=0,b=1)
    }
}
dev.off()

pdf("results_v2/hist.pdf", width = 8,height = 8)
par(mfrow=c(4,3),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:12){
    hist(vsd[,i],n=100,xlab=sample_data$name[i],main='')
}
dev.off()

#================================================
# make result
#================================================
# lfcs list
lfcs_list <- list(#NF2null = lfcs_NF2null,
                  WT=lfcs_wt)
                  #interaction=lfcs_interaction,
                  #NF2null_WT_control=lfcs_NF2null_WT_control,
                  #NF2null_WT_JQ1=lfcs_NF2null_WT_JQ1)

# DEseq2 normalized count
dat_nf2null <- counts(dds_wt, normalized=TRUE)
dat <- counts(dds_wt, normalized=TRUE)

dat <- cbind(dat,dat_nf2null)

pdf("results_v2/MA_wt.pdf", width = 4,height = 4)
#par(mfrow=c(2,2),cex.axis=1.2,cex.lab=1.2,font.lab=1.2)
for(i in 1:length(lfcs_list)){
    DESeq2::plotMA(lfcs_list[[i]], ylim=c(-5,5),main=names(lfcs_list)[i])
    tmp <- data.frame(lfcs_list[[i]][,c(1,2,4,5)],stringsAsFactors=F)
    colnames(tmp) <- paste0(names(lfcs_list)[i],c('.mean','.lfc','.pv','.fdr'))
    #dat <- tmp[match(rownames(dat),rownames(tmp)),]
    dat <- cbind(dat,tmp[match(rownames(dat),rownames(tmp)),])
}
dev.off()

# add gene names
dat$ID <- rownames(dat)
dat <- merge(id_gene,dat,by="ID")

write.table(dat,file='results_v2/wt_JQ1_result_normalized.txt',sep='\t',row.names=F,quote=F)

# idx <- genes$gene_type[match(rownames(dat),genes$gene_id)]=='protein_coding' ## only protein coding genes
# plot(log2(dat[idx,4]+1),log2(dat[idx,5]+1),pch=20,col=rgb(0,0,0,alpha=0.1))

#------------------------------------------------------------------------------
# volcano plots
#------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

dat <- read.table("results_v2/wt_JQ1_result_normalized.txt",sep='\t',header=T)


# one condition
vol_df <- dat[,c(2,15,16,18)]
colnames(vol_df) <- c("Symbol","mean","lfc","fdr")

# remove NA
vol_df <- vol_df[which(!is.na(vol_df$fdr) | !is.na(vol_df$fdr)),]
vol_df$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
vol_df$diffexpressed[vol_df$lfc > 0.8 & vol_df$fdr < 0.05 ] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
vol_df$diffexpressed[vol_df$lfc < -0.8 & vol_df$fdr < 0.05] <- "DOWN"

# get significant (lfc > 0.8, fdr < 0.05, mean > 50)
vol_df_sig <- vol_df[which(vol_df$fdr < 0.05 & abs(vol_df$lfc) > 0.8 & vol_df$mean > 50),]
# get top 50
vol_df_top <- head(vol_df_sig[order(-abs(vol_df_sig$lfc)),],50)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


volcano <- ggplot(data=vol_df,aes(x=lfc,y=-log10(fdr),
                       col = diffexpressed, label=Symbol))+
    geom_point(size=1)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size = 12))+
    geom_text_repel(data = vol_df_top,
                    max.overlaps = 50,size = 3)+
    geom_vline(xintercept=c(-0.8, 0.8), col="red",linetype='dashed') +
    geom_hline(yintercept=-log10(0.05), col="red",linetype='dashed')+
    scale_color_manual(values = mycolors)+
    ggtitle("JQ1 +/- in WT")+
    ylab("-log10(FDR)")+
    xlab("log2 fold change")

pdf(file = "results_v2/wt_volcano.pdf",width = 8, height = 7)
volcano
dev.off()


#------------------------------------------------------------------------------
# heatmap
#------------------------------------------------------------------------------
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendextend)

# one condition
temp_dat <- dat[,c(1:8,15,16,18)]
colnames(temp_dat)[3:11] <- c("X393_1","X393_2","X393_3","X393_4",
                              "X393_5","X393_6","mean","lfc","fdr")
# temp_dat <- dat[,c(1:2,
#                    6:8,
#                    12:14,
#                    33,34,36)]
#colnames(temp_dat)[9:11] <- c("mean","lfc","fdr")
sig_gene <- temp_dat[which(temp_dat$fdr < 0.05 &
                               abs(temp_dat$lfc)>0.8 &
                               temp_dat$mean> 50),]

# get top 50
sig_gene_top <- head(sig_gene[order(-abs(sig_gene$lfc)),],50)

# make heatmap matrix
sig_gene_plot_df <- sig_gene_top[,c(2:8)]
#sig_gene_plot_df <- sig_gene_top[,c(2:14)]
row.names(sig_gene_plot_df) <- sig_gene_plot_df[,1]

# make sample groups
sample_group <- data.frame(Conditions=rep(c("HSC_2L_control",
                                            "HSC_2L_JQ1"),
                                            #"HSC_NF2null_control",
                                            #"HSC_NF2null_JQ1"),
                                          c(3,3)))
row.names(sample_group) <- colnames(sig_gene_plot_df[,-1])
sample_group$Conditions <- factor(sample_group$Conditions,levels = c("HSC_2L_control",
                                                                     "HSC_2L_JQ1"))
                                                                     #"HSC_NF2null_control",
                                                                     #"HSC_NF2null_JQ1"))
group_color <- brewer.pal(4,"Dark2")
annotation_color <- list(Conditions = c(HSC_2L_control=group_color[1],
                                    HSC_2L_JQ1=group_color[2]))
                                    #HSC_NF2null_control=group_color[3],
                                    #HSC_NF2null_JQ1=group_color[4]))

heatmap<- pheatmap(sig_gene_plot_df[,-1],
         #color = rev(brewer.pal(9,"RdBu")),
         border_color = "black",
         show_rownames = T,
         #cluster_cols = F,
         #cellwidth = 15,
         #cellheight = 1,
         labels_row = sig_gene_plot_df$Symbol,
         annotation_colors = annotation_color,
         #gaps_col = 3,
         annotation_col = sample_group,
         scale="row",
         cellwidth=30,
         main = "JQ1 +/- in WT (Top 50)",
         treeheight_col = 0
         #fontsize_row = 0.4,
         #annotation_row = gene_col
         )

pdf(file = "results_v2/WT_heatmap.pdf",width = 8, height = 8)
heatmap
dev.off()

# vinn diagram
main_result <- read.table("results_v2/main_result_normalized.txt",sep="\t",header=T)

control_top_50 <- read.table("results_v2/NF2null_WT_control_top50",sep="\t",header=F)
jq_top_50 <- read.table("results_v2/NF2null_WT_JQ1_top50",sep="\t",header=F)

p1_plot_df <- list(Control = control_top_50$V1,
                   JQ1 = jq_top_50$V1)

p1 <- ggvenn(p1_plot_df)+
    ggtitle("NF2null vs. WT")

pdf(file = "results_v2/JQ1_control_vinn.pdf",width = 4, height = 4)
p1
dev.off()
