## salmon results_v2
setwd('~/Projects/JQ1/')


#------------------------------------------------------------------------------
# heatmap
#------------------------------------------------------------------------------
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendextend)

#------------------------------------------------------------------------------
# volcano plots
#------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

dat <- read.table("new_PTK_pathway/JQ1_NF2null.txt",sep='\t',header=T)

gene_list <- read.table("new_PTK_pathway/PI3K_Akt_pathway_gl",sep="\t",header=F)

dat <- dat[which(dat$Symbol %in% gene_list$V1),]

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

# target gene

vol_df_top <- vol_df_sig[which(vol_df_sig$Symbol %in% gene_list$V1),]

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

# pdf(file = "results_v2/wt_volcano.pdf",width = 8, height = 7)
# volcano
# dev.off()


#------------------------------------------------------------------------------
# heatmap
#------------------------------------------------------------------------------
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendextend)

dat <- read.table("new_PTK_pathway/JQ1_NF2null.txt",sep='\t',header=T)
gene_list <- read.table("new_PTK_pathway/PI3K_Akt_pathway_gl",sep="\t",header=F)


# one condition
temp_dat <- dat

sig_gene_top <- temp_dat[which(temp_dat$Symbol %in% gene_list$V1),]

write.table(sig_gene_top,"new_PTK_pathway/JQ1_PI3K_AKT_Pathway_NF2null.txt",sep="\t",quote = F,row.names = F)

# make heatmap matrix
sig_gene_plot_df <- sig_gene_top[,c(2:8)]
#sig_gene_plot_df <- sig_gene_top[,c(2:14)]
row.names(sig_gene_plot_df) <- sig_gene_plot_df[,1]

# remove all 0
sig_gene_plot_df <- sig_gene_plot_df[rowSums(sig_gene_plot_df[,-1])>0,]



# make sample groups
sample_group <- data.frame(Conditions=rep(c("WT_control",
                                            "WT_JQ1"),
                                          #"HSC_NF2null_control",
                                          #"HSC_NF2null_JQ1"),
                                          c(3,3)))
row.names(sample_group) <- colnames(sig_gene_plot_df[,-1])
sample_group$Conditions <- factor(sample_group$Conditions,levels = c("WT_control",
                                                                     "WT_JQ1"))
#"HSC_NF2null_control",
#"HSC_NF2null_JQ1"))
group_color <- brewer.pal(4,"Dark2")
annotation_color <- list(Conditions = c(WT_control=group_color[1],
                                        WT_JQ1=group_color[2]))
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
                   main = "JQ1 +/- in WT (Genes in PI3K-AKT Pathway)",
                   treeheight_col = 0,
                   fontsize_row = 5
                   #annotation_row = gene_col
)

pdf(file = "new_PTK_pathway/report/JQ1_WT_heatmap_2.pdf",width = 6, height = 18)
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
