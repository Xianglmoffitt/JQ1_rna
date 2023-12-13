# plot scripts

setwd("~/Projects/JQ1/new_PTK_pathway/GSEA/")

library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)
library(data.table)

gsea_pos <- read.table("JQ1_NF2null/JQ1_WT.GseaPreranked.1690380992232/gsea_report_for_na_pos_1690380992232.tsv",header=T,sep = "\t")
gsea_neg <- read.table("JQ1_NF2null/JQ1_WT.GseaPreranked.1690380992232/gsea_report_for_na_neg_1690380992232.tsv",header=T,sep = "\t")

pos_temp <- gsea_pos[,c(1,8)]
pos_temp$score <- -log10(pos_temp$FDR.q.val)
pos_temp$group <- "Up"

neg_temp <- gsea_neg[,c(1,8)]
neg_temp$score <- log10(neg_temp$FDR.q.val)
neg_temp$group <- "Down"

barplot_df <- rbind(pos_temp,neg_temp)

barplot_df$score[which(barplot_df$score=="-Inf")] <- -3

bar_plot <- ggplot(barplot_df,aes(x = reorder(NAME, score), y = score,fill=group)) +
  # set overall appearance of the plot
  theme_classic() +
  geom_bar(stat="identity") +
  # Set main and axis titles
  ggtitle("JQ1_NF2null") +
  xlab("PTK2 Related Pathway") +
  ylab("Negative Log10 FDR") +
  # Add a line showing the alpha = 0.05 level
  geom_hline(yintercept = -log10(0.05), size = 0.5, color = "black",linetype="dashed") +
  geom_hline(yintercept = log10(0.05), size = 0.5, color = "black",linetype="dashed") +
  # Flip the x and y axes
  coord_flip()+
  scale_fill_manual(values=c("blue","red"))+
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6))

pdf(file = "../report/JQ1_NF2null_PTK2_pathway.pdf",width = 5, height = 4)
bar_plot
dev.off()

