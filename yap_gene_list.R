load("~/Projects/JQ1/exp_190/deseq2.sub.rda")
yap <- res

yap_low_wt_low <- yap[,c("gene_id","gene_name","gene_type","YAP.KO.LOW.WT.LOW.mean",
                         "YAP.KO.LOW.WT.LOW.lfc","YAP.KO.LOW.WT.LOW.pv","YAP.KO.LOW.WT.LOW.fdr")]

write.table(yap_low_wt_low,"~/Projects/JQ1/yap_low_wt_low.txt",quote = F,row.names = F,sep='\t')
