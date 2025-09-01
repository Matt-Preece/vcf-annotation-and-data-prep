setwd("/data/home/qp241207/vcf_all/annotated/")

libraries <- c("dplyr","ReactomePA","enrichplot","clusterProfiler","ggplot2","org.Hs.eg.db")

for (lib in libraries) {
  if (require(package = lib, character.only = TRUE)) {
    successful <- "Successful"
  } else {
    installing <- "Installing"
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkgs = lib, suppressUpdates = T)
    library(lib, character.only = TRUE )
  }
}

samples <-read.delim("names.txt", header = F)
samples <- sort(samples[,1])

groups <- substr(samples,1,1)
groups <- unique(groups)

for(group in groups) {

keep <- grep(group,samples)
tmp_samp <- samples[keep]
tmp_snv <- data.frame()
  
  for(samp in tmp_samp) {
    tmp <- read.csv(paste0("./snv/",samp,"snv_annotated.csv"),header = T)
    tmp_snv <- rbind(tmp_snv,tmp)
  }

tmp_snv <- tmp_snv[,-8]
tmp_snv <- unique(tmp_snv)

gene_id_snv <- tmp_snv[,"annot.gene_id"]
gene_PA <- enrichPathway(as.character(gene_id_snv), readable = T, qvalueCutoff = 0.05)
gene_PA_df <- as.data.frame(gene_PA)

	if(nrow(gene_PA_df) > 0) {
		gene_PA_df <- gene_PA_df[,c("Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/Pathway/snv/group_",group,"_snv_path_enrich.csv"))
		png(file = paste0("./enrichment/Pathway/snv/group_",group,"_snv_PA.png"), width = 1000, height = 1000)
		p1 <- dotplot(gene_PA, showCategory=15) + ggtitle(paste("Reactome Pathway Enrichment Analysis for Group",group,"SNVs"))
		print(p1)
		dev.off()
	} else {
		print(paste("No significant pathway enrichment for group", group, "SNVs"))
	}
                         
gene_GO <- enrichGO(as.character(gene_id_snv), OrgDb = 'org.Hs.eg.db', ont = 'ALL', readable = TRUE, qvalueCutoff = 0.05)
gene_GO_df <- as.data.frame(gene_GO)

	if(nrow(gene_GO_df) > 0) {
		gene_GO_df <- gene_GO_df[,c("ONTOLOGY","Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_GO_df, file = paste0("./enrichment/GO/snv/group_",group,"_snv_ontology_enrich.csv"))
		png(file = paste0("./enrichment/GO/snv/group_",group,"_snv_GO.png"), width = 1500, height = 1000)
		p2 <- dotplot(gene_GO, showCategory=15) + ggtitle(paste("GO Enrichment Analysis for Group",group,"SNVs"))
		print(p2)
		dev.off()
	} else {
		print(paste("No significant GO enrichment for Group", group, "SNVs"))
	}

tmp_indel <- data.frame()
  
  for(samp in tmp_samp) {
    tmp <- read.csv(paste0("./indel/",samp,"indel_annotated.csv"),header = T)
    tmp_indel <- rbind(tmp_indel,tmp)
  }
  
tmp_indel <- tmp_indel[,-8]
tmp_indel <- unique(tmp_indel)
  
gene_id_indel <- tmp_indel[,"annot.gene_id"]
gene_PA <- enrichPathway(as.character(gene_id_indel), readable = T, qvalueCutoff = 0.05)
gene_PA_df <- as.data.frame(gene_PA)

	if(nrow(gene_PA_df) > 0) {
		gene_PA_df <- gene_PA_df[,c("Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/Pathway/indel/group_",group,"_indel_path_enrich.csv"))
		png(file = paste0("./enrichment/Pathway/indel/group_",group,"_indel_PA.png"), width = 1000, height = 1000)
		p1 <- dotplot(gene_PA, showCategory=15) + ggtitle(paste("Reactome Pathway Enrichment Analysis for Group",group,"Indels"))
		print(p1)
		dev.off()
	} else {
		print(paste("No significant pathway enrichment for group", group, "Indels"))
	}
                         
gene_GO <- enrichGO(as.character(gene_id_indel), OrgDb = 'org.Hs.eg.db', ont = 'ALL', readable = TRUE, qvalueCutoff = 0.05)
gene_GO_df <- as.data.frame(gene_GO)

	if(nrow(gene_GO_df) > 0) {
		gene_GO_df <- gene_GO_df[,c("ONTOLOGY","Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_GO_df, file = paste0("./enrichment/GO/indel/group_",group,"_indel_ontology_enrich.csv"))
		png(file = paste0("./enrichment/GO/indel/group_",group,"_indel_GO.png"), width = 1500, height = 1000)
		p2 <- dotplot(gene_GO, showCategory=15) + ggtitle(paste("GO Enrichment Analysis for Group",group,"Indels"))
		print(p2)
		dev.off()
	} else {
		print(paste("No significant GO enrichment for Group", group, "Indels"))
	}
  
}
