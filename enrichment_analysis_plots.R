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

for(samp in samples) {

tmp <- read.csv(paste0("./snv/", samp, "snv_annotated.csv"), header = T)
gene_id_snv <- tmp[,"annot.gene_id"]
gene_id_snv <- unique(gene_id_snv)
gene_snv <- tmp[,"annot.symbol"]
gene_snv <- unique(gene_snv)
mir_snv <- length(grep("MIR", gene_snv))
linc_snv <- length(grep("LINC|LOC", gene_snv))
as_snv <- length(grep("-AS|-DT", gene_snv))
coding_snv <- length(gene_snv) - mir_snv - linc_snv - as_snv
total_snv <- length(gene_snv)

title <- gsub("_|r","",samp)

gene_PA <- enrichPathway(as.character(gene_id_snv), readable = TRUE, qvalueCutoff = 0.05)
gene_PA_df <- as.data.frame(gene_PA)

	if(nrow(gene_PA_df) > 0) {
		gene_PA_df <- gene_PA_df[,c("Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/Pathway/snv/",title,"_snv_path_enrich.csv"))
		png(file = paste0("./enrichment/Pathway/snv/",title,"_snv_PA.png"), width = 1000, height = 1000)
		p1 <- dotplot(gene_PA, showCategory=15) + ggtitle(paste("Reactome Pathway Enrichment Analysis for",title,"SNVs"))
		print(p1)
		dev.off()
	} else {
		print(paste("No significant pathway enrichment for", title, "SNVs"))
	}

gene_GO <- enrichGO(as.character(gene_id_snv), OrgDb = 'org.Hs.eg.db', ont = 'ALL', readable = TRUE, qvalueCutoff = 0.05)
gene_GO_df <- as.data.frame(gene_GO)

	if(nrow(gene_GO_df) > 0) {
		gene_GO_df <- gene_GO_df[,c("ONTOLOGY","Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/GO/snv/",title,"_snv_ontology_enrich.csv"))
	
		png(file = paste0("./enrichment/GO/snv/",title,"_snv_GO.png"), width = 1500, height = 1000)
		p2 <- dotplot(gene_GO, showCategory=15) + ggtitle(paste("GO Enrichment Analysis for",title,"SNVs"))
		print(p2)
		dev.off()
	} else {
		print(paste("No significant GO enrichment for", title, "SNVs"))
	}

tmp <- read.csv(paste0("./indel/", samp, "indel_annotated.csv"), header = T)
gene_id_indel <- tmp[,"annot.gene_id"]
gene_id_indel <- unique(gene_id_indel)
gene_indel <- tmp[,"annot.symbol"]
gene_indel <- unique(gene_indel)
mir_indel <- length(grep("MIR", gene_indel))
linc_indel <- length(grep("LINC|LOC", gene_indel))
as_indel <- length(grep("-AS|-DT", gene_indel))
coding_indel <- length(gene_indel) - mir_indel - linc_indel - as_indel
total_indel <- length(gene_indel)

title <- gsub("_|r","",samp)

gene_PA <- enrichPathway(as.character(gene_id_indel), readable = TRUE, qvalueCutoff = 0.05)
gene_PA_df <- as.data.frame(gene_PA)

	if(nrow(gene_PA_df) > 0) {
		gene_PA_df <- gene_PA_df[,c("Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/Pathway/indel/",title,"_indel_path_enrich.csv"))
		
		png(file = paste0("./enrichment/Pathway/indel/",title,"_indel_PA.png"), width = 1000, height = 1000)
		p1 <- dotplot(gene_PA, showCategory=15) + ggtitle(paste("Reactome Pathway Enrichment Analysis for",title,"Indels"))
		print(p1)
		dev.off()
		
	} else {
		print(paste("No significant pathway enrichment for", title, "Indels"))
	}

gene_GO <- enrichGO(as.character(gene_id_indel), OrgDb = 'org.Hs.eg.db', ont = 'ALL', readable = TRUE, qvalueCutoff = 0.05)
gene_GO_df <- as.data.frame(gene_GO)

		if(nrow(gene_GO_df) > 0) {
		gene_GO_df <- gene_GO_df[,c("ONTOLOGY","Description","GeneRatio","BgRatio","qvalue","geneID")]
		write.csv(gene_PA_df, file = paste0("./enrichment/GO/indel/",title,"_indel_ontology_enrich.csv"))
		
		png(file = paste0("./enrichment/GO/indel/",title,"_indel_GO.png"), width = 1500, height = 1000)
		p2 <- dotplot(gene_GO, showCategory=15) + ggtitle(paste("GO Enrichment Analysis for",title,"Indels"))
		print(p2)
		dev.off()
		
	} else {
		print(paste("No significant GO enrichment for", title, "Indels"))
	}
}
