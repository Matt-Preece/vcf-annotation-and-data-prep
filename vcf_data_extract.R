setwd("/data/home/qp241207/vcf_all/")

libraries <- c("annotatr","org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg38.knownGene","dplyr","vcfR")

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

samples <- read.delim("names.txt", header = F)[,1]

annots <- "hg38_basicgenes"
annotations_gene <- build_annotations(genome = "hg38", annotations = annots)

for(samp in samples) {

	tmp  <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
	fixed_df <- getFIX(tmp)
	fixed_df <- as.data.frame(fixed_df)
	fixed_df <- fixed_df[,c("CHROM","POS","REF","ALT")]
	colnames(fixed_df) <- c("seqnames","start","ref","alt")
	fixed_df$end <- fixed_df[,"start"]
	fixed_df <- fixed_df[,c("seqnames","start","end","ref","alt")]
	
	fixed_grange <- as(fixed_df,Class = "GRanges")
	
	fixed_annot <- annotate_regions(fixed_grange,
									annotations = annotations_gene,
									ignore.strand = T,
									quiet = F)
	fixed_annot <- as.data.frame(fixed_annot)
	fixed_annot <- fixed_annot[,c("seqnames","start","end","ref","alt","annot.symbol","annot.gene_id","annot.type")]
	fixed_annot <- unique(fixed_annot)
	fixed_annot <- fixed_annot[!is.na(fixed_annot$annot.symbol),]
									
	write.csv(fixed_annot, file = paste("./annotated/snv/", samp, "snv_annotated.csv", sep = ""), row.names = F)
	
	tmp  <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
	fixed_df <- getFIX(tmp)
	fixed_df <- as.data.frame(fixed_df)
	fixed_df <- fixed_df[,c("CHROM","POS","REF","ALT")]
	colnames(fixed_df) <- c("seqnames","start","ref","alt")
	fixed_df$end <- fixed_df[,"start"]
	fixed_df <- fixed_df[,c("seqnames","start","end","ref","alt")]
	
	fixed_grange <- as(fixed_df,Class = "GRanges")
	
	fixed_annot <- annotate_regions(fixed_grange,
									annotations = annotations_gene,
									ignore.strand = T,
									quiet = F)
	fixed_annot <- as.data.frame(fixed_annot)
	fixed_annot <- fixed_annot[,c("seqnames","start","end","ref","alt","annot.symbol","annot.gene_id","annot.type")]
	fixed_annot <- unique(fixed_annot)
	fixed_annot <- fixed_annot[!is.na(fixed_annot$annot.symbol),]
									
	write.csv(fixed_annot, file = paste("./annotated/indel/", samp, "indel_annotated.csv", sep = ""), row.names = F)
}
