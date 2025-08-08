setwd("/data/home/qp241207/vcf_all/annotated/")

library(dplyr)

samples <- read.delim("names.txt", header = F)
samples <- sort(samples[,1])


all_stats_df <- data.frame(
  Feature = c(
    "Promoters SNV", "1to5kb SNV", "5UTRs SNV", "Exons SNV", "Introns SNV", "3UTRs SNV", "1",
    "Promoters Indel", "1to5kb Indel", "5UTRs Indel", "Exons Indel", "Introns Indel", "3UTRs Indel", "2",
    "Mean SNV per Gene", "Median SNV per Gene", "Max SNV per Gene", "Min SNV per Gene", "SD SNV per Gene", "3",
    "Mean Indels per Gene", "Median Indels per Gene", "Max Indels per Gene", "Min Indels per Gene", "SD Indels per Gene"
  )
)

for(samp in samples) {

tmp <- read.csv(paste0("./snv/",samp, "snv_annotated.csv"), header = T)
structure_snv <- tmp[,c("seqnames","start","end","ref","alt","annot.type")]
structure_snv <- unique(structure_snv)
remove <- duplicated(structure_snv$start)
remove <- structure_snv$start[remove]
remove <- which(structure_snv$start %in% remove)
structure_snv_unique <- structure_snv[-remove,]
structure_snv_count <- count(structure_snv_unique, annot.type)
structure_snv_count[,1] <- gsub("hg38_genes_", "",structure_snv_count[,1])
structure_snv_count[,1] <- paste0(toupper(substr(structure_snv_count[,1],1,1)), 
										  substr(structure_snv_count[,1],2, nchar(structure_snv_count[,1])))
structure_snv_count[,1] <- paste(structure_snv_count[,1], "SNV")
colnames(structure_snv_count) <- c("Feature","Count")


gene_snv <- tmp[,c("seqnames","start","end","ref","alt","annot.symbol")]
gene_snv <- unique(gene_snv)
gene_count_snv <- count(gene_snv, annot.symbol)
gene_mean_snv <- mean(gene_count_snv[,2])
gene_median_snv <- median(gene_count_snv[,2])
gene_max_snv <- max(gene_count_snv[,2])
gene_min_snv <- min(gene_count_snv[,2])
gene_sd_snv <- sd(gene_count_snv[,2])
gene_stats_snv <- data.frame(Feature = c("Mean SNV per Gene","Median SNV per Gene","Max SNV per Gene","Min SNV per Gene","SD SNV per Gene"), 
								Count = c(gene_mean_snv, gene_median_snv, gene_max_snv, gene_min_snv, gene_sd_snv))


tmp <- read.csv(paste0("./indel/", samp, "indel_annotated.csv"), header = T)
structure_indel <- tmp[,c("seqnames","start","end","ref","alt","annot.type")]
structure_indel <- unique(structure_indel)
remove <- duplicated(structure_indel$start)
remove <- structure_indel$start[remove]
remove <- which(structure_indel$start %in% remove)
structure_indel_unique <- structure_indel[-remove,]
structure_indel_count <- count(structure_indel_unique, annot.type)
structure_indel_count[,1] <- gsub("hg38_genes_", "",structure_indel_count[,1])
structure_indel_count[,1] <- paste0(toupper(substr(structure_indel_count[,1],1,1)), 
											substr(structure_indel_count[,1],2, nchar(structure_indel_count[,1])))
structure_indel_count[,1] <- paste(structure_indel_count[,1], "Indel")
colnames(structure_indel_count) <- c("Feature","Count")

gene_indel <- tmp[,c("seqnames","start","end","ref","alt","annot.symbol")]
gene_indel <- unique(gene_indel)
gene_count_indel <- count(gene_indel, annot.symbol)
gene_mean_indel <- mean(gene_count_indel[,2])
gene_median_indel <- median(gene_count_indel[,2])
gene_max_indel <- max(gene_count_indel[,2])
gene_min_indel <- min(gene_count_indel[,2])
gene_sd_indel <- sd(gene_count_indel[,2])
gene_stats_indel <- data.frame(Feature = c("Mean Indels per Gene","Median Indels per Gene","Max Indels per Gene","Min Indels per Gene","SD Indels per Gene"), 
									Count = c(gene_mean_indel, gene_median_indel, gene_max_indel, gene_min_indel, gene_sd_indel))

space_1 <- data.frame( Feature = "1", Count = " ")
space_2 <- data.frame( Feature = "2", Count = " ")
space_3 <- data.frame( Feature = "3", Count = " ")


tmp_df <- rbind(structure_snv_count, space_1,
					structure_indel_count,space_2,
					gene_stats_snv,space_3,
					gene_stats_indel)

samp <- gsub("_|r", "", samp)
colnames(tmp_df)[2] <- samp
all_stats_df <- left_join(all_stats_df, tmp_df, by = "Feature")

}

all_stats_df[is.na(all_stats_df)] <- 0
write.csv(all_stats_df, "feature_counts_all.csv", row.names = F)
