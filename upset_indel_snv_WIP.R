setwd("/data/home/qp241207/vcf_all/")

libraries <- c("dplyr","vcfR","UpSetR")

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

samples <- read.delim("names.txt", header = F)
samples <- sort(samples[,1])

groups <- substr(samples,1,1)
groups <- unique(groups)

for(group in groups[1:3]) {

keep <- grep(group, samples)
samp <- samples[keep]

one <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
one <- getFIX(one)
one <- as.data.frame(one)
one <- one[,c("CHROM","POS","REF","ALT")]
one$snv_id <- paste(one$CHROM, one$POS, one$REF, one$ALT, sep = ".")
one <- one$snv_id
one <- unique(one)

two <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
two <- getFIX(two)
two <- as.data.frame(two)
two <- two[,c("CHROM","POS","REF","ALT")]
two$snv_id <- paste(two$CHROM, two$POS, two$REF, two$ALT, sep = ".")
two <- two$snv_id
two <- unique(two)
	
three <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
three <- getFIX(three)
three <- as.data.frame(three)
three <- three[,c("CHROM","POS","REF","ALT")]
three$snv_id <- paste(three$CHROM, three$POS, three$REF, three$ALT, sep = ".")
three <- three$snv_id
three <- unique(three)

four <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
four <- getFIX(four)
four <- as.data.frame(four)
four <- four[,c("CHROM","POS","REF","ALT")]
four$snv_id <- paste(four$CHROM, four$POS, four$REF, four$ALT, sep = ".")
four <- four$snv_id
four <- unique(four)
	
five <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
five <- getFIX(five)
five <- as.data.frame(five)
five <- five[,c("CHROM","POS","REF","ALT")]
five$snv_id <- paste(five$CHROM, five$POS, five$REF, five$ALT, sep = ".")
five <- five$snv_id
five <- unique(five)

tmp <- gsub("_|r","", samp)
snvs <- list(one, two , three, four, five)
list_snvs <- setNames(snvs, tmp)
					
p1 <- upset(fromList(list_snvs),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "SNV Overlaps", sets.x.label = "Total SNVs per Sample",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/snv/", group, "_snv_loc_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

one <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
one <- getFIX(one)
one <- as.data.frame(one)
one <- one[,c("CHROM","POS","REF","ALT")]
one$snv_id <- paste(one$CHROM, one$POS, one$REF, one$ALT, sep = ".")
one <- one$snv_id
one <- unique(one)

two <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
two <- getFIX(two)
two <- as.data.frame(two)
two <- two[,c("CHROM","POS","REF","ALT")]
two$snv_id <- paste(two$CHROM, two$POS, two$REF, two$ALT, sep = ".")
two <- two$snv_id
two <- unique(two)
	
three <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
three <- getFIX(three)
three <- as.data.frame(three)
three <- three[,c("CHROM","POS","REF","ALT")]
three$snv_id <- paste(three$CHROM, three$POS, three$REF, three$ALT, sep = ".")
three <- three$snv_id
three <- unique(three)

four <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
four <- getFIX(four)
four <- as.data.frame(four)
four <- four[,c("CHROM","POS","REF","ALT")]
four$snv_id <- paste(four$CHROM, four$POS, four$REF, four$ALT, sep = ".")
four <- four$snv_id
four <- unique(four)
	
five <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
five <- getFIX(five)
five <- as.data.frame(five)
five <- five[,c("CHROM","POS","REF","ALT")]
five$snv_id <- paste(five$CHROM, five$POS, five$REF, five$ALT, sep = ".")
five <- five$snv_id
five <- unique(five)

tmp <- gsub("_|r","", samp)
indels <- list(one, two , three, four, five)
list_indels <- setNames(indels, tmp)
					
p1 <- upset(fromList(list_indels),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Indel Overlaps", sets.x.label = "Total Indels per Sample",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/indel/", group, "_indel_loc_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

}

group <- "D"

keep <- grep(group, samples)
samp <- samples[keep]

one <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
one <- getFIX(one)
one <- as.data.frame(one)
one <- one[,c("CHROM","POS","REF","ALT")]
one$snv_id <- paste(one$CHROM, one$POS, one$REF, one$ALT, sep = ".")
one <- one$snv_id
one <- unique(one)

two <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
two <- getFIX(two)
two <- as.data.frame(two)
two <- two[,c("CHROM","POS","REF","ALT")]
two$snv_id <- paste(two$CHROM, two$POS, two$REF, two$ALT, sep = ".")
two <- two$snv_id
two <- unique(two)
	
three <- read.vcfR(paste("./raw/caveman/", samp, "vs_E1_SNV_ASMD140_CLPM0.vcf.bgz", sep = ""))
three <- getFIX(three)
three <- as.data.frame(three)
three <- three[,c("CHROM","POS","REF","ALT")]
three$snv_id <- paste(three$CHROM, three$POS, three$REF, three$ALT, sep = ".")
three <- three$snv_id
three <- unique(three)


tmp <- gsub("_|r","", samp)
snvs <- list(one, two , three)
list_snvs <- setNames(snvs, tmp)
					
p1 <- upset(fromList(list_snvs),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "SNV Overlaps", sets.x.label = "Total SNVs per Sample",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/snv/", group, "_snv_loc_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

one <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
one <- getFIX(one)
one <- as.data.frame(one)
one <- one[,c("CHROM","POS","REF","ALT")]
one$snv_id <- paste(one$CHROM, one$POS, one$REF, one$ALT, sep = ".")
one <- one$snv_id
one <- unique(one)

two <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
two <- getFIX(two)
two <- as.data.frame(two)
two <- two[,c("CHROM","POS","REF","ALT")]
two$snv_id <- paste(two$CHROM, two$POS, two$REF, two$ALT, sep = ".")
two <- two$snv_id
two <- unique(two)
	
three <- read.vcfR(paste("./raw/pindel/", samp, "vs_E1_Indels_QUAL250_REP9.vcf.bgz", sep = ""))
three <- getFIX(three)
three <- as.data.frame(three)
three <- three[,c("CHROM","POS","REF","ALT")]
three$snv_id <- paste(three$CHROM, three$POS, three$REF, three$ALT, sep = ".")
three <- three$snv_id
three <- unique(three)


tmp <- gsub("_|r","", samp)
indels <- list(one, two , three)
list_indels <- setNames(indels, tmp)
					
p1 <- upset(fromList(list_indels),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Indel Overlaps", sets.x.label = "Total Indels per Sample",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/indel/", group, "_indel_loc_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()
