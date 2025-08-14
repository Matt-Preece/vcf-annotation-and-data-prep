library(UpSetR, dplyr)

setwd("/data/home/qp241207/vcf_all/annotated/")

samples <- read.delim("names.txt", header = F)
samples <- sort(samples[,1])

groups <- substr(samples,1,1)
groups <- unique(groups)

for(group in groups[1:3]) {

keep <- grep(group, samples)
samp <- samples[keep]

one <- read.csv(paste0("./snv/", samp[1], "snv_annotated.csv"), header = T)
one$snv_id <- paste(one$seqnames, one$start, one$ref, one$alt, sep = ".")
one <- one$snv_id
one <- unique(one)
two <- read.csv(paste0("./snv/", samp[2], "snv_annotated.csv"), header = T)
two$snv_id <- paste(two$seqnames, two$start, two$ref, two$alt, sep = ".")
two <- two$snv_id
two <- unique(two)
three <- read.csv(paste0("./snv/", samp[3], "snv_annotated.csv"), header = T)
three$snv_id <- paste(three$seqnames, three$start, three$ref, three$alt, sep = ".")
three <- three$snv_id
three <- unique(three)
four <- read.csv(paste0("./snv/", samp[4], "snv_annotated.csv"), header = T)
four$snv_id <- paste(four$seqnames, four$start, four$ref, four$alt, sep = ".")
four <- four$snv_id
four <- unique(four)
five <- read.csv(paste0("./snv/", samp[5], "snv_annotated.csv"), header = T)
five$snv_id <- paste(five$seqnames, five$start, five$ref, five$alt, sep = ".")
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

one <- read.csv(paste0("./indel/", samp[1], "indel_annotated.csv"), header = T)
one$indel_id <- paste(one$seqnames, one$start, one$ref, one$alt, sep = ".")
one <- one$indel_id
one <- unique(one)
two <- read.csv(paste0("./indel/", samp[2], "indel_annotated.csv"), header = T)
two$indel_id <- paste(two$seqnames, two$start, two$ref, two$alt, sep = ".")
two <- two$indel_id
two <- unique(two)
three <- read.csv(paste0("./indel/", samp[3], "indel_annotated.csv"), header = T)
three$indel_id <- paste(three$seqnames, three$start, three$ref, three$alt, sep = ".")
three <- three$indel_id
three <- unique(three)
four <- read.csv(paste0("./indel/", samp[4], "indel_annotated.csv"), header = T)
four$indel_id <- paste(four$seqnames, four$start, four$ref, four$alt, sep = ".")
four <- four$indel_id
four <- unique(four)
five <- read.csv(paste0("./indel/", samp[5], "indel_annotated.csv"), header = T)
five$indel_id <- paste(five$seqnames, five$start, five$ref, five$alt, sep = ".")
five <- five$indel_id
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

one <- read.csv(paste0("./snv/", samp[1], "snv_annotated.csv"), header = T)
one$snv_id <- paste(one$seqnames, one$start, one$ref, one$alt, sep = ".")
one <- one$snv_id
one <- unique(one)
two <- read.csv(paste0("./snv/", samp[2], "snv_annotated.csv"), header = T)
two$snv_id <- paste(two$seqnames, two$start, two$ref, two$alt, sep = ".")
two <- two$snv_id
two <- unique(two)
three <- read.csv(paste0("./snv/", samp[3], "snv_annotated.csv"), header = T)
three$snv_id <- paste(three$seqnames, three$start, three$ref, three$alt, sep = ".")
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

one <- read.csv(paste0("./indel/", samp[1], "indel_annotated.csv"), header = T)
one$indel_id <- paste(one$seqnames, one$start, one$ref, one$alt, sep = ".")
one <- one$indel_id
one <- unique(one)
two <- read.csv(paste0("./indel/", samp[2], "indel_annotated.csv"), header = T)
two$indel_id <- paste(two$seqnames, two$start, two$ref, two$alt, sep = ".")
two <- two$indel_id
two <- unique(two)
three <- read.csv(paste0("./indel/", samp[3], "indel_annotated.csv"), header = T)
three$indel_id <- paste(three$seqnames, three$start, three$ref, three$alt, sep = ".")
three <- three$indel_id
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
