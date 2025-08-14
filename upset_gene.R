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
one <- one$annot.symbol
one <- unique(one)
two <- read.csv(paste0("./snv/", samp[2], "snv_annotated.csv"), header = T)
two <- two$annot.symbol
two <- unique(two)
three <- read.csv(paste0("./snv/", samp[3], "snv_annotated.csv"), header = T)
three <- three$annot.symbol
three <- unique(three)
four <- read.csv(paste0("./snv/", samp[4], "snv_annotated.csv"), header = T)
four <- four$annot.symbol
four <- unique(four)
five <- read.csv(paste0("./snv/", samp[5], "snv_annotated.csv"), header = T)
five <- five$annot.symbol
five <- unique(five)

tmp <- gsub("_|r","", samp)
genes <- list(one, two , three, four, five)
list_genes <- setNames(genes, tmp)
					
p1 <- upset(fromList(list_genes),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Gene Overlaps - SNVs", sets.x.label = "Total Genes per Sample - SNVs",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/", group, "_snv_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

one <- read.csv(paste0("./indel/", samp[1], "indel_annotated.csv"), header = T)
one <- one$annot.symbol
one <- unique(one)
two <- read.csv(paste0("./indel/", samp[2], "indel_annotated.csv"), header = T)
two <- two$annot.symbol
two <- unique(two)
three <- read.csv(paste0("./indel/", samp[3], "indel_annotated.csv"), header = T)
three <- three$annot.symbol
three <- unique(three)
four <- read.csv(paste0("./indel/", samp[4], "indel_annotated.csv"), header = T)
four <- four$annot.symbol
four <- unique(four)
five <- read.csv(paste0("./indel/", samp[5], "indel_annotated.csv"), header = T)
five <- five$annot.symbol
five <- unique(five)

genes <- list(one, two , three, four, five)
list_genes <- setNames(genes, tmp)
					
p1 <- upset(fromList(list_genes),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Gene Overlaps - Indels", sets.x.label = "Total Genes per Sample - Indels",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/", group, "_indel_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

}

group <- "D"

keep <- grep(group, samples)
samp <- samples[keep]

one <- read.csv(paste0("./snv/", samp[1], "snv_annotated.csv"), header = T)
one <- one$annot.symbol
one <- unique(one)
two <- read.csv(paste0("./snv/", samp[2], "snv_annotated.csv"), header = T)
two <- two$annot.symbol
two <- unique(two)
three <- read.csv(paste0("./snv/", samp[3], "snv_annotated.csv"), header = T)
three <- three$annot.symbol
three <- unique(three)

tmp <- gsub("_|r","", samp)
genes <- list(one, two , three)
list_genes <- setNames(genes, tmp)
					
p1 <- upset(fromList(list_genes),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Gene Overlaps - SNVs", sets.x.label = "Total Genes per Sample - SNVs",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/", group, "_snv_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()

one <- read.csv(paste0("./indel/", samp[1], "indel_annotated.csv"), header = T)
one <- one$annot.symbol
one <- unique(one)
two <- read.csv(paste0("./indel/", samp[2], "indel_annotated.csv"), header = T)
two <- two$annot.symbol
two <- unique(two)
three <- read.csv(paste0("./indel/", samp[3], "indel_annotated.csv"), header = T)
three <- three$annot.symbol
three <- unique(three)


genes <- list(one, two , three)
list_genes <- setNames(genes, tmp)
					
p1 <- upset(fromList(list_genes),sets = tmp, keep.order = T, order.by = "freq",
			mainbar.y.label = "Gene Overlaps - Indels", sets.x.label = "Total Genes per Sample - Indels",
			point.size = 7, line.size = 3, text.scale = c(5, 5, 3, 3, 7, 2.5))

png(file = paste0("./upset/", group, "_indel_upset.png"), width = 1500, height = 1000)
print(p1)
dev.off()
