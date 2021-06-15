#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Texas A&M University, Department of Biomedical Engineering
#' Spring 2021
#' Description: This program annotates the genome
#' 
#load the related libraries
library(GenomicFeatures)
library(biomaRt)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(jamba)
library(splicejam)
library(biovizBase)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(stringr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
#color scheme: cadet blue, eggplant, Persian orange, Arylide yellow, yellow-green crayola 
color.scheme <- c('#51A3A3', '#75485E', '#CB904D', '#37D5A3', "#DFCC74", "#C3E991")
#source the functions script
source("/Users/annicenajafi/Downloads/sc/scripts/AnnotatingGenomeFunctions_Mod_April.R")
#Run the program
genome.hold <-"hg38"; blacklist.bed <- file.path("/Users/annicenajafi/Downloads/ENCFF356LFX.bed"); 
getBSgenome(genome.hold) -> BSgenome
#1. Make a TxDB object
txdb <- makeTxDbFromUCSC(genome.hold, "refGene")
#2. Get a list of gene names
gene_names <- TxToGene(txdb)
#3. Map gene names to transcripts
tx2gene  <- txdbTxToGene(txdb)
#4. Using the list above flatten transcripts to the genomic level for CDS, 3'UTR and 5'UTR and introns
cdsFlat <- CDSExonTxtoGene(txdb, tx2gene);threeFlat <- threeUTRExonTxtoGene(txdb, tx2gene);fiveFlat <- fiveUTRExonTxtoGene(txdb, tx2gene);flatIntronsByGene<- getIntrons(txdb, gene_names, tx2gene) 
#5. Get the remaining exons and label them as utr
rem_ex <- ExtractncRNAExons(cdsFlat, threeFlat, fiveFlat, txdb)
#7. Combine all annotations and check for overlaps and re-annotate and ADD entrezgene_ids
all <- CombineOverlap(txdb, flatIntronsByGene, cdsFlat, threeFlat, fiveFlat, rem_ex)
#8. Add flanking regions 
all <- flankAll(all)
#9. Add intergenic regions
all <- AddIntergenic(txdb, all)
#10. Transform to data table and reduce by gene
all.df <- as.data.table(all)
all.df <- all.df[, as.data.table(reduce(IRanges(start, end))), by = .(seqnames, strand,entrezgene_id, gene_name, feature_type, exon_name)]
#11. Check the number of ranges for each feature type
table(all.df$feature_type)
#12. Make a boxplot for the lengths of the ranges per feature type '#48D1CC'
ggplot(all.df, aes(x=feature_type, y=width)) + geom_boxplot(fill=color.scheme[1], color="black") +  coord_cartesian(ylim=c(0, 5000)) + ylab('Length(bp)') + ggtitle(genome.hold)
#12. Add three utr index
all <- AddThreeUTRIndex(GRanges(all.df)); all <- ConsiderBlacklist(all, blacklist.bed)
#13. Rank all features
all <- RankAll(all)
  #14. Calculate GC content <<CAUTION>> the following line takes a couple hours to run/// need to change it
  #all <- getGCPercentage(all, BSgenome)
#15. put labels for what range comes before and after every range
all <- FlankUpAndFlankDown(all)
all <- getProteinCoding(all, tx2gene)
all <- ExtractIntrons(all)
#16.convert to Granges
as.data.frame(all) -> all.df
newAll <- all.df[order(all.df$seqnames, all.df$start),]; GRanges(all.df) -> all
#newAll <- newAll[newAll$seqnames=="chr1",]
write.csv(newAll, "/Users/annicenajafi/Downloads/hg38_annotated_introns.csv")
#add RDS info
AddRDSInfo(all, genome.hold) -> all
#save the rds file
saveRDS(all, file = "/Users/annicenajafi/Desktop/hg38_annotated_introns_mar.rds")
#TTLL10, SLC35E2, 
hold.me <- export.bed(all, con='/Users/annicenajafi/Downloads/hg38_annice_final.bed')




