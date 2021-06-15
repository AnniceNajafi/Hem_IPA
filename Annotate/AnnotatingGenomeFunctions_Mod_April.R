#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Texas A&M University, Department of Biomedical Engineering
#' Spring 2021, Date: 12/26/2020
#' Description: This program annotates the  genome 
#' 
options(warn=-1)
#' txdbTxToGene 
#' Description: this function maps transcripts to genes
#' @return a dataframe with genes mapped to transcripts
txdbTxToGene <- function(txdb){
  #match the transcript id to gene id
  tx2gene<- suppressMessages(AnnotationDbi::select(txdb,  keys(txdb, "GENEID"),
                                                   columns=c("GENEID","TXNAME", "TXID"),keytype="GENEID"))
  tx2gene <- renameColumn(tx2gene,from=c("GENEID", "TXNAME", "TXTYPE"),
                          to=c("gene_id", "transcript_id", "transcript_type"))
  #step II*//uncomment to try
  #tx2gene$transcript_id <- genesList[match(tx2gene$TXNAME, genesList$ensembl_gene_id), 1]
  #entrez2gene<- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys(org.Hs.eg.db, "ENTREZID"),columns=c("ENTREZID","SYMBOL"),keytype="ENTREZID"))
  #get the gene names
  txdb <- txdb
  gene_names <- TxToGene(txdb)
  #relate all gene_ids to genes
  tx2gene$gene_name <- gene_names[as.character(tx2gene$gene_id)]
  tx2gene$biotype<-c(NA)
  tx2gene[substr(tx2gene$transcript_id, start=1, stop=2)=="NM",]$biotype<-"protein_coding"
  tx2gene[substr(tx2gene$transcript_id, start=1, stop=2)=="NR",]$biotype<-"non_protein_coding"
  return (tx2gene)
}

#' TxToGeneHg
#' Description: helper function that retrieves gene names given a txdb object
#' @return gene names
TxToGene <- function(txdb){
  gene_ids <- values(genes(txdb))$gene_id
  #gene_ids <- entrez2gene$ENTREZID
  #get all the gene names
  gene_names <- mget(gene_ids,org.Hs.egSYMBOL,ifnotfound=NA)
  ## Convert list to vector taking the first gene_name
  gene_names <- unlist(heads(S4Vectors::List(gene_names), 1))
  ## switch NA with LOC# format
  if (any(is.na(gene_names))) {
    gene_na <- which(is.na(gene_names))
    gene_names[gene_na] <- gsub("^([0-9]+)$", "LOC\\1",names(gene_names[gene_na]))
  }
  return (gene_names)
}

#' CDSExonTxtoGene
#' Description: This function receives a txdb object along with a dataframe with transcripts mapped to genes
#' and returns the cds exonic transcripts flattened to the genome level
#' @param txdb a TxDB object
#' @param tx2gene a dataframe with transcripts mapped to genes
#'
#' @return genomicrangeslist of cds regions
CDSExonTxtoGene <- function(txdb, tx2gene){
  #get CDS exons by transcript 
  cdsByTx <- cdsBy(txdb, by="tx",use.names=TRUE)
  #flatten CDS exon tx to genome level
  cdsFlat <- flattenExonsBy(exonsByTx=cdsByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  #label as cds exons
  GenomicRanges::values(cdsFlat@unlistData)[,"feature_type"] <- "cds_exon"
  return (cdsFlat)
  
}

#' threeUTRExonTxtoGene
#' Description: This function receives a txdb object along with a dataframe with transcripts mapped to genes
#' and returns the 3'UTR exonic transcripts flattened to the genome level
#' @param txdb a TxDB object
#' @param tx2gene a dataframe with transcripts mapped to genes
#'
#' @return genomicrangeslist of 3'UTR regions
threeUTRExonTxtoGene <- function(txdb, tx2gene){
  #get 3'UTR exons by transcript
  threeUTRByTx <- threeUTRsByTranscript(txdb, use.names=TRUE)
  #flatten 3'UTRexon tx to genome level
  threeFlat <- flattenExonsBy(exonsByTx=threeUTRByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  #label as 3'UTR exons
  GenomicRanges::values(threeFlat@unlistData)[,"feature_type"] <- "utr3"
  return (threeFlat)
  
}
#' fiveUTRExonTxtoGene
#' Description: This function receives a txdb object along with a dataframe with transcripts mapped to genes
#' and returns the 5'UTR exonic transcripts flattened to the genome level
#' @param txdb a TxDB object
#' @param tx2gene a dataframe with transcripts mapped to genes
#'
#' @return genomicrangeslist of 5'UTR regions
fiveUTRExonTxtoGene <- function(txdb, tx2gene){
  #get 5'UTR exons by transcript
  fiveUTRByTx <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  #flatten 5'UTRexon tx to genome level
  fiveFlat <- flattenExonsBy(exonsByTx=fiveUTRByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  #label as 5'UTR exons
  GenomicRanges::values(fiveFlat@unlistData)[,"feature_type"] <- "utr5"
  return (fiveFlat)
  
}

#' ExtractncRNAExons
#' Description: the following program finds the untranslated exons
#' @param cdsFlat 
#' @param threeFlat 
#' @param fiveFlat 
#' @param txdb 
#'
#' @return the utr regions on genome
ExtractncRNAExons <- function(cdsFlat, threeFlat, fiveFlat, txdb){
  exonsTx <- c(cdsFlat, threeFlat, fiveFlat)
  #get exons by transcript
  exonsByTx <- exonsBy(txdb,by="tx",use.names=TRUE)
  #get CDS exons by transcript 
  cdsByTx <- cdsBy(txdb, by="tx",use.names=TRUE)
  #At this point we have all the transcripts matched to gene for 3'UTR, 5'UTR and CDS exons 
  #We need to subtract these from the entire exons set to get the ncRNA exons
  exons_total <- exonsByTx
  #Flatten to genome level
  exons_total  <- flattenExonsBy(exonsByTx=exonsByTx, cdsByTx=cdsByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  #get the remaining exons by subtracting the 3'UTR, 5'UTR and cds regions from all exons and label as utr
  exons_total <- SubtractOverlap(exons_total, threeFlat, "utr")
  #do for 3'utr regions
  exons_total <- SubtractOverlap(exons_total, fiveFlat, "utr")
  #for cds regions
  exons_total <- SubtractOverlap(exons_total, cdsFlat, "utr")
  #store in a gr
  rem_ex <- exons_total
  #add column for gene symbol
  #GenomicRanges::values(rem_ex@unlistData)[,"gene_name"] <- rep(names(rem_ex), elementNROWS(rem_ex))
  #label as utr
  #GenomicRanges::values(rem_ex@unlistData)[,"feature_type"]  <- 'utr'
  return (rem_ex)
}

#' CheckUTRxthree
#' Description: this function finds the overlapping 3'UTR of a transcript annotated as utr with introns
#' @param flatIntronsByGene intronic regions flattened to genome
#' @param rem_ex utr regions
#' @param threeFlat three utr regions flattened to genome
#'
#' @return the overlapping 3'UTR region of a transcript annotated as utr with introns
CheckUTRxthree <- function(txdb, flatIntronsByGene, rem_ex, threeFlat){
  #get the regions by tx
  exonsByTx <- exonsBy(txdb,by="tx",use.names=TRUE)
  threeUTRByTx <- threeUTRsByTranscript(txdb, use.names=TRUE)
  intronsByTx <- intronsByTranscript(txdb, use.names=TRUE)
  #find the overlapping region between introns and utrs
  hold.ol.utrxthree <- FindOverlapByRanges(exonsByTx, FindOverlapByRanges(flatIntronsByGene, rem_ex))
  #get the tx names
  hold.tx.name <- unique(hold.ol.utrxthree$i.group_name)
  #get the gene names
  gene.names <- unique(hold.ol.utrxthree$gene_name)
  #convert intronstotx to datatable 
  intronsByTx.df <- as.data.table(intronsByTx)
  #check if tx has a 3utr region if so then it is noise get rid of it
  hold.ol.utrxthree<- hold.ol.utrxthree[!(hold.ol.utrxthree$i.group_name %in% names(threeUTRByTx)),]
  #check of tx has an intron and not just a utr floating around
  hold.ol.utrxthree<- hold.ol.utrxthree[(hold.ol.utrxthree$i.group_name %in% intronsByTx.df$group_name),]
  #check if the gene has an actual 3utr region
  hold.ol.utrxthree<- hold.ol.utrxthree[(hold.ol.utrxthree$gene_name %in% names(threeFlat)),]
  #filter for negative strand
  hold.ol.utrxthree.ne <- hold.ol.utrxthree[strand(hold.ol.utrxthree) == "-"]
  #filter for positive strand
  hold.ol.utrxthree <- hold.ol.utrxthree[strand(hold.ol.utrxthree) == "+"]
  #filter introns in negative  strand
  intronsByTx.ne <- intronsByTx[strand(intronsByTx) == "-"]
  # " " in positive strand
  intronsByTx <- intronsByTx[strand(intronsByTx) == "+"]
  #get the transcript names again after the step done above
  hold.tx.name <- unique(hold.ol.utrxthree$i.group_name)
  #for every transcript check if there exists an intron that is larger than the utr if so delete that transcript and do not consider it
  for(i in hold.tx.name){
    if(nrow(as.data.table(intronsByTx[i]))!=0 && length(hold.ol.utrxthree[hold.ol.utrxthree$i.group_name == i])!=0){
      if((GRanges(tail(as.data.table(intronsByTx[names(intronsByTx) == i]), 1))) > hold.ol.utrxthree[hold.ol.utrxthree$i.group_name == i]){
      
          hold.ol.utrxthree<- hold.ol.utrxthree[!(hold.ol.utrxthree$i.group_name == i),]
         
      }
      else if((GRanges(head(as.data.table(intronsByTx[names(intronsByTx) == i]), 1))) > hold.ol.utrxthree[hold.ol.utrxthree$i.group_name == i]){
        
        hold.ol.utrxthree<- hold.ol.utrxthree[!(hold.ol.utrxthree$i.group_name == i),]
        
      }
    }
  }
  
  hold.tx.name.ne <- unique(hold.ol.utrxthree.ne$i.group_name)
  gene.names.ne <- unique(hold.ol.utrxthree.ne$gene_name)
 #do the same for the other strand
  for(i in hold.tx.name.ne){
    if(nrow(as.data.table(intronsByTx.ne[i]))!=0 && length(hold.ol.utrxthree.ne[hold.ol.utrxthree.ne$i.group_name == i])!=0){
      if((GRanges(tail(as.data.table(intronsByTx.ne[names(intronsByTx.ne) == i]), 1))) < hold.ol.utrxthree.ne[hold.ol.utrxthree.ne$i.group_name == i]){
        
        hold.ol.utrxthree.ne<- hold.ol.utrxthree.ne[!(hold.ol.utrxthree.ne$i.group_name == i),]
        
      }
      else if((GRanges(head(as.data.table(intronsByTx.ne[names(intronsByTx.ne) == i]), 1))) < hold.ol.utrxthree.ne[hold.ol.utrxthree.ne$i.group_name == i]){
        
        hold.ol.utrxthree.ne<- hold.ol.utrxthree.ne[!(hold.ol.utrxthree.ne$i.group_name == i),]
        
      }
    }
  }
  #MIB2
  #SARS1
  #PAK3
  hold.ol.utrxthree <- c(hold.ol.utrxthree.ne, hold.ol.utrxthree)
  return (hold.ol.utrxthree)
}

#' getIntrons
#' Description: This function receives a txdb object along with a gene_names list and 
#' a list mapping transcript ids to gene symbols and flattens the intronic transcripts
#' to the genome level
#' @param txdb TxDB object
#' @param gene_names dataframe with gene_names
#' @param tx2gene dataframe with gene names mapped to transcripts
#'
#' @return intronic transcripts flattened to the genome level
getIntrons <- function(txdb, gene_names, tx2gene){
  #get the introns by transcript
  intronsByTx <- intronsByTranscript(txdb, use.names=TRUE)
  #flatten to genome
  flatIntronsByGene <- flattenIntronsBy(intronsByTx = intronsByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  #annotate
  GenomicRanges::values(flatIntronsByGene@unlistData)[,"subclass"] <- "intron"
  return (flatIntronsByGene)
}

#' getIntergenic
#' Description: Receives a txdb object and find the intergenic regions by calculating the gap between 
#' data entries
#' @param txdb TxDB object
#'
#' @return intergenic regions
AddIntergenic <- function(txdb, all){
  all <- c(all, getGeneOverlap(txdb))
  #Find gaps between entries
  intergenicRegions <- gaps((all))
  #label as intergenic
  GenomicRanges::values(intergenicRegions)[,"feature_type"] <- "intergenic"
  #combine with the rest
  intergeneAdded <- c(all, intergenicRegions)
  intergeneAdded <- GRanges(as.data.table(intergeneAdded)[!(as.data.table(intergeneAdded)$strand == "*"),])
  return (intergeneAdded)
}




#' CombineOverlap
#' Description: This function receives the introns and exons flattened to the genome 
#' and finds the intron.3utr region then 
#' @param flatIntronsByGene intronic regions flattened to genome
#' @param cdsFlat 5cds regions flattened to genome
#' @param threeFlat 3'UTR regions flattened to genome
#' @param fiveFlat 5'UTR regions flattened to genome
#' @param rem_ex utr regions flattened to genome
#'
#' @return the gr of all annotated regions
CombineOverlap <- function(txdb, flatIntronsByGene, cdsFlat, threeFlat, fiveFlat, rem_ex){
  utrIntronOverlap <- CheckUTRxthree(txdb, flatIntronsByGene, rem_ex, threeFlat)
  flatIntronsByGene <- unique(flatIntronsByGene@unlistData[-queryHits(findOverlaps(flatIntronsByGene@unlistData, threeFlat@unlistData, type="within")),])
  #Find overlapping region between 3UTR and introns
  threeIntronOverlap <- FindOverlapByRanges(flatIntronsByGene, threeFlat)
  #label as intron.3utr
  threeIntronOverlap$feature_type <- "intron.3utr"
  utrIntronOverlap$feature_type <- "intron.3utr"
  #switch group_name to gene_name to be consistent with the rest of the grs
  threeIntronOverlap$gene_name <- threeIntronOverlap$group_name; threeIntronOverlap$entrezgene_id <- NULL;threeIntronOverlap$group_name<- NULL
  utrIntronOverlap$gene_name <- utrIntronOverlap$group_name; utrIntronOverlap$entrezgene_id <- NULL;utrIntronOverlap$group_name<- NULL
  #remove the overlapping 3 utr region from introns
  flatIntronsByGene <-  SubtractOverlap(flatIntronsByGene, threeFlat, "intron")
  #remove the overlapping cds region from introns
  flatIntronsByGene <-  SubtractOverlap(flatIntronsByGene,  cdsFlat, "intron")
  #remove the overlapping 5 utr region from introns
  flatIntronsByGene <-  SubtractOverlap(flatIntronsByGene,  fiveFlat,  "intron")
  #remove the overlapping utr region from introns
  flatIntronsByGene <-  SubtractOverlap(flatIntronsByGene, rem_ex, "intron")
  #remove the overlapping region from introns
  flatIntronsByGene <- SubtractOverlap(flatIntronsByGene, getGeneOverlap(txdb), "intron")
  #reduce the intronic annotation by gene
  #flatIntronsByGene <- reduce(split(flatIntronsByGene,flatIntronsByGene$gene_name))
  GRanges(as.data.frame(flatIntronsByGene))
  #subtract the overlapping region from utrs
  rem_ex <- SubtractOverlap(rem_ex, threeIntronOverlap, "utr")
  rem_ex <- SubtractOverlap(rem_ex, utrIntronOverlap, "utr")
  rem_ex <- SubtractOverlap(rem_ex, fiveFlat, "utr")
  rem_ex <- SubtractOverlap(rem_ex, threeFlat, "utr")
  rem_ex <- SubtractOverlap(rem_ex, cdsFlat, "utr")
  #
  #remove 3 utr region from 5utr
  fiveFlat <- SubtractOverlap(fiveFlat, threeFlat, "utr5")
  #subtract the overlapping region from 5'UTRs
  fiveFlat <- SubtractOverlap(fiveFlat, threeIntronOverlap, "utr5")
  #subtract the cds region from 5'UTRs
  fiveFlat <- SubtractOverlap(fiveFlat, cdsFlat, "utr5")
  
  threeFlat <- SubtractOverlap(threeFlat, cdsFlat, "utr3")
  threeFlat <- SubtractOverlap(threeFlat, threeIntronOverlap, "utr3")
  threeFlat <- GRanges(as.data.frame(threeFlat))
  #subtract the overlapping region from cds
  cdsFlat <- SubtractOverlap(cdsFlat, threeIntronOverlap, "cds_exon")
  #combine the three exonic regions
  allExons <- GRanges(as.data.frame(c(threeFlat, fiveFlat, cdsFlat, rem_ex)))
  #combine intronic regions, exonic regions which are non overlapping with 3'UTR intronic regions
  all <- c(allExons, flatIntronsByGene, threeIntronOverlap, utrIntronOverlap)
  #transform into dataframe to remove unnecessary columns
  all.df <- as.data.frame(all)
  ##all.df %>% filter(gene_name=="LOC101927604")
  #remove some of the columns
  all.df$group_name<-NULL; all.df$subclass<-NULL; all.df$gene_nameExon<-NULL; all.df$group<-NULL;all.df$gene_nameIntron<-NULL;
  #add the entrezgene_id based on gene symbol
  all.df$entrezgene_id <- tx2gene[match(all.df$gene_name, tx2gene$gene_name), 1]
  #
  keys <- c("intron", "cds_exon", "utr5*", "utr3*", "utr5", "utr", "utr3", "intron.3utr")
  keyDF <- data.frame(key=keys,weight=1:length(keys))
  merged <- merge(all.df,keyDF,by.x='feature_type',by.y='key',all.x=T,all.y=F)
  all.df <- merged[order(merged$weight),c("seqnames", "start", "end", "width", "strand", "gene_name", "feature_type", "entrezgene_id")]
  all.df[!duplicated(all.df[,c("seqnames", "start", "end", "strand")]),]
  #convert to gr
  all <- GRanges(all.df)
  return(all)
}


#Note: Convert grangeslist to dataframe then look for psuedogenes by filtering for UTRs then 
#take the smallest start for every pseudogene and subtract 5000 and label as utr3* 
#take the largest end for every pseudogene and add 1000 and label as utr5*
#' AddPseudoFlank
#' Description: This function adds flanking regions for pseudogenes
#' @param all receives the entire gr
#'
#' @return the grange with flanking utr regions for pseudogenes added
AddPseudoFlank <- function(all){
  all <- as.data.table(all)
  #store entire gr as dataframe
  hold.all <- all
  #filter for utr regions
  #delete all rows where a gene that has a 3utr exists
  genesToDl <- hold.all %>% dplyr::filter(feature_type == "utr3")
  #sort the first df
  setkey(genesToDl, seqnames, start, end)
  #sort the second df
  setkey(hold.all, seqnames, start, end)
  
  hold.all <- hold.all[!(hold.all$gene_name %in% genesToDl$gene_name),]
  
  hold.all.utr.plus <- hold.all %>% dplyr::filter(feature_type == 'utr') %>% dplyr::filter(strand == "+")
  #order by start
  hold.utr.start.sorted.plus <- hold.all.utr.plus[order(hold.all.utr.plus$start),]
  #Take only the datapoint on top meaning the one with the smallest starting site
  flanked.five.plus <- hold.utr.start.sorted.plus[match(unique(hold.utr.start.sorted.plus$gene_name), hold.utr.start.sorted.plus$gene_name),]
  #annotate flanking 5'
  flanked.five.plus$end <- flanked.five.plus$start - 1
  #
  flanked.five.plus$start <- flanked.five.plus$start - 1000
  #
  flanked.five.plus$feature_type <- "utr5*"
  flanked.five.plus <- as.data.table(SubtractOverlap(GRanges(flanked.five.plus), all, "utr5*"))
  #sort this time in the opposite direction
  hold.utr.sorted.three.plus <- hold.all.utr.plus[order(hold.all.utr.plus$start, decreasing = TRUE),]
  #Take only the datapoint on top meaning the one with the smallest starting site
  flanked.three.plus <- hold.utr.sorted.three.plus[match(unique(hold.utr.sorted.three.plus$gene_name), hold.utr.sorted.three.plus$gene_name),]
  #find the region add 5000 to the start
  flanked.three.plus$start <- flanked.three.plus$end +1
  #
  flanked.three.plus$end <- flanked.three.plus$start + 4999
  #label as utr3*
  flanked.three.plus$feature_type <- "utr3*"
  flanked.three.plus <- as.data.table(SubtractOverlap(GRanges(flanked.three.plus), all, "utr3*"))
  #combine all flanking regions with gr
  allWithFlank.plus <- rbind(flanked.three.plus, flanked.five.plus)
  allWithFlank.plus$entrezgene_id <- tx2gene[match(allWithFlank.plus$gene_name, tx2gene$gene_name), 1]
  #Negative strands ...
  hold.all.utr.minus <- hold.all %>% dplyr::filter(feature_type == 'utr') %>% dplyr::filter(strand == "-")
  #order by start
  hold.utr.start.sorted.minus <- hold.all.utr.minus[order(hold.all.utr.minus$start),]
  #Take only the datapoint on top meaning the one with the smallest starting site
  flanked.five.minus <- hold.utr.start.sorted.minus[match(unique(hold.utr.start.sorted.minus$gene_name), hold.utr.start.sorted.minus$gene_name),]
  #annotate flanking 5'
  flanked.five.minus$end <- flanked.five.minus$start - 1
  #
  flanked.five.minus$start <- flanked.five.minus$start - 5000
  #
  flanked.five.minus$feature_type <- "utr3*"
  flanked.five.minus <- as.data.table(SubtractOverlap(GRanges(flanked.five.minus), all, "utr3*"))
  #sort this time in the opposite direction
  hold.utr.sorted.three.minus <- hold.all.utr.minus[order(hold.all.utr.minus$start, decreasing = TRUE),]
  #Take only the datapoint on top meaning the one with the smallest starting site
  flanked.three.minus <- hold.utr.sorted.three.minus[match(unique(hold.utr.sorted.three.minus$gene_name), hold.utr.sorted.three.minus$gene_name),]
  #find the region add 5000 to the start
  flanked.three.minus$start <- flanked.three.minus$end +1
  #
  flanked.three.minus$end <- flanked.three.minus$start + 999
  #label as utr3*
  flanked.three.minus$feature_type <- "utr5*"
  flanked.three.minus <- as.data.table(SubtractOverlap(GRanges(flanked.three.minus), all, "utr5*"))
  #combine all flanking regions with gr
  allWithFlank.minus <- rbind(flanked.three.minus, flanked.five.minus)
  allWithFlank.minus$entrezgene_id <- tx2gene[match(allWithFlank.minus$gene_name, tx2gene$gene_name), 1]
  allWithFlank <- rbind(allWithFlank.minus, allWithFlank.plus, all)
  #transform to gr
  allWithFlank <- GRanges(allWithFlank)
  return(allWithFlank)
}

#' AddGeneFlank
#' Description: The following function adds flanking UTR regions 
#' for gene
#' @param all receives the entire gr
#'
#' @return the grange with flanking utr regions for genes added
AddGeneFlank <- function(all){
  #Transform to data frame
  hold.all <- as.data.table(all)
  #dplyr::filter for 5' utr regions
  hold.all.exon.five.plus <- hold.all %>% dplyr::filter(feature_type == "utr5") %>% dplyr::filter(strand == "+")
  #sort the regions such that the one starting first appears on top
  hold.exon.five.sorted.plus <- hold.all.exon.five.plus[order(hold.all.exon.five.plus$start),]
  #
  flanked.five.plus <- hold.exon.five.sorted.plus [match(unique(hold.exon.five.sorted.plus$gene_name), hold.exon.five.sorted.plus$gene_name),]
  #
  flanked.five.plus$end <- flanked.five.plus$start - 1
  #
  flanked.five.plus$start <- flanked.five.plus$start - 1000
  #Label as utr5*
  flanked.five.plus$feature_type <- "utr5*"
  flanked.five.plus <- as.data.table(SubtractOverlap(GRanges(flanked.five.plus), all, "utr5*"))
  #now dplyr::filter for 3'UTr regions
  hold.all.exon.three.plus <- hold.all %>% dplyr::filter(feature_type == "utr3") %>% dplyr::filter(strand =="+")
  #sort this time in the opposite direction
  hold.exon.three.sorted.plus <- hold.all.exon.three.plus[order(hold.all.exon.three.plus$start, decreasing = TRUE),]
  #take the one entry on top
  flanked.three.plus <- hold.exon.three.sorted.plus [match(unique(hold.exon.three.sorted.plus$gene_name), hold.exon.three.sorted.plus$gene_name),]
  #annotate
  flanked.three.plus$start <- flanked.three.plus$end +1
  #
  flanked.three.plus$end <- flanked.three.plus$start + 4999
  #
  flanked.three.plus$feature_type <- "utr3*" 
  flanked.three.plus <- as.data.table(SubtractOverlap(GRanges(flanked.three.plus), all, "utr3*"))
  #combine the flanking regions with the rest of the annotated parts
  allWithFlank.plus <- rbind(flanked.three.plus, flanked.five.plus)
  #dplyr::filter for 3' utr regions on NEGATIVE STRAND
  hold.all.exon.five.minus <- hold.all %>% dplyr::filter(feature_type == "utr3") %>% dplyr::filter(strand == "-")
  #sort the regions such that the one starting first appears on top
  hold.exon.five.sorted.minus <- hold.all.exon.five.minus[order(hold.all.exon.five.minus$start),]
  #
  flanked.five.minus <- hold.exon.five.sorted.minus [match(unique(hold.exon.five.sorted.minus$gene_name), hold.exon.five.sorted.minus$gene_name),]
  #
  flanked.five.minus$end <- flanked.five.minus$start - 1
  #
  flanked.five.minus$start <- flanked.five.minus$start - 5000
  #Label as utr3*
  flanked.five.minus$feature_type <- "utr3*"
  flanked.five.minus <- as.data.table(SubtractOverlap(GRanges(flanked.five.minus), all, "utr3*"))
  #now dplyr::filter for 5'UTr regions
  hold.all.exon.three.minus <- hold.all %>% dplyr::filter(feature_type == "utr5") %>% dplyr::filter(strand == "-")
  #sort this time in the opposite direction
  hold.exon.three.sorted.minus <- hold.all.exon.three.minus[order(hold.all.exon.three.minus$start, decreasing = TRUE),]
  #take the one entry on top
  flanked.three.minus <- hold.exon.three.sorted.minus [match(unique(hold.exon.three.sorted.minus$gene_name), hold.exon.three.sorted.minus$gene_name),]
  #annotate
  flanked.three.minus$start <- flanked.three.minus$end +1
  #
  flanked.three.minus$end <- flanked.three.minus$start + 999
  #
  flanked.three.minus$feature_type <- "utr5*" 
  flanked.three.minus <- as.data.table(SubtractOverlap(GRanges(flanked.three.minus), all, "utr5*"))
  #combine the flanking regions with the rest of the annotated parts
  allWithFlank.minus <- rbind(flanked.three.minus, flanked.five.minus)
  allWithFlank.minus$entrezgene_id <- tx2gene[match(allWithFlank.minus$gene_name, tx2gene$gene_name), 1]
  allWithFlank.plus$entrezgene_id <- tx2gene[match(allWithFlank.plus$gene_name, tx2gene$gene_name), 1]
  #transform to gr
  #bind all and convert to gr
  allWithFlank <- rbind(allWithFlank.plus, allWithFlank.minus, hold.all)
  #
  allWithFlank <- GRanges(allWithFlank)
  return(allWithFlank)
}

##sort all rows by 
#' AddRestFlank
#' Description: The following function receives a dataframe with all regions annotated except some flanking regions
#' and annotates all remaining flanking regions
#' @param all Receives dataframe of all annotations
#'
#' @return the dataframe with all flanking regions annotated
AddRestFlank <- function(all){
  #Receive all annotations and convert to dataframe
  hold.all <- as.data.frame(all)
  #sort by starting point
  hold.all.negative <- hold.all[order(hold.all$start),]
  #dplyr::filter for negative strands
  hold.all.negative <- hold.all.negative %>% dplyr::filter(strand == "-")
  #take the 
  hold.all.negative <- hold.all.negative[!duplicated(hold.all.negative$gene_name),]
  #dplyr::filter for those who don't have the flanking utr3* region annotated
  hold.all.negative <- hold.all.negative %>% dplyr::filter(feature_type!="utr3*")
  #
  flanked.three.minus<- hold.all.negative
  #annotate
  flanked.three.minus$end <- flanked.three.minus$start - 1
  #
  flanked.three.minus$start <- flanked.three.minus$start - 5000
  #Label as utr3*
  flanked.three.minus$feature_type <- "utr3*"
  flanked.three.minus <- as.data.table(SubtractOverlap(GRanges(flanked.three.minus), all, "utr3*"))
  flanked.three.minus$width <- flanked.three.minus$end - flanked.three.minus$start + 1
  #change the width info
  #
  hold.all <- as.data.frame(all)
  #
  hold.all.negative <- hold.all[order(hold.all$start, decreasing = TRUE),]
  #
  hold.all.negative <- hold.all.negative %>% dplyr::filter(strand == "-")
  #
  hold.all.negative <- hold.all.negative[!duplicated(hold.all.negative$gene_name),]
  #
  hold.all.negative <- hold.all.negative %>% dplyr::filter(feature_type!="utr5*")
  #
  flanked.five.minus <- hold.all.negative
  #annotate
  flanked.five.minus$start <- flanked.five.minus$end +1
  #
  flanked.five.minus$end <- flanked.five.minus$start + 999
  #Label as utr5*
  #
  flanked.five.minus$feature_type <- "utr5*"
  flanked.five.minus <- as.data.table(SubtractOverlap(GRanges(flanked.five.minus), all, "utr5*"))
  flanked.five.minus$width <- flanked.five.minus$end - flanked.five.minus$start +1
  #
  hold.all <- as.data.frame(all)
  #
  hold.all.positive <- hold.all[order(hold.all$start),]
  #
  hold.all.positive <- hold.all.positive %>% dplyr::filter(strand == "+")
  #
  hold.all.positive <- hold.all.positive[!duplicated(hold.all.positive$gene_name),]
  #
  hold.all.positive <- hold.all.positive %>% dplyr::filter(feature_type!="utr5*")
  #
  flanked.five.plus <- hold.all.positive
  #annotate
  flanked.five.plus$end <- flanked.five.plus$start - 1
  #
  flanked.five.plus$start <- flanked.five.plus$start - 1000
  #Label as utr5*
  #
  flanked.five.plus$feature_type <- "utr5*"
  flanked.five.plus <- as.data.table(SubtractOverlap(GRanges(flanked.five.plus), all, "utr5*"))
  flanked.five.plus$width <- flanked.five.plus$end - flanked.five.plus$start +1
  #
  hold.all <- as.data.frame(all)
  #
  hold.all.positive <- hold.all[order(hold.all$start, decreasing = TRUE),]
  #
  hold.all.positive <- hold.all.positive %>% dplyr::filter(strand == "+")
  #
  hold.all.positive <- hold.all.positive[!duplicated(hold.all.positive$gene_name),]
  #
  hold.all.positive <- hold.all.positive %>% dplyr::filter(feature_type!="utr3*")
  #
  flanked.three.plus <- hold.all.positive
  #annotate
  flanked.three.plus$start <- flanked.three.plus$end +1
  #
  flanked.three.plus$end <- flanked.three.plus$start + 4999
  flanked.three.plus <- as.data.table(SubtractOverlap(GRanges(flanked.three.plus), all, "utr5*"))
  flanked.three.plus$width <- flanked.three.plus$end - flanked.three.plus$start + 1
  #
  #
  flanked.three.plus$feature_type <- "utr3*" 
  #
  allWithFlank <- rbind(flanked.five.minus, flanked.three.plus, flanked.five.plus, flanked.three.minus)
  allWithFlank$entrezgene_id <- tx2gene[match(allWithFlank$gene_name, tx2gene$gene_name), 1]
  allWithFlank <- rbind(allWithFlank, hold.all)
  #
  allWithFlank <- GRanges(allWithFlank)
  return(allWithFlank)
}

#' FindOverlapsByRanges
#' Description: this function receives two granges and finds the overlapping region between the two 
#' @param flat1 the first grange
#' @param flat2 the second grange
#'
#' @return the overlapping region between them
FindOverlapByRanges <- function(flat1, flat2){
  #Transform as data table
  flat1D <- as.data.table(flat1)
  #Transform as data table
  flat2D <- as.data.table(flat2)
  flat1D.pos <- flat1D %>% dplyr::filter(strand == "+")
  #Transform as data table
  flat2D.pos <- flat2D %>% dplyr::filter(strand == "+")
  #sort the first df
  setkey(flat1D.pos, seqnames, start, end)
  #sort the second df
  setkey(flat2D.pos, seqnames, start, end)
  #find the overlapping intervals using foverlaps
  resultTable.pos <- foverlaps(flat1D.pos, flat2D.pos, nomatch = 0)
  #take the maximum start
  resultTable.pos[, start := pmax(start, i.start)]
  #take thew minimum end
  resultTable.pos[, end := pmin(end, i.end)]
  #remove unnecessary columns
  resultTable.pos[, `:=`(i.start = NULL, i.end = NULL)]
  #transform into df
  resultTable.pos <- as.data.frame(resultTable.pos)
  #remove more unnecessary columns
  resultTable.pos$i.width <- NULL; resultTable.pos$i.strand <- NULL; resultTable.pos$i.gene_name <- NULL; resultTable.pos$i.feature_type <- NULL;
  resultTable.pos$i.group <- NULL; resultTable.pos$i.group_name;
  #add entrezgene_id
  if("group_name" %in% colnames(resultTable.pos)){
    resultTable.pos$group_name -> resultTable.pos$gene_name
  }
  resultTable.pos$entrezgene_id <- tx2gene[match(resultTable.pos$gene_name, tx2gene$gene_name), 1]
  flat1D.ne <- flat1D %>% dplyr::filter(strand == "-")
  #Transform as data table
  flat2D.ne <- flat2D %>% dplyr::filter(strand == "-")
  #sort the first df
  setkey(flat1D.ne, seqnames, start, end)
  #sort the second df
  setkey(flat2D.ne, seqnames, start, end)
  #find the overlapping intervals using foverlaps
  resultTable.ne <- foverlaps(flat1D.ne, flat2D.ne, nomatch = 0)
  #take the maximum start
  resultTable.ne[, start := pmax(start, i.start)]
  #take thew minimum end
  resultTable.ne[, end := pmin(end, i.end)]
  #remove unnecessary columns
  resultTable.ne[, `:=`(i.start = NULL, i.end = NULL)]
  #transform into df
  resultTable.ne <- as.data.frame(resultTable.ne)
  #remove more unnecessary columns
  resultTable.ne$i.width <- NULL; resultTable.ne$i.strand <- NULL; resultTable.ne$i.gene_name <- NULL; resultTable.ne$i.feature_type <- NULL;
  resultTable.ne$i.group <- NULL; resultTable.ne$i.group_name;
  #add entrezgene_id
  if("group_name" %in% colnames(resultTable.ne)){
    resultTable.ne$group_name -> resultTable.ne$gene_name
  }
  resultTable.ne$entrezgene_id <- tx2gene[match(resultTable.ne$gene_name, tx2gene$gene_name), 1]
  #convert to gr
  resultTable <- rbind(resultTable.ne, resultTable.pos)
  resultTable <- GRanges(resultTable)
  return(resultTable)
}

#' SubtractOverlap
#' Description: This function receives two granges and subtracts the intervals
#' of one grange from the other 
#' @param flatBig The grange to subtract from
#' @param flatSmall The grange that will be subtracted
#' @param featureName the feature name assigned 
#'
#' @return a grange containing the difference between the two granges intervals
SubtractOverlap <- function(flatBig, flatSmall, featureName){
  #Transform the small grangelist to grange / in case it is a grl
  flatSmall.unlisted <- GRanges(as.data.frame(flatSmall, keep.extra.columns=TRUE))
  #Transform the other grangelist to grange / in case it is a grl
  flatBig.unlisted <- GRanges(as.data.frame(flatBig, keep.extra.columns=TRUE))
  #sort the small grangelist **important for matching with gene symbols later
  flatSmall.unlisted <- sort(flatSmall.unlisted)
  #sort the other grange
  flatBig.unlisted <- sort(flatBig.unlisted)
  #use setdiff to find the difference in intervals
  #flatDiff <- setdiff(flatBig.unlisted,FindOverlapByRanges(flatBig.unlisted, flatSmall.unlisted))
  flatDiff <- setdiff(GRanges(as.data.frame(flatBig.unlisted)),FindOverlapByRanges(flatBig.unlisted, flatSmall.unlisted))
  #sort the resulting grange
  flatDiff <- sort(flatDiff)
  #transform to dataframe
  flatBig.df <- as.data.frame(flatBig.unlisted, keep.extra.columns=TRUE)
  #match the meta data using findOverlaps
  flatsComb <- c(flatBig.unlisted, flatSmall.unlisted)
  #subjectHits <- unique(subjectHits(findOverlaps(flatDiff, flatsComb, type="any")))
  findOverlaps(flatDiff, flatsComb, type="any", select='first') -> idx
  #store the index
  #idx <- sort(subjectHits)
  #retrieve the gene_name
  flatDiff.df.new <- GRanges(flatBig.df)[idx]
  #store the feature type from original grange
  flatDiff.df<- as.data.table(flatDiff)
  flatDiff.df$gene_name <- as.data.table(flatDiff.df.new)$gene_name
  
  flatDiff.df$feature_type <- featureName
  #transform to grange
  flatDiff.df.res <- GRanges(flatDiff.df)
  return (flatDiff.df.res)
}

#' AddRDSInfo
#' Description: Receives a Granges object and adds the seqinfo to it
#' @param all granges object
#' @param genome.hold the name fo the genome
#'
#' @return the granges object along with seqinfo
AddRDSInfo <- function(all, genome.hold){
  #change the genome info
  genome(all) <- genome.hold
  #change the circular info
  isCircular.hold <- c(NA)
  #
  isCircular.hold[1:length(isCircular(all))] <- FALSE
  #
  isCircular(all) <- isCircular.hold
  #retrieve the chromosome names
  unique(levels(seqnames(all))) -> chr.vals
  #
  res <- c(NA)
  #get the length of each chromosome
  for(i in 1:(length(chr.vals))){
    hold.here <- GRanges(as.data.frame(all) %>% dplyr::filter(seqnames == chr.vals[i]))
    if(length(hold.here) == 0){
      res[i] <- 0
      next()
    }
    #store in res
    res[i]<- sum(width(hold.here))
  }
  #
  seqlengths(all)<- res
  return(all)
}

#' AddThreeUTRIndex
#' Description: The following function adds the 3 utr index to the granges object
#' @param all granges object with all annotations
#'
#' @return granges object with the 3 utr index added
AddThreeUTRIndex <- function(all){
  #dplyr::filter for the rows where intron.3utr and 3 utr exon is present
  all.df <- as.data.table(all)
  #set the index of all equal to 0 initially
  all.df$three.utr.index<- rep(0, nrow(all.df))
  #dplyr::filter for when intron.3utr or 3 utr is present
  all.df.hold <- all.df %>% dplyr::filter(feature_type %in% c("intron.3utr" ,"utr3"))
  #dplyr::filter for positive strand
  all.df.hold.pos <- all.df.hold %>% dplyr::filter(strand == "+")
  #sort in increasing order
  all.df.hold.pos  <- all.df.hold.pos[order(all.df.hold.pos$start, decreasing = TRUE),]
  #give index to the ones with intron.3utrs and on the + strand
  all.df.hold.pos <- all.df.hold.pos  %>% arrange(gene_name, start) %>% group_by(gene_name) %>% mutate(three.utr.index = rank(start, ties.method = "first"))
  #dplyr::filter for + strand again
  all.df.hold.pos.rest <- all.df %>% dplyr::filter(feature_type == "utr3*") %>% dplyr::filter(strand == "+")
  #take the entry on top
  matchToPos <- all.df.hold.pos[!duplicated(all.df.hold.pos$gene_name, fromLast = T),]
  #take the genes with gene name in the variable above
  all.df.hold.pos.rest <- all.df.hold.pos.rest[(all.df.hold.pos.rest$gene_name %in% matchToPos$gene_name)]
  #match the gene names and retrieve the 3utr index to add to the 3utr* 
  all.df.hold.pos.rest$three.utr.index <- matchToPos[match(all.df.hold.pos.rest$gene_name, matchToPos$gene_name), 9]
  #Do the same for the negative strand
  all.df.hold.ne <- all.df.hold %>% dplyr::filter(strand == "-")
  #this time sort in decreasing order
  all.df.hold.ne  <- all.df.hold.ne[order(all.df.hold.ne$start, decreasing = TRUE),]
  #put the ranks in the opposite way
  all.df.hold.ne <- all.df.hold.ne  %>% arrange(gene_name, start) %>% group_by(gene_name) %>% mutate(three.utr.index = rank(-start, ties.method = "first"))
  #dplyr::filter for utr3*
  all.df.hold.ne.rest <- all.df %>% dplyr::filter(feature_type == "utr3*") %>% dplyr::filter(strand == "-")
  #take the entry on top
  matchToNe <- all.df.hold.ne[!duplicated(all.df.hold.ne$gene_name),]
  #match the index for utr3 to utr3*
  all.df.hold.ne.rest <- all.df.hold.ne.rest[(all.df.hold.ne.rest$gene_name %in% matchToNe$gene_name)]
  #
  all.df.hold.ne.rest$three.utr.index <- matchToNe[match(all.df.hold.ne.rest$gene_name, matchToNe$gene_name),9]
  #bind all together
  allWithIndI <- rbind(all.df.hold.pos, all.df.hold.ne);allWithIndII <- rbind(all.df.hold.ne.rest, all.df.hold.pos.rest)
  #
  allWithInd <- rbind(as.data.frame(allWithIndI), as.data.frame(allWithIndII)) 
  #for all other utr3* regions set index to 1
  all.df.hold.rest <- all.df[!(all.df$gene_name %in% allWithInd$gene_name),]
  #
  all.df.hold.rest <- all.df.hold.rest  %>% dplyr::filter(feature_type == "utr3*")
  #
  all.df.hold.rest$three.utr.index <- 1 
  #bind all
  allIndexed <- rbind(allWithInd, all.df.hold.rest)
  #
  IndexMin <- all.df[!(all.df$start %in% allIndexed$start),]
  #bind and convert to granges
  result <- GRanges(rbind(IndexMin, allIndexed))
  return (result)
}
flankAll <- function(all){
  all<- AddPseudoFlank(all)
  #
  all <- AddGeneFlank(all)
  #
  all <- AddRestFlank(all)
  return (all)
}
#' ConsiderBlacklist
#' Description: The following fucntion finds the blacklisted regions and label them as 1 in a new blacklist column
#' @param all the granges without blacklist
#' @param filePathTobed the path to the blacklist bed file
#'
#' @return the entire granges
ConsiderBlacklist <- function(all, filePathTobed){
  blacklist <- import.bed(filePathTobed)
  blacklisted.regions <- subsetByOverlaps(blacklist, all)
  all$blacklist <- 0
  all[all %over% blacklisted.regions]$blacklist <- 1
  return (all)
}

#' getGeneOverlap
#' Description: the fucntion finds the overlap regions where there are a couple lonely genes overlapping in between
#' @param txdb a txdb object
#'
#' @return will return the overlap regions
getGeneOverlap <- function(txdb){
  exonsByTx <- exonsBy(txdb,by="tx",use.names=TRUE)
  #get CDS exons by transcript 
  cdsByTx <- cdsBy(txdb, by="tx",use.names=TRUE)
  exons_total <- exonsByTx
  #Flatten to genome level
  exons_total  <- flattenExonsBy(exonsByTx=exonsByTx, cdsByTx=cdsByTx, by="gene",genes=gene_names, tx2geneDF=tx2gene,verbose=FALSE)
  myGeneOL <- unique(exonsByTx@unlistData[-queryHits(findOverlaps(exonsByTx@unlistData, exons_total@unlistData, type="any")),])
  myGeneOL$exon_rank <- NULL; myGeneOL$exon_id <- NULL; myGeneOL$exon_rank <- NULL
  myGeneOL$feature_type <- "overlap"
  return(myGeneOL)
}


#' RankAll
#' Description: Rank every single feature type for every gene in the granges
#' @param all entire granges
#'
#' @return updated granges
RankAll <- function(all){
  #convert to dt
  all.df <- as.data.table(all)
  #set the index of all equal to 0 initially
  all.df$rank.all <- 0
  #dplyr::filter for positive strand
  all.df.hold.pos <- all.df %>% dplyr::filter(strand == "+")
  #sort in increasing order
  all.df.hold.pos  <- all.df.hold.pos[order(all.df.hold.pos$start, decreasing = TRUE),]
  #Rank
  all.df.hold.pos <- all.df.hold.pos  %>% arrange(gene_name, start) %>% group_by(gene_name, feature_type) %>% mutate(rank.all = rank(start, ties.method = "first"))
  #dplyr::filter for positive strand
  all.df.hold.ne <- all.df %>% dplyr::filter(strand == "-")
  #sort in increasing order
  all.df.hold.ne <- all.df.hold.ne[order(all.df.hold.ne$start, decreasing = TRUE),]
  #Rank
  all.df.hold.ne <- all.df.hold.ne  %>% arrange(gene_name, start) %>% group_by(gene_name, feature_type) %>% mutate(rank.all = rank(-start, ties.method = "first"))
  result <- rbind(all.df.hold.ne, all.df.hold.pos)
  GRanges(result)-> result
  return (result)
}

#' getGCPercentage
#'
#' @param all the entire granges
#' @param BSgenome 
#'
#' @return the entire granges with gc content percentage calculated
getGCPercentage <- function(all,BSgenome){
  all.df <- as.data.table(all)
  all.df$gc.content <- 0
  all.df <- all.df[all.df$seqnames %in% seqnames(BSgenome)]
  for(i in 1:length(all.df$seqnames)){
    if(all.df$seqnames[i] %in% seqnames(BSgenome)){
    hold.me.seq <- getSeq(BSgenome, names = all.df$seqnames[i], start=all.df$start[i], end=all.df$end[i], strand=all.df$strand[i])
    all.df$gc.content[i] <- (str_count(hold.me.seq, "C") + str_count(hold.me.seq, "G"))/length(hold.me.seq)
    }
  }
  return(Granges(all.df))

}
#' FlankUpAndFlankDown
#' Description: The following function finds the feature type before and after every range.
#' For the negative strand looks at the range after the range. for the positive range it finds the range before it stores it in -> before
#' @param all granges
#'
#' @return all granges
FlankUpAndFlankDown <- function(all){
  as.data.table(all)-> all.df
  all.df %>% dplyr::filter(strand == "-") -> all.dt.ne
  all.df %>% dplyr::filter(strand == "+") -> all.dt.pos
  all.dt.pos[order(all.dt.pos$seqnames, all.dt.pos$start),] -> all.dt.pos
  all.dt.ne[order(all.dt.ne$seqnames, all.dt.ne$start),] -> all.dt.ne
  mutate(all.dt.ne, before = lead(feature_type)) -> all.dt.ne
  mutate(all.dt.pos, before = lag(feature_type)) -> all.dt.pos
  mutate(all.dt.ne, after = lag(feature_type)) -> all.dt.ne
  mutate(all.dt.pos, after = lead(feature_type)) -> all.dt.pos
  rbind(all.dt.ne, all.dt.pos)-> all       
  return(GRanges(all))
}

getProteinCoding<- function(all, tx2gene){
  
  all$biotype<-c(NA)
  protein_coding <- unique(tx2gene[tx2gene$biotype == "protein_coding",]$gene_id)
  all[all$entrezgene_id %in% protein_coding]$biotype="protein_coding" 
  all[!all$entrezgene_id %in% protein_coding]$biotype="non_protein_coding"
  return(all)
}
# ExtractRetainedIntrons <- function(){
#   myMart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#   library(EnsDb.Hsapiens.v86)
#   edb <- EnsDb.Hsapiens.v86       
#   edb.df <- as.data.frame(AnnotationDbi::select(edb, keys=gene_names, keytype = "GENENAME",columns = c("GENEID", "GENENAME", "TXID","TXBIOTYPE", "SEQNAME", "TXSEQSTART","TXSEQEND")))
#   protein_coding<- unique(edb.df[!edb.df$TXBIOTYPE %in% c("processed_transcript","miRNA", "snoRNA", "scRNA", "scaRNA", "lincRNA", "LRG_gene", "snRNA"),]$GENENAME)
#   edb.df[edb.df$TXBIOTYPE=="retained_intron",]
#   
# }
#'

#   ````````````````````````````````````````````   
#   /dsoooooooooooooooooooooooooooooooooooooydd-   
#   /mdo///////////////////////////////////sh+d-   
#   /mmmmmmmmmmmmmmmmmmmmmo:hmmmmmmmmmmmmmmm.`d-   
#   /mmmmmmm+/++++++++hmmmo oy+++++++++smmmm.`d.   
#   +mmmmmmm.        `ymmmo os`        /mmmm+.d.   
#   /hhhhhhh.        `ymmmo os`        :hhhhhyh.   
#   ````````         `ymmmo os`        ````````    
#                    `ymmmo os`                    
#      -ooooo/`      `ymmmo os` `+ooo+`    `/oooo. 
#      :hmmmdo`      `ymmmo os` .sdmmm+`   /dmmdy. 
#      .ymhmd/`      `ymmms os`  `ymddd+` :dmhmd.  
#     `smd:omd:      `ymmms os`  `ymy+md/:dmsomd.  
#    `ommo//hmh-     `ymmms os`  `ymy.omdhmy.omd.  
#    /mmdddhddmy.    `ymmms os`  `ymy``smmy.`omd.  
# ``:dms-....:dms``  `ymmms os`  .ymy. .sh- `omd.  
# :ydmmh.   `+dmmys` `ymmms +s` .sdmds. .- `ohmmy- 
# -+++++.   `:++++/` `ymmms os` .+ooo+`    `/oooo- 
#                    `ymmms os`                    
#                ````.hmmms os.````                
#                `+hhhhdmmms odhhyy/                
#                `ommmmmmmmhymmmy.y+                
#                `smmmmmmmmmmmmms.y+                
#                `oddddddddddddddhh/                
#                ``````````````````-+   
#Spring 2021#####################################################################################      
###############################################################################################
SubtractOverlapIgnoreStrand <- function(flatBig, flatSmall, featureName){
  #Transform the small grangelist to grange / in case it is a grl
  flatSmall.unlisted <- GRanges(as.data.frame(flatSmall, keep.extra.columns=TRUE))
  #Transform the other grangelist to grange / in case it is a grl
  flatBig.unlisted <- GRanges(as.data.frame(flatBig, keep.extra.columns=TRUE))
  #sort the small grangelist **important for matching with gene symbols later
  flatSmall.unlisted <- sort(flatSmall.unlisted)
  #sort the other grange
  flatBig.unlisted <- sort(flatBig.unlisted)
  #use setdiff to find the difference in intervals
  #flatDiff <- setdiff(flatBig.unlisted,FindOverlapByRanges(flatBig.unlisted, flatSmall.unlisted))
  as.data.table(flatSmall.unlisted) %>% dplyr::filter(strand=="+")->flatSmall.unlisted.pos
  flatSmall.unlisted.pos$strand<-"-"
  as.data.table(flatSmall.unlisted) %>% dplyr::filter(strand=="-")->flatSmall.unlisted.ne
  flatSmall.unlisted.pos$strand<-"+"
  GRanges(rbind(flatSmall.unlisted.pos,flatSmall.unlisted.ne))->flatSmall.unlisted
  
  flatDiff <- setdiff(GRanges(as.data.frame(flatBig.unlisted)),FindOverlapByRanges(flatBig.unlisted, flatSmall.unlisted))
  #sort the resulting grange
  flatDiff <- sort(flatDiff)
  #transform to dataframe
  flatBig.df <- as.data.frame(flatBig.unlisted, keep.extra.columns=TRUE)
  #match the meta data using findOverlaps
  flatsComb <- c(flatBig.unlisted, flatSmall.unlisted)
  #subjectHits <- unique(subjectHits(findOverlaps(flatDiff, flatsComb, type="any")))
  findOverlaps(flatDiff,flatsComb,  type="any", select='first', ignore.strand=TRUE) -> idx
  #store the index
  #idx <- sort(subjectHits)
  #retrieve the gene_name
  flatDiff.df.new <- GRanges(flatBig.df)[idx]
  #store the feature type from original grange
  flatDiff.df<- as.data.table(flatDiff)
  flatDiff.df$gene_name <- as.data.table(flatDiff.df.new)$gene_name
  
  flatDiff.df$feature_type <- featureName
  #transform to grange
  flatDiff.df.res <- GRanges(flatDiff.df)
  return (flatDiff.df.res)
}
FindOverlapByRangesIgnoreStrand <- function(flat1, flat2){
  #Transform as data table
  flat1D <- as.data.table(flat1)
  #Transform as data table
  flat2D <- as.data.table(flat2)
  #sort the first df
  setkey(flat1D, seqnames, start, end)
  #sort the second df
  setkey(flat2D, seqnames, start, end)
  #find the overlapping intervals using foverlaps
  resultTable <- foverlaps(flat1D, flat2D, nomatch = 0)
  #take the maximum start
  resultTable[, start := pmax(start, i.start)]
  #take thew minimum end
  resultTable[, end := pmin(end, i.end)]
  #remove unnecessary columns
  resultTable[, `:=`(i.start = NULL, i.end = NULL)]
  #transform into df
  resultTable <- as.data.frame(resultTable)
  #remove more unnecessary columns
  resultTable$i.width <- NULL; resultTable$i.strand <- NULL; resultTable$i.gene_name <- NULL; resultTable$i.feature_type <- NULL;
  resultTable$i.group <- NULL; resultTable$i.group_name;
  #add entrezgene_id
  if("group_name" %in% colnames(resultTable)){
    resultTable$group_name -> resultTable$gene_name
  }
  resultTable$entrezgene_id <- tx2gene[match(resultTable$gene_name, tx2gene$gene_name), 1]
  resultTable <- GRanges(resultTable)
  return(resultTable)
}
ExtractIntrons<- function(ref){
  ref <- as.data.frame(ref)
  ref %>% dplyr::filter(strand=="+")-> ref.df.pos
  ref %>% dplyr::filter(strand=="-")-> ref.df.ne
  mutate(ref.df.pos, NextExonEnd = lead(end))-> ref.df.pos
  mutate(ref.df.ne, NextExonEnd = lag(start))-> ref.df.ne
  rbind(ref.df.ne, ref.df.pos)-> ref
  #dplyr::filter for introns
  ref %>% dplyr::filter(feature_type == "intron") -> ref1
  #add the intron.3utrs
  ref %>% dplyr::filter(feature_type == "intron.3utr") -> ref2
  #bind the two
  rbind(ref1, ref2)-> ref
  ref %>% dplyr::filter(sub(".*chr", "", ref$seqnames) %in% c(1:22, 'X', 'Y')) -> ref
  #dplyr::filter for regions that don't fall in blacklists
  ref <- ref %>% dplyr::filter(blacklist==0)
  #dplyr::filter to not go over these twice
  #ref[!(ref$feature_type=="intron" && ref$before=="intron.3utr"),]-> ref
  ref[!duplicated(ref),]->ref
  ref <- filterStrandNoise(ref, all)
  ref$exon_name<-NULL
  GRanges(as.data.frame(ref) %>% dplyr::filter(biotype=="protein_coding"))->ref
  return(ref)
}



filterStrandNoise<- function(IntronicRegions, rdsSource){
  #filter the annotation source/ remove introns
  as.data.table(rdsSource) %>% dplyr::filter(feature_type!="intron")-> anno1
  #remove intron.3utr
  as.data.table(anno1) %>% dplyr::filter(feature_type!="intron.3utr")-> anno2
  #remove intergenic
  as.data.table(anno2) %>% dplyr::filter(feature_type!="intergenic")-> anno
  #now get the dplyr::filtered intronic regions
  as.data.table(IntronicRegions)-> IntronicRegions
  #subtract the strand overlapping coverage with introns of the other strand
  findOverlaps(GRanges(anno), GRanges(ref), ignore.strand=TRUE)->idx
  GRanges(ref)[subjectHits(idx)]-> forbiddenIntrons
  IntronicRegions %>% dplyr::filter(seqnames %in% as.data.frame(forbiddenIntrons)$seqnames) %>% dplyr::filter(!start %in% as.data.frame(forbiddenIntrons)$start)->filtered_intronicRegions
  return(filtered_intronicRegions)
}

#'
#'
#'
#'
#'
#'
#'
#'
#'
#Functions below this line were taken from other packages###################################### 
#####Function was taken from the Splicejam/jamba package and modified for introns
flattenIntronsBy <- function(intronsByTx,
                             tx2geneDF,
                             by=c("gene", "tx"),
                             detectedTx=NULL,
                             genes=NULL,
                             txColname="transcript_id",
                             geneColname="gene_name",
                             cdsByTx=NULL,
                             cdsByGene=NULL,
                             verbose=FALSE)
{
  ##
  by <- match.arg(by);
  # take the union of the two grangeslist
  iTxs <- unique(c(names(intronsByTx)))
  ## Optionally subset tx2geneDF by genes //// UNIMPORTANT STEP
  if (length(genes) > 0) {
    tx2geneDF <- subset(tx2geneDF,
                        tx2geneDF[[geneColname]] %in% genes);
    if (nrow(tx2geneDF) == 0) {
      stop("tx2geneDF[[geneColname]] contains no values matching the supplied genes.");
    }
  }
  ## Validate iTxs in tx2geneDF //// UNIMPORTANT STEP just validating that transcripts have associated genes
  iTxs <- intersect(iTxs,tx2geneDF[[txColname]]);
  tx2geneDF <- subset(tx2geneDF,tx2geneDF[[txColname]] %in% iTxs);
  if (length(intronsByTx) > 0) {
    intronsByTx <- intronsByTx[names(intronsByTx) %in% tx2geneDF[[txColname]]];
  }
  if (length(cdsByTx) > 0) {
    cdsByTx <- cdsByTx[names(cdsByTx) %in% tx2geneDF[[txColname]]];
  }
  
  ## Subset intronsByTx and add gene annotations
  iTxintronsGRL <- intronsByTx[iTxs];
  iTxMatch <- match(names(iTxintronsGRL),
                    tx2geneDF[[txColname]]);
  GenomicRanges::values(iTxintronsGRL@unlistData)[,geneColname] <- rep(
    as.character(tx2geneDF[iTxMatch,geneColname]),
    S4Vectors::elementNROWS(iTxintronsGRL));
  
  ## split introns by gene
  
  if ("gene" %in% by) {
    intronsByGene <- GenomicRanges::GRangesList(
      GenomicRanges::split(
        iTxintronsGRL@unlistData,
        GenomicRanges::values(iTxintronsGRL@unlistData)[[geneColname]])
    );
  } else {
    GenomicRanges::values(iTxintronsGRL@unlistData)[,txColname] <- rep(
      names(iTxintronsGRL),
      S4Vectors::elementNROWS(iTxintronsGRL));
    intronsByGene <- iTxintronsGRL[,c(txColname,geneColname)];
  }
  ## Disjoin introns within each gene GRL // EXCTRACT the non-overlapping intron parts
  if ("gene" %in% by) {
    
    iGeneintronsDisGRL <- GenomicRanges::disjoin(intronsByGene);
  } else {
    iGeneintronsDisGRL <- intronsByGene;
  }
  ## Add gene annotation to each entry
  if ("gene" %in% by) {
    GenomicRanges::values(iGeneintronsDisGRL@unlistData)[,geneColname] <- rep(
      names(iGeneintronsDisGRL),
      elementNROWS(iGeneintronsDisGRL));
  }
  ## Assign intron names and numbers
  if ("gene" %in% by) {
    iGeneintronsDisGRL <- assignGRLintronNames(iGeneintronsDisGRL,
                                               geneSymbolColname=geneColname,
                                               verbose=FALSE);
    GenomicRanges::values(iGeneintronsDisGRL)[,geneColname] <- "intron"
    #names(iGeneintronsDisGRL);
  } else {
    iGeneintronsDisGRL <- assignGRLintronNames(iGeneintronsDisGRL,
                                               geneSymbolColname=txColname,
                                               verbose=FALSE);
    GenomicRanges::values(iGeneintronsDisGRL)[,txColname] <- names(iGeneintronsDisGRL);
    GenomicRanges::values(iGeneintronsDisGRL@unlistData)[,txColname] <- rep(
      names(iGeneintronsDisGRL),
      S4Vectors::elementNROWS(iGeneintronsDisGRL)
    )
    txMatch <- match(names(iGeneintronsDisGRL),
                     tx2geneDF[[txColname]]);
    GenomicRanges::values(iGeneintronsDisGRL)[,geneColname] <- tx2geneDF[txMatch, geneColname];
    GenomicRanges::values(iGeneintronsDisGRL@unlistData)[,geneColname] <- rep(
      GenomicRanges::values(iGeneintronsDisGRL)[,geneColname],
      S4Vectors::elementNROWS(iGeneintronsDisGRL)
    )
  }
  GenomicRanges::values(iGeneintronsDisGRL@unlistData)[,"feature_type"] <- "intron";
  return(iGeneintronsDisGRL);
}




assignGRLintronNames <- function(GRL,
                                 geneSymbolColname="geneSymbol",
                                 intronNameColname=paste0(geneSymbolColname, "Intron"),
                                 suffix="_intron",
                                 renameOnes=FALSE,
                                 filterTwoStrand=TRUE,
                                 checkDisjoin=c("warn","none","stop"),
                                 assignGRLnames=TRUE,
                                 verbose=FALSE,
                                 ...)
{
  ## Purpose is to assign intron names using numbers
  ## to represent contiguous segments, and lowercase
  ## letters to represent subsections of each intron.
  ##
  ## filterTwoStrand=TRUE will remove entries which have two strands
  ## for the same GRL entry
  ##
  ## This function is a light wrapper for renumberGRanges()
  ##
  ## checkDisjoin="stop" will check to make sure introns in each set
  ## of GRanges are disjoint, otherwise the numbering can be
  ## problematic.
  ## The most common symptom is negative strand embedded introns receive
  ## sub-numbers in opposite order.
  ##
  checkDisjoin <- match.arg(checkDisjoin);
  
  ## First verify that incoming data is valid per assumptions
  ## that introns for a transcript would all appear only on one strand
  GRLstrandL <- unique(GenomicRanges::strand(GRL));
  if (filterTwoStrand && any(S4Vectors::elementNROWS(GRLstrandL) > 1)) {
    iRemove <- which(S4Vectors::elementNROWS(GRLstrandL) > 1);
    GRL <- GRL[-iRemove];
  } 
  ## check disjoint GRanges
  if (checkDisjoin %in% c("warn","stop")) {
    
    GRLdis <- GenomicRanges::disjoin(GRL);
    if (!all(S4Vectors::elementNROWS(GRLdis) == S4Vectors::elementNROWS(GRL))) {
      if (checkDisjoin %in% "stop") {
        stop("assignGRLintronNames() detected overlapping GRanges, stopping.");
      } else {
        printDebug("assignGRLintronNames(): ",
                   "detected overlapping GRanges, continuing.",
                   fgText=c("red","orange"));
      }
    }
  }
  
  ## Reduce entries
  
  GRLred <- GenomicRanges::reduce(GRL);
  
  ## Add geneSymbolColname if it does not already exist
  if (!geneSymbolColname %in% colnames(GenomicRanges::values(GRLred@unlistData))) {
    GenomicRanges::values(GRLred@unlistData)[,geneSymbolColname] <- rep(names(GRLred),
                                                                        S4Vectors::elementNROWS(GRLred));
  }
  
  GRLredStrand <- unlist(unique(strand(GRLred)));
  GRLredStrandP <- which(GRLredStrand %in% "+");
  GRLredStrandN <- which(GRLredStrand %in% "-");
  
  ## Stranded intron numbering
  GenomicRanges::values(GRLred@unlistData)[,intronNameColname] <- "";
  
  if (length(GRLredStrandP) > 0) {
    GenomicRanges::values(GRLred[GRLredStrandP]@unlistData)[,intronNameColname] <- jamba::makeNames(
      GenomicRanges::values(GRLred[GRLredStrandP]@unlistData)[,geneSymbolColname],
      suffix=suffix,
      renameOnes=TRUE);
  }
  if (length(GRLredStrandN) > 0) {
    GenomicRanges::values(GRLred[GRLredStrandN]@unlistData)[,intronNameColname] <- rev(jamba::makeNames(
      GenomicRanges::values(GRLred[rev(GRLredStrandN)]@unlistData)[,geneSymbolColname],
      suffix=suffix,
      renameOnes=TRUE));
  }
  
  
  ## Add lowercase letter suffix
  GRLcolnames <- unvigrep(
    paste0(intronNameColname
           ,"(_v[0-9]|)$"
    ),
    colnames(GenomicRanges::values(GRL@unlistData)));
  
  GRLnew <- annotateGRLfromGRL(GRL1=GRL[,GRLcolnames],
                               GRL2=GRLred[,intronNameColname],
                               verbose=verbose);
  
  GRLnewStrand <- unlist(unique(strand(GRLnew)));
  GRLnewStrandP <- which(GRLnewStrand %in% "+");
  GRLnewStrandN <- which(GRLnewStrand %in% "-");
  GRLnewStrandNn <- names(GRLnew[GRLnewStrandN]@unlistData);
  subFeatureNumberStyle <- "letters";
  subFeatureSuffix <- "";
  intronNameColname1 <- paste0(intronNameColname, "1");
  
  GenomicRanges::values(GRLnew@unlistData)[,intronNameColname] <-
    GenomicRanges::values(GRLnew@unlistData)[,intronNameColname];
  GenomicRanges::values(GRLnew[GRLnewStrandP]@unlistData)[,intronNameColname] <- (
    jamba::makeNames(
      GenomicRanges::values(GRLnew[GRLnewStrandP]@unlistData)[,intronNameColname],
      numberStyle=subFeatureNumberStyle,
      suffix=subFeatureSuffix,
      renameOnes=renameOnes));
  GenomicRanges::values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,intronNameColname] <- (
    jamba::makeNames(
      GenomicRanges::values(GRLnew@unlistData[rev(GRLnewStrandNn),])[,intronNameColname],
      numberStyle=subFeatureNumberStyle,
      suffix=subFeatureSuffix,
      renameOnes=renameOnes));
  
  if (assignGRLnames) {
    names(GRLnew@unlistData) <- jamba::makeNames(
      GenomicRanges::values(GRLnew@unlistData)[,intronNameColname]);
  }
  
  
  return(GRLnew);
  
}


                 
                 
               
               
               
