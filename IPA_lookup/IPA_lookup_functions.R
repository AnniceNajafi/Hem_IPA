#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Texas A&M University, Department of Biomedical Engineering
#' Spring 2021
#' Description: This program detects IPA events using the KMMD method

#load the related libraries
library(Rsamtools)
library(edgeR)
library(GenomicAlignments)
library(GenomicFeatures)
library(tidyverse)
library(dplyr)
library(data.table)
library(rtracklayer)
library(foreach)
library(doParallel)
library(rslurm)
#libraries from previous script needed
library(org.Hs.eg.db)
library(jamba)

#' runFilter
#'
#' @param input.data.path 
#'
#' @return
#' description: The following function receives the 
runChangePoint<- function(input.data.path, wd){
  coverage.dir <- "coverage_data"
  if(!dir.exists(coverage.dir)){
  #make a new directory to store the coverage data
  dir.create(coverage.dir)
  }
  if(!dir.exists('final_result')){
             dir.create('final_result')
   }
  if(!dir.exists('slurm_submission')){
  	    dir.create('slurm_submission')
  }
  #Read the excel file which contains the sample name and the location of the bam file
  data.input <- read.delim(input.data.path, sep="\t", header=T)
  #save the sample name
  sampleNames <- as.character(data.input$sampleName)
  #For each sample:
  sapply(sampleNames, function(sample){
    #Retrieve the annotation
    anno.loc <- as.character(subset(data.input, sampleName==sample)$rdsSource)
    #If the annotation file is absent show an error message
    if(!file.exists(anno.loc)){
      print("annotation file is missing!")
    }
    #otherwise read the related annotation file
    rdsSource <- readRDS(anno.loc)
    #store the path for bam file
    bam.loc <- as.character(subset(data.input, sampleName==sample)$bamPath)
    #store the bam file directory
    bam.dir <- dirname(bam.loc)
    #in case bam file is missing show an error message
    if(!file.exists(bam.loc)){
      print("ERROR! Bam file does not exist!")
    }
    #unmapped reads bam file loc
    unmapped.bam.loc <- paste0(bam.dir, "/",sample, "_unmapped.bam")
    #uniquely mapped reads bam file loc
    uniq.bam.loc <- paste0(bam.dir, "/",sample, "_uniq.bam")
    #if the unmapped bam file does not exist
    if(!file.exists(unmapped.bam.loc)){
    #filter for unmapped reads
    system(paste0("samtools view -b -f 4 ", bam.loc, ">", bam.dir, "/",sample, "_unmapped.bam"))
    #create the index file for the bam file
    system(paste0("samtools index ", bam.dir, "/",sample,"_unmapped.bam ",bam.dir,"/", sample, "_unmapped.bam.bai"))
    }
    #if the uniquely mapped bam file does not exist
    if(!file.exists(uniq.bam.loc)){
    #filter for uniquely mapped reads from bam file
    system(paste0("samtools view -q 10 -b ", bam.loc, ">", bam.dir, "/", sample, "_uniq.bam"))
    #create the index file for uniquely mapped reads
    system(paste0("samtools index ", bam.dir, "/",sample, "_uniq.bam ",bam.dir, "/", sample, "_uniq.bam.bai"))
    }
    #print bam file made if the above steps pass
    print("New bam files were made.")
    #make a directory within coverage with the sample's name to store the coverage object
    if(!dir.exists(paste0(file.path("coverage_data/"), sample))){
      #the directory should start with the sample name
      dir.create(paste0("coverage_data/",sample))
    }
    ref.anno <- as.character(subset(data.input, sampleName==sample)$referenceAnnotation)
    ref.anno<- readRDS(ref.anno)
    readGAlignments(uniq.bam.loc)->bam.reads
    findTPM(ref.anno, bam.reads)->tpm.ls
    if(!dir.exists(paste0("annotation/",sample))){
      dir.create(paste0("annotation/",sample))
    }
    saveRDS(tpm.ls, paste0("annotation/",sample, "/tpm.RDS"))
    
    ExtractJunctions(bam.reads)->junc.df
    saveRDS(junc.df, paste0("annotation/", sample,"/junc.RDS"))
    #store the annotation file in a variable to be passed to function
    as.data.table(rdsSource)-> hold
    #combine with the bam file address and sample name for submission to function
    ref<- hold
    #coverage column contains the path of the coverage file
    ref$coverage<- file.path(paste0(wd, "/coverage_data/",sample,"/",  "cvg.Rda"))
    if(length(list.files(file.path(paste0(wd, "/coverage_data/",sample,"/",  "cvg.Rda"))))==0){
     coverage(bam.reads)->cov
     saveRDS(cov, file.path(paste0(wd, "/coverage_data/",sample,"/",  "cvg.Rda")))
    }
    #submit to slurm
    #break the job over chromosomes
    chrs <- unique(ref$seqnames)
    #for every chromosome apply the function
    save("goOver", file="slurm_submission/goOver.Rdata")
    if(!dir.exists(paste0("final_result/", sample))){
          dir.create(paste0("final_result/", sample))
      }
    if(!dir.exists(file.path(paste0("slurm_submission/", sample)))){
          dir.create(paste0("slurm_submission/", sample))
      }
    if(!dir.exists(file.path('logs'))){
          dir.create('logs')
    }
    num_groups<-30
    ref <- ref %>% 
   group_by((row_number()-1) %/% (n()/num_groups)) %>%
   nest %>% pull(data)
    for (i in 1:length(ref)){
      #Get the chromosome names from GG.locs
      #GG.locs <- ref %>% filter(seqnames==chrs[i])
      GG.locs <- ref[[i]]
      GG.locs$output <- file.path(paste0("final_result", '/', sample))
      #make a new folder for submission
      results.sample.dir <- file.path(paste0("slurm_submission/", sample, "/"))
      script.name <- file.path(paste0(results.sample.dir,'/',i, "_goOverRun.R"))
      print('ready to make GG.locs Rdata objects')
      #save GG.locs as rdata object
      save(GG.locs, file=paste0(results.sample.dir, sample,"_", i,"GGlocs.Rdata"))
      #create a new script
      sink(file=script.name)
      cat("
        \nlibrary(Rsamtools)
        \nlibrary(edgeR)
        \nlibrary(GenomicAlignments)
        \nlibrary(GenomicFeatures)
        \nlibrary(tidyverse)
        \nlibrary(dplyr)
        \nlibrary(data.table)
        \nlibrary(rtracklayer)
        \nlibrary(foreach)
        \nlibrary(doParallel)
        \nlibrary(org.Hs.eg.db)
        \nlibrary(jamba)
        \nlibrary(DESeq)
        \nlibrary(zoo)
        \nlibrary(kernlab)
        \nlibrary(dgof)
        \nlibrary(DESeq2)")
      #paste in goOver
      cat(paste0("\nname<-",i))
      cat(paste0("\nsample<-'",sample,"'"))
      cat(paste0("\nload (\'slurm_submission/goOver.Rdata\')"))
      #load the previously saved GG.locs object
      cat(paste0("\nload (\'", results.sample.dir, sample, "_", i, "GGlocs.Rdata\')"))
      #Read the tpm values store in variable
      cat(paste0('\ntpm.res <- readRDS(file.path(paste0("annotation/",sample,"/tpm.RDS")))'))
      #Read the junctions file
      cat(paste0('\njunc.df <- readRDS(file.path(paste0("annotation/",sample, "/junc.RDS")))'))
      #call the function
      cat(paste0("\ngoOver(GG.locs)"))
      #close
      sink()
      #create the bash script
      bash.file.location <- file.path(paste0(results.sample.dir, i, 'submit.sh'))
      slurm.jobname <- sprintf("#SBATCH --job-name=%s_IPA_Final", sample)
      slurm.time <- sprintf("#SBATCH --time=24:00:00")
      slurm.mem <- sprintf("#SBATCH --mem=12G")
      slurm.tasks <- sprintf("#SBATCH --ntasks=8")
      slurm.output <- sprintf("#SBATCH --output=%s%s_IPA_Final_%s",'logs/', sample, i)
      sbatch.line<- sprintf("/scratch/group/isinghlab/usr/local/lib64/R/bin/Rscript --vanilla %s", script.name)
      file.conn <- file(bash.file.location)
      writeLines(c("#!/bin/bash", slurm.jobname, slurm.time, slurm.mem, slurm.tasks, slurm.output, sbatch.line), file.conn)
      close(file.conn)
      system(paste0("sbatch ", bash.file.location))
      
    }
    
    
    
    print(paste0("Running algorithm... ", "job for ", chrs[i], " submitted."))
    print("job submitted")
    #stop for a minute...
    Sys.sleep(30)
    #First check if previous job has been completed if not postpone for fifteen minutes
    while(grepl("IPA", system(paste0("squeue -u ", user), intern = TRUE)[2])){
       Sys.sleep(900)
     }
     list.files(file.path(paste0('final_result/', sample)))-> final.file.paths
     if(is.null(final.file.paths)){
         #if not show an error message
         print("ERROR! No files were found from previous step in directory... Function DID NOT Work...")
       }else{
         res <- do.call("rbind", lapply(final.file.paths, readRDS))
         write.csv(res, paste0(wd, "/final_result/finals.csv"))
         saveRDS(res, paste0(wd, "/final_result/finals.RDS"))
       }

  })
  
  
}


#' Description: This function goes over introns and intron.3utrs and finds the maximum difference in change in mean
#' 
#' @param cvg 
#' @param GG.locs 
#'
#' @return the dataframe with changes in means
goOver <- function(GG.locs){
  #define a shift function 
  shift_seq <- function(v, s) {c(tail(v, s), head(v, -s))}
  #print statement ensure the function is running
  print("running function")
  #filter sequences with length exceeding 502 nts
  GG.locs %>% filter(width >800) -> GG.locs
  #get the path to the coverage file
  GG.locs$coverage[1]-> path
  #save the path 
  file.path(path)->path
  #read coverage from rds file 
  cvg <- readRDS(path)
  #convert to datatable
  as.data.table(GG.locs)-> GG.locs
  #filter for negative strand
  GG.locs %>% filter(strand == "-") -> ref38.dt.ne
  #filter for positive strand
  GG.locs %>% filter(strand == "+") -> ref38.dt.pos
  #sort the positive and negative strands
  as.data.table(ref38.dt.pos)[order(seqnames, start),] -> ref38.dt.pos
  #
  as.data.table(ref38.dt.ne)[order(seqnames, start),] -> ref38.dt.ne
  #add the starting point for intron.utr
  mutate(ref38.dt.ne, lagStart = lead(end)) -> ref38.dt.ne
  #
  mutate(ref38.dt.pos, lagStart = lag(start)) -> ref38.dt.pos

  mutate(ref38.dt.ne, lagEnd = lag(start)) -> ref38.dt.ne
  mutate(ref38.dt.pos, lagEnd = lead(end)) -> ref38.dt.pos

  #combine the two
  rbind(ref38.dt.pos, ref38.dt.ne) -> GG.locs
  #for every intron or intron.3utr region
  tpm.res %>% filter(TPM>1)->tpm.res
  GG.locs %>% filter(gene_name %in% tpm.res$gene_name)->GG.locs
  
   checkCigarJunc<- function(chr.hold, x.loc , strand.hold, junc.df){
     
     junc.df %>% distinct(seqnames, start, end, .keep_all=TRUE)->junc.hold
     junc.hold %>% filter(seqnames==chr.hold)-> junc.hold
     junc.hold %>% filter(strand=="+")->junc.hold.positive
     junc.hold %>% filter(strand=="-")->junc.hold.negative
     
     junc.hold.positive -> junc.hold.positive.start
     junc.hold.positive -> junc.hold.positive.end
     #store start and add 200 nt to start
     junc.hold.positive.start$end<-junc.hold.positive.start$start+100
     #store end and subtract 200nt
     junc.hold.positive.end$start<-junc.hold.positive.end$end-100
     #
     junc.hold.negative -> junc.hold.negative.start
     junc.hold.negative -> junc.hold.negative.end
     #store
     junc.hold.negative.start$end<-junc.hold.negative.start$start
     junc.hold.negative.start$start<- junc.hold.negative.start$start-100
     junc.hold.negative.end$start <- junc.hold.negative.end$end
     junc.hold.negative.end$end <- junc.hold.negative.end$end + 100
     rbind(junc.hold.negative.end, junc.hold.negative.start, junc.hold.positive.start, junc.hold.positive.end)->junc.all
     
     GRanges(Rle(chr.hold), IRanges(as.numeric(x.loc-100), as.numeric(x.loc+100)), Rle(strand.hold))->GR
     
     if(length(findOverlaps(GR, GRanges(junc.all), type="any", ignore.strand=TRUE))!=0){
       isOverlapping<- TRUE
     }else{
       isOverlapping<- FALSE
     }
     
     return(isOverlapping)
     
   }
   shift_seq <- function(v, s) {c(tail(v, s), head(v, -s))}
     result<-list()
     KMDtest<-list()
     ChangePt.stat<-list()
     xx<-list()
     foreach (i = 1:length(GG.locs$start)) %do% {
       print(i)
       #if intron.3utr after an intron then add the preceeding intronic region
       if(GG.locs$feature_type[i] == "intron.3utr" && GG.locs$before[i] == "intron" && GG.locs$after[i] == "intron"){
         #check the strand and find the start and end points according to the strand
         hold.me <- NULL
         if(GG.locs$strand[i] == "+"){
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$lagStart[i];end.hold<-GG.locs$lagEnd[i];
         }
         else{
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$lagEnd[i];end.hold<-GG.locs$lagStart[i];
         }
         
         
       }
       else if(GG.locs$feature_type[i] == "intron.3utr" && GG.locs$before[i] == "intron"){
         #check the strand and find the start and end points according to the strand
         hold.me <- NULL
         if(GG.locs$strand[i] == "+"){
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$lagStart[i];end.hold<-GG.locs$end[i];
         }
         else{
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$start[i];end.hold<-GG.locs$lagStart[i];
         }
         
         
       }else if(GG.locs$feature_type[i] == "intron.3utr" && GG.locs$after[i] == "intron"){
         #check the strand and find the start and end points according to the strand
         hold.me <- NULL
         if(GG.locs$strand[i] == "+"){
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$start[i];end.hold<-GG.locs$lagEnd[i];
         }
         else{
           chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$lagEnd[i];end.hold<-GG.locs$end[i];
         }
         
       }
       
       else{
         #instantiate a dataframe
         hold.me <- NULL
         #
         chr.hold<-GG.locs$seqnames[i];start.hold<-GG.locs$start[i];end.hold<-GG.locs$end[i];
         
       }
       if(!is.na(start.hold) && !is.na(end.hold)){
       #retrieve the read counts for the range
       cvg.hold <- as.numeric(cvg[GRanges(Rle(chr.hold), IRanges(as.numeric(start.hold), as.numeric(end.hold)))][[1]])
       
       }
       #have to unfactor to avoid errors
       chr.hold<-unfactor(chr.hold)
       #find the span of contiguous coverage exceeding a certain number of reads
       vals<- (with(rle(cvg.hold >=16), lengths[values]))
       #find the location where this is happening
       num <- which.max(((with(rle(cvg.hold >=16), lengths[values]))))
       #if contiguous reads do not exist say the span is 0 at location 1
       if(length(vals)==0){
         vals<-0
         num<-1
       }
       #find where this contiguity stops
       which(data.table(cvg.hold)<16)->datdat
       #store the location of contiguity in a variable called locloc
       datdat[num]->locloc
       #in case entire range is contiguous store the length of the entire region in variable
       if(length(datdat)==0){
         
       locloc<-length(cvg.hold)
       
       }
       #In case contiguity falls to close to the nearby exons ignore it
       if(locloc<100){
       locloc<-NA
       }
       #In case if it falls after that region but within 300 nts of the start then add 200 nts so that we consider the region before it in calculations
       if(!is.na(locloc) && locloc<300){
       locloc<-locloc + 200
       } 
       #Make sure it is not set to NA
	if(!is.na(locloc)){
	     #instantiate variable to store cleavage location
	     Cleavage.loc<-'None'
	     #instantiate variable to store MMD statistic
       mmd.s<-'None'
       #if in the vicinity of a junction ignore otherwise preceed
       if(!checkCigarJunc(chr.hold, locloc , GG.locs$strand[i], junc.df)){
         #Make sure the span exceeds threshold and the length of coverage exceeds 200 nts 
         if(max(vals)>12 && (length(cvg.hold)>200)){
           #ignore the first 100 and last 100 nts
           cvg.hold[100:(length(cvg.hold)-100)]->cvg.hold.me
           #consider a range of after first 100 nts and before 100 nts from the end
           coolRange<- seq(101, (length(cvg.hold.me)-100))
           #if the location of contiguous reads is less than 100 nts from the end only consider 200 nts before the contiguous reads and 800 nts after the contiguous reads
           #Case I: location of contiguous reads is less than 100 nts from the end
           #Case II: The entire range is contiguous
           #Case III: contiguous reads happen within 1000 nts from the end
           if(locloc<(length(cvg.hold.me)-1000)){
             
           coolRange<-seq(locloc-200,(locloc+800))
           
           }else  if (locloc >= (length(cvg.hold) -200)){
           coolRange<- seq(201, (length(cvg.hold.me)-200))
           
           }else{
             coolRange<-seq(locloc-200,length(cvg.hold.me))
           }
           system.time({
             kernel="splinedot"
             coolList<- lapply(coolRange, function(x) kmmd(as.matrix(cvg.hold.me[(x-200):(x)]), as.matrix(cvg.hold.me[(x):(x+200)]), kernel=kernel))
             Hy <- lapply(coolList, function(x) H0(x))
             mmd.stats <- lapply(coolList, function(x) mmdstats(x)[1])
             data.frame(as.matrix(Hy), as.matrix(mmd.stats))->coolStats
             rownames(coolStats)->coolStats$loc
             as.data.table(coolStats)->coolStats
             
             which.max((coolStats %>% filter(as.matrix.Hy.=='TRUE'))$as.matrix.mmd.stats.)->Cleavage.loc
             max(unlist((coolStats %>% filter(as.matrix.Hy.=='TRUE'))$as.matrix.mmd.stats.))->mmd.s
             remove(coolStats)
             remove(coolList)
           })
           
         }
         
       result <- append(result, paste0(chr.hold, ",", start.hold, ",", end.hold,",",GG.locs$feature_type[i],",", GG.locs$gene_name[i], ",",locloc,",", max(vals), ",", GG.locs$width[i],",", Cleavage.loc, ",", mmd.s))
       
       }
       }
     }
    
     
     data.table::fread(paste(result, collapse = "\n"))->result
     saveRDS(result, paste0(GG.locs$output[1], "/result_continu_final_",name , ".RDS"))
     return(GG.locs)
}




#' getBamCvg
#' Description: The function retrieves the coverage for every nucleotide
#' @param myBamFile path to bam file
#'
#' @return the coverage per nt
getBamCvg <- function(myBamFile){
  #get the read coverage
  BamReads.CLL11 <- readGAlignments(myBamFile)
  #get the coverage
  cvg <- coverage(BamReads.CLL11)
  return(cvg)
}

#' Title
#'
#' @param chainPath 
#'
#' @return
#' @export
#'
#' @examples
hg19Tohg38Lift <- function(chainPath){
  #liftover using rtracklayer
  path <- file.path("hg19ToHg38.over.chain")
  chain <- import.chain(path)
  #convert to granges for conversion
  locs.conv <- NULL
  #
  locs.conv$seqnames <- ipa_locs$seqnames
  locs.conv <- as.data.table(locs.conv)
  locs.conv$start <- ipa_locs$start
  locs.conv$end <- ipa_locs$end
  locs.conv$strand <- ipa_locs$strand
  locs.conv$id <- ipa_locs$id
  locs.conv <- GRanges(locs.conv)
  #liftover with rtracklayer
  locations_lifted <- liftOver(locs.conv, chain)
  locations_lifted <- as.data.table(unlist(locations_lifted))
  locations_lifted <- locations_lifted[!duplicated(locations_lifted$id),]
  data.hold <- ipa_locs[match(ipa_locs$id, locations_lifted$id),]
  data.hold$start <- locations_lifted$start
  data.hold$end <- locations_lifted$end
  return (data.hold)
}

#' Title
#'
#' @param myBamFile 
#' @param rdsSource 
#' @param sampleName 
#'
#' @return
#' @export
#'
#' @examples
ExtractIntrons<- function(myBamFile, rdsSource, sampleName){
  ref<- as.data.frame(rdsSource)
  #get the read alignemnts
  BamReads <- readGAlignments(myBamFile)
  cvg<- coverage(BamReads)
  
  saveRDS(cvg, file=paste0("coverage_data/",sampleName,"/", "cvg.Rda"))
  
  ref %>% filter(strand=="+")-> ref.df.pos
  ref %>% filter(strand=="-")-> ref.df.ne
  mutate(ref.df.pos, NextExonEnd = lead(end))-> ref.df.pos
  mutate(ref.df.ne, NextExonEnd = lag(start))-> ref.df.ne
  rbind(ref.df.ne, ref.df.pos)-> ref
  #filter for introns
  ref %>% filter(feature_type == "intron") -> ref1
  #add the intron.3utrs
  ref %>% filter(feature_type == "intron.3utr") -> ref2
  #bind the two
  rbind(ref1, ref2)-> ref
  ref %>% filter(sub(".*chr", "", ref$seqnames) %in% c(1:22, 'X', 'Y')) -> ref
  #filter for regions that don't fall in blacklists
  ref <- ref %>% filter(blacklist==0)
  #store the coverage for every range
  ref$cvg <- countOverlaps(GRanges(ref),BamReads)
  #filter for regions where coverage is non zero
  ref <- ref %>% filter(!cvg ==0)
  #filter to not go over these twice
  ref[!(ref$feature_type=="intron" && ref$before=="intron.3utr"),]-> ref
  ref <- filterStrandNoise(ref, rdsSource)
  ref <- as.data.table(ref) %>% filter(!grepl("LINC", gene_name))
  ref <- as.data.table(ref) %>% filter(!grepl("MIR", gene_name))
  ref <- as.data.table(ref) %>% filter(!grepl("SNORD", gene_name))
  
  return(ref)

}

filterStrandNoise<- function(IntronicRegions, rdsSource){
#filter the annotation source/ remove introns
as.data.table(rdsSource) %>% filter(feature_type!="intron")-> anno1
#remove intron.3utr
as.data.table(anno1) %>% filter(feature_type!="intron.3utr")-> anno2
#remove intergenic
as.data.table(anno2) %>% filter(feature_type!="intergenic")-> anno
#now get the filtered intronic regions
as.data.table(IntronicRegions)-> IntronicRegions
#subtract the strand overlapping coverage with introns of the other strand
SubtractOverlapIgnoreStrand(GRanges(IntronicRegions), GRanges(anno), "intron")-> filtered_intronicRegions

return(filtered_intronicRegions)
}

getUntemplatedAReads<- function(bamPath){
  bamFile <- BamFile(bamPath)
  aln <- scanBam(bamFile)
  aln <- aln[[1]]
  #filter for reads more than 21 nts
  aln[width(aln$seq)>21]-> aln2
  as.data.frame(aln2$seq)-> aln2
  #extract the last four nts
  aln2$x->aln2
  substr(aln2,(nchar(aln2)+1)-4,nchar(aln2)) -> lastFourNT
  
    
}
ExtractJunctions<- function(bam.reads){
  as.data.frame(junctions(bam.reads))->junc.df
  junc.df %>% distinct(seqnames, start, end, strand, .keep_all=TRUE)->junc.df
  return(junc.df)
}


findTPM<- function(anno, bam.reads){
  
  GRanges(as.data.frame(anno) %>% filter(sub(".*chr", "", as.data.frame(anno)$seqnames) %in% c(1:22, 'X', 'Y')))->anno
  countOverlaps(anno, bam.reads)->anno$cvg
  as.data.frame(anno) %>% group_by(gene_name) %>% summarize(sumCVG=sum(cvg))
  cvg.gene<- as.data.frame(anno) %>% group_by(gene_name) %>% summarize(sum.cvg=sum(cvg))
  width.gene <- as.data.frame(anno) %>% group_by(gene_name) %>% summarize(sum.width=sum(width))
  data.frame(cvg.gene, width.gene)->gene.associated;gene.associated$gene_name.1<-NULL
  gene.associated$RPK<- 1000*gene.associated$sum.cvg/(gene.associated$sum.width)
  sum(gene.associated$RPK)/1000000->perMillionSF
  gene.associated$RPK/perMillionSF->gene.associated$TPM
  gene.associated %>% filter(TPM>1)->res
  return(res)
}



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
##Taken from the previous script for annotation
#####goal is to overlap introns stored in an object with the rest of the ranges excluding introns regardless of their strands

