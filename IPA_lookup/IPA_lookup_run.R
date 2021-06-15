#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Texas A&M University, Department of Biomedical Engineering
#' Spring 2021

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

#set working directory
wd <- "/scratch/user/annicenajafi/042821"; user<- "annicenajafi"; 
setwd(wd)

source("scripts/IPA_lookup_functions.R")



#color scheme: aquamarine, cadet blue, eggplant, Persian orange, Arylide yellow, yellow-green crayola 
color.scheme <- c('#48D1CC', '#51A3A3', '#75485E', '#CB904D', '#37D5A3', "#DFCC74", "#C3E991")
#read the bam file
myBamFile <- file.path("/scratch/group/isinghlab/projects/Hem_IPA/bams/rnaseq/CLL/RNA_CLL11.bam"); 
#read the rds file
runChangePoint("data_tables/sample_input.txt", wd)
print("job submitted")


