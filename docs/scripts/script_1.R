#LIBRARIES
#------------------------------libraries----------------------------------------
if (!require("pacman")) install.packages("pacman") 
#install required libraries (if needed) and load them
pacman::p_load(dplyr,gread,rtracklayer,GenomicRanges,HelloRanges,seqinr,stringr,GenomicAlignments,tximport,tximportData,edgeR,DESeq2,readr,
               pheatmap,DEGreport,vsn,data.table,reshape,gplots,VennDiagram,limma,WGCNA,sm,network,tidyr,
               viridis,hrbrthemes,spgs,ggplot2,tidyverse,strex,ncRNAtools,tRNA,ggrepel,ggh4x,nptest,ggpubr,rstatix,datarium,fdrtool,RColorBrewer)
library(ggplot2)
library(RRNA)
library(dplyr)
library(stringr)
library("xlsx")
library(rvest)
library(plotrix)
library(ggh4x)
#tRNA_fasta=read.fasta("F:/PARN_ELAC_silencing/smallRNA/tRNA/pretRNA_from_chr_and_mtDNA_multi_dedup.fasta")
# tRNA_fasta=read.fasta("F:/PARN_ELAC_silencing/smallRNA/tRNA/pretRNA_from_chr_and_mtDNA_multi_dedup.fasta")
# getName(tRNA_fasta)
# getLength(tRNA_fasta)
# tRNA_trailer_pos=data.frame(ref=getName(tRNA_fasta),
#                             ref_tRNA_len=getLength(tRNA_fasta))
# #tRNA_trailer_pos$strand=substr(tRNA_trailer_pos$ref,nchar(tRNA_trailer_pos$ref)-1,nchar(tRNA_trailer_pos$ref)-1)
# tRNA_trailer_pos$three_prime_end_start=tRNA_trailer_pos$ref_tRNA_len-19
# 
# duplicates=read.table("F:/PARN_ELAC_silencing/smallRNA/smallRNAwithAdapters/miRNA/corrected_duplicate_clusters.tsv",header = TRUE)

#loading all genome and transcriptome
(trans_list=list.files("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_transcriptome/mapped_seq_with_strand/", pattern=NULL, all.files=FALSE,
                       full.names=FALSE))
(genome_list=list.files("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_genome/mapped_seq_with_strand/", pattern=NULL, all.files=FALSE,
                        full.names=FALSE))
(sample_list=list.files("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_RNAcentral/mapped_seq_with_strand_new/", pattern=NULL, all.files=FALSE,
                        full.names=FALSE))

#duplicates=read.table("F:/PARN_ELAC_silencing/smallRNA/smallRNAwithAdapters/miRNA/corrected_duplicate_clusters.tsv",header = TRUE)
# set1=read.table(paste("F:/smallRNAwithAdapters/miRNA/mapped_corrected_seq/",sample_list[1],sep=""))
# set1=read.table("F:/smallRNAwithAdapters/miRNA/mapped_seq/test.txt")

# colnames(duplicates) # "RetainedRef"  "DuplicateRef"
# for (j in 1:length(trans_list)) {
#   print(paste("Proceeding file",j))
#   trans_seq=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_transcriptome/mapped_seq_with_strand/",trans_list[j],
#                              sep=""),col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
#   genome_seq=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_genome/mapped_seq_with_strand/",genome_list[j],
#                               sep=""),col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
#   new_genome_seq=subset(genome_seq,!(genome_seq$read %in% trans_seq$read))
#   mapped_seq=rbind(trans_seq,new_genome_seq)
#   mapped_seq$RNA_type=ifelse(startsWith(mapped_seq$ref,"dd_Smed"),"mRNA fragments","genome")
#   mapped_seq$cigar_interpr=cigarOpTable(mapped_seq$cigar)
#   mapped_seq$cigar_del=mapped_seq$cigar_interpr[,3]
#   mapped_seq$cigar_ins=mapped_seq$cigar_interpr[,2]
#   mapped_seq$cigar_match=mapped_seq$cigar_interpr[,1]
#   mapped_seq$cigar_mismatch=(mapped_seq$cigar_interpr[,2]+mapped_seq$cigar_interpr[,3]+mapped_seq$cigar_interpr[,4]+mapped_seq$cigar_interpr[,5]+mapped_seq$cigar_interpr[,6])
#   mapped_seq$orig_seq=ifelse(mapped_seq$strand==0,mapped_seq$seq,reverseComplement(mapped_seq$seq,case="upper"))
#   #mapped_seq$polyA_start=ifelse(startsWith(mapped_seq$orig_seq,"AAAA"),"polyA_start","no_polyA_start")
#   #mapped_seq$polyT_start=ifelse(startsWith(mapped_seq$orig_seq,"TTTT"),"polyT_start","no_polyT_start")
#   #mapped_seq$polyA_end=ifelse(endsWith(mapped_seq$orig_seq,"AAAA"),"polyA_end","no_polyA_end")
#   #mapped_seq$polyT_end=ifelse(endsWith(mapped_seq$orig_seq,"TTTT"),"polyT_end","no_polyT_end")
#   assign(paste("trans_genome_RNA",j,sep=""),mapped_seq)
# }



# save(trans_genome_RNA1,trans_genome_RNA2,trans_genome_RNA3,trans_genome_RNA4,trans_genome_RNA5,trans_genome_RNA6,trans_genome_RNA7,
#      trans_genome_RNA8,trans_genome_RNA9,trans_genome_RNA10,trans_genome_RNA11,trans_genome_RNA12, trans_genome_RNA13,
#      trans_genome_RNA14,trans_genome_RNA15,trans_genome_RNA16,trans_genome_RNA17,trans_genome_RNA18,trans_genome_RNA19,
#      trans_genome_RNA20,trans_genome_RNA21,trans_genome_RNA22,trans_genome_RNA23,trans_genome_RNA24,file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_trans_genome.RData') #transcriptome, genome data
# 


#load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_trans_genome.RData')
#colnames(trans_genome_RNA1)
#head(trans_genome_RNA1)
# [1] "read"           "strand"         "ref"            "position"       "qual"           "cigar"          "seq"            "XM"             "MD"            
#"RNA_type"       "cigar_interpr"  "cigar_del"      "cigar_ins"      "cigar_match"    "cigar_mismatch"
# [16] "orig_seq"
#loading all scnRNA
(sample_list=list.files("E:/Illumina/PARN_ELAC2_silencing/smallRNA/smallRNAwithAdapters/miRNA/bowtie2_mapped_SM_tRNA_seq_with_strand/", pattern=NULL, all.files=FALSE,
                        full.names=FALSE))
# (sample_list2=list.files("F:/PARN_ELAC_silencing/smallRNA/tRNA/mapped_to_tRNA/mapped_seq/", pattern=NULL, all.files=FALSE,
#                         full.names=FALSE))
#/smallRNAwithAdapters/miRNA/bowtie2_mapped_SM_tRNA_seq_with_strand
# Define a function to extract substring before 6th "_"
# extract_substring <- function(x) {
#   output <- gsub("^(([^_]*_){5}[^_]*)_.*", "\\1", x)
#   return(output)
# }
original_set=c("PARN13S","ELAC23S","GFP33S","WT15S","GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S",
               "ELAC15S","PARN15S","ELAC25S","GFP35S","WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")

for (j in 1:length(sample_list)) {
  mapped_seq=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_RNAcentral/mapped_seq_with_strand_new/",sample_list[j],sep=""),
                        col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  mapped_seq_tRNA=read.table(paste("E:/Illumina/PARN_ELAC2_silencing/smallRNA/smallRNAwithAdapters/miRNA/bowtie2_mapped_SM_tRNA_seq_with_strand/",sample_list[j],sep=""),
                             col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  mapped_seq_rRNA=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_rRNA_and_genome/rRNA_seq_with_strand/",sample_list[j],sep=""),
                        col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  #mapped_seq=subset(mapped_seq,!(mapped_seq$read %in% mapped_seq_tRNA$read))
  print(paste("Proceeding file",j))
  mapped_seq_tRNA$RNA_type="tRNA fragments"
  mapped_seq_rRNA$RNA_type="rRNA fragments"
  #change annotation, remove "_smed_chr2_pos_10496669-10496742(-)" 
  #mapped_seq_tRNA$ref=ifelse(startsWith(mapped_seq_tRNA$ref,"dd_Smed_g"),lapply(mapped_seq_tRNA$ref, extract_substring),mapped_seq_tRNA$ref)
  
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("ribosomal", ignore_case = TRUE)),"rRNA fragments","other")
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("rRNA", ignore_case = TRUE)),"rRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("tRNA", ignore_case = TRUE)),"tRNA fragments",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("lncRNA", ignore_case = TRUE)),"lncRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("piRNA", ignore_case = TRUE)),"piRNA",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("RNase_P_RNA", ignore_case = TRUE)),"RNase_P_RNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("snoRNA", ignore_case = TRUE)),"snoRNA fragments",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("siRNA", ignore_case = TRUE)),"siRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("miRNA", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("pre_miRNA", ignore_case = TRUE)),"pre_miRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("snRNA", ignore_case = TRUE)),"snRNA fragments",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("Sme-", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)  
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("spliceosomal-", ignore_case = TRUE)),"snRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("spliceosomal", ignore_case = TRUE)),"snRNA fragments",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("lnc", ignore_case = TRUE)),"lncRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("microRNA", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  #Mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ID, regex("long_non-coding", ignore_case = TRUE)),"lncRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("7SK", ignore_case = TRUE)),"snRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("nucleolar", ignore_case = TRUE)),"snoRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("nuclear", ignore_case = TRUE)),"snRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("miR", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("transfer", ignore_case = TRUE)),"tRNA fragments",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("sme-lin", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("sme-let", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  mapped_seq$RNA_type=ifelse(str_detect(mapped_seq$ref, regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",mapped_seq$RNA_type)
  
  new_mapped_seq=subset(mapped_seq,!(mapped_seq$read %in% mapped_seq_tRNA$read))
  mapped_seq_all=rbind(new_mapped_seq,mapped_seq_tRNA)
  mapped_seq_rRNA=subset(mapped_seq_rRNA,!(mapped_seq_rRNA$read %in% mapped_seq_all$read))
  mapped_seq_all=rbind(mapped_seq_all,mapped_seq_rRNA)
  trans_seq=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_transcriptome/mapped_seq_with_strand/",trans_list[j],
                             sep=""),col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  # genome_seq=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_genome/mapped_seq_with_strand/",genome_list[j],
  #                             sep=""),col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  new_trans_seq=subset(trans_seq,!(trans_seq$read %in% mapped_seq_all$read))
  new_trans_seq$RNA_type="mRNA fragments"
  mapped_seq_all=rbind(mapped_seq_all,new_trans_seq)
  #mapped_seq_all=rbind(mapped_seq_all,mapped_seq_rRNA)
  mapped_seq_all$cigar_interpr=cigarOpTable(mapped_seq_all$cigar)
  
  mapped_seq_all$cigar_del=mapped_seq_all$cigar_interpr[,3]
  mapped_seq_all$cigar_ins=mapped_seq_all$cigar_interpr[,2]
  mapped_seq_all$cigar_match=mapped_seq_all$cigar_interpr[,1]
  mapped_seq_all$cigar_mismatch=(mapped_seq_all$cigar_interpr[,2]+mapped_seq_all$cigar_interpr[,3]+mapped_seq_all$cigar_interpr[,4]+mapped_seq_all$cigar_interpr[,5]+mapped_seq_all$cigar_interpr[,6])
  mapped_seq_all$orig_seq=ifelse(mapped_seq_all$strand==0,mapped_seq_all$seq,reverseComplement(mapped_seq_all$seq,case="upper"))
 
  mapped_seq_all$new_XM=as.numeric(substr(mapped_seq_all$XM,6,nchar(mapped_seq_all$XM)))
  mapped_seq_all$all_mm=as.numeric(mapped_seq_all$new_XM+mapped_seq_all$cigar_mismatch)
  #new_set=subset(tableset,tableset$XM=="XM:i:0"|tableset$XM=="XM:i:1"|tableset$XM=="XM:i:2"|tableset$XM=="XM:i:3")
  #new_set$filter=ifelse(new_set$all_mm>3, "low_quality_no_tail", "good_quality_or_has_tail")
  mapped_seq_all=subset(mapped_seq_all,mapped_seq_all$all_mm<4)
  mapped_seq_all=subset(mapped_seq_all,mapped_seq_all$strand==0)
  #new_set$DuplicateRef=new_set$ref
  #new_set_merged=merge(new_set,duplicates,by="DuplicateRef",all.x = TRUE)
  #new_set_merged$RetainedRef=ifelse(is.na(new_set_merged$RetainedRef),new_set_merged$DuplicateRef,new_set_merged$RetainedRef)
  assign(paste(original_set[j]),mapped_seq_all)
  #assign(paste("smallRNA",j,sep=""),mapped_seq_all)
}

#check if there are any reads dublicates
nrow(PARN13S)
length(unique(PARN13S$read))

nrow(ELAC23S)
length(unique(ELAC23S$read))
#save all sets
save(PARN13S,ELAC23S,GFP33S,WT15S,GFP25S,WT35S,PARN23S,ELAC33S,GFP15S,WT25S,PARN33S,ELAC15S,PARN15S,ELAC25S,GFP35S,
     WT13S,PARN25S,GFP23S,PARN35S,ELAC13S,ELAC35S,GFP13S,WT23S,WT33S, file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_no_genome_bad_anno.RData') #transcriptome and scnRNA (both nuclear and mtDNA tRNA) data

load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_no_genome_bad_anno.RData')
ELAC13S$set="ELAC3S"
ELAC23S$set="ELAC3S"
ELAC33S$set="ELAC3S"

ELAC15S$set="ELAC5S"
ELAC25S$set="ELAC5S"
ELAC35S$set="ELAC5S"

GFP13S$set="GFP3S"
GFP23S$set="GFP3S"
GFP33S$set="GFP3S"

GFP15S$set="GFP5S"
GFP25S$set="GFP5S"
GFP35S$set="GFP5S"

WT13S$set="WT3S"
WT23S$set="WT3S"
WT33S$set="WT3S"

WT15S$set="WT5S"
WT25S$set="WT5S"
WT35S$set="WT5S"

all_sets_table=rbind(ELAC23S,GFP33S,WT15S,GFP25S,WT35S,
                     ELAC33S,GFP15S,WT25S,
                     ELAC15S,ELAC25S,GFP35S, WT13S,
                     GFP23S,ELAC13S,ELAC35S,
                     GFP13S,WT23S,WT33S)
#check how many are reverse strand
#table(all_sets_table$strand)
#890809/(890809+20931385) 0.04% are on reverse strand
#all_sets_table=subset(all_sets_table,all_sets_table$strand==0)
nrow(all_sets_table) #27935936
length(unique(all_sets_table$orig_seq)) #3101937
#-------------------------EXTRACT_UNIQUE_SEQUNCES_WITH_ANNOTATION---------------
library(dplyr)
library(stringr)

# Create a new column by concatenating orig_seq and annotation columns
all_sets_table$combined <- paste(all_sets_table$position,all_sets_table$orig_seq, all_sets_table$ref, sep = "_")
#add position to annotation
all_sets_table$pos_ref=paste(all_sets_table$position, all_sets_table$ref, sep = "_")
#unique_combinations <- unique(all_sets_table[!duplicates, c('orig_seq', 'pos_ref')])
unique_combinations <- unique(all_sets_table[, c('orig_seq', 'pos_ref')])
unique_combinations$pos_ref=as.character(unique_combinations$pos_ref)

#check if there are sequences both annotated as rRNA and tRNA
# duplicates_in_unique_combinations <- duplicated(unique_combinations$orig_seq)
# diff_RNA_type=unique_combinations[duplicates_in_unique_combinations,]
#diff_RNA_type=diff_RNA_type[sort(diff_RNA_type$orig_seq),]
sorted_diff_RNA_type <- unique_combinations[order(unique_combinations$orig_seq),]
colnames(sorted_diff_RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("ribosomal", ignore_case = TRUE)),"rRNA fragments","other")
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("rRNA", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("tRNA", ignore_case = TRUE)),"tRNA fragments",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("lncRNA", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("piRNA", ignore_case = TRUE)),"piRNA",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("RNase_P_RNA", ignore_case = TRUE)),"RNase_P_RNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("snoRNA", ignore_case = TRUE)),"snoRNA fragments",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("siRNA", ignore_case = TRUE)),"siRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("miRNA", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("pre_miRNA", ignore_case = TRUE)),"pre_miRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("snRNA", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("Sme-", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)  
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("spliceosomal-", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("spliceosomal", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("lnc", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("microRNA", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
#Mapped_seq$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$ID, regex("long_non-coding", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("7SK", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("nucleolar", ignore_case = TRUE)),"snoRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("nuclear", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("miR", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("transfer", ignore_case = TRUE)),"tRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("sme-lin", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("sme-let", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("dd_Smed_v6", ignore_case = TRUE)),"mRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("ITS1", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("ITS2", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("SpacerA", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("28S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("Schmed_cloneH735c", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("12S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("16S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("5.8S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("long_non-coding", ignore_case = TRUE)),"lncRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type$RNA_type=ifelse(str_detect(sorted_diff_RNA_type$pos_ref, regex("lnc", ignore_case = TRUE)),"lncRNA fragments",sorted_diff_RNA_type$RNA_type)
sorted_diff_RNA_type_df=as.data.frame(table(sorted_diff_RNA_type$orig_seq,sorted_diff_RNA_type$RNA_type))
colnames(sorted_diff_RNA_type_df)=c("seq","RNA","Freq")
sorted_diff_RNA_type_wide_df=reshape(sorted_diff_RNA_type_df, idvar = "seq", timevar = "RNA", direction = "wide")
nrow(sorted_diff_RNA_type_wide_df) #3101937
# Check for rows with more than 1 non-zero column
rows_with_multiple_nonzero <- apply(sorted_diff_RNA_type_wide_df[,-1], 1, function(x) sum(x != 0) > 1)
# Print rows with more than 1 non-zero column
seq_dif_type=sorted_diff_RNA_type_wide_df[rows_with_multiple_nonzero, ]
nrow(seq_dif_type) #18918 sequenceds with not clear annotation --> 0.00609877% of all sequences
18918/3101937
nrow(sorted_diff_RNA_type)
# not_clear_RNA_type=unique_combinations[unique_combinations$orig_seq %in% seq_dif_type$seq,]
# colnames(not_clear_RNA_type)
# not_clear_RNA_type=not_clear_RNA_type[sort(not_clear_RNA_type$orig_seq),]
# seq_dif_tRNA_type=seq_dif_type[seq_dif_type$`Freq.tRNA fragments`>0,]
# 
# head(sorted_diff_RNA_type)
# 
# ?head
# seq_to_test=seq_dif_tRNA_type$seq[1:10]
# test_df=sorted_diff_RNA_type[sorted_diff_RNA_type$orig_seq %in% seq_to_test,]



# sorted_diff_RNA_type
# 
# 
# test_df_corrected <- test_df %>%
#   group_by(orig_seq) %>%
#   mutate(num_rna_types = n_distinct(RNA_type)) %>%
#   mutate(pos_ref = case_when(
#     num_rna_types == 1 ~ sample(pos_ref, 1),
#     num_rna_types > 1 & any(RNA_type == "tRNA fragments")  ~ 
#       pos_ref[which(RNA_type == "tRNA fragments")][1],
#     num_rna_types > 1 & any(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE))) ~ 
#       pos_ref[which(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE)))][1],
#     num_rna_types > 1 & any(RNA_type == "other") & !any(RNA_type == "tRNA fragments") ~ 
#       pos_ref[which(RNA_type == "other")][1],
#     num_rna_types > 1 & all(RNA_type != "other") ~ "multiple_hits"
#   )) %>%
#   ungroup() %>%
#   select(-num_rna_types)

sorted_diff_RNA_type_corrected <- sorted_diff_RNA_type %>%
  group_by(orig_seq) %>%
  mutate(num_rna_types = n_distinct(RNA_type)) %>%
  mutate(pos_ref = case_when(
    num_rna_types == 1 ~ sample(pos_ref, 1),
    num_rna_types > 1 & any(RNA_type == "tRNA fragments")  ~ 
      pos_ref[which(RNA_type == "tRNA fragments")][1],
    num_rna_types > 1 & any(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE))) ~ 
      pos_ref[which(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE)))][1],
    num_rna_types > 1 & any(RNA_type == "other") & !any(RNA_type == "tRNA fragments") ~ 
      pos_ref[which(RNA_type == "other")][1],
    num_rna_types > 1 & all(RNA_type != "other") ~ "multiple_hits"
  )) %>%
  ungroup() %>%
  select(-num_rna_types)

# sorted_diff_RNA_type_corrected <- sorted_diff_RNA_type %>%
#   group_by(orig_seq) %>%
#   mutate(num_rna_types = n_distinct(RNA_type)) %>%
#   mutate(pos_ref = case_when(
#     num_rna_types == 1 ~ sample(pos_ref, 1),
#     num_rna_types > 1 & any(RNA_type == "tRNA fragments")  ~ 
#       pos_ref[which(RNA_type == "tRNA fragments")][1],
#     num_rna_types > 1 & any(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE))) ~ 
#       pos_ref[which(str_detect(pos_ref, regex("Schmidtea", ignore_case = TRUE)))][1],
#     num_rna_types > 1 & any(RNA_type == "other") & !any(RNA_type == "tRNA fragments") ~ 
#       pos_ref[which(RNA_type == "other")][1],
#     num_rna_types > 1 & all(RNA_type != "other") ~ "multiple_hits"
#   )) %>%
#   ungroup() %>%
#   select(-num_rna_types)


sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("ribosomal", ignore_case = TRUE)),"rRNA fragments","other")
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("rRNA", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("tRNA", ignore_case = TRUE)),"tRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("lncRNA", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("piRNA", ignore_case = TRUE)),"piRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("RNase_P_RNA", ignore_case = TRUE)),"RNase_P_RNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("snoRNA", ignore_case = TRUE)),"snoRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("siRNA", ignore_case = TRUE)),"siRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("miRNA", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("pre_miRNA", ignore_case = TRUE)),"pre_miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("snRNA", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("Sme-", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)  
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("spliceosomal-", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("spliceosomal", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("lnc", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("microRNA", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
#Mapped_seq$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$ID, regex("long_non-coding", ignore_case = TRUE)),"lncRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("7SK", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("nucleolar", ignore_case = TRUE)),"snoRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("nuclear", ignore_case = TRUE)),"snRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("miR", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("transfer", ignore_case = TRUE)),"tRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("sme-lin", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("sme-let", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("dd_Smed_v6", ignore_case = TRUE)),"mRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("ITS1", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("ITS2", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("SpacerA", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("28S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("Schmed_cloneH735c", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("12S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("16S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("5.8S", ignore_case = TRUE)),"rRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("multiple_hits", ignore_case = TRUE)),"multiple_hits",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("lnc", ignore_case = TRUE)),"lncRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)
sorted_diff_RNA_type_corrected$RNA_type_corrected=ifelse(str_detect(sorted_diff_RNA_type_corrected$pos_ref, regex("long_non-coding", ignore_case = TRUE)),"lncRNA fragments",sorted_diff_RNA_type_corrected$RNA_type_corrected)


#table(scnRNA_anno_uniq_corrected$RNA_type_corrected)

#print(new_test_df,n=100)
unique(sorted_diff_RNA_type_corrected[sorted_diff_RNA_type_corrected$RNA_type_corrected=="multiple_hits",])
sorted_diff_RNA_type_corrected=sorted_diff_RNA_type_corrected[,-3]
unique_sorted_diff_RNA_type_corrected=unique(sorted_diff_RNA_type_corrected)
nrow(unique_sorted_diff_RNA_type_corrected)
colnames(unique_sorted_diff_RNA_type_corrected)=c("orig_seq","correct_annotation","correct_RNA_type")
table(unique_sorted_diff_RNA_type_corrected$correct_RNA_type)
# 
# test_df=head(sorted_diff_RNA_type_wide_df)
# rownames(test_df)=test_df$seq
# test_df=test_df[,-1]
# 
# # Create example matrix
# m <- matrix(c(1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1), ncol = 3, byrow = TRUE)
# 
# # Check for rows with more than 1 non-zero column
# ?apply
# rows_with_multiple_nonzero <- apply(test_df, 1, function(x) sum(x != 0) > 1)
# 
# # Print rows with more than 1 non-zero column
# m[rows_with_multiple_nonzero, ]
# 

# 
# # Group by sncRNA and select the most frequent annotation
# scnRNA_anno_uniq_corrected <- unique_combinations %>%
#   group_by(orig_seq) %>%
#   summarise(pos_ref = names(which.max(table(pos_ref))))

# save(scnRNA_anno_uniq_corrected,file="F:/PARN_ELAC_silencing/smallRNA/scnRNA_anno_uniq_corrected_no_genome.RData")
# sorted_diff_RNA_type_corrected
save(unique_sorted_diff_RNA_type_corrected,file="F:/PARN_ELAC_silencing/smallRNA/sorted_diff_RNA_type_corrected_new.RData")
load("F:/PARN_ELAC_silencing/smallRNA/sorted_diff_RNA_type_corrected_new.RData")
unique_sorted_diff_RNA_type_corrected[unique_sorted_diff_RNA_type_corrected$correct_RNA_type=="other",]
library(stringr)
str_view(unique_sorted_diff_RNA_type_corrected$correct_annotation, "79327")
#load("F:/PARN_ELAC_silencing/smallRNA/sorted_diff_RNA_type_corrected.RData")
#colnames(all_sets_table)
#colnames(sorted_diff_RNA_type_corrected)=c("orig_seq","correct_annotation","correct_RNA_type")

#----------------------CORRECT_ANNOTATION----------------------------------------
#load("F:/PARN_ELAC_silencing/smallRNA/scnRNA_anno_uniq_corrected.RData")
#scnRNA_all_sets_corrected=merge(all_sets_table,sorted_diff_RNA_type_corrected,by="orig_seq",all.x = TRUE)

#scnRNA_all_sets_corrected=subset(scnRNA_all_sets_corrected,scnRNA_all_sets_corrected$strand==0)
#scnRNA_all_sets_corrected=subset(scnRNA_all_sets_corrected,scnRNA_all_sets_corrected$all_mm<4)
#save(scnRNA_all_sets_corrected,file="F:/PARN_ELAC_silencing/smallRNA/scnRNA_all_sets_corrected_no_genome.RData")
#load("F:/PARN_ELAC_silencing/smallRNA/scnRNA_all_sets_corrected_no_genome.RData")

#subset(scnRNA_all_sets_corrected,scnRNA_all_sets_corrected$orig_seq=="AAAACCCTTAGTTGAGC")
################################################################################
#plot pie chart for RNA types
# load("F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_analysis_filtered.RData")
# colnames(PARN13S_filtered)
#PARN13S_filtered$RNA_pos
# PARN13S$set="PARN3S"
# PARN23S$set="PARN3S"
# PARN33S$set="PARN3S"
# 
# PARN15S$set="PARN5S"
# PARN25S$set="PARN5S"
# PARN35S$set="PARN5S"

#unique(WT35S$RNA_type)


#for PARN
# PARN_sets_table=rbind(PARN13S,PARN23S,PARN33S,PARN13S,PARN25S,PARN35S,
#                      GFP13S,GFP23S,GFP33S,GFP13S,GFP25S,GFP35S)
# save(PARN_sets_table, file = "E:/project/smallRNA/calculated_data_bowtie2_end_to_end/PARN_sets_table.RData")
# #PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("spliceosomal", ignore_case = TRUE)),"snRNA",PARN_sets_table$RNA_type)
# #PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("snRNA", ignore_case = TRUE)),"snRNA fragments",PARN_sets_table$RNA_type)
# #PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("snoRNA", ignore_case = TRUE)),"snoRNA fragments",PARN_sets_table$RNA_type)
# #PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("tRNA", ignore_case = TRUE)),"tRNA fragments",PARN_sets_table$RNA_type)
# #PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("rRNA", ignore_case = TRUE)),"rRNA fragments",PARN_sets_table$RNA_type)
# PARN_sets_table=subset(PARN_sets_table,PARN_sets_table$RNA_type!="genome")
# PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$ref, regex("lnc", ignore_case = TRUE)),"lncRNA fragments",PARN_sets_table$RNA_type)
# PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$ref, regex("long", ignore_case = TRUE)),"lncRNA fragments",PARN_sets_table$RNA_type)
# PARN_sets_table$RNA_type=ifelse(str_detect(PARN_sets_table$RNA_type, regex("other", ignore_case = TRUE)),"other fragments",PARN_sets_table$RNA_type)
# table(PARN_sets_table$RNA_type)
# 
# PARN_sets_table$length=nchar(PARN_sets_table$seq)
# #nchar(PARN_sets_table$seq[1])
# 
# #length distribution
# ?geom_histogram
# 
# PARN_sets_table=subset(PARN_sets_table,PARN_sets_table$length<85)
# (p <- PARN_sets_table %>%
#     mutate(RNA_type = fct_reorder(RNA_type, length)) %>%
#     ggplot( aes(x=length, color=RNA_type, fill=RNA_type)) +
#     geom_histogram(alpha=0.6, binwidth = 1) +
#     scale_fill_viridis(discrete=TRUE) +
#     scale_color_viridis(discrete=TRUE) +
#     theme_ipsum() +
#     xlab("length") +
#     ylab("number of RNA species ") +
#     facet_wrap(~RNA_type, scales = "free")+
#     theme(
#       legend.position="none",
#       panel.spacing = unit(0.1, "lines"),
#       strip.text.x = element_text(size = 20,hjust=0.5),
#       axis.title.x = element_text(size = 15,hjust=0.5,face = "plain"),
#       axis.title.y = element_text(size = 15,hjust=0.5,face = "plain"))+
#     scale_x_continuous(breaks = seq(15, 115, by = 10)))
# ?element_text
# 
# 
# PARN_sets_table_df=as.data.frame(table(PARN_sets_table$set,PARN_sets_table$RNA_type))
# colnames(PARN_sets_table_df)=c("Comparison","Expression","Frequency")
# PARN_sets_table_sum=aggregate(Frequency~Comparison, PARN_sets_table_df,sum)
# colnames(PARN_sets_table_sum)=c("Comparison","sum")
# PARN_sets_table_df=merge(all_sets_table_df,PARN_sets_table_sum,by="Comparison",all.x = TRUE)
# all_sets_table_df$perc=round((all_sets_table_df$Frequency/all_sets_table_df$sum)*100, digits = 2)
# colnames(all_sets_table_df)=c("Set","RNA_type","Frequency","sum","perc")
# #write.csv(all_sets_table_df,"E:/project/smallRNA/RNA_type_perc.csv")
# stest$RNA_type <- factor(stest$RNA_type, levels=c("genome","piRNA","miRNA","mRNA","other",
#                                                   "rRNA","snoRNA","tRNA","snRNA"))
# ggplot(all_sets_table_df, aes(x = 3, y = perc, fill = RNA_type)) +
#   geom_col(color = "white") +
#   geom_text_repel(aes(label = label),
#                   position = position_stack(vjust = 0.5)) +
#   coord_polar(theta = "y") +
#   scale_fill_brewer(palette = "Set3") +
#   xlim(c(0.2, 3.5)) +
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank())+
#   facet_wrap(.~ Set,nrow=2) +theme_void()+theme(
#     plot.title = element_text(size=20, face="bold",hjust = 0.5),
#     text = element_text(size=25))


#####################################################################################
#remove duplicates once more
#####################################################################################
# all_sets_table
# # filtered_all_sets_table=rbind(PARN13S,ELAC23S,GFP33S,WT15S,GFP25S,WT35S,PARN23S,ELAC33S,GFP15S,WT25S,PARN33S,ELAC15S,PARN15S,ELAC25S,GFP35S,
# #                               WT13S,PARN25S,GFP23S,PARN35S,ELAC13S,ELAC35S,GFP13S,WT23S,WT33S)
# filtered_all_sets_table=subset(all_sets_table,all_sets_table$strand==0)
# colnames(filtered_all_sets_table)
# #unique(filtered_all_sets_table$strand)
# # [1] "DuplicateRef"   "strand"         "ref"            "position"       "qual"           "cigar"          "seq"            "XM"             "MD"             "RNA_type"      
# # [11] "cigar_interpr"  "cigar_match"    "cigar_mismatch" "orig_seq"       "polyA_start"    "polyT_start"    "polyA_end"      "polyT_end"      "new_XM"         "all_mm"        
# # [21] "filter"         "RetainedRef"  
# uniq_filtered_all_sets_table=unique(filtered_all_sets_table[c("orig_seq","RetainedRef","position","XM","MD","cigar")])
# nrow(uniq_filtered_all_sets_table) #10211610
# length(unique(uniq_filtered_all_sets_table$orig_seq)) #7207365
# seq_freq=as.data.frame(table(uniq_filtered_all_sets_table$orig_seq))
# nrow(seq_freq) #7207365
# dupl_seq=subset(seq_freq,seq_freq$Freq>1)
# nrow(dupl_seq) #301478
# dupl_seq_set=uniq_filtered_all_sets_table[uniq_filtered_all_sets_table$orig_seq %in% dupl_seq$Var1,]
# nrow(dupl_seq_set) #2957004
# 
# #write.table(dupl_seq_set,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/new_dupl_seq_end_to_end.txt",row.names=FALSE,quote=FALSE)
# #dupl_seq_set_data=read.table("F:/smallRNAwithAdapters/miRNA/new_dupl_seq.txt")
# #dupl_seq_set[sort(dupl_seq_set$orig_seq),]
# sorted_dupl_seq_set_data=dupl_seq_set[order(dupl_seq_set$orig_seq),]
# length(sorted_dupl_seq_set_data$orig_seq) #2957004
# length(unique(sorted_dupl_seq_set_data$orig_seq) ) #301478
# #colnames(sorted_dupl_seq_set_data)=c("orig_seq","RNA","pos")
# sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl=paste(sorted_dupl_seq_set_data$RetainedRef,sorted_dupl_seq_set_data$position,
#                                                              sorted_dupl_seq_set_data$XM,sorted_dupl_seq_set_data$MD,sorted_dupl_seq_set_data$cigar, sep=" " )
# new_sorted_dupl_seq_set_data=aggregate(data=sorted_dupl_seq_set_data,RNA_pos_uniq_XM_MD_cigar_dupl~orig_seq,FUN=paste)
# test_sorted=new_sorted_dupl_seq_set_data
# for (i in 1:nrow(test_sorted)) {
#   new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl[[i]]=new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl[[i]][1]
# }
# 
# 
# 
# dupl_uniq_seq=new_sorted_dupl_seq_set_data
# colnames(dupl_uniq_seq)
# colnames(test_sorted)
# new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl=unlist(new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl) 
# #write.table(new_sorted_dupl_seq_set_data,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/new_dupl_uniq_seq_end_to_end.txt",row.names=FALSE,quote=FALSE)
# colnames(new_sorted_dupl_seq_set_data) #orig_seq RNA_pos_uniq_XM_MD_cigar_dupl
# #new_sorted_dupl_seq_set_data=read.table("F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/new_dupl_uniq_seq_end_to_end.txt")
# # new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl=paste(new_sorted_dupl_seq_set_data$V2,new_sorted_dupl_seq_set_data$V3,
# #                                                                  new_sorted_dupl_seq_set_data$V4,new_sorted_dupl_seq_set_data$V5,new_sorted_dupl_seq_set_data$V6, sep=" ")
# # new_sorted_dupl_seq_set_data=new_sorted_dupl_seq_set_data[c(1,7)]
# colnames(new_sorted_dupl_seq_set_data)=c("orig_seq", "RNA_pos_uniq_XM_MD_cigar_dupl")
# nrow(new_sorted_dupl_seq_set_data)
# length(new_sorted_dupl_seq_set_data$orig_seq)
# length(new_sorted_dupl_seq_set_data$RNA_pos_uniq_XM_MD_cigar_dupl)
# head(new_sorted_dupl_seq_set_data)
# save(new_sorted_dupl_seq_set_data, file = 'F:/PARN_ELAC_silencing/smallRNA/new_sorted_dupl_seq_set_data.RData') #transcriptome, genome and scnRNA (both nuclear and mtDNA tRNA) data
# save(dupl_uniq_seq, file = 'F:/PARN_ELAC_silencing/smallRNA/dupl_uniq_seq.RData')
#################################################################################### 
# load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_new_new_all.RData')
filtered_sets=list(ELAC23S,GFP33S,WT15S,GFP25S,WT35S,ELAC33S,GFP15S,WT25S,ELAC15S,ELAC25S,GFP35S,
                   WT13S,GFP23S,ELAC13S,ELAC35S,GFP13S,WT23S,WT33S)
original_set=c("ELAC23S","GFP33S","WT15S","GFP25S","WT35S","ELAC33S","GFP15S","WT25S","ELAC15S","ELAC25S","GFP35S",
               "WT13S","GFP23S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")
# 
# 
# original_set=c("PARN13S","ELAC23S","GFP33S","WT15S","GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S",
#                "ELAC15S","PARN15S","ELAC25S","GFP35S","WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")
# 
# filtered_sets=list(PARN13S,ELAC23S,GFP33S,WT15S,GFP25S,WT35S,PARN23S,ELAC33S,GFP15S,WT25S,PARN33S,ELAC15S,PARN15S,ELAC25S,GFP35S,
#                    WT13S,PARN25S,GFP23S,PARN35S,ELAC13S,ELAC35S,GFP13S,WT23S,WT33S)




#scnRNA_all_sets_corrected=subset(scnRNA_all_sets_corrected,scnRNA_all_sets_corrected$strand==0)
#scnRNA_all_sets_corrected=subset(scnRNA_all_sets_corrected,scnRNA_all_sets_corrected$all_mm<4)


table(unique_sorted_diff_RNA_type_corrected$correct_RNA_type)
for (i in 1:length(filtered_sets)){
  print(paste("Proceeding file",i))
  tableset=filtered_sets[[i]]
  #tableset=subset(tableset,tableset$strand==0)
  #tableset=subset(tableset,tableset$all_mm<4)
  new_set_merged=merge(tableset,unique_sorted_diff_RNA_type_corrected,by="orig_seq",all.x = TRUE)
  #tableset$RNA_pos_uniq_XM_MD_cigar=paste(tableset$RetainedRef,tableset$position,tableset$XM,tableset$MD,tableset$cigar, sep=" " )
  #tableset=tableset[c("orig_seq","RNA_pos_uniq_XM_MD_cigar")]
  #new_set_merged=merge(tableset,new_sorted_dupl_seq_set_data,by="orig_seq",all.x = TRUE)
  #new_set_merged$RNA_pos=ifelse(is.na(new_set_merged$RNA_pos_uniq_XM_MD_cigar_dupl),new_set_merged$RNA_pos_uniq_XM_MD_cigar,new_set_merged$RNA_pos_uniq_XM_MD_cigar_dupl)
  
  assign(paste(original_set[i],"filtered",sep="_"),new_set_merged)
}
# colnames(PARN13S)
# PARN13S[c("seq","RetainedRef","position")]
save(ELAC23S_filtered,GFP33S_filtered,WT15S_filtered,GFP25S_filtered,WT35S_filtered,ELAC33S_filtered,
     GFP15S_filtered,WT25S_filtered,ELAC15S_filtered,ELAC25S_filtered,GFP35S_filtered,WT13S_filtered,
     GFP23S_filtered,ELAC13S_filtered,ELAC35S_filtered,GFP13S_filtered,WT23S_filtered,WT33S_filtered, 
     file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_ELAC_filtered_no_genome_good_anno.RData') #transcriptome and scnRNA (both nuclear and mtDNA tRNA) data
table(ELAC23S_filtered$correct_RNA_type)

load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_ELAC_filtered_no_genome_good_anno.RData')
################################################################################
#load genome data
original_set=c("PARN13S","ELAC23S","GFP33S","WT15S","GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S",
               "ELAC15S","PARN15S","ELAC25S","GFP35S","WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")
for (j in 1:length(sample_list)) {
  #mapped_seq=subset(mapped_seq,!(mapped_seq$read %in% mapped_seq_tRNA$read))
  print(paste("Proceeding file",j))
  mapped_seq_all=read.table(paste("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_genome/mapped_seq_with_strand/",genome_list[j],
                              sep=""),col.names = c("read","strand","ref","position","qual","cigar","seq","XM","MD"))
  mapped_seq_all$RNA_type="genome"
  mapped_seq_all$cigar_interpr=cigarOpTable(mapped_seq_all$cigar)
  mapped_seq_all$cigar_del=mapped_seq_all$cigar_interpr[,3]
  mapped_seq_all$cigar_ins=mapped_seq_all$cigar_interpr[,2]
  mapped_seq_all$cigar_match=mapped_seq_all$cigar_interpr[,1]
  mapped_seq_all$cigar_mismatch=(mapped_seq_all$cigar_interpr[,2]+mapped_seq_all$cigar_interpr[,3]+mapped_seq_all$cigar_interpr[,4]+mapped_seq_all$cigar_interpr[,5]+mapped_seq_all$cigar_interpr[,6])
  mapped_seq_all$orig_seq=ifelse(mapped_seq_all$strand==0,mapped_seq_all$seq,reverseComplement(mapped_seq_all$seq,case="upper"))
  mapped_seq_all$new_XM=as.numeric(substr(mapped_seq_all$XM,6,nchar(mapped_seq_all$XM)))
  mapped_seq_all$all_mm=as.numeric(mapped_seq_all$new_XM+mapped_seq_all$cigar_mismatch)
  #new_set=subset(tableset,tableset$XM=="XM:i:0"|tableset$XM=="XM:i:1"|tableset$XM=="XM:i:2"|tableset$XM=="XM:i:3")
  #new_set$filter=ifelse(new_set$all_mm>3, "low_quality_no_tail", "good_quality_or_has_tail")
  mapped_seq_all=subset(mapped_seq_all,mapped_seq_all$all_mm<4)
  mapped_seq_all=subset(mapped_seq_all,mapped_seq_all$strand==0)
  #new_set$DuplicateRef=new_set$ref
  #new_set_merged=merge(new_set,duplicates,by="DuplicateRef",all.x = TRUE)
  #new_set_merged$RetainedRef=ifelse(is.na(new_set_merged$RetainedRef),new_set_merged$DuplicateRef,new_set_merged$RetainedRef)
  #assign(paste(original_set[j]),mapped_seq_all)
  assign(paste(original_set[j],"genome",sep="_"),mapped_seq_all)
  #assign(paste("smallRNA",j,sep=""),mapped_seq_all)
}
save(ELAC23S_genome,GFP33S_genome,WT15S_genome,GFP25S_genome,WT35S_genome,ELAC33S_genome,
     GFP15S_genome,WT25S_genome,ELAC15S_genome,ELAC25S_genome,GFP35S_genome,WT13S_genome,
     GFP23S_genome,ELAC13S_genome,ELAC35S_genome,GFP13S_genome,WT23S_genome,WT33S_genome,
     PARN13S_genome,PARN23S_genome,PARN33S_genome,PARN15S_genome,PARN25S_genome,PARN35S_genome,
     file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_all_just_genome.RData')
load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_all_just_genome.RData')
colnames(ELAC23S_filtered)
colnames(ELAC23S_genome)
filtered_sets_no_genome=list(ELAC23S_filtered,GFP33S_filtered,WT15S_filtered,GFP25S_filtered,WT35S_filtered,ELAC33S_filtered,
                             GFP15S_filtered,WT25S_filtered,ELAC15S_filtered,ELAC25S_filtered,GFP35S_filtered,WT13S_filtered,
                             GFP23S_filtered,ELAC13S_filtered,ELAC35S_filtered,GFP13S_filtered,WT23S_filtered,WT33S_filtered)
filtered_sets_genome=list(ELAC23S_genome,GFP33S_genome,WT15S_genome,GFP25S_genome,WT35S_genome,ELAC33S_genome,
                          GFP15S_genome,WT25S_genome,ELAC15S_genome,ELAC25S_genome,GFP35S_genome,WT13S_genome,
                          GFP23S_genome,ELAC13S_genome,ELAC35S_genome,GFP13S_genome,WT23S_genome,WT33S_genome)

original_set=c("ELAC23S","GFP33S","WT15S","GFP25S","WT35S","ELAC33S","GFP15S","WT25S","ELAC15S","ELAC25S","GFP35S",
               "WT13S","GFP23S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")

for (i in 1:length(filtered_sets_genome)){
  print(paste("Proceeding file",i))
  tableset_no_genome=filtered_sets_no_genome[[i]]
  tableset_genome=filtered_sets_genome[[i]]
  tableset_genome=subset(tableset_genome,!(tableset_genome$read %in% tableset_no_genome$read))
  tableset_genome$set=original_set[i]
  tableset_genome$correct_annotation=tableset_genome$ref
  tableset_genome$correct_RNA_type="genome"
  tableset_all=rbind(tableset_no_genome,tableset_genome)
  #tableset=subset(tableset,tableset$strand==0)
  #tableset=subset(tableset,tableset$all_mm<4)
  #new_set_merged=merge(tableset,unique_sorted_diff_RNA_type_corrected,by="orig_seq",all.x = TRUE)
  #tableset$RNA_pos_uniq_XM_MD_cigar=paste(tableset$RetainedRef,tableset$position,tableset$XM,tableset$MD,tableset$cigar, sep=" " )
  #tableset=tableset[c("orig_seq","RNA_pos_uniq_XM_MD_cigar")]
  #new_set_merged=merge(tableset,new_sorted_dupl_seq_set_data,by="orig_seq",all.x = TRUE)
  #new_set_merged$RNA_pos=ifelse(is.na(new_set_merged$RNA_pos_uniq_XM_MD_cigar_dupl),new_set_merged$RNA_pos_uniq_XM_MD_cigar,new_set_merged$RNA_pos_uniq_XM_MD_cigar_dupl)
  
  assign(paste(original_set[i],"filtered_with_genome",sep="_"),tableset_all)
}
save(ELAC23S_filtered_with_genome,GFP33S_filtered_with_genome,WT15S_filtered_with_genome,GFP25S_filtered_with_genome,WT35S_filtered_with_genome,ELAC33S_filtered_with_genome,
     GFP15S_filtered_with_genome,WT25S_filtered_with_genome,ELAC15S_filtered_with_genome,ELAC25S_filtered_with_genome,GFP35S_filtered_with_genome,WT13S_filtered_with_genome,
     GFP23S_filtered_with_genome,ELAC13S_filtered_with_genome,ELAC35S_filtered_with_genome,GFP13S_filtered_with_genome,WT23S_filtered_with_genome,WT33S_filtered_with_genome,
     file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_ELAC_filtered_with_genome_good_anno.RData') #transcriptome,genome and scnRNA (both nuclear and mtDNA tRNA) data



# save(PARN13S_filtered,ELAC23S_filtered,GFP33S_filtered,WT15S_filtered,GFP25S_filtered,WT35S_filtered,PARN23S_filtered,ELAC33S_filtered,
#      GFP15S_filtered,WT25S_filtered,PARN33S_filtered,ELAC15S_filtered,PARN15S_filtered,ELAC25S_filtered,GFP35S_filtered,WT13S_filtered,
#      PARN25S_filtered,GFP23S_filtered,PARN35S_filtered,ELAC13S_filtered,ELAC35S_filtered,GFP13S_filtered,WT23S_filtered,WT33S_filtered, 
#      file = 'F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_new_new_all_filtered.RData') #transcriptome and scnRNA (both nuclear and mtDNA tRNA) data
ELAC13S_filtered_with_genome$set="ELAC3S"
ELAC23S_filtered_with_genome$set="ELAC3S"
ELAC33S_filtered_with_genome$set="ELAC3S"

ELAC15S_filtered_with_genome$set="ELAC5S"
ELAC25S_filtered_with_genome$set="ELAC5S"
ELAC35S_filtered_with_genome$set="ELAC5S"

GFP13S_filtered_with_genome$set="GFP3S"
GFP23S_filtered_with_genome$set="GFP3S"
GFP33S_filtered_with_genome$set="GFP3S"

GFP15S_filtered_with_genome$set="GFP5S"
GFP25S_filtered_with_genome$set="GFP5S"
GFP35S_filtered_with_genome$set="GFP5S"

WT13S_filtered_with_genome$set="WT3S"
WT23S_filtered_with_genome$set="WT3S"
WT33S_filtered_with_genome$set="WT3S"

WT15S_filtered_with_genome$set="WT5S"
WT25S_filtered_with_genome$set="WT5S"
WT35S_filtered_with_genome$set="WT5S"
all_sets_table=rbind(ELAC23S_filtered_with_genome,GFP33S_filtered_with_genome,WT15S_filtered_with_genome,GFP25S_filtered_with_genome,WT35S_filtered_with_genome,ELAC33S_filtered_with_genome,
                     GFP15S_filtered_with_genome,WT25S_filtered_with_genome,ELAC15S_filtered_with_genome,ELAC25S_filtered_with_genome,GFP35S_filtered_with_genome,WT13S_filtered_with_genome,
                     GFP23S_filtered_with_genome,ELAC13S_filtered_with_genome,ELAC35S_filtered_with_genome,GFP13S_filtered_with_genome,WT23S_filtered_with_genome,WT33S_filtered_with_genome)
save(all_sets_table, file = "F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_ELAC_filtered_sets_table.RData")
# all_sets_table=rbind(ELAC23S,GFP33S,WT15S,GFP25S,WT35S,
#                      ELAC33S,GFP15S,WT25S,
#                      ELAC15S,ELAC25S,GFP35S, WT13S,
#                      GFP23S,ELAC13S,ELAC35S,
#                      GFP13S,WT23S,WT33S)
#check how many are reverse strand
#table(all_sets_table$strand)
#890809/(890809+20931385) 0.04% are on reverse strand
#all_sets_table=subset(all_sets_table,all_sets_table$strand==0)
################################################################################

# #for ELAC
# PARN_sets_table=rbind(ELAC13S,ELAC23S,ELAC33S,ELAC13S,ELAC25S,ELAC35S,
#                       GFP13S,GFP23S,GFP33S,GFP13S,GFP25S,GFP35S)

#


#save(all_sets_table, file = "F:/PARN_ELAC_silencing/smallRNA/all_sets_table.RData")
#load("F:/PARN_ELAC_silencing/smallRNA/all_sets_table.RData")
#table(all_sets_table$RNA_type)
#unique(all_sets_table$ref[all_sets_table$RNA_type=="other"])
colnames(all_sets_table)
#all_sets_table$RNA_type=ifelse(str_detect(all_sets_table$ref, regex("spliceosomal", ignore_case = TRUE)),"snRNA",all_sets_table$RNA_type)
all_sets_table_df=as.data.frame(table(all_sets_table$set,all_sets_table$correct_RNA_type))
colnames(all_sets_table_df)=c("Comparison","Expression","Frequency")
all_sets_table_sum=aggregate(Frequency~Comparison, all_sets_table_df,sum)
colnames(all_sets_table_sum)=c("Comparison","sum")
all_sets_table_df=merge(all_sets_table_df,all_sets_table_sum,by="Comparison",all.x = TRUE)
all_sets_table_df$perc=round((all_sets_table_df$Frequency/all_sets_table_df$sum)*100, digits = 2)
colnames(all_sets_table_df)=c("Set","RNA_type","Frequency","sum","perc")
#write.csv(all_sets_table_df,"F:/PARN_ELAC_silencing/smallRNA/RNA_type_perc.csv")
unique(all_sets_table_df$RNA_type)
#levels(all_sets_table_df$RNA_type)
# all_sets_table_df_copy=all_sets_table_df
# all_sets_table_df=all_sets_table_df_copy
#all_sets_table_df$RNA_type=as.factor(all_sets_table_df$RNA_type)
all_sets_table_df$RNA_type <- as.character(all_sets_table_df$RNA_type)
#levels(all_sets_table_df$RNA_type)
#all_sets_table_df$RNA_type=ifelse(all_sets_table_df$RNA_type=="mRNA","mRNA fragments",all_sets_table_df$RNA_type)

#all_sets_table_df$RNA_type[all_sets_table_df$RNA_type == 'mRNA'] <- 'mRNA fragments'
# all_sets_table_df$RNA_type <- ifelse(all_sets_table_df$RNA_type == "mRNA", "mRNA fragments", all_sets_table_df$RNA_type)
all_sets_table_df$RNA_type <- factor(all_sets_table_df$RNA_type, levels=c("genome","piRNA","miRNA","mRNA fragments","other",
                                                                          "rRNA fragments","snoRNA fragments","tRNA fragments","snRNA fragments","lncRNA fragments","multiple_hits"))
# all_sets_table_df$RNA_type <- factor(all_sets_table_df$RNA_type, levels=c("genome","piRNA","miRNA","mRNA fragments","other",
#                                                                           "rRNA fragments","snoRNA fragments","tRNA fragments","snRNA fragments"))
all_sets_table_df$Set=factor(all_sets_table_df$Set, levels=c("ELAC3S","GFP3S","WT3S","ELAC5S","GFP5S","WT5S"))
ggplot(all_sets_table_df, aes(x = 3, y = perc, fill = RNA_type)) +
  geom_col(color = "white") +
  geom_text_repel(aes(label = paste (perc, "%")),
                  position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  xlim(c(0.2, 3.5)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  facet_wrap(.~ Set,nrow=2) +theme_void()+theme(
    plot.title = element_text(size=20, face="bold",hjust = 0.5),
    text = element_text(size=25))



head(all_sets_table)
all_sets_table$length=nchar(all_sets_table$seq)
#nchar(PARN_sets_table$seq[1])

#length distribution
?geom_histogram
all_sets_table$correct_RNA_type <- as.character(all_sets_table$correct_RNA_type)
#all_sets_table$RNA_type[all_sets_table$RNA_type == 'mRNA'] <- 'mRNA fragments'
#PARN_sets_table=subset(PARN_sets_table,PARN_sets_table$length<85)
(p <- all_sets_table %>%
    mutate(correct_RNA_type = fct_reorder(correct_RNA_type, length)) %>%
    ggplot( aes(x=length, color=correct_RNA_type, fill=correct_RNA_type)) +
    geom_histogram(alpha=0.6, binwidth = 1) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    xlab("length") +
    ylab("number of RNA species ") +
    scale_y_continuous(labels = scales::scientific)+
    facet_wrap(~correct_RNA_type, scales = "free")+
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 20,hjust=0.5),
      axis.title.x = element_text(size = 15,hjust=0.5,face = "plain"),
      axis.title.y = element_text(size = 15,hjust=0.5,face = "plain"))+
    scale_x_continuous(breaks = seq(15, 115, by = 10)))
?element_text

# all_sets_table_trna=subset(all_sets_table,all_sets_table$RNA_type=="tRNA fragments")
# length(unique(all_sets_table$orig_seq))
#PARN13S_filtered

#load("F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_analysis_filtered.RData")
#load("F:/smallRNA/calculated_data_bowtie2_end_to_end/all_sets_table.RData")
# load('F:/PARN_ELAC_silencing/smallRNA/sncRNAlist_filtered_new_new_all_filtered.RData')
# load("F:/PARN_ELAC_silencing/smallRNA/scnRNA_anno_uniq_corrected.RData")

# filtered_sets_for_table=list(PARN13S_filtered,ELAC23S_filtered,GFP33S_filtered,WT15S_filtered,GFP25S_filtered,WT35S_filtered,PARN23S_filtered,ELAC33S_filtered,
#                              GFP15S_filtered,WT25S_filtered,PARN33S_filtered,ELAC15S_filtered,PARN15S_filtered,ELAC25S_filtered,GFP35S_filtered,WT13S_filtered,
#                              PARN25S_filtered,GFP23S_filtered,PARN35S_filtered,ELAC13S_filtered,ELAC35S_filtered,GFP13S_filtered,WT23S_filtered,WT33S_filtered)
#colnames(PARN13S_filtered)
# original_set=c("PARN13S","ELAC23S","GFP33S","WT15S","GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S",
#                "ELAC15S","PARN15S","ELAC25S","GFP35S","WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S"
# )

filtered_sets_for_table=list(ELAC23S_filtered,GFP33S_filtered,WT15S_filtered,GFP25S_filtered,WT35S_filtered,ELAC33S_filtered,GFP15S_filtered,WT25S_filtered,ELAC15S_filtered,ELAC25S_filtered,GFP35S_filtered,
                             WT13S_filtered,GFP23S_filtered,ELAC13S_filtered,ELAC35S_filtered,GFP13S_filtered,WT23S_filtered,WT33S_filtered)
colnames(ELAC23S_filtered)
tail(ELAC23S_filtered)
for (i in 1:length(filtered_sets_for_table)) {
  print(paste("Proceeding file",i))
  tableset=filtered_sets_for_table[[i]]
  #tableset=merge(tableset,scnRNA_anno_uniq_corrected,by="orig_seq",all.x = TRUE)
  tableset$seq_anno=paste(tableset$orig_seq,tableset$correct_annotation,sep = " ")
  set_table=as.data.frame(table(tableset$seq_anno))
  colnames(set_table)=c("Var1",paste(original_set[i]))
  assign(paste(original_set[i],"table",sep="_"),set_table)
}
# all_tables=list(PARN13S_table,ELAC23S_table,GFP33S_table,WT15S_table,GFP25S_table,WT35S_table,PARN23S_table,ELAC33S_table,GFP15S_table,WT25S_table,
#                 PARN33S_table,ELAC15S_table,PARN15S_table,ELAC25S_table,GFP35S_table,WT13S_table,PARN25S_table,GFP23S_table,
#                 PARN35S_table,ELAC13S_table,ELAC35S_table,GFP13S_table,WT23S_table,WT33S_table)

all_tables=list(ELAC23S_table,GFP33S_table,WT15S_table,GFP25S_table,WT35S_table,ELAC33S_table,GFP15S_table,WT25S_table,
                ELAC15S_table,ELAC25S_table,GFP35S_table,WT13S_table,GFP23S_table,
                ELAC13S_table,ELAC35S_table,GFP13S_table,WT23S_table,WT33S_table)


MyMerge <- function(x, y){
  df<- merge(x, y, by= "Var1", all=TRUE)
  return(df)
}
dat <- Reduce(MyMerge, all_tables)
#sum(dat$PARN13S)
#sum(dat$ELAC23S)
#?write.csv()
dat[is.na(dat)] <- 0
#write.csv(dat,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/corrected_merged_raw_counts_with_tags_end_to_end.csv",row.names=FALSE,quote=FALSE)
#dat=read.csv("F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/corrected_merged_raw_counts_with_tags_end_to_end.csv")
#nrow(PARN13S_table) #186560
nrow(dat) #7207365 new-3101937
length(unique(dat$Var1))#3101937
#dat_table=as.data.frame(table(dat$Var1))
#dat_table_dupl=subset(dat_table,dat_table$Freq>1)
#replace NA with zeros


#new_dat$Var1
# new_dat=dat
# new_dat[c('Sequence', 'Annotation', 'Position')] <- str_split_fixed(new_dat$Var1, ' ', 3)
# colnames(new_dat)
# #[1] "Var1"       "PARN13S"    "ELAC23S"    "GFP33S"     "WT15S"      "GFP25S"     "WT35S"      "PARN23S"    "ELAC33S"    "GFP15S"     "WT25S"      "PARN33S"    "ELAC15S"    "PARN15S"    "ELAC25S"    "GFP35S"    
# #[17] "WT13S"      "PARN25S"    "GFP23S"     "PARN35S"    "ELAC13S"    "ELAC35S"    "GFP13S"     "WT23S"      "WT33S"      "Sequence"   "Annotation" "Position"  
# new_dat=new_dat[c("Sequence","Position","Annotation","PARN13S","ELAC23S","GFP33S","WT15S" ,"GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S","ELAC15S","PARN15S","ELAC25S","GFP35S",
#                   "WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")]
##################################################################################
#DGE analysis

#count matrix data
typeof(dat)
class(dat)
ncol(dat)
head(dat)
smallRNA_counts=apply(as.matrix.noquote(dat[c(2:19)]),2,as.numeric) #create a numeric matrix
rownames(smallRNA_counts)=dat$Var1
colnames(smallRNA_counts)
#save(smallRNA_counts,file='F:/PARN_ELAC_silencing/smallRNA/smallRNA_counts_new.RData')
################################################################################
#EdgeR
myCPM_new=cpm(smallRNA_counts)
rownames(myCPM_new)
#load("F:/smallRNA/calculated_data_bowtie2_end_to_end/raw_cpm.RData")
#save(myCPM, file = "F:/smallRNA/calculated_data_bowtie2_end_to_end/raw_cpm.RData")



# #load old cpm table
load("G:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/raw_cpm.RData")
myCPM
thresh <- myCPM > 10
#thresh <- myCPM > 1
table(rowSums(thresh)) ## There are 2872 RNAs that have TRUEs in all 18 samples.
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep_new <- rowSums(thresh) >= 2
#keep <- rowSums(thresh) >= 3
#counts.keep_new <- smallRNA_counts[keep_new,]
cpm.keep_new <- myCPM[keep_new,]
# thresh_old <- myCPM > 10
# #thresh <- myCPM > 1
# table(rowSums(thresh_old))
# table(rowSums(thresh))## There are 2872 RNAs that have TRUEs in all 18 samples.
# # we would like to keep genes that have at least 3 TRUES in each row of thresh
# keep <- rowSums(thresh_old) >= 2
# #keep <- rowSums(thresh) >= 3
# counts.keep <- smallRNA_counts[keep,]
# cpm.keep <- myCPM[keep,]


#save(cpm.keep,file='F:/PARN_ELAC_silencing/smallRNA/cpm_keep.RData')
# plot(cpm.keep[,1],counts.keep[,1])
# plot(cpm.keep[,2],counts.keep[,2])
# plot(cpm.keep[,3],counts.keep[,3])
#write.csv(cpm.keep,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/cpm_filtered.csv",row.names=TRUE,quote=FALSE)
#write.csv(cpm.keep,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/cpm_filtered_with_tags_end_to_end.csv",row.names=TRUE,quote=FALSE)
#write.csv(cpm.keep,"F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/cpm_filtered_with_tags_end_to_end_cpm10.csv",row.names=TRUE,quote=FALSE)
#cpm.keep=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/cpm_filtered_with_tags_end_to_end.csv")
#cpm(cpm.keep)
# rownames(cpm.keep)
# cpm.keep[1,]
# colnames(cpm.keep)
# cpm.keep[,1]
################################################################################
#---------------------------RNA accumulation------------------------------------
RNA_accum=as.data.frame(cpm.keep_new)
rownames(RNA_accum)

RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("ribosomal", ignore_case = TRUE)),"rRNA fragments","other fragments")
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("rRNA", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("tRNA", ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("lncRNA", ignore_case = TRUE)),"lncRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("piRNA", ignore_case = TRUE)),"piRNA",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("RNase_P_RNA", ignore_case = TRUE)),"RNase_P_RNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("snoRNA", ignore_case = TRUE)),"snoRNA fragments",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("siRNA", ignore_case = TRUE)),"siRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("miRNA", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("pre_miRNA", ignore_case = TRUE)),"pre_miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("snRNA", ignore_case = TRUE)),"snRNA fragments",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",RNA_accum$RNA_type)
#RNA_accum$RNA_type=ifelse(str_detect(RNA_accum$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("Sme-", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)  
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("spliceosomal-", ignore_case = TRUE)),"snRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("spliceosomal", ignore_case = TRUE)),"snRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("lnc", ignore_case = TRUE)),"lncRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("microRNA", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("long", ignore_case = TRUE)),"lncRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("7SK", ignore_case = TRUE)),"snRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("nucleolar", ignore_case = TRUE)),"snoRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("nuclear", ignore_case = TRUE)),"snRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("miR", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("transfer", ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("dd_Smed_v6", ignore_case = TRUE)),"mRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("dd_Smes", ignore_case = TRUE)),"genome",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("sme-lin", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("sme-let", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("multiple", ignore_case = TRUE)),"multiple_hits",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("ITS1", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("ITS2", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("SpacerA", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("28S", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("Schmed_cloneH735c", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("12S", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("16S", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex("5.8S", ignore_case = TRUE)),"rRNA fragments",RNA_accum$RNA_type)

###############################################################################
#select miRNA reference:
miRNA_set=subset(RNA_accum,RNA_accum$RNA_type=="miRNA")
nrow(miRNA_set)
colnames(miRNA_set)
#miRNA_set$var=var(miRNA_set[,1:24])
#(mad(as.numeric(DGE_set_with_norm_counts[1,ELAC3]))/median(as.numeric(DGE_set_with_norm_counts[1,ELAC3]))) * 100
for (i in 1: nrow(miRNA_set)) {
  miRNA_set$var_coef[i]=(mad(as.numeric(miRNA_set[i,1:24]))/median(as.numeric(miRNA_set[i,1:24]))) * 100
}
summary(miRNA_set$var_coef)

constant_miRNA_set=subset(miRNA_set,miRNA_set$var_coef<20)

constant_miRNA_set_low=subset(constant_miRNA_set,constant_miRNA_set$PARN13S<20)
constant_miRNA_set_low_med=subset(constant_miRNA_set,constant_miRNA_set$PARN13S<300 & constant_miRNA_set$PARN13S>20)
constant_miRNA_set_med=subset(constant_miRNA_set,constant_miRNA_set$PARN13S<1000  & constant_miRNA_set$PARN13S>300)
constant_miRNA_set_high=subset(miRNA_set,miRNA_set$PARN13S>10000)










###############################################################################
subset(RNA_accum,RNA_accum$RNA_type=="other fragments")
subset(RNA_accum,RNA_accum$RNA_type=="multiple_hits")
SM_other_to_identify=RNA_accum[endsWith(rownames(RNA_accum),"79327"),]
(urs <- sub(".*?(URS\\w+_\\d+).*", "\\1", rownames(SM_other_to_identify)))
uniq_urs=unique(urs)
uniq_urs_df=as.data.frame(uniq_urs)
colnames(uniq_urs_df)
for(i in 1:length(uniq_urs)){
  rnaCentralEntry <- rnaCentralRetrieveEntry(uniq_urs[i])
  uniq_urs_df$correct_annotation[i]=paste(rnaCentralEntry$description)
  #print(rnaCentralEntry)
}
uniq_urs_df$RNA_type=ifelse(str_detect(uniq_urs_df$correct_annotation, regex("tRNA", ignore_case = TRUE)),"tRNA fragments","other")
uniq_urs_df$RNA_type=ifelse(str_detect(uniq_urs_df$correct_annotation, regex("mir", ignore_case = TRUE)),"miRNA",uniq_urs_df$RNA_type)
#RNA_accum[(str_detect(rownames(RNA_accum), "79327"))& RNA_accum$RNA_type=="other fragments",]
?rnaCentralRetrieveEntry
SM_tRNA_weird=uniq_urs_df$uniq_urs[uniq_urs_df$RNA_type=="tRNA fragments"]
SM_miRNA_weird=uniq_urs_df$uniq_urs[uniq_urs_df$RNA_type=="miRNA"]
length(SM_tRNA_weird)
SM_tRNA_weird[1]
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_miRNA_weird, ignore_case = TRUE)),"miRNA",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[1], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[2], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[3], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[4], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[5], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[6], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[7], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[8], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[9], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[10], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[11], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[12], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)
RNA_accum$RNA_type=ifelse(str_detect(rownames(RNA_accum), regex(SM_tRNA_weird[13], ignore_case = TRUE)),"tRNA fragments",RNA_accum$RNA_type)


table(RNA_accum$RNA_type)












#RNA_accum=subset(RNA_accum,RNA_accum$RNA_type!="genome")
subset(ELAC23S_filtered,ELAC23S_filtered$orig_seq=="AATGTAGTAGTATAACGAGTATATTTAGTTTACATCTAAAAGGTATTGAT")


subset(RNA_accum,RNA_accum$RNA_type=="mRNA fragments")
#---------------------------------------check miRNA length----------------------
miRNA_cpm=subset(RNA_accum,RNA_accum$RNA_type=="miRNA")
nrow(miRNA_cpm)
colnames(miRNA_cpm)
#miRNA_cpm$PARN_3dpa=(miRNA_cpm$PARN13S+miRNA_cpm$PARN23S+miRNA_cpm$PARN33S)
#miRNA_cpm$PARN_5dpa=(miRNA_cpm$PARN15S+miRNA_cpm$PARN25S+miRNA_cpm$PARN35S)
miRNA_cpm$GFP_3dpa=(miRNA_cpm$GFP13S+miRNA_cpm$GFP23S+miRNA_cpm$GFP33S)
miRNA_cpm$GFP_5dpa=(miRNA_cpm$GFP15S+miRNA_cpm$GFP25S+miRNA_cpm$GFP35S)
miRNA_cpm$ELAC_3dpa=(miRNA_cpm$ELAC13S+miRNA_cpm$ELAC23S+miRNA_cpm$ELAC33S)
miRNA_cpm$ELAC_5dpa=(miRNA_cpm$ELAC15S+miRNA_cpm$ELAC25S+miRNA_cpm$ELAC35S)
miRNA_cpm$WT_3dpa=(miRNA_cpm$WT13S+miRNA_cpm$WT23S+miRNA_cpm$WT33S)
miRNA_cpm$WT_5dpa=(miRNA_cpm$WT15S+miRNA_cpm$WT25S+miRNA_cpm$WT35S)
miRNA_cpm=miRNA_cpm[c("GFP_3dpa","GFP_5dpa","ELAC_3dpa","ELAC_5dpa","WT_3dpa","WT_5dpa","RNA_type")]
miRNA_cpm[c('Sequence', 'Annotation', 'Position', 'Mismatch', 'Cigar', 'Match')] <- str_split_fixed(rownames(miRNA_cpm), ' ', 6)
miRNA_cpm$length=nchar(miRNA_cpm$Sequence)
miRNA_cpm_long <- gather(miRNA_cpm, condition, cpm, GFP_3dpa:WT_5dpa, factor_key=TRUE)
aggregate(cpm~condition,miRNA_cpm_long,sum)
miRNA_cpm_long=transform(miRNA_cpm_long, percent = ave(cpm, condition, FUN = prop.table))
miRNA_cpm_long$condition=as.character(miRNA_cpm_long$condition)
miRNA_cpm_long$dpa=substr(miRNA_cpm_long$condition,nchar(miRNA_cpm_long$condition)-3,nchar(miRNA_cpm_long$condition))
miRNA_cpm_long=subset(miRNA_cpm_long,miRNA_cpm_long$length<27)
miRNA_cpm_long$gene=substr(miRNA_cpm_long$condition,1,nchar(miRNA_cpm_long$condition)-5)
# miRNA_cpm_long_dpa3=subset(miRNA_cpm_long,miRNA_cpm_long$dpa=="3dpa" & miRNA_cpm_long$length<27)
# miRNA_cpm_long_dpa5=subset(miRNA_cpm_long,miRNA_cpm_long$dpa=="5dpa" & miRNA_cpm_long$length<27)
miRNA_cpm_long$mononucleotide_tail=substr(miRNA_cpm_long$Sequence,nchar(miRNA_cpm_long$Sequence),nchar(miRNA_cpm_long$Sequence))
miRNA_cpm_long$dinucleotide_tail=substr(miRNA_cpm_long$Sequence,nchar(miRNA_cpm_long$Sequence)-1,nchar(miRNA_cpm_long$Sequence))
miRNA_cpm_long_sum=aggregate(percent~gene+dpa+length,miRNA_cpm_long,sum)

ggplot(miRNA_cpm_long_sum, aes(fill=gene, y=percent, x=length)) + 
  geom_bar(position="dodge", stat="identity")+theme_ipsum()+
  theme(text = element_text(size = 15),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  ggtitle("miRNA length distribution")+facet_grid(.~dpa)


miRNA_cpm_long_sum_lastnucl=aggregate(percent~gene+dpa+mononucleotide_tail,miRNA_cpm_long,sum)

ggplot(miRNA_cpm_long_sum_lastnucl, aes(fill=gene, y=percent, x=mononucleotide_tail)) + 
  geom_bar(position="dodge", stat="identity")+theme_ipsum()+
  theme(text = element_text(size = 15),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  ggtitle("miRNA mononucleotide tail")+facet_grid(.~dpa)
miRNA_cpm_long_sum_lastdinucl=aggregate(percent~gene+dpa+dinucleotide_tail,miRNA_cpm_long,sum)

ggplot(miRNA_cpm_long_sum_lastdinucl, aes(fill=gene, y=percent, x=dinucleotide_tail)) + 
  geom_bar(position="dodge", stat="identity")+theme_ipsum()+
  theme(text = element_text(size = 15),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  ggtitle("miRNA dinucleotide tail")+facet_grid(.~dpa)

#write.csv(miRNA_cpm_long,"F:/smallRNA/calculated_data_bowtie2_end_to_end/miRNA_cpm_filtered.csv",row.names=FALSE,quote=FALSE)
miRNA_filtered=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/miRNA_cpm_filtered.csv")


################################################################################
#accumulation of RNA species for ELAC
ELAC_accum=RNA_accum[c("ELAC13S","ELAC23S","ELAC33S","ELAC15S","ELAC25S","ELAC35S",
                       "GFP13S","GFP23S","GFP33S","GFP15S","GFP25S","GFP35S",
                       "WT13S","WT23S","WT33S","WT15S","WT25S","WT35S",
                       "RNA_type")]
ELAC_accum$ELAC_3dpa=(ELAC_accum$ELAC13S+ELAC_accum$ELAC23S+ELAC_accum$ELAC33S)
ELAC_accum$ELAC_5dpa=(ELAC_accum$ELAC15S+ELAC_accum$ELAC25S+ELAC_accum$ELAC35S)
ELAC_accum$GFP_3dpa=(ELAC_accum$GFP13S+ELAC_accum$GFP23S+ELAC_accum$GFP33S)
ELAC_accum$GFP_5dpa=(ELAC_accum$GFP15S+ELAC_accum$GFP25S+ELAC_accum$GFP35S)
ELAC_accum$WT_3dpa=(ELAC_accum$WT13S+ELAC_accum$WT23S+ELAC_accum$WT33S)
ELAC_accum$WT_5dpa=(ELAC_accum$WT15S+ELAC_accum$WT25S+ELAC_accum$WT35S)
ELAC_accum=ELAC_accum[c("ELAC_3dpa","ELAC_5dpa","GFP_3dpa","GFP_5dpa","WT_3dpa","WT_5dpa","RNA_type")]
colnames(ELAC_accum)=c("ELAC2 3dpa","ELAC2 5dpa","GFP 3dpa","GFP 5dpa","WT 3dpa","WT 5dpa","RNA type")
#ELAC_accum$`RNA type`=ifelse(ELAC_accum$`RNA type`=="other","other fragments",ELAC_accum$`RNA type`)

ELAC_long <- melt(ELAC_accum, id = "RNA type", variable_name = "variable")
colnames(ELAC_long)[3] ="cpm"
colnames(ELAC_long)[2] ="Gene"
ELAC_accum_sum=aggregate(cpm~`RNA type`+Gene, ELAC_long,sum)
ELAC_gene_sum=aggregate(cpm~Gene, ELAC_accum_sum,sum)
colnames(ELAC_gene_sum)=c("Gene","sum")
ELAC_table_df=merge(ELAC_accum_sum,ELAC_gene_sum,by="Gene",all.x = TRUE)
ELAC_table_df$perc=round((ELAC_table_df$cpm/ELAC_table_df$sum)*100, digits = 2)
ELAC_table_df$label=paste(ELAC_table_df$perc,"%",sep=" ")
ELAC_table_df$`RNA type`
unique(ELAC_table_df$`RNA type`)
ELAC_table_df$`RNA type` <- factor(ELAC_table_df$`RNA type`,
                                   levels = c("other fragments","mRNA fragments","lncRNA fragments",
                                              "tRNA fragments","snoRNA fragments","rRNA fragments",
                                              "snRNA fragments","miRNA","piRNA"
                                              ), ordered = TRUE)
ELAC_table_df
ELAC_table_df$Gene=as.character(ELAC_table_df$Gene)
ELAC_table_df$dpa=substr(ELAC_table_df$Gene,nchar(ELAC_table_df$Gene)-3,nchar(ELAC_table_df$Gene))
ELAC_table_df$gene=substr(ELAC_table_df$Gene,1,nchar(ELAC_table_df$Gene)-5)
new_ELAC_table=ELAC_table_df %>% arrange(Gene,desc(`RNA type`)) %>%
  group_by(Gene) %>% 
  mutate(text_y = cumsum(perc) - perc/2)
unique(new_ELAC_table$`RNA type`)
save(new_ELAC_table,file="D:/Illumina_silencing/final_results/new_ELAC_table_snRNA.RData")
save(ELAC_table_df,file="D:/Illumina_silencing/final_results/ELAC_table_snRNA.RData")
# new_ELAC_table$RNA_type <- factor(new_ELAC_table$RNA_type,
#                              levels = c("tRNA fragments","rRNA fragments" ,"mRNA fragments","miRNA","other fragments",
#                                         "snRNA fragments","snoRNA fragments","piRNA","lncRNA fragments"), ordered = TRUE)
#unique(new_ELAC_table$RNA_type)

# ?geom_label
# ggplot(data = new_ELAC_table, aes(x = "", y = perc, fill = `RNA type`))+ylab("total normalized read count") + 
#   geom_bar(stat = "identity") +
#   geom_label(aes(label = paste0(perc, "%"),y = text_y),show.legend = FALSE,label.size = 0.5,label.padding = unit(0.25, "lines")) +
#   facet_grid(.~Gene)+
#   scale_fill_brewer(palette = "Pastel1")+theme_void()+theme(
#     plot.title = element_text(size=20, face="bold",hjust = 0.5),
#     text = element_text(size=20))+labs(fill="RNA species")
ggplot(data = new_ELAC_table, aes(x = "", y = perc, fill = `RNA type`))+ylab("total normalized read count") +
  geom_bar(stat = "identity") +
  geom_label(aes(label = paste0(perc, "%"),y = text_y),show.legend = FALSE,label.size = 0.5,label.padding = unit(0.25, "lines")) +
  facet_nested(.~dpa+gene)+
  scale_fill_brewer(palette = "Pastel1")+theme_bw()+theme(
    plot.title = element_text(size=20, face="bold",hjust = 0.5),
    text = element_text(size=25))+labs(fill="RNA species",x="")

#for PARN and ELAC
# PARN_ELAC_accum=RNA_accum[c("PARN13S","PARN23S","PARN33S","PARN15S","PARN25S","PARN35S","ELAC13S","ELAC23S","ELAC33S","ELAC15S","ELAC25S","ELAC35S",
#                             "GFP13S","GFP23S","GFP33S","GFP15S","GFP25S","GFP35S","RNA_type")]
# 
# 
# PARN_ELAC_accum$ELAC_3dpa=(PARN_ELAC_accum$ELAC13S+PARN_ELAC_accum$ELAC23S+PARN_ELAC_accum$ELAC33S)
# PARN_ELAC_accum$ELAC_5dpa=(PARN_ELAC_accum$ELAC15S+PARN_ELAC_accum$ELAC25S+PARN_ELAC_accum$ELAC35S)
# PARN_ELAC_accum$PARN_3dpa=(PARN_ELAC_accum$PARN13S+PARN_ELAC_accum$PARN23S+PARN_ELAC_accum$PARN33S)
# PARN_ELAC_accum$PARN_5dpa=(PARN_ELAC_accum$PARN15S+PARN_ELAC_accum$PARN25S+PARN_ELAC_accum$PARN35S)
# PARN_ELAC_accum$GFP_3dpa=(PARN_ELAC_accum$GFP13S+PARN_ELAC_accum$GFP23S+PARN_ELAC_accum$GFP33S)
# PARN_ELAC_accum$GFP_5dpa=(PARN_ELAC_accum$GFP15S+PARN_ELAC_accum$GFP25S+PARN_ELAC_accum$GFP35S)
# PARN_ELAC_accum=PARN_ELAC_accum[c("ELAC_3dpa","ELAC_5dpa","PARN_3dpa","PARN_5dpa","GFP_3dpa","GFP_5dpa","RNA_type")]
# colnames(PARN_ELAC_accum)=c("ELAC2 3dpa","ELAC2 5dpa","PARN 3dpa","PARN 5dpa","GFP 3dpa","GFP 5dpa","RNA_type")
# #PARN_ELAC_accum$`RNA type`=ifelse(PARN_ELAC_accum$RNA_type=="other","other fragments",PARN_ELAC_accum$RNA_type)
# 
# ELAC_long <- melt(PARN_ELAC_accum, id = "RNA_type", variable_name = "Gene")
# colnames(ELAC_long)[3] ="cpm"
# ELAC_accum_sum=aggregate(cpm~RNA_type+Gene, ELAC_long,sum)
# ELAC_gene_sum=aggregate(cpm~Gene, ELAC_accum_sum,sum)
# colnames(ELAC_gene_sum)=c("Gene","sum")
# ELAC_table_df=merge(ELAC_accum_sum,ELAC_gene_sum,by="Gene",all.x = TRUE)
# ELAC_table_df$perc=round((ELAC_table_df$cpm/ELAC_table_df$sum)*100, digits = 2)
# ELAC_table_df$label=paste(ELAC_table_df$perc,"%",sep=" ")
# ELAC_table_df$RNA_type <- factor(ELAC_table_df$RNA_type,
#                                  levels = c("lncRNA fragments","tRNA fragments","other fragments","rRNA fragments" ,"snRNA fragments","mRNA fragments",
#                                             "snoRNA fragments","miRNA","piRNA"), ordered = TRUE)
# new_ELAC_table=ELAC_table_df %>% arrange(Gene,desc(RNA_type)) %>%
#   group_by(Gene) %>% 
#   mutate(text_y = cumsum(perc) - perc/2)
# unique(new_ELAC_table$RNA_type)
# # new_ELAC_table$RNA_type <- factor(new_ELAC_table$RNA_type,
# #                              levels = c("tRNA fragments","rRNA fragments" ,"mRNA fragments","miRNA","other fragments",
# #                                         "snRNA fragments","snoRNA fragments","piRNA","lncRNA fragments"), ordered = TRUE)
# unique(new_ELAC_table$RNA_type)
# 
# ?geom_label
# ggplot(data = new_ELAC_table, aes(x = "", y = perc, fill = RNA_type))+ylab("total normalized read count") + 
#   geom_bar(stat = "identity") +
#   geom_label(aes(label = paste0(perc, "%"),y = text_y),show.legend = FALSE,label.size = 0.5,label.padding = unit(0.25, "lines")) +
#   facet_grid(.~Gene)+
#   scale_fill_brewer(palette = "Pastel1")+theme_void()+theme(
#     plot.title = element_text(size=20, face="bold",hjust = 0.5),
#     text = element_text(size=20))+labs(fill="RNA species")
# colnames(new_ELAC_table)
# new_ELAC_table$Gene=as.character(new_ELAC_table$Gene)
# new_ELAC_table$gene=substr(new_ELAC_table$Gene,1,nchar(new_ELAC_table$Gene)-5)
# new_ELAC_table$dpa=substr(new_ELAC_table$Gene,nchar(new_ELAC_table$Gene)-3,nchar(new_ELAC_table$Gene))
# ggplot(data = new_ELAC_table, aes(x = "", y = perc, fill = RNA_type))+ylab("total normalized read count") + 
#   geom_bar(stat = "identity") +
#   geom_label(aes(label = paste0(perc, "%"),y = text_y),show.legend = FALSE,label.size = 0.5,label.padding = unit(0.25, "lines")) +
#   facet_nested(.~dpa+gene)+
#   scale_fill_brewer(palette = "Pastel1")+theme_bw()+theme(
#     plot.title = element_text(size=20, face="bold",hjust = 0.5),
#     text = element_text(size=25))+labs(fill="RNA species",x="")
# 
# ?facet_nested

################################################################################
#DESEq
#sample=c(sample_list,sample_list)

# myCPM_new=cpm(smallRNA_counts)
# rownames(myCPM_new)
# #load("F:/smallRNA/calculated_data_bowtie2_end_to_end/raw_cpm.RData")
# #save(myCPM, file = "F:/smallRNA/calculated_data_bowtie2_end_to_end/raw_cpm.RData")
# thresh <- myCPM_new > 10
# #thresh <- myCPM > 1
# table(rowSums(thresh)) ## There are 2872 RNAs that have TRUEs in all 18 samples.
# # we would like to keep genes that have at least 2 TRUES in each row of thresh
# keep_new <- rowSums(thresh) >= 2
# #keep <- rowSums(thresh) >= 3
# counts.keep_new <- smallRNA_counts[keep_new,]
# cpm.keep_new <- myCPM_new[keep_new,]

# keep <- (rowSums( counts(dds[,c(ELAC3)], normalized=TRUE) >= 10 ) >= 2) | (rowSums( counts(dds[,c(ELAC5)], normalized=TRUE) >= 10 ) >= 2) |
#   (rowSums( counts(dds[,c(WT3)], normalized=TRUE) >= 10 ) >= 2) | (rowSums( counts(dds[,c(WT5)], normalized=TRUE) >= 10 ) >= 2) |
#   (rowSums( counts(dds[,c(GFP3)], normalized=TRUE) >= 10 ) >= 2) | (rowSums( counts(dds[,c(GFP5)], normalized=TRUE) >= 10 ) >= 2) 
# 
# dds <- dds[keep,]



counts.keep_new
salmon_samples=as.data.frame(colnames(smallRNA_counts))
colnames(salmon_samples)="run"
#?str_extract
salmon_samples$gene=substr(salmon_samples$run,1,(nchar(salmon_samples$run)-3))
salmon_samples$dpa=paste("dpa",substr(salmon_samples$run,(nchar(salmon_samples$run)-1),(nchar(salmon_samples$run)-1)),sep="")
salmon_samples$replicate=paste("rep",substr(salmon_samples$run,(nchar(salmon_samples$run)-2),(nchar(salmon_samples$run)-2)),sep="")
salmon_samples$condition=paste(salmon_samples$gene,salmon_samples$dpa,sep="_")
rownames(salmon_samples)=salmon_samples$run

head(smallRNA_counts,2)
all(rownames(salmon_samples) %in% colnames(smallRNA_counts)) #TRUE
all(rownames(salmon_samples) == colnames(smallRNA_counts))
#load all raw counts
dds <- DESeqDataSetFromMatrix(countData = smallRNA_counts,
                              colData = salmon_samples,
                              design = ~ replicate+ condition)
#normalize
dds <- estimateSizeFactors(dds)
#select only genes that have at least 2 TRUES in each row of thresh
dds <- dds[keep_new,]

?DESeq
deseqddsColl_boot500=DESeq(dds,betaPrior=FALSE, minRep=Inf)
plotDispEsts(deseqddsColl_boot500)
#with all samples
rld <- rlog(deseqddsColl_boot500, blind=TRUE)
#DESeq2::plotPCA(rlogTransformation(deseqddsColl_boot500), intgroup=c("condition"))
(plotPCA_gene_old=plotPCA(rld, intgroup="condition",ntop=10000)+ geom_point(size = 4)+ggtitle("PCA plot for all samples.\n Low accumulated sncRNA were removed"))
#?DESeq2::plotPCA
#normalized counts
# normalized_counts=as.data.frame(counts(deseqddsColl_boot500, normalized=TRUE))
# colnames(normalized_counts)
# rownames(normalized_counts)
# class(normalized_counts)
# deseqddsColl_boot500@assays@data$counts
# normalized_counts[c('Sequence', 'Annotation', 'Position')] <- str_split_fixed(rownames(normalized_counts), ' ', 3)
# colnames(normalized_counts)
# new_normalized_counts=normalized_counts[c("Sequence","Position","Annotation")]
# new_normalized_counts$GFP3=mean(normalized_counts$GFP33S+normalized_counts$GFP13S+normalized_counts$GFP23S)
# new_normalized_counts$GFP5=(normalized_counts$GFP35S+normalized_counts$GFP15S+normalized_counts$GFP25S)
# new_normalized_counts$WT3=(normalized_counts$WT33S+normalized_counts$WT13S+normalized_counts$WT23S)
# new_normalized_counts$WT5=(normalized_counts$WT35S+normalized_counts$WT15S+normalized_counts$WT25S)
# new_normalized_counts$PARN3=(normalized_counts$PARN33S+normalized_counts$PARN13S+normalized_counts$PARN23S)
# new_normalized_counts$PARN5=(normalized_counts$PARN35S+normalized_counts$PARN15S+normalized_counts$PARN25S)
# new_normalized_counts$ELAC3=(normalized_counts$ELAC33S+normalized_counts$ELAC13S+normalized_counts$ELAC23S)
# new_normalized_counts$ELAC5=(normalized_counts$ELAC35S+normalized_counts$ELAC15S+normalized_counts$ELAC25S)
# # new_dat=new_dat[c("Sequence","Position","Annotation","PARN13S","ELAC23S","GFP33S","WT15S" ,"GFP25S","WT35S","PARN23S","ELAC33S","GFP15S","WT25S","PARN33S","ELAC15S","PARN15S","ELAC25S","GFP35S",
# #                   "WT13S","PARN25S","GFP23S","PARN35S","ELAC13S","ELAC35S","GFP13S","WT23S","WT33S")]
# ?write.csv
# write.csv(new_normalized_counts,"F:/smallRNAwithAdapters/miRNA/normalized_counts.csv",row.names=FALSE,quote=FALSE)
# write.csv(normalized_counts,"F:/smallRNAwithAdapters/miRNA/normalized_counts_with_repl.csv",row.names=FALSE,quote=FALSE)
# write.table(normalized_counts, file="F:/smallRNAwithAdapters/miRNA/normalized_counts_with_repl.txt", sep="\t", quote=F, row.names=FALSE)
# subset(normalized_counts,normalized_counts$Sequence=="CAGTCGGTAGAGCATCAGAC")
resultsNames(deseqddsColl_boot500)
#comparison_to_change=c("GFP_vs_WT_dpa3","Elac_vs_GFP_dpa3","PARN_vs_WT_dpa3","Elac_vs_WT_dpa3","GFP_vs_WT_dpa5","PARN_vs_GFP_dpa5","Elac_vs_GFP_dpa5","PARN_vs_WT_dpa5","Elac_vs_WT_dpa5")
GFP_vs_WT_dpa3 <- results(deseqddsColl_boot500, contrast=c("condition","GFP_dpa3","WT_dpa3"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
GFP_vs_WT_dpa3 <- GFP_vs_WT_dpa3[ !is.na(GFP_vs_WT_dpa3$padj), ]
GFP_vs_WT_dpa3 <- GFP_vs_WT_dpa3[ !is.na(GFP_vs_WT_dpa3$pvalue), ]
# PARN_vs_GFP_dpa3 <- results(deseqddsColl_boot500, contrast=c("condition", "PARN_dpa3","GFP_dpa3"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
# PARN_vs_GFP_dpa3 <- PARN_vs_GFP_dpa3[ !is.na(PARN_vs_GFP_dpa3$padj), ]
# PARN_vs_GFP_dpa3 <- PARN_vs_GFP_dpa3[ !is.na(PARN_vs_GFP_dpa3$pvalue), ]
Elac_vs_GFP_dpa3 <- results(deseqddsColl_boot500, contrast=c("condition", "ELAC_dpa3","GFP_dpa3"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
Elac_vs_GFP_dpa3 <- Elac_vs_GFP_dpa3[ !is.na(Elac_vs_GFP_dpa3$padj), ]
Elac_vs_GFP_dpa3 <- Elac_vs_GFP_dpa3[ !is.na(Elac_vs_GFP_dpa3$pvalue), ]
# PARN_vs_WT_dpa3 <- results(deseqddsColl_boot500, contrast=c("condition", "PARN_dpa3","WT_dpa3"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
# PARN_vs_WT_dpa3 <- PARN_vs_WT_dpa3[ !is.na(PARN_vs_WT_dpa3$padj), ]
# PARN_vs_WT_dpa3 <- PARN_vs_WT_dpa3[ !is.na(PARN_vs_WT_dpa3$pvalue), ]
Elac_vs_WT_dpa3 <- results(deseqddsColl_boot500, contrast=c("condition","ELAC_dpa3","WT_dpa3"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
Elac_vs_WT_dpa3 <- Elac_vs_WT_dpa3[ !is.na(Elac_vs_WT_dpa3$padj), ]
Elac_vs_WT_dpa3 <- Elac_vs_WT_dpa3[ !is.na(Elac_vs_WT_dpa3$pvalue), ]

GFP_vs_WT_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","GFP_dpa5","WT_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
GFP_vs_WT_dpa5 <- GFP_vs_WT_dpa5[ !is.na(GFP_vs_WT_dpa5$padj), ]
GFP_vs_WT_dpa5 <- GFP_vs_WT_dpa5[ !is.na(GFP_vs_WT_dpa5$pvalue), ]
# PARN_vs_GFP_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition", "PARN_dpa5","GFP_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
# PARN_vs_GFP_dpa5 <- PARN_vs_GFP_dpa5[ !is.na(PARN_vs_GFP_dpa5$padj), ]
# PARN_vs_GFP_dpa5 <- PARN_vs_GFP_dpa5[ !is.na(PARN_vs_GFP_dpa5$pvalue), ]
Elac_vs_GFP_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition", "ELAC_dpa5","GFP_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
Elac_vs_GFP_dpa5 <- Elac_vs_GFP_dpa5[ !is.na(Elac_vs_GFP_dpa5$padj), ]
Elac_vs_GFP_dpa5 <- Elac_vs_GFP_dpa5[ !is.na(Elac_vs_GFP_dpa5$pvalue), ]
# PARN_vs_WT_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition", "PARN_dpa5","WT_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
# PARN_vs_WT_dpa5 <- PARN_vs_WT_dpa5[ !is.na(PARN_vs_WT_dpa5$padj), ]
# PARN_vs_WT_dpa5 <- PARN_vs_WT_dpa5[ !is.na(PARN_vs_WT_dpa5$pvalue), ]
Elac_vs_WT_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","ELAC_dpa5","WT_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
Elac_vs_WT_dpa5 <- Elac_vs_WT_dpa5[ !is.na(Elac_vs_WT_dpa5$padj), ]
Elac_vs_WT_dpa5 <- Elac_vs_WT_dpa5[ !is.na(Elac_vs_WT_dpa5$pvalue), ]


WT_dpa3_vs_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","WT_dpa3", "WT_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
WT_dpa3_vs_dpa5 <- WT_dpa3_vs_dpa5[ !is.na(WT_dpa3_vs_dpa5$padj), ]
WT_dpa3_vs_dpa5 <- WT_dpa3_vs_dpa5[ !is.na(WT_dpa3_vs_dpa5$pvalue), ]
Elac_dpa3_vs_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","ELAC_dpa3", "ELAC_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
Elac_dpa3_vs_dpa5 <- Elac_dpa3_vs_dpa5[ !is.na(Elac_dpa3_vs_dpa5$padj), ]
Elac_dpa3_vs_dpa5 <- Elac_dpa3_vs_dpa5[ !is.na(Elac_dpa3_vs_dpa5$pvalue), ]
# PARN_dpa3_vs_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","PARN_dpa3", "PARN_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
# PARN_dpa3_vs_dpa5 <- PARN_dpa3_vs_dpa5[ !is.na(PARN_dpa3_vs_dpa5$padj), ]
# PARN_dpa3_vs_dpa5 <- PARN_dpa3_vs_dpa5[ !is.na(PARN_dpa3_vs_dpa5$pvalue), ]
GFP_dpa3_vs_dpa5 <- results(deseqddsColl_boot500, contrast=c("condition","GFP_dpa3", "GFP_dpa5"), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=FALSE)
GFP_dpa3_vs_dpa5 <- GFP_dpa3_vs_dpa5[ !is.na(GFP_dpa3_vs_dpa5$padj), ]
GFP_dpa3_vs_dpa5 <- GFP_dpa3_vs_dpa5[ !is.na(GFP_dpa3_vs_dpa5$pvalue), ]


unique(GFP_dpa3_vs_dpa5$padj)

log2cutoff <- 2
qvaluecutoff <- 0.05
GFP_vs_WT_dpa3_filtered=as.data.frame(subset(GFP_vs_WT_dpa3, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#PARN_vs_GFP_dpa3_filtered=as.data.frame(subset(PARN_vs_GFP_dpa3, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
Elac_vs_GFP_dpa3_filtered=as.data.frame(subset(Elac_vs_GFP_dpa3, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#PARN_vs_WT_dpa3_filtered=as.data.frame(subset(PARN_vs_WT_dpa3, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
Elac_vs_WT_dpa3_filtered=as.data.frame(subset(Elac_vs_WT_dpa3, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
GFP_vs_WT_dpa5_filtered=as.data.frame(subset(GFP_vs_WT_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#PARN_vs_GFP_dpa5_filtered=as.data.frame(subset(PARN_vs_GFP_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
Elac_vs_GFP_dpa5_filtered=as.data.frame(subset(Elac_vs_GFP_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#PARN_vs_WT_dpa5_filtered=as.data.frame(subset(PARN_vs_WT_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
Elac_vs_WT_dpa5_filtered=as.data.frame(subset(Elac_vs_WT_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))


WT_dpa3_vs_dpa5_filtered=as.data.frame(subset(WT_dpa3_vs_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
GFP_dpa3_vs_dpa5_filtered=as.data.frame(subset(GFP_dpa3_vs_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
Elac_dpa3_vs_dpa5_filtered=as.data.frame(subset(Elac_dpa3_vs_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))
#PARN_dpa3_vs_dpa5_filtered=as.data.frame(subset(PARN_dpa3_vs_dpa5, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff))

GFP_vs_WT_dpa3_filtered$set="GFP_vs_WT_dpa3"
#PARN_vs_GFP_dpa3_filtered$set="PARN_vs_GFP_dpa3"
Elac_vs_GFP_dpa3_filtered$set="Elac_vs_GFP_dpa3"
#PARN_vs_WT_dpa3_filtered$set="PARN_vs_WT_dpa3"
Elac_vs_WT_dpa3_filtered$set="Elac_vs_WT_dpa3"
GFP_vs_WT_dpa5_filtered$set="GFP_vs_WT_dpa5"
#PARN_vs_GFP_dpa5_filtered$set="PARN_vs_GFP_dpa5"
Elac_vs_GFP_dpa5_filtered$set="Elac_vs_GFP_dpa5"
#PARN_vs_WT_dpa5_filtered$set="PARN_vs_WT_dpa5"
Elac_vs_WT_dpa5_filtered$set="Elac_vs_WT_dpa5"
WT_dpa3_vs_dpa5_filtered$set="WT_dpa3_vs_dpa5"
GFP_dpa3_vs_dpa5_filtered$set="GFP_dpa3_vs_dpa5" #no RNAs
Elac_dpa3_vs_dpa5_filtered$set="Elac_dpa3_vs_dpa5" #no RNAs
#PARN_dpa3_vs_dpa5_filtered$set="PARN_dpa3_vs_dpa5" #no RNAs

GFP_vs_WT_dpa3_filtered$RNA=rownames(GFP_vs_WT_dpa3_filtered)
#PARN_vs_GFP_dpa3_filtered$RNA=rownames(PARN_vs_GFP_dpa3_filtered)
Elac_vs_GFP_dpa3_filtered$RNA=rownames(Elac_vs_GFP_dpa3_filtered)
#PARN_vs_WT_dpa3_filtered$RNA=rownames(PARN_vs_WT_dpa3_filtered)
Elac_vs_WT_dpa3_filtered$RNA=rownames(Elac_vs_WT_dpa3_filtered)
GFP_vs_WT_dpa5_filtered$RNA=rownames(GFP_vs_WT_dpa5_filtered)
#PARN_vs_GFP_dpa5_filtered$RNA=rownames(PARN_vs_GFP_dpa5_filtered)
Elac_vs_GFP_dpa5_filtered$RNA=rownames(Elac_vs_GFP_dpa5_filtered)
#PARN_vs_WT_dpa5_filtered$RNA=rownames(PARN_vs_WT_dpa5_filtered)
Elac_vs_WT_dpa5_filtered$RNA=rownames(Elac_vs_WT_dpa5_filtered)
WT_dpa3_vs_dpa5_filtered$RNA=rownames(WT_dpa3_vs_dpa5_filtered)
GFP_dpa3_vs_dpa5_filtered$RNA=rownames(GFP_dpa3_vs_dpa5_filtered)
Elac_dpa3_vs_dpa5_filtered$RNA=rownames(Elac_dpa3_vs_dpa5_filtered)
#PARN_dpa3_vs_dpa5_filtered$RNA=rownames(PARN_dpa3_vs_dpa5_filtered)


# 
# nrow(t1)+nrow(t2)+nrow(t3)+nrow(t4)+nrow(t5)+nrow(t6)+nrow(t7)+nrow(t8)+nrow(t9)+nrow(t10)+nrow(t11)+nrow(t12)+nrow(t13)+nrow(t14)+nrow(t15)+nrow(t16)
# t4[which(rownames(t4)=="dd_Smed_v6_8992_0_11"),]
# DGE_set=rbind(GFP_vs_WT_dpa3_filtered,PARN_vs_GFP_dpa3_filtered,Elac_vs_GFP_dpa3_filtered,PARN_vs_WT_dpa3_filtered,Elac_vs_WT_dpa3_filtered,GFP_vs_WT_dpa5_filtered,PARN_vs_GFP_dpa5_filtered,Elac_vs_GFP_dpa5_filtered,
#               PARN_vs_WT_dpa5_filtered, Elac_vs_WT_dpa5_filtered, WT_dpa3_vs_dpa5_filtered,GFP_dpa3_vs_dpa5_filtered,Elac_dpa3_vs_dpa5_filtered,PARN_dpa3_vs_dpa5_filtered)
DGE_set=rbind(GFP_vs_WT_dpa3_filtered,Elac_vs_GFP_dpa3_filtered,Elac_vs_WT_dpa3_filtered,GFP_vs_WT_dpa5_filtered,
              Elac_vs_GFP_dpa5_filtered,Elac_vs_WT_dpa5_filtered,WT_dpa3_vs_dpa5_filtered,GFP_dpa3_vs_dpa5_filtered,Elac_dpa3_vs_dpa5_filtered)

nrow(DGE_set) #1541
#length()

colnames(DGE_set)



DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("ribosomal", ignore_case = TRUE)),"rRNA_fragments","other")
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("rRNA", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("tRNA", ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("lncRNA", ignore_case = TRUE)),"lncRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("piRNA", ignore_case = TRUE)),"piRNA",DGE_set$RNA_type)
#DGE_set$RNA_type=ifelse(str_detect(DGE_set$ID, regex("RNase_P_RNA", ignore_case = TRUE)),"RNase_P_RNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("snoRNA", ignore_case = TRUE)),"snoRNA_fragments",DGE_set$RNA_type)
#DGE_set$RNA_type=ifelse(str_detect(DGE_set$ID, regex("siRNA", ignore_case = TRUE)),"siRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("miRNA", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
#DGE_set$RNA_type=ifelse(str_detect(DGE_set$ID, regex("pre_miRNA", ignore_case = TRUE)),"pre_miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("snRNA", ignore_case = TRUE)),"snRNA_fragments",DGE_set$RNA_type)
#DGE_set$RNA_type=ifelse(str_detect(DGE_set$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",DGE_set$RNA_type)
#DGE_set$RNA_type=ifelse(str_detect(DGE_set$ID, regex("guide_RNA", ignore_case = TRUE)),"guide_RNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("Sme-", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)  
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("spliceosomal-", ignore_case = TRUE)),"snRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("lnc", ignore_case = TRUE)),"lncRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("microRNA", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("long_non-coding", ignore_case = TRUE)),"lncRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("7SK", ignore_case = TRUE)),"snRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("nucleolar", ignore_case = TRUE)),"snoRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("nuclear", ignore_case = TRUE)),"snRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("miR", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("transfer", ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("synthetic", ignore_case = TRUE)),"synthetic",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("Cloning", ignore_case = TRUE)),"Cloning_vector",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("dd_Smed_v6", ignore_case = TRUE)),"mRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("sme-lin", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("sme-let", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("Sme-Bantam", ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("multiple_hits", ignore_case = TRUE)),"multiple_hits",DGE_set$RNA_type)

DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("ITS1", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("ITS2", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("SpacerA", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("28S", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("Schmed_cloneH735c", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("12S", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("16S", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("5.8S", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_miRNA_weird, ignore_case = TRUE)),"miRNA",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[1], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[2], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[3], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[4], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[5], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[6], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[7], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[8], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[9], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[10], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[11], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[12], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex(SM_tRNA_weird[13], ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("5S", ignore_case = TRUE)),"rRNA_fragments",DGE_set$RNA_type)
DGE_set$RNA_type=ifelse(str_detect(DGE_set$RNA, regex("UNDET", ignore_case = TRUE)),"tRNA_fragments",DGE_set$RNA_type)


as.data.frame(table(DGE_set$RNA_type))

DGE_set[DGE_set$RNA_type=="other",]
RNA_type_to_remove=c("Cloning_vector","synthetic","multiple_hits")
DGE_set_weird=subset(DGE_set,DGE_set$RNA_type %in% RNA_type_to_remove)
pattern <- ".*?_(.*?)_.*"
DGE_set_weird$URS=sub(pattern, "\\1", DGE_set_weird$RNA)
unique(DGE_set_weird$URS)

ELAC_mRNA_fr=which(str_detect(DGE_set$RNA, regex("dd_Smed_v6_8992_0_1", ignore_case = TRUE)))
DGE_set=DGE_set[-ELAC_mRNA_fr,]
DGE_set
#before pattern
sub(" xxx.*", "", x) 
#after pattern
sub(".*xxx ", "", x) 
DGE_set$Sequence=sub(" .*", "", DGE_set$RNA) 
pos_anno=sub(".* ", "", DGE_set$RNA)
DGE_set$Position=sub("_.*", "", pos_anno)
DGE_set$Seq_length=nchar(DGE_set$Sequence)
DGE_set$Annotation=sub("^[^_]+_", "", DGE_set$RNA)
DGE_set$Annotation=ifelse(DGE_set$Annotation=="hits","multiple_hits",DGE_set$Annotation)
colnames(DGE_set)


DGE_set$expression=ifelse(DGE_set$log2FoldChange>0,"accumulation","loss")
gene_norm_count=as.data.frame(counts(deseqddsColl_boot500, normalized=TRUE))
colnames(gene_norm_count)
gene_norm_count=gene_norm_count[ , order(names(gene_norm_count))]
gene_norm_count$RNA=rownames(gene_norm_count)
DGE_set_with_norm_counts=merge(DGE_set,gene_norm_count,by="RNA",all.x=TRUE)



#check for a consistency within replicates
#coefficient of dispersion
unique(DGE_set_with_norm_counts$set)
unique_sets=c("GFP_vs_WT_dpa3","Elac_vs_WT_dpa5","Elac_vs_WT_dpa3","Elac_vs_GFP_dpa5","GFP_vs_WT_dpa5","Elac_vs_GFP_dpa3")
#unique_sets=unique(DGE_set_with_norm_counts$set)
# "Elac_dpa3_vs_dpa5" "GFP_dpa3_vs_dpa5"  "GFP_vs_WT_dpa5"    "Elac_vs_GFP_dpa5"  "WT_dpa3_vs_dpa5"   "Elac_vs_GFP_dpa3"  "Elac_vs_WT_dpa3"   "Elac_vs_WT_dpa5"   "GFP_vs_WT_dpa3" 


#which columns correspond to each sets
#substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==3
ELAC3=which((startsWith(colnames(DGE_set_with_norm_counts),"ELAC")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==3))
ELAC5=which((startsWith(colnames(DGE_set_with_norm_counts),"ELAC")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==5))
WT3=which((startsWith(colnames(DGE_set_with_norm_counts),"WT")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==3))
WT5=which((startsWith(colnames(DGE_set_with_norm_counts),"WT")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==5))
GFP3=which((startsWith(colnames(DGE_set_with_norm_counts),"GFP")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==3))
GFP5=which((startsWith(colnames(DGE_set_with_norm_counts),"GFP")) & (substr(colnames(DGE_set_with_norm_counts),nchar(colnames(DGE_set_with_norm_counts))-1,nchar(colnames(DGE_set_with_norm_counts))-1)==5))
#corresponding to unique sets columns
unique_sets_list=list(list(GFP3,WT3),#"GFP_vs_WT_dpa3"
                      list(ELAC5,WT5), #"Elac_vs_WT_dpa5"
                      list(ELAC3,WT3), #"Elac_vs_WT_dpa3"
                      list(ELAC5,GFP5), #"Elac_vs_GFP_dpa5"
                      list(GFP5,WT5), #"GFP_vs_WT_dpa5"
                      list(ELAC3,GFP3)) #"Elac_vs_GFP_dpa3"





#cv_df=as.data.frame(counts(dds, normalized=TRUE))
as.numeric(DGE_set_with_norm_counts[1,ELAC3])
(mad(as.numeric(DGE_set_with_norm_counts[1,ELAC3]))/median(as.numeric(DGE_set_with_norm_counts[1,ELAC3]))) * 100
for(i in 1:nrow(DGE_set_with_norm_counts)){
  DGE_set_with_norm_counts$ELAC3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,ELAC3]))/median(as.numeric(DGE_set_with_norm_counts[i,ELAC3]))) * 100
  DGE_set_with_norm_counts$ELAC5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,ELAC5]))/median(as.numeric(DGE_set_with_norm_counts[i,ELAC5]))) * 100
  
  DGE_set_with_norm_counts$GFP3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,GFP3]))/median(as.numeric(DGE_set_with_norm_counts[i,GFP3]))) * 100
  DGE_set_with_norm_counts$GFP5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,GFP5]))/median(as.numeric(DGE_set_with_norm_counts[i,GFP5]))) * 100
  
  DGE_set_with_norm_counts$WT3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,WT3]))/median(as.numeric(DGE_set_with_norm_counts[i,WT3]))) * 100
  DGE_set_with_norm_counts$WT5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,WT5]))/median(as.numeric(DGE_set_with_norm_counts[i,WT5]))) * 100
  
}

#(median(as.numeric(DGE_set_with_norm_counts[1,ELAC3])))+(4*(mad(as.numeric(DGE_set_with_norm_counts[1,ELAC3]))))
DGE_set_with_norm_counts=replace(DGE_set_with_norm_counts, is.na(DGE_set_with_norm_counts), 0)
#The MAD Ratio can be interpreted as a measure of relative variability. 
#A small MAD Ratio indicates that the data is tightly clustered around the median, indicating a high degree of consistency within the dataset. 
#A large MAD Ratio indicates that the data is widely dispersed, indicating a lower degree of consistency within the dataset.
#relative median absolute deviation
cv_df=DGE_set_with_norm_counts[,c("ELAC3_cv","ELAC5_cv","GFP3_cv","GFP5_cv","WT3_cv","WT5_cv")]
summary(cv_df$ELAC3_cv)
cv_df_long <- gather(cv_df, condition, cv, ELAC3_cv:WT5_cv, factor_key=TRUE)
p<-ggplot(cv_df_long, aes(x=condition, y=cv, fill=condition)) +geom_boxplot()
p+scale_fill_brewer(palette="Dark2")

for(i in 1:nrow(DGE_set_with_norm_counts)){
  DGE_set_with_norm_counts$ELAC3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,ELAC3]))/median(as.numeric(DGE_set_with_norm_counts[i,ELAC3]))) * 100
  DGE_set_with_norm_counts$ELAC5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,ELAC5]))/median(as.numeric(DGE_set_with_norm_counts[i,ELAC5]))) * 100
  
  DGE_set_with_norm_counts$GFP3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,GFP3]))/median(as.numeric(DGE_set_with_norm_counts[i,GFP3]))) * 100
  DGE_set_with_norm_counts$GFP5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,GFP5]))/median(as.numeric(DGE_set_with_norm_counts[i,GFP5]))) * 100
  
  DGE_set_with_norm_counts$WT3_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,WT3]))/median(as.numeric(DGE_set_with_norm_counts[i,WT3]))) * 100
  DGE_set_with_norm_counts$WT5_cv[i]=(mad(as.numeric(DGE_set_with_norm_counts[i,WT5]))/median(as.numeric(DGE_set_with_norm_counts[i,WT5]))) * 100
  
}





for(i in 1:nrow(DGE_set_with_norm_counts)) {
  for(j in 1: length(unique_sets)){
    if(DGE_set_with_norm_counts$set[i]==unique_sets[j]){
      
      
      DGE_set_with_norm_counts$consistency_strict[i]=ifelse(sum(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))])>median(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))]))-4*mad(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))])))==3&
                                                              sum(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))])<median(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))]))+4*mad(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][1]))])))==3&
                                                              sum(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))])>median(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))]))-4*mad(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))])))==3&
                                                              sum(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))])<median(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))]))+4*mad(as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list[[j]][2]))])))==3,
                                                            "reliable","not_reliable")
    }
  }
}






unique_sets_list
which(colnames(DGE_set_with_norm_counts)=="GFP3_cv")

unique_sets_list_cv=list(list(which(colnames(DGE_set_with_norm_counts)=="GFP3_cv"),which(colnames(DGE_set_with_norm_counts)=="WT3_cv")),#"GFP_vs_WT_dpa3"
                         list(which(colnames(DGE_set_with_norm_counts)=="ELAC5_cv"),which(colnames(DGE_set_with_norm_counts)=="WT5_cv")), #"Elac_vs_WT_dpa5"
                         list(which(colnames(DGE_set_with_norm_counts)=="ELAC3_cv"),which(colnames(DGE_set_with_norm_counts)=="WT3_cv")), #"Elac_vs_WT_dpa3"
                         list(which(colnames(DGE_set_with_norm_counts)=="ELAC5_cv"),which(colnames(DGE_set_with_norm_counts)=="GFP5_cv")), #"Elac_vs_GFP_dpa5"
                         list(which(colnames(DGE_set_with_norm_counts)=="GFP5_cv"),which(colnames(DGE_set_with_norm_counts)=="WT5_cv")), #"GFP_vs_WT_dpa5"
                         list(which(colnames(DGE_set_with_norm_counts)=="ELAC3_cv"),which(colnames(DGE_set_with_norm_counts)=="GFP3_cv"))) #"Elac_vs_GFP_dpa3"

for(i in 1:nrow(DGE_set_with_norm_counts)) {
  for(j in 1: length(unique_sets)){
    if(DGE_set_with_norm_counts$set[i]==unique_sets[j]){
      DGE_set_with_norm_counts$consistency_mad_rmad[i]=ifelse((as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list_cv[[j]][1]))]) <= 2 * mad(DGE_set_with_norm_counts[,c(unlist(unique_sets_list_cv[[j]][1]))]))&
                                                                (as.numeric(DGE_set_with_norm_counts[i,c(unlist(unique_sets_list_cv[[j]][2]))]) <= 2 * mad(DGE_set_with_norm_counts[,c(unlist(unique_sets_list_cv[[j]][2]))])),
                                                              "reliable","not_reliable")
    }
  }
}


table(DGE_set_with_norm_counts$consistency_strict)
table(DGE_set_with_norm_counts$consistency_mad_rmad)

consistency_strict=subset(DGE_set_with_norm_counts,DGE_set_with_norm_counts$consistency_strict=="reliable")
consistency_mad_rmad=subset(DGE_set_with_norm_counts,DGE_set_with_norm_counts$consistency_mad_rmad=="reliable")

table(consistency_strict$RNA_type)
table(consistency_mad_rmad$RNA_type)
save(DGE_set_with_norm_counts,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_norm_counts_corrected_final.RData")
# (as.numeric(DGE_set_with_norm_counts[1,c(unlist(unique_sets_list_cv[[1]][1]))]) <= 2 * mad(DGE_set_with_norm_counts[,c(unlist(unique_sets_list_cv[[1]][1]))]))&
#   (as.numeric(DGE_set_with_norm_counts[1,c(unlist(unique_sets_list_cv[[1]][2]))]) <= 2 * mad(DGE_set_with_norm_counts[,c(unlist(unique_sets_list_cv[[1]][2]))]))

load("G:/Elac2/smallRNA/DGE_set_with_norm_counts_corrected_final.RData")
#save(DGE_set_with_norm_counts,file="F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_new_ELAC.RData")
library("xlsx")
write.xlsx(DGE_set_with_norm_counts, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_norm_counts_ELAC_new.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
# 
# 
# DGE_set$RNA_type_without_underscore=DGE_set$RNA_type
# DGE_set$RNA_type_without_underscore=ifelse(DGE_set$RNA_type_without_underscore=="snoRNA_fragments","snoRNA fragments",DGE_set$RNA_type_without_underscore)
# DGE_set$RNA_type_without_underscore=ifelse(DGE_set$RNA_type_without_underscore=="tRNA_fragments","tRNA fragments",DGE_set$RNA_type_without_underscore)
# DGE_set$RNA_type_without_underscore=ifelse(DGE_set$RNA_type_without_underscore=="rRNA_fragments","rRNA fragments",DGE_set$RNA_type_without_underscore)
################################################################################
#change names
# unique(DGE_set$set)
# colnames(DGE_set)
# comparison_to_be_changed=c("WT_vs_GFP_dpa3","GFP_vs_Elac_dpa3","WT_vs_PARN_dpa3","WT_vs_Elac_dpa3","WT_vs_GFP_dpa5","GFP_vs_PARN_dpa5","GFP_vs_Elac_dpa5","WT_vs_PARN_dpa5","WT_vs_Elac_dpa5")
# comparison_to_change=c("GFP_vs_WT_dpa3","Elac_vs_GFP_dpa3","PARN_vs_WT_dpa3","Elac_vs_WT_dpa3","GFP_vs_WT_dpa5","PARN_vs_GFP_dpa5","Elac_vs_GFP_dpa5","PARN_vs_WT_dpa5","Elac_vs_WT_dpa5")
# 
# DGE_set$new_comparison=DGE_set$set
# #DGE_set$new_exp=all_GO_elim$exp
# for (i in 1:length(comparison_to_be_changed)) {
#   DGE_set$new_comparison=ifelse(all_GO_elim$new_comparison==comparison_to_be_changed[i],comparison_to_change[i],all_GO_elim$new_comparison)
#   # all_GO_elim$new_exp=ifelse(all_GO_elim$new_comparison==comparison_to_change[i] & all_GO_elim$new_exp=="under","over",all_GO_elim$new_exp)
#   # all_GO_elim$new_exp=ifelse(all_GO_elim$new_comparison==comparison_to_change[i] & all_GO_elim$new_exp=="over","under",all_GO_elim$new_exp)
# }
################################################################################

#DGE_set$set
#DGE_set$Sequence
# DGE_set_up=subset(DGE_set,DGE_set$log2FoldChange>0)
# DGE_set_down=subset(DGE_set,DGE_set$log2FoldChange<0)
# DGE_set_up_table=as.data.frame(table(DGE_set_up$set,DGE_set_up$RNA_type))
# colnames(DGE_set_up_table)=c("Comparison","RNA_type","Frequency")
# DGE_set_up_table=subset(DGE_set_up_table,DGE_set_up_table$Frequency>0)
# unique(DGE_set_up_table$Comparison)
# DGE_set_up_table=subset(DGE_set_up_table,DGE_set_up_table$Comparison!="Elac_vs_GFP_dpa3" & DGE_set_up_table$Comparison!="Elac_vs_GFP_dpa5" & 
#                           DGE_set_up_table$Comparison!="PARN_vs_GFP_dpa5")
# # DGE_set_up_table=subset(DGE_set_up_table,DGE_set_up_table$Comparison!="Elac_vs_GFP_dpa3" & DGE_set_up_table$Comparison!="Elac_vs_GFP_dpa5" & 
# #                           DGE_set_up_table$Comparison!="PARN_vs_GFP_dpa5")
# 
# DGE_set_down_table=as.data.frame(table(DGE_set_down$set,DGE_set_down$RNA_type))
# colnames(DGE_set_down_table)=c("Comparison","RNA_type","Frequency")
# DGE_set_down_table=subset(DGE_set_down_table,DGE_set_down_table$Frequency>0)
# unique(DGE_set_down_table$Comparison)
# DGE_set_down_table=subset(DGE_set_down_table,DGE_set_down_table$Comparison!="Elac_vs_GFP_dpa3" & DGE_set_down_table$Comparison!="Elac_vs_GFP_dpa5" & 
#                             DGE_set_down_table$Comparison!="PARN_vs_GFP_dpa5")
# max(DGE_set_down_table$Frequency)
# typeof(DGE_set_up_table$RNA_type)
# colnames(DGE_set)
# DGE_set_table=as.data.frame(table(DGE_set$set,DGE_set$RNA_type_without_underscore,DGE_set$expression))
# colnames(DGE_set_table)=c("Comparison","RNA_type","Expression","Frequency")
# DGE_set_table=subset(DGE_set_table,DGE_set_table$Frequency>0)
# unique(DGE_set_table$Comparison)
# DGE_set_table=subset(DGE_set_table,DGE_set_table$Comparison!="Elac_vs_GFP_dpa3" & DGE_set_table$Comparison!="Elac_vs_GFP_dpa5" &
#                        DGE_set_table$Comparison!="PARN_vs_GFP_dpa5" & DGE_set_table$Comparison!="PARN_vs_GFP_dpa3")
# # DGE_set_table=subset(DGE_set_table,DGE_set_table$Comparison!="Elac_vs_WT_dpa3" & DGE_set_table$Comparison!="Elac_vs_WT_dpa5" & 
# #                        DGE_set_table$Comparison!="PARN_vs_WT_dpa5" & DGE_set_table$Comparison!="GFP_dpa3_vs_dpa5" & DGE_set_table$Comparison!="GFP_vs_WT_dpa3" &
# #                        DGE_set_table$Comparison!="GFP_vs_WT_dpa5")
# # DGE_set_up_table$RNA_type=as.factor(DGE_set_up_table$RNA_typ)
# # DGE_set_down_table$RNA_type=as.factor(DGE_set_down_table$RNA_typ)
# # DGE_set_down_table$RNA_type=ifelse(DGE_set_down_table$RNA_type=="rRNA","rRF",DGE_set_down_table$RNA_type)
# # DGE_set_up_table$RNA_type=ifelse(DGE_set_up_table$RNA_type=="rRNA","rRF",DGE_set_up_table$RNA_type)
# # DGE_set_down_table$RNA_type=ifelse(DGE_set_down_table$RNA_type=="tRNA","tRF",DGE_set_down_table$RNA_type)
# # DGE_set_up_table$RNA_type=ifelse(DGE_set_up_table$RNA_type=="tRNA","tRF",DGE_set_up_table$RNA_type)
# # class(DGE_set_up_table$Comparison)
# # sort(unique(DGE_set_up_table$Frequency))
# #DGE_set_up_table$Frequency=DGE_set_up_table$Frequency*10
# upplot=ggplot(DGE_set_up_table) +
#   geom_point(aes(x = RNA_type, y = Comparison, color = RNA_type, size=Frequency),alpha = .5) +
#   theme(text = element_text(size = 20)) +
#   xlab("RNA type") +
#   ylab("Comparison") +
#   ggtitle("Upregulated sncRNAs")+
#   scale_size_continuous(breaks = seq(0, 100, by = 20),
#                         range = c(1, 10))
# downplot=ggplot(DGE_set_down_table) +
#   geom_point(aes(x = RNA_type, y = Comparison, color = RNA_type, size=Frequency),alpha = .5) +
#   theme(text = element_text(size = 20)) +
#   xlab("RNA type") +
#   ylab("Comparison") +
#   ggtitle("Downregulated sncRNAs")+
#   scale_size_continuous(breaks = seq(0, 100, by = 20),
#                         range = c(1, 10))
# facet_grid(~RNA_type)
# ggarrange(upplot,downplot)
# write.csv(DGE_set, "F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation.csv",
#           row.names = FALSE, quote = FALSE)
# 
# 
# ggplot(DGE_set_table) +theme_bw()+
#   geom_point(aes(x = RNA_type, y = Comparison, color = RNA_type, size=Frequency),alpha = .5)+
#   theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap(~Expression, ncol=2, scales = "fixed")+
#   labs(title = "RNA fragments", x = "RNA type", y = "Comparison", color = "RNA type", size="Frequency")
# ?substr
#not_to_include=c("PARN_dpa3_vs_dpa5","WT_dpa3_vs_dpa5","GFP_dpa3_vs_dpa5","Elac_vs_GFP_dpa3","Elac_vs_GFP_dpa5","PARN_vs_GFP_dpa5")
DGE_set_to_analyse=DGE_set_with_norm_counts[DGE_set_with_norm_counts$consistency_strict=="reliable",]
unique(DGE_set_to_analyse$set)
DGE_set_GFP=subset(DGE_set_to_analyse,startsWith(DGE_set_to_analyse$set,"GFP"))
DGE_set_GFP_loss=subset(DGE_set_GFP,DGE_set_GFP$expression=="loss")
DGE_set_GFP_accumulation=subset(DGE_set_GFP,DGE_set_GFP$expression=="accumulation")
DGE_set_PARN_ELAC=subset(DGE_set_to_analyse,DGE_set_to_analyse$set=="Elac_vs_WT_dpa3"|DGE_set_to_analyse$set=="Elac_vs_WT_dpa5")
DGE_set_PARN_ELAC_loss=subset(DGE_set_PARN_ELAC,DGE_set_PARN_ELAC$expression=="loss")
DGE_set_PARN_ELAC_accumulation=subset(DGE_set_PARN_ELAC,DGE_set_PARN_ELAC$expression=="accumulation")
GFP_ELAC=subset(DGE_set_PARN_ELAC,DGE_set_PARN_ELAC$set=="Elac_vs_GFP_dpa5"|DGE_set_PARN_ELAC$set=="Elac_vs_GFP_dpa3")
GFP_ELAC_loss=subset(GFP_ELAC,GFP_ELAC$expression=="loss")
GFP_ELAC_accumulation=subset(GFP_ELAC,GFP_ELAC$expression=="accumulation")

unique_for_GFP_loss=subset(DGE_set_GFP_loss,!(DGE_set_GFP_loss$Sequence%in%GFP_ELAC_loss$Sequence))
unique_for_GFP_accumulation=subset(DGE_set_GFP_accumulation,!(DGE_set_GFP_accumulation$Sequence%in%GFP_ELAC_accumulation$Sequence))

DGE_set_PARN_ELAC_loss=subset(DGE_set_PARN_ELAC_loss,!(DGE_set_PARN_ELAC_loss$Sequence%in%unique_for_GFP_loss$Sequence))
DGE_set_PARN_ELAC_accumulation=subset(DGE_set_PARN_ELAC_accumulation,!(DGE_set_PARN_ELAC_accumulation$Sequence%in%unique_for_GFP_accumulation$Sequence))
DGE_set_all=rbind(DGE_set_PARN_ELAC_loss,DGE_set_PARN_ELAC_accumulation)
DGE_set_all=subset(DGE_set_all,DGE_set_all$RNA_type!="multiple_hits")
unique(DGE_set_all$RNA_type)
table(DGE_set_all$RNA_type)

DGE_set_table=as.data.frame(table(DGE_set_all$set,DGE_set_all$RNA_type,DGE_set_all$expression))
colnames(DGE_set_table)=c("Comparison","RNA_type","Expression","Frequency")
DGE_set_table=subset(DGE_set_table,DGE_set_table$Frequency>0)
unique(DGE_set_table$Comparison)
# DGE_set_table=subset(DGE_set_table,DGE_set_table$Comparison!="Elac_vs_GFP_dpa3" & DGE_set_table$Comparison!="Elac_vs_GFP_dpa5" &
#                        DGE_set_table$Comparison!="PARN_vs_GFP_dpa5" & DGE_set_table$Comparison!="PARN_vs_GFP_dpa3")




(DGE_set_table$Set=substr(DGE_set_table$Comparison,1,nchar(as.character(DGE_set_table$Comparison))-5))
DGE_set_table$Set=ifelse(DGE_set_table$Set=="Elac_vs_WT","ELAC2","PARN")
(DGE_set_table$dpa=substr(DGE_set_table$Comparison,nchar(as.character(DGE_set_table$Comparison))-3,nchar(as.character(DGE_set_table$Comparison))))
(p=ggplot(DGE_set_table) +
    geom_point(aes(x = RNA_type, y = Set, color = RNA_type, size=Frequency),alpha = .5)+
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(dpa ~Expression,scale="free",space="free")+
    labs(title = "RNA fragments accumulation", x = "RNA type", y = "Comparison", color = "RNA type", size="Number of RNA species")+
    theme(ggh4x.facet.nestline = element_line(linetype = 3),plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size = 25))+
    geom_rect(aes(fill = dpa),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.01)) + guides(dpa = "none")

DGE_set_table_all=as.data.frame(table(DGE_set_to_analyse$set,DGE_set_to_analyse$RNA_type,DGE_set_to_analyse$expression))
colnames(DGE_set_table_all)=c("Comparison","RNA_type","Expression","Frequency")
DGE_set_table_all=subset(DGE_set_table_all,DGE_set_table_all$Frequency>0)
(DGE_set_table_all$Set=substr(DGE_set_table_all$Comparison,1,nchar(as.character(DGE_set_table_all$Comparison))-5))
DGE_set_table_all$Set=ifelse(DGE_set_table_all$Set=="Elac_vs_WT","ELAC2","GFP")
(DGE_set_table_all$dpa=substr(DGE_set_table_all$Comparison,nchar(as.character(DGE_set_table_all$Comparison))-3,nchar(as.character(DGE_set_table_all$Comparison))))
(p=ggplot(DGE_set_table_all) +
    geom_point(aes(x = RNA_type, y = Set, color = RNA_type, size=Frequency),alpha = .5)+
    theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(dpa ~Expression,scale="free",space="free")+
    labs(title = "RNA fragments accumulation", x = "RNA type", y = "Comparison", color = "RNA type", size="Number of RNA species")+
    theme(ggh4x.facet.nestline = element_line(linetype = 3),plot.title = element_text(hjust = 0.5,face = "bold"),text = element_text(size = 25))+
    geom_rect(aes(fill = dpa),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.01)) + guides(dpa = "none")

################################################################################
#check length of original sequences


# save(DGE_set_with_anno, file = "F:/silencing_with_batch/new_DGE_set_with_annotation.RData")
#DGE_set=read.csv("E:/project/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation_new.csv",header=TRUE,sep=";")
#cpm.keep=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/cpm_filtered_with_tags_end_to_end.csv")
# cpm.keep[1,]
# colnames(cpm.keep)
# rownames(cpm.keep)=cpm.keep$X
# nrow(DGE_set_all)
# cpm.keep_test=as.data.frame(cpm.keep)
# cpm.keep=as.data.frame(cpm.keep)
# cpm.keep[c('Sequence', 'Annotation', 'Position', 'Mismatch', 'Cigar', 'Match')] <- str_split_fixed(rownames(cpm.keep), ' ', 6)
# cpm.keep_test
# colnames(cpm.keep_test)
# #get cpm for samples
# colnames(DGE_set_all)
# DGE_set_all$first_condition_cpm_sum="first_condition"
# DGE_set_all$second_condition_cpm_sum="second_condition"
# DGE_set_all$Sequence
# unique(DGE_set_all$set)
# 
# 
# for (i in 1:nrow(DGE_set_all)){
#   if(DGE_set_all$set[i]=="GFP_vs_WT_dpa3"){
#     DGE_set_all$first_condition_cpm_sum[[i]]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"3S"))])
#   }
#   if(DGE_set_all$set[i]=="GFP_vs_WT_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"5S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="PARN_vs_WT_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"PARN")&endsWith(colnames(cpm.keep),"5S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="Elac_vs_WT_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"ELAC")&endsWith(colnames(cpm.keep),"5S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="PARN_vs_WT_dpa3"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"PARN")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"3S"))])
#   }
#   if(DGE_set_all$set[i]=="Elac_vs_WT_dpa3"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"ELAC")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"3S"))])
#   }
#   if(DGE_set_all$set[i]=="PARN_vs_GFP_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"PARN")&endsWith(colnames(cpm.keep),"5S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="Elac_vs_GFP_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"ELAC")&endsWith(colnames(cpm.keep),"5S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="Elac_vs_GFP_dpa3"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"ELAC")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"3S"))])
#   }
#   if(DGE_set_all$set[i]=="PARN_dpa3_vs_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"PARN")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"PARN")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="GFP_dpa3_vs_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"GFP")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   if(DGE_set_all$set[i]=="WT_dpa3_vs_dpa5"){
#     DGE_set_all$first_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"3S"))])
#     DGE_set_all$second_condition_cpm_sum[i]=sum(cpm.keep[cpm.keep$Sequence==DGE_set_all$Sequence[i],which(startsWith(colnames(cpm.keep),"WT")&endsWith(colnames(cpm.keep),"5S"))])
#   }
#   
# }

# write.csv(DGE_set, "F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation_with_cpm.csv",
#           row.names = FALSE, quote = FALSE)
# write.table(DGE_set, "F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation_with_cpm_new.csv",sep=";",
#           row.names = FALSE, quote = FALSE)
# write.csv(DGE_set[,c("Annotation","first_condition_cpm_sum","second_condition_cpm_sum")], "F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation_with_cpm_just_cpm.csv",
#           row.names = FALSE, quote = FALSE)
tail(DGE_set_all)
DGE_set_all$Annotation
colnames(DGE_set_all)
DGE_set_all=subset(DGE_set_all,DGE_set_all$RNA_type!="mRNA_fragments")
DGE_set_all_SM_tRNA=subset(DGE_set_all,startsWith(DGE_set_all$Annotation,"dd_Smed"))
DGE_set_all_SM_rRNA=subset(DGE_set_all,(!startsWith(DGE_set_all$Annotation,"URS"))&DGE_set_all$RNA_type=="rRNA_fragments")
DGE_set_all=subset(DGE_set_all,startsWith(DGE_set_all$Annotation,"URS"))
#get ID and extract len and seq
DGE_set_all$rnacentral_ID=str_before_nth(DGE_set_all$Annotation, "_", 2)
DGE_set_all[which(is.na(DGE_set_all$rnacentral_ID)),]
DGE_set_all$rnacentral_ID=ifelse(is.na(DGE_set_all$rnacentral_ID),DGE_set_all$Annotation,DGE_set_all$rnacentral_ID)
#DGE_set_all=subset(DGE_set_all,DGE_set_all$rnacentral_ID!="dd_Smed")
for (i in 1:nrow(DGE_set_all)){
  DGE_set_all$rnacentral_Seq[i]=rnaCentralRetrieveEntry(DGE_set_all$rnacentral_ID[i])$sequence
  DGE_set_all$rnacentral_len[i]=rnaCentralRetrieveEntry(DGE_set_all$rnacentral_ID[i])$sequenceLength
}


#read SM tRNA
tRNA_fasta=read.fasta("F:/PARN_ELAC_silencing/smallRNA/tRNA/SM_tRNA.fasta")
tRNA_SM_df=data.frame(Annotation=getName(tRNA_fasta),
                      rnacentral_len=getLength(tRNA_fasta),
                      rnacentral_Seq=unlist(getSequence(tRNA_fasta,as.string = TRUE)))
DGE_set_all_SM_tRNA=merge(DGE_set_all_SM_tRNA,tRNA_SM_df,by="Annotation",all.x=TRUE)
DGE_set_all_SM_tRNA$rnacentral_ID="no_ID"
#read SM rRNA
rRNA_fasta=read.fasta("F:/genome/rRNA/FINAL_rRNA_Schmed.fa")
rRNA_SM_df=data.frame(Annotation=getName(rRNA_fasta),
                      rnacentral_len=getLength(rRNA_fasta),
                      rnacentral_Seq=unlist(getSequence(rRNA_fasta,as.string = TRUE)))
DGE_set_all_SM_rRNA=merge(DGE_set_all_SM_rRNA,rRNA_SM_df,by="Annotation",all.x=TRUE)
DGE_set_all_SM_rRNA$rnacentral_ID="no_ID"

#merge all


colnames(DGE_set_all)
colnames(DGE_set_all_SM_tRNA)
DGE_set_all=rbind(DGE_set_all,DGE_set_all_SM_tRNA)
DGE_set_all=rbind(DGE_set_all,DGE_set_all_SM_rRNA)
# DGE_set$expression=ifelse(DGE_set$log2FoldChange>0,"Upregulated","Downregulated")
# write.csv(DGE_set, "F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation.csv",
#           row.names = FALSE, quote = FALSE)
#save(DGE_set_all,file="F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_filtered_ELAC_new.RData")
#load("G:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_filtered.RData")

colnames(DGE_set_all)
DGE_set_all$seq_start=as.numeric(DGE_set_all$Position)
DGE_set_all$seq_end=DGE_set_all$seq_start+DGE_set_all$Seq_length
DGE_set_all$cleavage="both ends cleavage"
DGE_set_all$cleavage=ifelse(DGE_set_all$Seq_length==DGE_set_all$rnacentral_len,"not cleaved",DGE_set_all$cleavage)
DGE_set_all$cleavage=ifelse(DGE_set_all$seq_start==1,"5'end retained",DGE_set_all$cleavage)
DGE_set_all$cleavage=ifelse(DGE_set_all$seq_end==DGE_set_all$rnacentral_len,"3'end retained",DGE_set_all$cleavage)
save(DGE_set_all,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad.RData")
write.xlsx(DGE_set_all, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_norm_counts_ELAC_filtered.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

################################################################################

################################################################################
#-------------------------------------mRNA--------------------------------------
# DGE_set=read.csv("E:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_annotation_new.csv",
#                  header=TRUE,sep=";")
# mRNA_set=subset(DGE_set,DGE_set$RNA_type=="mRNA")
# nrow(mRNA_set)
# colnames(mRNA_set)
# mRNA_anno=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/swiss_pfam_planmine_rnacentral_annotation.csv",
#                    header=TRUE)
# nrow(mRNA_anno) #31008
# colnames(mRNA_anno)
# 
# mRNA_anno=mRNA_anno[c("Transcript","UniProt_protein","UniProt_gene","UniProt_GO_BP","UniProt_GO_CC",
#                       "UniProt_GO_MF","Pfam_GO","Planmine_GO_MF","Planmine_GO_BP")]
# mRNA_anno$Annotation=mRNA_anno$Transcript
# mRNA_annotated=merge(mRNA_set,mRNA_anno,by="Annotation",all.x=TRUE)
# write.table(mRNA_annotated, "G:/smallRNA/calculated_data_bowtie2_end_to_end/mRNA_annotated.csv",
#             row.names = FALSE, quote = FALSE,sep=";")
################################################################################
#-------------------------------------rRNA--------------------------------------
rRNA_set=subset(DGE_set_all,DGE_set_all$RNA_type=="rRNA_fragments")
rRNA_set[1,]
rRNA_set$rRNA_type="rRNA_fragments"
nrow(rRNA_set)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("LSU", ignore_case = TRUE)),"LSU",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("SSU", ignore_case = TRUE)),"SSU",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("large", ignore_case = TRUE)),"LSU",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("28S", ignore_case = TRUE)),"28S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("16S", ignore_case = TRUE)),"16S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("18S", ignore_case = TRUE)),"18S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("23S", ignore_case = TRUE)),"23S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("5.8S", ignore_case = TRUE)),"5.8S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("5S", ignore_case = TRUE)),"5S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("small", ignore_case = TRUE)),"SSU",rRNA_set$rRNA_type)

LSU_to_check=rRNA_set[rRNA_set$rRNA_type=="LSU",]
LSU_to_check=unique(LSU_to_check$rnacentral_ID)
SSU_to_check=rRNA_set[rRNA_set$rRNA_type=="SSU",]
SSU_to_check=unique(SSU_to_check$rnacentral_ID)
rnaCentralRetrieveEntry("URS00021748B9_256318")

write.xlsx(c(LSU_to_check,SSU_to_check), "F:/PARN_ELAC_silencing/smallRNA/rRNA_to_check.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(c(LSU_to_check,SSU_to_check), "F:/PARN_ELAC_silencing/smallRNA/rRNA_to_check.txt", 
           col.names = FALSE, row.names = FALSE, append = FALSE,quote = FALSE)
# for (i in 1:nrow(DGE_set_all)){
#   DGE_set_all$rnacentral_Seq[i]=rnaCentralRetrieveEntry(DGE_set_all$rnacentral_ID[i])$sequence
#   DGE_set_all$rnacentral_len[i]=rnaCentralRetrieveEntry(DGE_set_all$rnacentral_ID[i])$sequenceLength
# }


#install.packages("rvest")
library(rvest)
rRNA_to_check=c(LSU_to_check,SSU_to_check)
rRNA_to_check_df=as.data.frame(rRNA_to_check)
for (i in 1:nrow(rRNA_to_check_df)) {
  url <- paste0("https://rnacentral.org/rna/", gsub("_", "/", rRNA_to_check[i]))
  webpage <- read_html(url)
  rna_type_elements <- html_nodes(webpage, "ol.breadcrumb a")
  rna_types <- html_text(rna_type_elements)
  rna_type <- paste(rna_types, collapse = " ")
  rRNA_to_check_df$rRNA_type[i]=rna_type
  #cat("ID:", id, "\n")
  #cat("RNA Type:", rna_type, "\n\n")
}
rRNA_to_check_df$correct_RNA=substr(rRNA_to_check_df$rRNA_type,nchar(rRNA_to_check_df$rRNA_type)-7,nchar(rRNA_to_check_df$rRNA_type))
rRNA_to_check_df$correct_RNA=sub(" ", "_", rRNA_to_check_df$correct_RNA)
unique(rRNA_to_check_df$correct_RNA)


rRNA_set$rRNA_type=ifelse(rRNA_set$rnacentral_ID %in% rRNA_to_check_df$rRNA_to_check[rRNA_to_check_df$correct_RNA=="28S_rRNA"],"28S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(rRNA_set$rnacentral_ID %in% rRNA_to_check_df$rRNA_to_check[rRNA_to_check_df$correct_RNA=="18S_rRNA"],"18S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(rRNA_set$rnacentral_ID %in% rRNA_to_check_df$rRNA_to_check[rRNA_to_check_df$correct_RNA=="23S_rRNA"],"23S_rRNA",rRNA_set$rRNA_type)
rRNA_set$rRNA_type=ifelse(rRNA_set$rnacentral_ID %in% rRNA_to_check_df$rRNA_to_check[rRNA_to_check_df$correct_RNA=="16S_rRNA"],"16S_rRNA",rRNA_set$rRNA_type)
table(rRNA_set$rRNA_type)
rRNA_set[rRNA_set$rRNA_type=="rRNA_fragments",]
rRNA_set$rRNA_type=ifelse(str_detect(rRNA_set$Annotation, regex("Putative_SpacerA", ignore_case = TRUE)),"28S_rRNA",rRNA_set$rRNA_type)

#[1] "28S_rRNA" ""         "LSU_rRNA" "23S_rRNA" "18S_rRNA" "16S_rRNA"








#URS0000C24971_6313, URS00003C9343_50054,URS0000C80C4F_1005395 - LSU, 
head(rRNA_set$rnacentral_ID)
unique(rRNA_set$rRNA_type)
#subset(rRNA_set,rRNA_set$rRNA_type=="rRNA_fragments")
#rnaCentralRetrieveEntry("URS00003278A1_297250")

# rRNA=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_rRNA_new.csv",header=TRUE)
# table(rRNA$rRNA_type) #28S_rRNA twice????
# 
# which(rRNA[rRNA$rRNA_type=="28S_rRNA"])
# rRNA$rRNA_type=="28S_rRNA"
# subset_normal=subset(rRNA,rRNA$rRNA_type=="28S_rRNA" | rRNA$rRNA_type=="16S_rRNA" | rRNA$rRNA_type=="18S_rRNA" | 
#                        rRNA$rRNA_type=="23S_rRNA" | rRNA$rRNA_type=="5.8S_rRNA" | rRNA$rRNA_type=="5S_rRNA"  | 
#                        rRNA$rRNA_type=="LSU rRNA" | rRNA$rRNA_type=="SSU_rRNA")
# subset_strange=subset(rRNA,rRNA$rRNA_type!="28S_rRNA" & rRNA$rRNA_type!="16S_rRNA" & rRNA$rRNA_type!="18S_rRNA" & 
#                         rRNA$rRNA_type!="23S_rRNA" & rRNA$rRNA_type!="5.8S_rRNA" & rRNA$rRNA_type!="5S_rRNA"  & 
#                         rRNA$rRNA_type!="LSU rRNA" & rRNA$rRNA_type!="SSU_rRNA")
# subset_strange$rRNA_type="28S_rRNA"
# rRNA_set_normal=rbind(subset_normal,subset_strange)
# nrow(rRNA_set_normal)
# table(rRNA_set_normal$rRNA_type)

rRNA_df=as.data.frame(table(rRNA_set$rRNA_type,rRNA_set$set,rRNA_set$expression))

#typeof(rRNA$rRNA_type)
#unique(rRNA$rRNA_type)
colnames(rRNA_df)=c("rRNA_type","comparison","Expression","Frequency")
#trna_df=subset(trna_df,trna_df$Frequency>0)
#trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
unique(rRNA_df$comparison)
#Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
# Elac_vs_GFP_dpa3_df=rRNA_df[rRNA_df$comparison=="Elac_vs_GFP_dpa3",]
# Elac_vs_GFP_dpa5_df=rRNA_df[rRNA_df$comparison=="Elac_vs_GFP_dpa5",]
Elac_vs_WT_dpa3_df=rRNA_df[rRNA_df$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df=rRNA_df[rRNA_df$comparison=="Elac_vs_WT_dpa5",]
# GFP_vs_WT_dpa3_df=rRNA_df[rRNA_df$comparison=="GFP_vs_WT_dpa3",]
# GFP_vs_WT_dpa5_df=rRNA_df[rRNA_df$comparison=="GFP_vs_WT_dpa5",]
# PARN_vs_GFP_dpa5_df=rRNA_df[rRNA_df$comparison=="PARN_vs_GFP_dpa5",]
# PARN_vs_WT_dpa3_df=rRNA_df[rRNA_df$comparison=="PARN_vs_WT_dpa3",]
# PARN_vs_WT_dpa5_df=rRNA_df[rRNA_df$comparison=="PARN_vs_WT_dpa5",]
# WT_dpa3_vs_dpa5_df=rRNA_df[rRNA_df$comparison=="WT_dpa3_vs_dpa5",]

# Elac_vs_GFP_dpa3_reshape=reshape(Elac_vs_GFP_dpa3_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_GFP_dpa3_reshape[is.na(Elac_vs_GFP_dpa3_reshape)] <- 0
# Elac_vs_GFP_dpa5_reshape=reshape(Elac_vs_GFP_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_GFP_dpa5_reshape[is.na(Elac_vs_GFP_dpa5_reshape)] <- 0
Elac_vs_WT_dpa3_reshape=reshape(Elac_vs_WT_dpa3_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape[is.na(Elac_vs_WT_dpa3_reshape)] <- 0
Elac_vs_WT_dpa5_reshape=reshape(Elac_vs_WT_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape[is.na(Elac_vs_WT_dpa5_reshape)] <- 0
# GFP_vs_WT_dpa3_reshape=reshape(GFP_vs_WT_dpa3_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa3_reshape[is.na(GFP_vs_WT_dpa3_reshape)] <- 0
# GFP_vs_WT_dpa5_reshape=reshape(GFP_vs_WT_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa5_reshape[is.na(GFP_vs_WT_dpa5_reshape)] <- 0
# PARN_vs_WT_dpa3_reshape=reshape(PARN_vs_WT_dpa3_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa3_reshape[is.na(PARN_vs_WT_dpa3_reshape)] <- 0
# PARN_vs_WT_dpa5_reshape=reshape(PARN_vs_WT_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa5_reshape[is.na(PARN_vs_WT_dpa5_reshape)] <- 0
# WT_dpa3_vs_dpa5_reshape=reshape(WT_dpa3_vs_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# WT_dpa3_vs_dpa5_reshape[is.na(WT_dpa3_vs_dpa5_reshape)] <- 0
# PARN_vs_GFP_dpa5_reshape=reshape(PARN_vs_GFP_dpa5_df, idvar = "rRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_GFP_dpa5_reshape[is.na(PARN_vs_GFP_dpa5_reshape)] <- 0

# elac_test=trna_df[trna_df$comparison=="WT_vs_Elac_dpa3",]
# par(mar=pyramid.plot(elac_test$Frequency[elac_test$Expression=="Upregulated",],elac_test$Frequency[elac_test$Expression=="Upregulated",],labels=agelabels,
#                      main="Australian population pyramid 2002",lxcol=mcol,rxcol=fcol,
#                      gap=0.5,show.values=TRUE))
#install.packages("plotrix")
library(plotrix)
#display.brewer.all()
#cols <- brewer.pal(3, "BuGn")
par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
# pyramid.plot(Elac_vs_GFP_dpa3_reshape$Frequency.Downregulated, Elac_vs_GFP_dpa3_reshape$Frequency.Upregulated,labels= Elac_vs_GFP_dpa3_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=1.3, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(Elac_vs_GFP_dpa3_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_GFP_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in Elac_vs_GFP_dpa3")
# pyramid.plot(Elac_vs_GFP_dpa5_reshape$Frequency.Downregulated, Elac_vs_GFP_dpa5_reshape$Frequency.Upregulated,labels= Elac_vs_GFP_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=0.5, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(Elac_vs_GFP_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_GFP_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in Elac_vs_GFP_dpa5")
pyramid.plot(Elac_vs_WT_dpa3_reshape$Frequency.loss, Elac_vs_WT_dpa3_reshape$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
             gap=11, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.accumulation)),main="Differentially accumulated rRFs at 3 dpa")
pyramid.plot(Elac_vs_WT_dpa5_reshape$Frequency.loss, Elac_vs_WT_dpa5_reshape$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
             gap=6, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.accumulation)),main="Differentially accumulated rRFs at 5 dpa")
# 
# pyramid.plot(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated, PARN_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa3_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=2.25, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in PARN_vs_WT_dpa3")
# pyramid.plot(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated, PARN_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=6.35, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in PARN_vs_WT_dpa5")

# pyramid.plot(WT_dpa3_vs_dpa5_reshape$Frequency.Downregulated, WT_dpa3_vs_dpa5_reshape$Frequency.Upregulated,labels= WT_dpa3_vs_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=0.25, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(WT_dpa3_vs_dpa5_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(WT_dpa3_vs_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in WT_dpa3_vs_dpa5")

# pyramid.plot(PARN_vs_GFP_dpa5_reshape$Frequency.Downregulated, PARN_vs_GFP_dpa5_reshape$Frequency.Upregulated,labels= PARN_vs_GFP_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=0.25, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(PARN_vs_GFP_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_vs_GFP_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in PARN_vs_GFP_dpa5")
# pyramid.plot(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated, GFP_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa3_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=3.2, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in GFP_vs_WT_dpa3")
# pyramid.plot(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated, GFP_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa5_reshape$rRNA_type,lxcol="lightgreen", rxcol="tomato1",unit = "Frequency",
#              gap=1.25, space=0.15, top.labels = c("Downregulated", "rRF type","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated rRFs in GFP_vs_WT_dpa5")
# 
rRNA_set
write.xlsx(rRNA_set, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_rRNA_ELAC.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
################################################################################
#-------------------------------------tRNA--------------------------------------

DGE_set_with_len_trna=subset(DGE_set_all,DGE_set_all$RNA_type=="tRNA_fragments")
DGE_set_with_len_trna$tRNA_type_sub=paste("tRNA",sub(".*tRNA", "", DGE_set_with_len_trna$Annotation),sep="")
DGE_set_with_len_trna$rnacentral_Seq=print(gsub("T", "U", DGE_set_with_len_trna$rnacentral_Seq) )
DGE_set_with_len_trna$Sequence=print(gsub("T", "U", DGE_set_with_len_trna$Sequence) )
colnames(DGE_set_with_len_trna)
length(unique(DGE_set_with_len_trna$rnacentral_Seq))
nrow(unique(DGE_set_with_len_trna[,c("rnacentral_Seq","Annotation")]))
anno_seq=unique(DGE_set_with_len_trna[,c("rnacentral_Seq","Annotation")])
as.data.frame(table(anno_seq$rnacentral_Seq))
anno_seq[anno_seq$rnacentral_Seq=="GCCCAGAUGGUGAAAUUGGUAGACACGCCAGCUUCAGGUGCUGGUGACCUUACGGUCGUGGAAGUUCGAGUCUUCUUCUGGGCACCA",]
unique_tRNA=unique(unique(DGE_set_with_len_trna[, c('rnacentral_ID', 'rnacentral_Seq')]))
unique_tRNA$rnacentral_ID=ifelse(unique_tRNA$rnacentral_ID=="no_ID",paste("tRNA",sample(10),sep=""),unique_tRNA$rnacentral_ID)
nrow(unique_tRNA)
write.table(unique_tRNA$rnacentral_ID, "F:/PARN_ELAC_silencing/smallRNA/tRNA_ELAC2_DEG_unique_ID.txt", 
           col.names = FALSE, row.names = FALSE, append = FALSE,quote = FALSE)
write.table(unique_tRNA$rnacentral_Seq, "F:/PARN_ELAC_silencing/smallRNA/tRNA_ELAC2_DEG_unique_seq.txt", 
            col.names = FALSE, row.names = FALSE, append = FALSE,quote = FALSE)
unique_tRNA_tRNAScan=read.xlsx("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/unique_trna_seq_ELAC.xlsx",1)
colnames(unique_tRNA_tRNAScan)=c("rnacentral_Seq","tRNA_Scan_type","tRNA-AAA","Structure","tRNA_number")
unique_tRNA_tRNAScan_merge=merge(unique_tRNA,unique_tRNA_tRNAScan,by="rnacentral_Seq",all.x=TRUE)
unique_tRNA_tRNAScan_merge_NA=unique_tRNA_tRNAScan_merge[is.na(unique_tRNA_tRNAScan_merge$Structure),]
write.table(unique_tRNA_tRNAScan_merge_NA$rnacentral_ID, "F:/PARN_ELAC_silencing/smallRNA/tRNA_ELAC2_DEG_unique_ID.txt", 
            col.names = FALSE, row.names = FALSE, append = FALSE,quote = FALSE)
write.table(unique_tRNA_tRNAScan_merge_NA$rnacentral_Seq, "F:/PARN_ELAC_silencing/smallRNA/tRNA_ELAC2_DEG_unique_seq.txt", 
            col.names = FALSE, row.names = FALSE, append = FALSE,quote = FALSE)
# write.xlsx(DGE_set_with_len_trna, "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_trna_ELAC.xlsx", 
#            col.names = TRUE, row.names = FALSE, append = FALSE)
# unique(DGE_set_with_len_trna$rnacentral_Seq)
# write.xlsx(unique(DGE_set_with_len_trna$rnacentral_Seq), "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/unique_trna_seq_ELAC.xlsx", 
#            col.names = TRUE, row.names = FALSE, append = FALSE)
# ?read.xlsx

unique_tRNA_tRNAScan=read.xlsx("F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/unique_trna_seq_ELAC.xlsx",1)
colnames(unique_tRNA_tRNAScan)=c("rnacentral_Seq","tRNA_Scan_type","tRNA-AAA","Structure","tRNA_number")
unique_tRNA=unique(unique(DGE_set_with_len_trna[, c('rnacentral_ID', 'rnacentral_Seq')]))
nrow(unique_tRNA)
unique_tRNA$rnacentral_ID=ifelse(unique_tRNA$rnacentral_ID=="no_ID",paste("tRNA",sample(10),sep=""),unique_tRNA$rnacentral_ID)

unique(DGE_set_with_len_trna$tRNA_type_sub)
colnames(DGE_set_with_len_trna)

DGE_set_with_len_trna=merge(DGE_set_with_len_trna,unique_tRNA_tRNAScan,by="rnacentral_Seq",all.x = TRUE)
DGE_set_with_len_trna=DGE_set_with_len_trna[DGE_set_with_len_trna$Structure != "no_structure_predicted",]
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNAURS00001E3937_7091_Bombyx_mori_(domestic_silkworm)_transfer_RNA-Gly","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Gly-GCC-9-1)","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA_His_GUG","tRNA-His",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Gly_for_anticodon_GCC","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Gly-GCC-5-1)","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Arg-ACG-1-1)","tRNA-Arg",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA_Gly_GCC","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Asn-GTT-2-1)","tRNA-Asn",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA_Gly","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# unique(DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNAURS0000418301_6210_Echinococcus_granulosus_transfer_RNA-Pro","tRNA-Pro",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA_Leu_CAG","tRNA-Leu",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA,_Gly","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNAURS000032F6AC_10090_Mus_musculus_(house_mouse)_transfer_RNA-Lys","tRNA-Lys",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Asp-GTC-1_1_to_3)","tRNA-Asp",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Arg-CCT-1-1)","tRNA-Arg",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Gly-GCC-11-1)","tRNA-Gly",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA-Pro_for_anticodon_CGG","tRNA-Pro",DGE_set_with_len_trna$tRNA_type_sub)
# DGE_set_with_len_trna$tRNA_type_sub=ifelse(DGE_set_with_len_trna$tRNA_type_sub=="tRNA_Lys_CUU","tRNA-Lys",DGE_set_with_len_trna$tRNA_type_sub)
# 
# not_anno_DGE_set_with_len_trna=subset(DGE_set_with_len_trna,!startsWith(DGE_set_with_len_trna$tRNA_type_sub,"tRNA-"))
# nrow(not_anno_DGE_set_with_len_trna) #35
# write.csv(DGE_set_with_len_trna, "F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_trna.csv",
#           row.names = FALSE, quote = FALSE)
# DGE_set_with_len_trna=read.csv("F:/smallRNA/calculated_data_bowtie2_end_to_end/DGE_tRNA_new.csv")
# DGE_set_with_len_trna_structure=read.csv("E:/project/smallRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_trna_structure.csv",header=TRUE,sep=";")
# DGE_set_with_len_trna_structure$structure_len=nchar(DGE_set_with_len_trna_structure$Parent_seq_structure)
# DGE_set_with_len_trna_structure$seq_end_corrected=ifelse(DGE_set_with_len_trna_structure$structure_len<DGE_set_with_len_trna_structure$seq_end,DGE_set_with_len_trna_structure$structure_len,DGE_set_with_len_trna_structure$seq_end)
# DGE_set_with_len_trna_structure$tRF_structure=substr(DGE_set_with_len_trna_structure$Parent_seq_structure, DGE_set_with_len_trna_structure$seq_start,DGE_set_with_len_trna_structure$seq_end_corrected)
# DGE_set_with_len_trna$Sequence
# 
# write.table(DGE_set_with_len_trna_structure, "G:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_trna_with_structure.csv",
#             row.names = FALSE, quote = FALSE,sep=";")
# DGE_set_with_len_trna_structure=read.csv("F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_set_with_len_trna_with_structure.csv",header=TRUE,sep=";")
# write.table(DGE_set_with_len_trna$Sequence, "F:/smallRNAwithAdapters/miRNA/calculated_data_bowtie2_end_to_end/DGE_trna_list.txt",
#             row.names = FALSE, quote = FALSE)
#DGE_set_with_len_trna$expr=ifelse(DGE_set_with_len_trna$log2FoldChange>0,"Upregulated","Downregulated")

#data <- as.matrix(mtcars)
#plot for tRNA_codon_type
DGE_set_with_len_trna=subset(DGE_set_with_len_trna,DGE_set_with_len_trna$tRNA_Scan_type!="not_tRNA")

# DGE_set_with_len_trna_structure$tRNAScan_type
# DGE_set_with_len_trna$tRNA_type_sub
# DGE_set_with_len_trna_structure$tRNA_type_sub=paste("tRNA-",DGE_set_with_len_trna_structure$tRNAScan_type,sep="")
# DGE_set_with_len_trna_structure$tRNA_type_sub=ifelse(DGE_set_with_len_trna_structure$tRNA_type_sub=="tRNA-Ser/Supressor","tRNA-Ser",DGE_set_with_len_trna_structure$tRNA_type_sub)
trna_df=as.data.frame(table(DGE_set_with_len_trna$set,DGE_set_with_len_trna$tRNA_Scan_type,DGE_set_with_len_trna$expression))
colnames(trna_df)=c("comparison","tRNA_type","Expression","Frequency")
#trna_df=subset(trna_df,trna_df$Frequency>0)
# trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
# unique(trna_df$comparison)
#Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
Elac_vs_WT_dpa3_df=trna_df[trna_df$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df=trna_df[trna_df$comparison=="Elac_vs_WT_dpa5",]
# GFP_vs_WT_dpa3_df=trna_df[trna_df$comparison=="GFP_vs_WT_dpa3",]
# GFP_vs_WT_dpa5_df=trna_df[trna_df$comparison=="GFP_vs_WT_dpa5",]
# PARN_vs_WT_dpa3_df=trna_df[trna_df$comparison=="PARN_vs_WT_dpa3",]
# PARN_vs_WT_dpa5_df=trna_df[trna_df$comparison=="PARN_vs_WT_dpa5",]
# PARN_dpa3_vs_dpa5_df=trna_df[trna_df$comparison=="PARN_dpa3_vs_dpa5",]

Elac_vs_WT_dpa3_reshape=reshape(Elac_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape[is.na(Elac_vs_WT_dpa3_reshape)] <- 0
Elac_vs_WT_dpa5_reshape=reshape(Elac_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape[is.na(Elac_vs_WT_dpa5_reshape)] <- 0
# GFP_vs_WT_dpa3_reshape=reshape(GFP_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa3_reshape[is.na(GFP_vs_WT_dpa3_reshape)] <- 0
# GFP_vs_WT_dpa5_reshape=reshape(GFP_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa5_reshape[is.na(GFP_vs_WT_dpa5_reshape)] <- 0
# PARN_vs_WT_dpa3_reshape=reshape(PARN_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa3_reshape[is.na(PARN_vs_WT_dpa3_reshape)] <- 0
# PARN_vs_WT_dpa5_reshape=reshape(PARN_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa5_reshape[is.na(PARN_vs_WT_dpa5_reshape)] <- 0
# PARN_dpa3_vs_dpa5_reshape=reshape(PARN_dpa3_vs_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_dpa3_vs_dpa5_reshape[is.na(PARN_dpa3_vs_dpa5_reshape)] <- 0

# elac_test=trna_df[trna_df$comparison=="WT_vs_Elac_dpa3",]
# par(mar=pyramid.plot(elac_test$Frequency[elac_test$Expression=="Upregulated",],elac_test$Frequency[elac_test$Expression=="Upregulated",],labels=agelabels,
#                      main="Australian population pyramid 2002",lxcol=mcol,rxcol=fcol,
#                      gap=0.5,show.values=TRUE))
#install.packages("plotrix")

#display.brewer.all()
#cols <- brewer.pal(3, "BuGn")
par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
pyramid.plot(Elac_vs_WT_dpa3_reshape$Frequency.loss, Elac_vs_WT_dpa3_reshape$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
             gap=10, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.accumulation)),main="Differentially accumulated tRF at 3 dpa")+theme(
               plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))

pyramid.plot(Elac_vs_WT_dpa5_reshape$Frequency.loss, Elac_vs_WT_dpa5_reshape$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
             gap=4.5, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.accumulation)),main="Differentially accumulated tRF at 5 dpa")

#get tRF structure type
colnames(DGE_set_with_len_trna)
DGE_set_with_len_trna$rnacentral_Seq=str_to_upper(DGE_set_with_len_trna$rnacentral_Seq, locale = "en")
DGE_set_with_len_trna$New_str=print(gsub(">", "(", DGE_set_with_len_trna$Structure) )
DGE_set_with_len_trna$New_str=print(gsub("<", ")", DGE_set_with_len_trna$New_str) )
DGE_set_with_len_trna$parent_len=nchar(DGE_set_with_len_trna$rnacentral_Seq)
DGE_set_with_len_trna$parent_str_len=nchar(DGE_set_with_len_trna$Structure)
DGE_set_with_len_trna$str_seq_same_len=ifelse(DGE_set_with_len_trna$parent_len==DGE_set_with_len_trna$parent_str_len,"same_len","structure_len_diff") 
DGE_set_with_len_trna$parent_seq=ifelse(DGE_set_with_len_trna$str_seq_same_len=="same_len",DGE_set_with_len_trna$rnacentral_Seq,"to_be_trimmed")
table(DGE_set_with_len_trna$str_seq_same_len)
DGE_set_with_len_trna$parent_seq=ifelse((DGE_set_with_len_trna$parent_seq=="to_be_trimmed" & endsWith(DGE_set_with_len_trna$rnacentral_Seq,"CCA")),
                                        substr(DGE_set_with_len_trna$rnacentral_Seq,1,nchar(DGE_set_with_len_trna$rnacentral_Seq)-3),DGE_set_with_len_trna$parent_seq)
DGE_set_with_len_trna$parent_seq=ifelse((DGE_set_with_len_trna$parent_seq=="to_be_trimmed" & endsWith(DGE_set_with_len_trna$rnacentral_Seq,"CC")),
                                        substr(DGE_set_with_len_trna$rnacentral_Seq,1,nchar(DGE_set_with_len_trna$rnacentral_Seq)-2),DGE_set_with_len_trna$parent_seq)
                                       
DGE_set_with_len_trna[DGE_set_with_len_trna$str_seq_same_len=="structure_len_diff",]  
DGE_set_with_len_trna$parent_len_dif=DGE_set_with_len_trna$parent_len-DGE_set_with_len_trna$parent_str_len  
DGE_set_with_len_trna[DGE_set_with_len_trna$parent_len_dif>3,] 
#validate manualy GAUGCGGAUCAGUGGUAGAAUGCUCGCCUGCCACGCGGGCGGCCCGGGUUCGGUUCCCGGCCGAUGCACCA      and GGCCGUGAUCGUCUAGUGGUUAGGACCCCACGUUGUGGCCGUGGUAACCCAGGUUCGAAUCCUGGUCACGGCACCA 
DGE_set_with_len_trna$parent_seq=ifelse(DGE_set_with_len_trna$rnacentral_Seq=="GAUGCGGAUCAGUGGUAGAAUGCUCGCCUGCCACGCGGGCGGCCCGGGUUCGGUUCCCGGCCGAUGCACCA",
                                        "GCGGAUCAGUGGUAGAAUGCUCGCCUGCCACGCGGGCGGCCCGGGUUCGGUUCCCGGCCGAUGCA",DGE_set_with_len_trna$parent_seq)
DGE_set_with_len_trna$parent_seq=ifelse(DGE_set_with_len_trna$rnacentral_Seq=="GGCCGUGAUCGUCUAGUGGUUAGGACCCCACGUUGUGGCCGUGGUAACCCAGGUUCGAAUCCUGGUCACGGCACCA",
                                        "GCCGUGAUCGUCUAGUGGUUAGGACCCCACGUUGUGGCCGUGGUAACCCAGGUUCGAAUCCUGGUCACGGCA",DGE_set_with_len_trna$parent_seq)
save(DGE_set_with_len_trna,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_filtered_mad.RData")
DGE_set_with_len_trna$parent_seq=print(gsub("T", "U", DGE_set_with_len_trna$parent_seq) )


DGE_set_with_len_trna_unique=unique(DGE_set_with_len_trna[,c("Sequence","Position","Seq_length","Annotation","parent_seq","Structure")])
nrow(DGE_set_with_len_trna_unique)
DGE_set_with_len_trna_unique$tRF_ID=sprintf("tRF%d",seq(1:nrow(DGE_set_with_len_trna_unique)))

write.xlsx(DGE_set_with_len_trna_unique, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_unique.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_filtered_mad.RData")
DGE_set_with_len_trna
################################################################################
#-------------------------------classify tRNA fragments type--------------------
DGE_set_with_len_trna_unique_structure=unique(DGE_set_with_len_trna_unique$Structure)
DGE_set_with_len_trna_unique$start=DGE_set_with_len_trna_unique$Position
DGE_set_with_len_trna_unique$end=as.numeric(DGE_set_with_len_trna_unique$start)+as.numeric(DGE_set_with_len_trna_unique$Seq_length)-1
DGE_set_with_len_trna_unique_tRF3=DGE_set_with_len_trna_unique[as.numeric(DGE_set_with_len_trna_unique$start)<3,]
DGE_set_with_len_trna_unique_tRF3$start=1

output_dir ="F:/PARN_ELAC_silencing/smallRNA/plots/ELAC2_final/tRNA_fragments/"
setwd(output_dir)
### Create a CT file from bracket notation
# ct=makeCt("(((...(((...)))...(((...)))...)))","AAAUUUCCCAAAGGGUUUAAAGGGUUUCCCUUU")
# coord=ct2coord(ct)
# RNAPlot(coord,hl=c("GGGUUU","AAAUUU"),seqcols=c(2,4),labTF=TRUE)
for (i in 1:nrow(DGE_set_with_len_trna_unique_tRF3)){
  new_str=print(gsub(">", "(", DGE_set_with_len_trna_unique_tRF3$Structure[i]) )
  new_str=print(gsub("<", ")", new_str) )
  ct=makeCt(new_str,DGE_set_with_len_trna_unique_tRF3$parent_seq[i])
  coord=ct2coord(ct)
  ranges=data.frame(min=DGE_set_with_len_trna_unique_tRF3$start[i],max=DGE_set_with_len_trna_unique_tRF3$end[i],
                    col=2,desc=paste("tRF",DGE_set_with_len_trna_unique_tRF3$Sequence[i]))
  pdf(file= paste0(output_dir, DGE_set_with_len_trna_unique_tRF3$tRF_ID[i], ".pdf"))
  RNAPlot(coord,nt=TRUE,ranges = ranges,
          labTF=TRUE,main=DGE_set_with_len_trna_unique_tRF3$tRF_ID[i])
  # call this function to save the file 
  dev.off()
}


#draw separately tRNA

tRNA_seq="GCTGTTCAGTGGTAGAATGCTCGCCTGCCACGTCGGGCGACCCGGGTTCGATTCCGGGCCGAGCA"
tRNA_str=">>>>>>>.......<<<<.>>>>>.......<.<<<<....>>.>>.......<<.<<...<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=6,max=26,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF65", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCGAUAGUGGUCCAACGGCUAUGAUUGGUGCUUCCCAAGCACCUGACUCGGGUUCGACUCCCGGCUAUCGCA"
tRNA_str=">>>>>>>..>>>..........<<<.>>>>>.......<<<<<....>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=49,max=75,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF63", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCCGCGGUGGCGGAACUGGCAGACGCAAGGGACUUAAAAUCCCUCGGAUAGAAAUAUCCGUACCGGUUCGAUUCCGGUUCGCGGCA"
tRNA_str=">>>>>>>..>>>...........<<<.>>>>>.......<<<<<.>>>>>....<<<<<..>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=74,max=89,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF58", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCCGACAUGGUGGAAUUGGUAGACACGCUAUCUUGAGGGGGUAGUGGCCCCAGGCUGUGCGAGUUCGAGUCUCGCUGUCGGCA"
tRNA_str=">>>>>>>..>>>...........<<<.>>>>>.......<<<<<.>>>>...<<<<..>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=49,max=86,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF57", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCCCAGAUGGUGAAAUUGGUAGACACGCCAGCUUCAGGUGCUGGUGACCUUACGGUCGUGGAAGUUCGAGUCUUCUUCUGGGCA"
tRNA_str=">>>>>>>..>>>...........<<<.>>>>>.......<<<<<.>>>>....<<<<..>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=22,max=87,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF55", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCCCAGAUGGUGAAAUUGGUAGACACGCCAGCUUCAGGUGCUGGUGACCUUACGGUCGUGGAAGUUCGAGUCUUCUUCUGGGCA"
tRNA_str=">>>>>>>..>>>...........<<<.>>>>>.......<<<<<.>>>>....<<<<..>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=49,max=87,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF54", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCCCAGAUAGCUCAGUCGGUAGAGCAGAGGAUUGAAAAUCCUCGUGUCGGCGGUUCGACUCCGUCUCUGGGCA"
tRNA_str=">>>>>>>..>>>>........<<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=55,max=75,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF53", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GCGGAUCAGUGGUAGAAUGCUCGCCUGCCACGCGGGCGGCCCGGGUUCGGUUCCCGGCCGAUGCA"
tRNA_str=">>.>.>>.......<<.<.>>>>>.......<<<<<....>>>>>.......<<<<<.....<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=57,max=71,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF49", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()

tRNA_seq="GAUACGAUGGCCGAGUGGUUAAGGCGAAGGAUGCAGGUUCCUUUGGGCAUUGCCCGCGCAGGUUCGAACCCUGCUCGUGUCG"
tRNA_str=">>>>>>>..>>>..........<<<.>>>>>.......<<<<<.>>>>...<<<<..>>>>>.......<<<<<<<<<<<<."
new_str=print(gsub(">", "(", tRNA_str) )
new_str=print(gsub("<", ")", new_str) )
ct=makeCt(new_str,tRNA_seq)
coord=ct2coord(ct)
ranges=data.frame(min=70,max=85,col=2,desc=paste("tRF"))
pdf(file= paste0(output_dir, "tRF48", ".pdf"))
RNAPlot(coord,nt=TRUE,ranges=ranges,labTF=TRUE)
dev.off()
################################################################################
unique_tRNA_fr_type=read.xlsx("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_unique.xlsx",1)
colnames(unique_tRNA_fr_type)
colnames(DGE_set_with_len_trna)
unique_tRNA_fr_type=unique_tRNA_fr_type[,c("Sequence","tRF_type")]

DGE_set_with_len_trna_with_type=merge(DGE_set_with_len_trna,unique_tRNA_fr_type,by="Sequence",all.x = TRUE)
DGE_set_with_len_trna_with_type$tRF_type
trna_df_type=as.data.frame(table(DGE_set_with_len_trna_with_type$set,DGE_set_with_len_trna_with_type$tRF_type,DGE_set_with_len_trna_with_type$expression))
colnames(trna_df_type)=c("comparison","tRNA_type","Expression","Frequency")
#trna_df=subset(trna_df,trna_df$Frequency>0)
# trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
# unique(trna_df$comparison)
#Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
Elac_vs_WT_dpa3_df_type=trna_df_type[trna_df_type$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df_type=trna_df_type[trna_df_type$comparison=="Elac_vs_WT_dpa5",]

Elac_vs_WT_dpa3_reshape_type=reshape(Elac_vs_WT_dpa3_df_type, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape_type[is.na(Elac_vs_WT_dpa3_reshape_type)] <- 0
Elac_vs_WT_dpa5_reshape_type=reshape(Elac_vs_WT_dpa5_df_type, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape_type[is.na(Elac_vs_WT_dpa5_reshape_type)] <- 0


#display.brewer.all()
#cols <- brewer.pal(3, "BuGn")
par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
pyramid.plot(Elac_vs_WT_dpa3_reshape_type$Frequency.loss, Elac_vs_WT_dpa3_reshape_type$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape_type$tRNA_type,lxcol="salmon", rxcol="navy",unit = "Frequency",
             gap=10, space=0.15, top.labels = c("Downregulated", "tRF-type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_type$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_type$Frequency.accumulation)),main="Differentially accumulated tRF at 3 dpa")+theme(
               plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))

pyramid.plot(Elac_vs_WT_dpa5_reshape_type$Frequency.loss, Elac_vs_WT_dpa5_reshape_type$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape_type$tRNA_type,lxcol="salmon", rxcol="navy",unit = "Frequency",
             gap=6, space=0.15, top.labels = c("Downregulated", "type-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_type$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_type$Frequency.accumulation)),main="Differentially accumulated tRF at 5 dpa")



#--------------------------------------------check tRNA trailer:
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_filtered_mad.RData")
DGE_set_with_len_trna
tRNA_fasta=read.fasta("F:/PARN_ELAC_silencing/smallRNA/tRNA/pretRNA_from_chr_and_mtDNA_multi_dedup.fasta")
getName(tRNA_fasta)
getLength(tRNA_fasta)
tRNA_trailer_pos=data.frame(ref=getName(tRNA_fasta),
                            ref_tRNA_len=getLength(tRNA_fasta))
#tRNA_trailer_pos$strand=substr(tRNA_trailer_pos$ref,nchar(tRNA_trailer_pos$ref)-1,nchar(tRNA_trailer_pos$ref)-1)
tRNA_trailer_pos$three_prime_end_start=tRNA_trailer_pos$ref_tRNA_len-19

unique(tRNA_trailer_pos$ref)
unique(DGE_set_with_len_trna$Annotation)
new_mapped_seq=merge(new_mapped_seq,tRNA_trailer_pos,by="ref",all.x=TRUE)




for (i in 1:nrow(new_mapped_seq)) {new_mapped_seq$maturetRNA_retained[i]=ifelse(((length(intersect(seq(new_mapped_seq$start[i],new_mapped_seq$end[i]),
                                                                                                   seq(21,new_mapped_seq$three_prime_end_start[i]))))/new_mapped_seq$mapped_len[i]>0.5),"tRF","just_tRF_end")}
#



#----------------------------combine all type info------------------------------
DGE_set_with_len_trna_with_type$tRF_type_AA_anticodon=paste(DGE_set_with_len_trna_with_type$tRF_type,
                                                            substr(DGE_set_with_len_trna_with_type$`tRNA-AAA`,6,nchar(DGE_set_with_len_trna_with_type$`tRNA-AAA`)),sep="-")
trna_df_type_AA_anticodon=as.data.frame(table(DGE_set_with_len_trna_with_type$set,DGE_set_with_len_trna_with_type$tRF_type_AA_anticodon,DGE_set_with_len_trna_with_type$expression))
colnames(trna_df_type_AA_anticodon)=c("comparison","tRNA_type","Expression","Frequency")
#trna_df=subset(trna_df,trna_df$Frequency>0)
# trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
# unique(trna_df$comparison)
#Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
Elac_vs_WT_dpa3_df_type_AA_anticodon=trna_df_type_AA_anticodon[trna_df_type_AA_anticodon$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df_type_AA_anticodon=trna_df_type_AA_anticodon[trna_df_type_AA_anticodon$comparison=="Elac_vs_WT_dpa5",]

Elac_vs_WT_dpa3_reshape_type_AA_anticodon=reshape(Elac_vs_WT_dpa3_df_type_AA_anticodon, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape_type_AA_anticodon[is.na(Elac_vs_WT_dpa3_reshape_type_AA_anticodon)] <- 0
Elac_vs_WT_dpa5_reshape_type_AA_anticodon=reshape(Elac_vs_WT_dpa5_df_type_AA_anticodon, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape_type_AA_anticodon[is.na(Elac_vs_WT_dpa5_reshape_type_AA_anticodon)] <- 0


#display.brewer.all()
#cols <- brewer.pal(3, "BuGn")
par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
pyramid.plot(Elac_vs_WT_dpa3_reshape_type_AA_anticodon$Frequency.loss, Elac_vs_WT_dpa3_reshape_type_AA_anticodon$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape_type_AA_anticodon$tRNA_type,lxcol="brown", rxcol="darkcyan",unit = "Frequency",
             gap=10, space=0.15, top.labels = c("Downregulated", "tRF-type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_type_AA_anticodon$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_type_AA_anticodon$Frequency.accumulation)),main="Differentially accumulated tRF at 3 dpa")+theme(
               plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))

pyramid.plot(Elac_vs_WT_dpa5_reshape_type_AA_anticodon$Frequency.loss, Elac_vs_WT_dpa5_reshape_type_AA_anticodon$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape_type_AA_anticodon$tRNA_type,lxcol="brown", rxcol="darkcyan",unit = "Frequency",
             gap=6, space=0.15, top.labels = c("Downregulated", "type-type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_type_AA_anticodon$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_type_AA_anticodon$Frequency.accumulation)),main="Differentially accumulated tRF at 5 dpa")




DGE_set_with_len_trna_with_type$Seq_length
DGE_set_with_len_trna_with_type[DGE_set_with_len_trna_with_type$tRF_type_AA_anticodon=="5'-tiRNA-Gly-GCC",]

save(DGE_set_with_len_trna_with_type,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_filtered_mad_all_types.RData")

write.xlsx(DGE_set_with_len_trna_with_type, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_all_type.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
DGE_set_with_len_trna_with_type



################################################################################
#-------------------------------miRNA-------------------------------------------
################################################################################
#load all sets
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad.RData")
colnames(DGE_set_all)
DGE_set_all[DGE_set_all$RNA_type=="other",]
table(DGE_set_all$set)
table(DGE_set_all$RNA_type)
DGE_set_with_len_mirna=subset(DGE_set_all,DGE_set_all$RNA_type=="miRNA")

DGE_set_with_len_mirna$miRNA=ifelse(str_detect(DGE_set_with_len_mirna$Annotation, regex("Schmidtea", ignore_case = TRUE)),
                                    str_extract(DGE_set_with_len_mirna$Annotation, "(?i)sme.*"),DGE_set_with_len_mirna$Annotation)

write.xlsx(DGE_set_with_len_mirna, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_mirna_ELAC_to_verify.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
?read.xlsx
correct_miRNA=read.xlsx("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_mirna_ELAC_correct.xlsx",1)
correct_miRNA=correct_miRNA[,-5]
colnames(DGE_set_with_len_mirna)
colnames(correct_miRNA)
incorrect_col <- c("Position","Seq_length","Annotation","rnacentral_len","seq_start",         
                  "seq_end" ,"cleavage" , "miRNA","cleavage_in_detail")
DGE_set_with_len_mirna_reduced = DGE_set_with_len_mirna[,!(names(DGE_set_with_len_mirna) %in% incorrect_col)]
DGE_set_with_len_mirna_corrected=merge(DGE_set_with_len_mirna_reduced,correct_miRNA,by="Sequence",all.x=TRUE)
#-----------------------------------------check cleavage for all sets-------------------------
colnames(DGE_set_all)

colnames(DGE_set_with_len_mirna_corrected)
DGE_set_all_not_miRNA=subset(DGE_set_all,DGE_set_all$RNA_type!="miRNA")
DGE_set_all_corrected=rbind(DGE_set_all_not_miRNA,DGE_set_with_len_mirna_corrected[,!(names(DGE_set_with_len_mirna_corrected) %in% c("miRNA","cleavage_in_detail"))])
save(DGE_set_all_corrected,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad_all_corrected.RData")
write.xlsx(DGE_set_all_corrected, "F:/PARN_ELAC_silencing/smallRNA/DGE_set_all_corrected_ELAC.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
#DGE_set_all_miRNA=subset(DGE_set_all,DGE_set_all$RNA_type=="miRNA")
# DGE_set$cleavage=ifelse(DGE_set$Seq_length==DGE_set$rnacentral_len,"not cleaved",DGE_set$cleavage)
# DGE_set$cleavage=ifelse(DGE_set$cleavage=="3'end cleavage","3'end retained",DGE_set$cleavage)
# DGE_set$cleavage=ifelse(DGE_set$cleavage=="5'end cleavage","5'end retained",DGE_set$cleavage)
DGE_set_cleavage=as.data.frame(table(DGE_set_all_corrected$RNA_type,DGE_set_all_corrected$cleavage))
colnames(DGE_set_cleavage)=c("RNA type","RNA fragment end","frequency")
#DGE_set_cleavage$`RNA type`
#DGE_set_cleavage$`RNA fragment end`
DGE_set_cleavage_table_sum=aggregate(frequency~`RNA type`, DGE_set_cleavage,sum)
colnames(DGE_set_cleavage_table_sum)=c("RNA type","sum")
DGE_set_cleavage=merge(DGE_set_cleavage,DGE_set_cleavage_table_sum,by="RNA type",all.x = TRUE)
DGE_set_cleavage$percentage=(DGE_set_cleavage$frequency/DGE_set_cleavage$sum)*100
#?ggplot
ggplot(data=DGE_set_cleavage,aes(`RNA fragment end`,percentage, fill = `RNA fragment end`)) +
  geom_col(position = position_dodge(width = 1)) +
  facet_grid(~`RNA type`) +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=25))

################################################################################
#---------------------------miRNA continuation----------------------------------
DGE_set_with_len_mirna_corrected
save(DGE_set_with_len_mirna_corrected,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad_miRNA_corrected.RData")
write.xlsx(DGE_set_with_len_mirna_corrected, "F:/PARN_ELAC_silencing/smallRNA/DGE_miRNA_corrected_ELAC.xlsx", 
           col.names = TRUE, row.names = FALSE, append = FALSE)

length(unique(DGE_set_with_len_mirna_corrected$miRNA))

miRNA_df=as.data.frame(table(DGE_set_with_len_mirna_corrected$set,DGE_set_with_len_mirna_corrected$miRNA,DGE_set_with_len_mirna_corrected$expression))
colnames(miRNA_df)=c("comparison","miRNA","Expression","Frequency")
#trna_df=subset(trna_df,trna_df$Frequency>0)
# trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
# unique(trna_df$comparison)
#Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
Elac_vs_WT_dpa3_df_miRNA=miRNA_df[miRNA_df$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df_miRNA=miRNA_df[miRNA_df$comparison=="Elac_vs_WT_dpa5",]

Elac_vs_WT_dpa3_reshape_miRNA=reshape(Elac_vs_WT_dpa3_df_miRNA, idvar = "miRNA", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape_miRNA[is.na(Elac_vs_WT_dpa3_reshape_miRNA)] <- 0
Elac_vs_WT_dpa5_reshape_miRNA=reshape(Elac_vs_WT_dpa5_df_miRNA, idvar = "miRNA", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape_miRNA[is.na(Elac_vs_WT_dpa5_reshape_miRNA)] <- 0


#display.brewer.all()
#cols <- brewer.pal(3, "BuGn")
par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
pyramid.plot(Elac_vs_WT_dpa3_reshape_miRNA$Frequency.loss, Elac_vs_WT_dpa3_reshape_miRNA$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape_miRNA$miRNA,lxcol="maroon4", rxcol="olivedrab",unit = "Frequency",
             gap=2, space=0.15, top.labels = c("Downregulated", "miRNA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_miRNA$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_miRNA$Frequency.accumulation)),main="Differentially accumulated miRNA at 3 dpa")+theme(
               plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))

pyramid.plot(Elac_vs_WT_dpa5_reshape_miRNA$Frequency.loss, Elac_vs_WT_dpa5_reshape_miRNA$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape_miRNA$miRNA,lxcol="maroon4", rxcol="olivedrab",unit = "Frequency",
             gap=1, space=0.15, top.labels = c("Downregulated", "miRNA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_miRNA$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_miRNA$Frequency.accumulation)),main="Differentially accumulated miRNA at 5 dpa")





################################################################################
#---------------------------------PLOT ALL TYPE WITH LFC------------------------
#load all RData
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad.RData")
DGE_set_all
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad_miRNA_corrected.RData")
DGE_set_with_len_mirna_corrected
load("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_trna_ELAC_filtered_mad_all_types.RData")
DGE_set_with_len_trna_with_type
?read.xlsx
rRNA_set=read.xlsx("F:/PARN_ELAC_silencing/smallRNA/DGE_set_with_len_rRNA_ELAC.xlsx",1)
#save(rRNA_set,file="F:/PARN_ELAC_silencing/smallRNA/DGE_set_just_ELAC_filtered_mad_rRNA_corrected.RData")


length(unique(DGE_set_with_len_trna_with_type$Sequence))

DGE_set_with_len_trna_with_type_tiRNA=DGE_set_with_len_trna_with_type[DGE_set_with_len_trna_with_type$tRF_type_AA_anticodon=="5'-tiRNA-Gly-GCC",]
length(unique(DGE_set_with_len_trna_with_type_tiRNA$Sequence))


data(diff_express)







#--------------------------PLOT_FOR_tRNA_AAA____________________________________
trna_df_AAA=as.data.frame(table(DGE_set_with_len_trna$set,DGE_set_with_len_trna$`tRNA-AAA`,DGE_set_with_len_trna$expression))
colnames(trna_df_AAA)=c("comparison","tRNA_type","Expression","Frequency")

Elac_vs_WT_dpa3_df_AAA=trna_df_AAA[trna_df_AAA$comparison=="Elac_vs_WT_dpa3",]
Elac_vs_WT_dpa5_df_AAA=trna_df_AAA[trna_df_AAA$comparison=="Elac_vs_WT_dpa5",]


Elac_vs_WT_dpa3_reshape_AAA=reshape(Elac_vs_WT_dpa3_df_AAA, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa3_reshape_AAA[is.na(Elac_vs_WT_dpa3_reshape_AAA)] <- 0
Elac_vs_WT_dpa5_reshape_AAA=reshape(Elac_vs_WT_dpa5_df_AAA, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
Elac_vs_WT_dpa5_reshape_AAA[is.na(Elac_vs_WT_dpa5_reshape_AAA)] <- 0

par( mfrow= c(1,2), mai = c(1, 0.1, 0.1, 0.1) )
pyramid.plot(Elac_vs_WT_dpa3_reshape_AAA$Frequency.loss, Elac_vs_WT_dpa3_reshape_AAA$Frequency.accumulation,labels= Elac_vs_WT_dpa3_reshape_AAA$tRNA_type,lxcol="blue", rxcol="yellow",unit = "Frequency",
             gap=14.5, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_AAA$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_AAA$Frequency.accumulation)),main="Differentially accumulated tRF at 3 dpa")+theme(
               plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))

pyramid.plot(Elac_vs_WT_dpa5_reshape_AAA$Frequency.loss, Elac_vs_WT_dpa5_reshape_AAA$Frequency.accumulation,labels= Elac_vs_WT_dpa5_reshape_AAA$tRNA_type,lxcol="blue", rxcol="yellow",unit = "Frequency",
             gap=7, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_AAA$Frequency.loss)), 
             raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape_AAA$Frequency.accumulation)),main="Differentially accumulated tRF at 5 dpa")
#################################################################################

# pyramid.plot(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated, PARN_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa3_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=0.25, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa3")
# pyramid.plot(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated, PARN_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa5_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=3.75, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa5")
# 
# pyramid.plot(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated, GFP_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa3_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=0.25, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa3")
# pyramid.plot(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated, GFP_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa5_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=0.5, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa5")
# 
# Elac_vs_WT_dpa3_reshape_selected=subset(Elac_vs_WT_dpa3_reshape,Elac_vs_WT_dpa3_reshape$tRNA_type=="tRNA-Gly")
# pyramid.plot(Elac_vs_WT_dpa3_reshape_selected$Frequency.Downregulated, Elac_vs_WT_dpa3_reshape_selected$Frequency.Upregulated,labels= "tRNA-Gly",lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=3.65, space=0.15, top.labels = c("Downregulated", "tRF-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_selected$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape_selected$Frequency.Upregulated)),main="Differentially accumulated tRF in Elac_vs_WT_dpa3")+theme(
#                plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size=20))
# #just for ELAC (downregulated)

# pyramid.plot(PARN_dpa3_vs_dpa5_reshape$Frequency.Downregulated, PARN_dpa3_vs_dpa5_reshape$Frequency.Upregulated,labels= PARN_dpa3_vs_dpa5_reshape$tRNA_type,lxcol="blueviolet", rxcol="hotpink",unit = "Frequency",
#              gap=0.35, space=0.15, top.labels = c("Downregulated", "tRF type","Upregulated"),laxlab=seq(0,max(PARN_dpa3_vs_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_dpa3_vs_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_dpa3_vs_dpa5")

################################################################################
#plot for tRNA_cleavage
# DGE_set_with_len_trna_structure
# DGE_set_with_len_trna_structure$tRNAScan_type=ifelse(DGE_set_with_len_trna_structure$tRNAScan_type=="Ser/Supressor","Ser",DGE_set_with_len_trna_structure$tRNAScan_type)
# DGE_set_with_len_trna_structure$tRNA_type_sub=paste("tRNA-",DGE_set_with_len_trna_structure$tRNAScan_type,sep="")
# DGE_set_with_len_trna_structure$tRNA_type_codon=paste(DGE_set_with_len_trna_structure$tRNA_type_sub,DGE_set_with_len_trna_structure$Anticodon, sep="-")
# trna_df=as.data.frame(table(DGE_set_with_len_trna$set,DGE_set_with_len_trna$tRNA_subtype,DGE_set_with_len_trna$expression))
# # trna_df_codon=as.data.frame(table(DGE_set_with_len_trna_structure$set,DGE_set_with_len_trna_structure$tRNA_type_codon,
# #                                   DGE_set_with_len_trna_structure$tRNA_subtype,DGE_set_with_len_trna_structure$expression))
# colnames(trna_df)=c("comparison","tRNA_type","Expression","Frequency")
# trna_df=subset(trna_df,trna_df$tRNA_type!="not_known")
# colnames(trna_df_codon)=c("comparison","tRNA_type","Expression","Frequency")
# #trna_df=subset(trna_df,trna_df$Frequency>0)
# unique(trna_df$comparison)
# #Levels: Elac_vs_WT_dpa3 Elac_vs_WT_dpa5 GFP_vs_WT_dpa3 GFP_vs_WT_dpa5 PARN_vs_WT_dpa3 PARN_vs_WT_dpa5
# Elac_vs_WT_dpa3_df=trna_df[trna_df$comparison=="Elac_vs_WT_dpa3",]
# Elac_vs_WT_dpa5_df=trna_df[trna_df$comparison=="Elac_vs_WT_dpa5",]
# GFP_vs_WT_dpa3_df=trna_df[trna_df$comparison=="GFP_vs_WT_dpa3",]
# GFP_vs_WT_dpa5_df=trna_df[trna_df$comparison=="GFP_vs_WT_dpa5",]
# PARN_vs_WT_dpa3_df=trna_df[trna_df$comparison=="PARN_vs_WT_dpa3",]
# PARN_vs_WT_dpa5_df=trna_df[trna_df$comparison=="PARN_vs_WT_dpa5",]
# 
# Elac_vs_WT_dpa3_reshape=reshape(Elac_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_WT_dpa3_reshape[is.na(Elac_vs_WT_dpa3_reshape)] <- 0
# Elac_vs_WT_dpa5_reshape=reshape(Elac_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_WT_dpa5_reshape[is.na(Elac_vs_WT_dpa5_reshape)] <- 0
# GFP_vs_WT_dpa3_reshape=reshape(GFP_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa3_reshape[is.na(GFP_vs_WT_dpa3_reshape)] <- 0
# GFP_vs_WT_dpa5_reshape=reshape(GFP_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa5_reshape[is.na(GFP_vs_WT_dpa5_reshape)] <- 0
# PARN_vs_WT_dpa3_reshape=reshape(PARN_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa3_reshape[is.na(PARN_vs_WT_dpa3_reshape)] <- 0
# PARN_vs_WT_dpa5_reshape=reshape(PARN_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa5_reshape[is.na(PARN_vs_WT_dpa5_reshape)] <- 0
# 
# # elac_test=trna_df[trna_df$comparison=="WT_vs_Elac_dpa3",]
# # par(mar=pyramid.plot(elac_test$Frequency[elac_test$Expression=="Upregulated",],elac_test$Frequency[elac_test$Expression=="Upregulated",],labels=agelabels,
# #                      main="Australian population pyramid 2002",lxcol=mcol,rxcol=fcol,
# #                      gap=0.5,show.values=TRUE))
# 
# 
# 
# par( mfrow= c(3,2), mai = c(1, 0.1, 0.1, 0.1) )
# pyramid.plot(Elac_vs_WT_dpa3_reshape$Frequency.Downregulated, Elac_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= Elac_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=4, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in Elac_vs_WT_dpa3")
# 
# pyramid.plot(Elac_vs_WT_dpa5_reshape$Frequency.Downregulated, Elac_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= Elac_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=1.25, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in Elac_vs_WT_dpa5")
# 
# pyramid.plot(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated, PARN_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=0.35, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa3")
# pyramid.plot(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated, PARN_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=6.5, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa5")
# pyramid.plot(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated, GFP_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=0.45, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa3")
# pyramid.plot(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated, GFP_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=0.55, space=0.15, top.labels = c("Downregulated", "tRNA type","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa5")
# 
# ?pyramid.plot
# ################################################################################
# #--------------------------------plots for ELAC2 poster-------------------------------
# ################################################################################
# DGE_set_with_len_trna_structure
# new_trna_df=DGE_set_with_len_trna_structure
# new_trna_df$tRNA_whole=paste(new_trna_df$tRNA_subtype,new_trna_df$tRNAScan_type,new_trna_df$Anticodon,sep="-")
# trna_df=as.data.frame(table(new_trna_df$set,new_trna_df$tRNA_whole,new_trna_df$expression))
# colnames(trna_df)=c("comparison","tRNA","expression","freq")
# dpa3_sets=c("")
# down_dpa3=subset(trna_df,(trna_df$comparison=="Elac_vs_WT_dpa3" | trna_df$comparison=="GFP_vs_WT_dpa3") & trna_df$expression=="Downregulated" & trna_df$tRNA=="5'-tiRNA-Gly-GCC")
# trna_df_notAA=as.data.frame(table(new_trna_df$set,new_trna_df$tRNA_subtype,new_trna_df$expression))
# colnames(trna_df_notAA)=c("comparison","tRNA","expression","freq")
# up_dpa3=subset(trna_df_notAA,(trna_df_notAA$comparison=="Elac_vs_WT_dpa3" | trna_df_notAA$comparison=="GFP_vs_WT_dpa3") 
#                & trna_df_notAA$expression=="Upregulated" & (trna_df_notAA$tRNA=="5'-tRF" | trna_df_notAA$tRNA=="3'-tRF"))
# up_dpa5=subset(trna_df_notAA,(trna_df_notAA$comparison=="Elac_vs_WT_dpa5" | trna_df_notAA$comparison=="GFP_vs_WT_dpa5") 
#                & trna_df_notAA$expression=="Upregulated" & (trna_df_notAA$tRNA=="5'-tRF" | trna_df_notAA$tRNA=="3'-tRF"))
# #convert to wide
# down_dpa3_reshape=reshape(down_dpa3, idvar = "tRNA", timevar = "comparison", direction = "wide")
# up_dpa3_reshape=reshape(up_dpa3, idvar = "tRNA", timevar = "comparison", direction = "wide")
# up_dpa5_reshape=reshape(up_dpa5, idvar = "tRNA", timevar = "comparison", direction = "wide")
# par( mfrow= c(3,1), mai = c(1, 0.1, 0.1, 0.1) )
# 
# 
# pyramid.plot(up_dpa3_reshape$freq.Elac_vs_WT_dpa3, up_dpa3_reshape$freq.GFP_vs_WT_dpa3,labels=up_dpa3_reshape$tRNA,lxcol="#67A9CF", rxcol="#EF8A62",unit = "number of RNA species",
#              gap=1.5, space=0.15, top.labels = c("ELAC2", "tRNA type","GFP"),laxlab=seq(0,max(up_dpa3_reshape$freq.Elac_vs_WT_dpa3)),
#              raxlab=seq(0,max(up_dpa3_reshape$freq.GFP_vs_WT_dpa3)),main="Accumulation of tRF in dpa3")
# 
# pyramid.plot(up_dpa5_reshape$freq.Elac_vs_WT_dpa5, up_dpa5_reshape$freq.GFP_vs_WT_dpa5,labels=up_dpa5_reshape$tRNA,lxcol="#67A9CF", rxcol="#EF8A62",unit = "number of RNA species",
#              gap=1.5, space=0.15, top.labels = c("ELAC2", "tRNA type","GFP"),laxlab=seq(0,max(up_dpa5_reshape$freq.Elac_vs_WT_dpa5)),
#              raxlab=seq(0,max(up_dpa5_reshape$freq.GFP_vs_WT_dpa5)),main="Accumulation of tRF in dpa5")
# pyramid.plot(down_dpa3_reshape$freq.Elac_vs_WT_dpa3, down_dpa3_reshape$freq.GFP_vs_WT_dpa3,labels="5'-tiRNA-Gly-GCC",lxcol="#67A9CF", rxcol="#EF8A62",unit = "number of RNA species",
#              gap=1.5, space=0.15, top.labels = c("ELAC2", "tRNA type","GFP"),laxlab=seq(0,max(down_dpa3_reshape$freq.Elac_vs_WT_dpa3)),
#              raxlab=seq(0,max(down_dpa3_reshape$freq.GFP_vs_WT_dpa3)),main="Loss of tRF in dpa3")
# 
# up_dpa3_reshape
# 
# ggplot(up_dpa3, aes(x = tRNA, y = freq, fill = comparison)) +
#   geom_bar(stat = "identity") + 
#   facet_share(~comparison, dir = "h", scales = "free", reverse_num = TRUE) +
#   coord_flip() +
#   theme_minimal()
# 
# # ggplot(nigeria, aes(x = Age, y = Population, fill = Gender)) + 
# #   geom_bar(subset = .(Gender == "Female"), stat = "identity") + 
# #   geom_bar(subset = .(Gender == "Male"), stat = "identity") + 
# #   scale_y_continuous(breaks = seq(-15000000, 15000000, 5000000), 
# #                      labels = paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")) + 
# #   coord_flip() + 
# #   scale_fill_brewer(palette = "Set1") + 
# #   theme_bw()
# 
# 
# 
# 
# down_dpa3$dpa="3dpa"
# up_dpa3$dpa="3dpa"
# up_dpa5$dpa="5dpa"
# poster_elac_trna=rbind(down_dpa3,up_dpa3,up_dpa5)
# poster_elac_trna$comparison=as.character(poster_elac_trna$comparison)
# poster_elac_trna$comparison=substr(poster_elac_trna$comparison,1,nchar(poster_elac_trna$comparison)-5)
# 
# ggplot(poster_elac_trna, aes(x = tRNA, y = freq, fill = comparison)) + 
#   geom_bar(subset = .(comparison == "Elac_vs_WT"), stat = "identity") + 
#   geom_bar(subset = .(comparison == "GFP_vs_WT"), stat = "identity") + 
#   coord_flip() + 
#   scale_fill_brewer(palette = "Set1") + 
#   facet_grid(.~ dpa)+
#   theme_bw()
# #facet_wrap(.~ Comparison)
# install.packages("apyramid")
# library(apyramid)
# ?pyramid.plot
# 
# 
# 
# apyramid::age_pyramid(
#   data = poster_elac_trna,
#   age_group = "tRNA",
#   split_by = "comparison",              # show percents, not counts
#   show_midpoint = FALSE, 
#   proportional = FALSE
#   # remove bar mid-point line
#   #pal = c("orange", "purple")      # can specify alt. colors here (but not labels)
# )+                 
#   
#   # additional ggplot commands
#   theme_minimal()+                               # simplfy background
#   scale_fill_manual(                             # specify colors AND labels
#     values = c("orange", "purple"))+
#   labs(y = "Percent of all cases",              # note x and y labs are switched
#        x = "Age categories",                          
#        fill = "Gender", 
#        caption = "My data source and caption here",
#        title = "Title of my plot",
#        subtitle = "Subtitle with \n a second line...")+
#   theme(
#     legend.position = "bottom",                          # legend to bottom
#     axis.text = element_text(size = 10, face = "bold"),  # fonts/sizes
#     axis.title = element_text(size = 12, face = "bold"))
# 
# install.packages("conmat", repos = "https://njtierney.r-universe.dev")
# library(conmat)
# library(tidyverse)
# abs_age_lga("Hobart (C)")
# two_abs_age_lga <- function(lga_1, lga_2){
#   bind_rows(
#     abs_age_lga(lga_1),
#     abs_age_lga(lga_2)
#   )
# }
# melb_syd <- two_abs_age_lga("Melbourne (C)", "Sydney (C)")
# melb_syd
# tail(melb_syd)
# 
# melb_syd_pyramid <- melb_syd %>% 
#   mutate(
#     population = case_when(
#       lga == "Sydney (C)" ~ -population,
#       TRUE ~ population
#     ),
#     lower.age.limit = as_factor(lower.age.limit)
#   )
# 
# 
# ggplot(melb_syd_pyramid,
#        aes(x = population,
#            y = lower.age.limit,
#            fill = lga)) +
#   geom_col() 
# 
poster_elac_trna$tRNA=as.factor(poster_elac_trna$tRNA)
# ggplot(poster_elac_trna,
#        aes(x = freq,
#            y = tRNA,
#            fill = comparison)) +
#   geom_col() 
poster_elac_trna_pyramid <- poster_elac_trna %>%
  mutate(
    freq = case_when(
      comparison == "GFP_vs_WT" ~ -freq,
      TRUE ~ freq
    ),
    tRNA = as_factor(tRNA)
  )
# ggplot(poster_elac_trna_pyramid,
#        aes(x = freq,
#            y = tRNA,
#            fill = comparison)) +
#   geom_col()

pop_range <- range(poster_elac_trna_pyramid$freq)
pop_range_seq <- seq(pop_range[1], pop_range[2], by = 4)

# pop_range_breaks <- pretty(pop_range, n = 7)
poster_elac_trna_pyramid$expression=ifelse(poster_elac_trna_pyramid$expression=="Upregulated","accumulation","loss")
poster_elac_trna_pyramid$expression <- factor(poster_elac_trna_pyramid$expression,
                                              levels = c("loss","accumulation"), ordered = TRUE)
poster_elac_trna_pyramid$comparison=ifelse(poster_elac_trna_pyramid$comparison=="Elac_vs_WT","ELAC2","GFP")

# poster_elac_trna_pyramid
# ggplot(poster_elac_trna_pyramid,
#        aes(x = freq,
#            y = tRNA,
#            fill = comparison)) +
#   geom_col() +
#   scale_x_continuous(breaks  = pop_range_seq,
#                      labels = abs(pop_range_seq))+
#   ggtitle("tRNA-derived fragments accumulation")+xlab("number of RNA species")+ylab("tRF type")  +
#   theme(legend.position = "top",plot.title = element_text(size=20, face="bold",hjust = 0.5),
#         strip.text.x = element_text(size = 30)) +facet_grid(.~ dpa+expression,scales = "free")+theme_bw()+scale_fill_brewer(palette = "Set2")+
#   theme(strip.text.x = element_text(size = 15),plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size = 15))

ggplot(poster_elac_trna_pyramid,
       aes(x = freq,
           y = tRNA,
           fill = comparison)) +
  geom_col() +
  scale_x_continuous(breaks  = pop_range_seq,
                     labels = abs(pop_range_seq))+
  ggtitle("tRNA-derived fragments accumulation")+xlab("number of RNA species")+ylab("tRF type")  +
  theme(legend.position = "top",plot.title = element_text(size=20, face="bold",hjust = 0.5),
        strip.text.x = element_text(size = 30)) +facet_grid(.~ dpa+expression)+theme_bw()+scale_fill_brewer(palette = "Set2")+
  theme(strip.text.x = element_text(size = 15),plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size = 15))+
  guides(fill = guide_legend(title = "gene"))


ggplot(poster_elac_trna_pyramid,
       aes(x = freq,
           y = tRNA,
           fill = comparison)) +
  geom_col() +
  scale_x_continuous(breaks  = pop_range_seq,
                     labels = abs(pop_range_seq))+
  ggtitle("tRNA-derived fragments accumulation")+xlab("number of RNA species")+ylab("tRF type")  +
  theme(legend.position = "top",plot.title = element_text(size=20, face="bold",hjust = 0.5),
        strip.text.x = element_text(size = 30)) +facet_grid(.~ dpa+expression)+theme_bw()+scale_fill_brewer(palette = "Set2")+
  theme(strip.text.x = element_text(size = 15),plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size = 15))+
  guides(fill = guide_legend(title = "gene"))




#alternative
poster_elac_trna=rbind(down_dpa3,up_dpa3,up_dpa5)
poster_elac_trna$comparison=as.character(poster_elac_trna$comparison)
poster_elac_trna$comparison=substr(poster_elac_trna$comparison,1,nchar(poster_elac_trna$comparison)-5)
poster_elac_trna$tRNA=as.factor(poster_elac_trna$tRNA)

poster_elac_trna_pyramid <- poster_elac_trna %>%
  mutate(
    freq = case_when(
      expression == "Downregulated" ~ -freq,
      TRUE ~ freq
    ),
    tRNA = as_factor(tRNA)
  )


pop_range <- range(poster_elac_trna_pyramid$freq)
pop_range_seq <- seq(pop_range[1], pop_range[2], by = 2)
pop_range_seq <- seq(-30, 30, by = 5)
# pop_range_breaks <- pretty(pop_range, n = 7)
poster_elac_trna_pyramid$expression=ifelse(poster_elac_trna_pyramid$expression=="Upregulated","accumulation","loss")
poster_elac_trna_pyramid$expression <- factor(poster_elac_trna_pyramid$expression,
                                              levels = c("loss","accumulation"), ordered = TRUE)
poster_elac_trna_pyramid$comparison=ifelse(poster_elac_trna_pyramid$comparison=="Elac_vs_WT","ELAC2","GFP")

#change scale or each plot
count <- 0
breaks_fun <- function(x) {
  count <<- count + 1L
  switch(
    count,
    pop_range_seq,
    seq(1,5,by=1),
    seq(1,12,by=2)
  )
}



my_breaks <- function(x) { if (max(x) < 21) seq(0, 12, 2) else pop_range_seq }



pop_range_seq
ggplot(poster_elac_trna_pyramid[apply(poster_elac_trna_pyramid!=0, 1, all),],
       aes(x = freq,
           y = tRNA,
           fill = expression)) +
  geom_col() +
  scale_x_continuous(breaks  = my_breaks)+
  ggtitle("tRNA-derived fragments accumulation")+xlab("number of RNA species")+ylab("tRF type")  +
  theme(legend.position = "top",plot.title = element_text(size=20, face="bold",hjust = 0.5),
        strip.text.x = element_text(size = 30)) +facet_grid(.~ dpa+comparison,scale="free",space="free_x")+theme_bw()+scale_fill_manual(values=c("red4","darkgreen"))+
  theme(strip.text.x = element_text(size = 15),plot.title = element_text(size=20, face="bold",hjust = 0.5),text = element_text(size = 15))+
  guides(fill = guide_legend(title = ""))
# + facetted_pos_scales(y = list(
#       dpa == "3dpa" &  comparison=="GFP" ~ scale_x_continuous(breaks = seq(1,5,by=1)),
#       dpa == "3dpa" &  comparison=="GFP" ~ scale_x_continuous(breaks = seq(1,12,by=2))
#     )
#   )
?scale_fill_manual
poster_elac_trna_pyramid[apply(poster_elac_trna_pyramid!=0, 1, all),]



count <- 0
breaks_fun <- function(x) {
  count <<- count + 1L
  switch(
    count,
    c(1, 3, 5, 7, 9),
    c(45, 55),
    c(0, 50, 100),
    seq(0, 8, 0.2)
  )
}

install.packages("ggh4x")
library(ggh4x)
ggplot(df) + 
  geom_point(aes(x, y)) + 
  facet_wrap(~id, scales = 'free_x') + 
  theme_bw() + 
  scale_x_continuous(breaks = breaks_fun, limits = c(0, NA)) + 
  labs(
    title = "Custom x-axis breaks for each facet", 
    subtitle = "Based upon plot index"
  )

plot <- ggplot(iris, aes(Sepal.Width, Sepal.Length)) +
  geom_point(aes(colour = Species)) +
  facet_wrap(Species ~ ., scales = "free_y")

# Reversing the y-axis in the second panel. When providing a list of scales,
# NULL indicates to use the default, global scale
plot +
  facetted_pos_scales(
    y = list(NULL, scale_y_continuous(trans = "reverse"))
  )

# Alternative for specifying scales with formula lists. The LHS can access
# columns in the plot's layout.
plot +
  facetted_pos_scales(
    y = list(
      Species == "virginica" ~ scale_y_continuous(breaks = c(6, 7)),
      Species == "versicolor" ~ scale_y_reverse()
    )
  )

?facetted_pos_scales
# load ggplot2
library("ggplot2")

# Data from the facet plot
x1 <- rnorm(100)
x2 <- rnorm(100)+x1
grp <- rbinom(100, 1, 0.1)
x1[grp == 1] <- x1[grp == 1] * 5
x2[grp == 1] <- x2[grp == 1] * 5

# Data from the facet plot
gfg <- data.frame(x1, x2, grp)

# facet plot with facet_wrap
gfg_plot <- ggplot(gfg, aes(x1, x2)) +
  geom_point() + facet_wrap(~ grp)

# Draw plot with free x-axis scales
gfg_plot + facet_wrap(~ grp, scales="free_x")









# 
# guides(fill = guide_legend(title = "RNA type"))+ theme(panel.spacing = unit(2, "lines"))+
#   theme(strip.text.x = element_text(size = 22))
# 
# scale_fill_brewer(palette = "Pastel1")+theme_void()+theme(
#   plot.title = element_text(size=20, face="bold",hjust = 0.5),
#   text = element_text(size=20))+labs(fill="RNA species")




########------------anticodon
# Elac_vs_WT_dpa3_df=trna_df_codon[trna_df_codon$comparison=="Elac_vs_WT_dpa3",]
# Elac_vs_WT_dpa5_df=trna_df_codon[trna_df_codon$comparison=="Elac_vs_WT_dpa5",]
# GFP_vs_WT_dpa3_df=trna_df_codon[trna_df_codon$comparison=="GFP_vs_WT_dpa3",]
# GFP_vs_WT_dpa5_df=trna_df_codon[trna_df_codon$comparison=="GFP_vs_WT_dpa5",]
# PARN_vs_WT_dpa3_df=trna_df_codon[trna_df_codon$comparison=="PARN_vs_WT_dpa3",]
# PARN_vs_WT_dpa5_df=trna_df_codon[trna_df_codon$comparison=="PARN_vs_WT_dpa5",]
# 
# Elac_vs_WT_dpa3_reshape=reshape(Elac_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_WT_dpa3_reshape[is.na(Elac_vs_WT_dpa3_reshape)] <- 0
# Elac_vs_WT_dpa5_reshape=reshape(Elac_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# Elac_vs_WT_dpa5_reshape[is.na(Elac_vs_WT_dpa5_reshape)] <- 0
# GFP_vs_WT_dpa3_reshape=reshape(GFP_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa3_reshape[is.na(GFP_vs_WT_dpa3_reshape)] <- 0
# GFP_vs_WT_dpa5_reshape=reshape(GFP_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# GFP_vs_WT_dpa5_reshape[is.na(GFP_vs_WT_dpa5_reshape)] <- 0
# PARN_vs_WT_dpa3_reshape=reshape(PARN_vs_WT_dpa3_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa3_reshape[is.na(PARN_vs_WT_dpa3_reshape)] <- 0
# PARN_vs_WT_dpa5_reshape=reshape(PARN_vs_WT_dpa5_df, idvar = "tRNA_type", timevar = "Expression", direction = "wide")
# PARN_vs_WT_dpa5_reshape[is.na(PARN_vs_WT_dpa5_reshape)] <- 0
# 
# # elac_test=trna_df[trna_df$comparison=="WT_vs_Elac_dpa3",]
# # par(mar=pyramid.plot(elac_test$Frequency[elac_test$Expression=="Upregulated",],elac_test$Frequency[elac_test$Expression=="Upregulated",],labels=agelabels,
# #                      main="Australian population pyramid 2002",lxcol=mcol,rxcol=fcol,
# #                      gap=0.5,show.values=TRUE))
# 
# 
# 
# par( mfrow= c(2,2), mai = c(1, 0.1, 0.1, 0.1) )
# pyramid.plot(Elac_vs_WT_dpa3_reshape$Frequency.Downregulated, Elac_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= Elac_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=5, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in Elac_vs_WT_dpa3")
# 
# pyramid.plot(Elac_vs_WT_dpa5_reshape$Frequency.Downregulated, Elac_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= Elac_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=0.7, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(Elac_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in Elac_vs_WT_dpa5")
# # pyramid.plot(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated, GFP_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
# #              gap=0.5, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Downregulated)),
# #              raxlab=seq(0,max(GFP_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa3")
# # pyramid.plot(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated, GFP_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= GFP_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
# #              gap=0.65, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
# #              raxlab=seq(0,max(GFP_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in GFP_vs_WT_dpa5")
# pyramid.plot(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated, PARN_vs_WT_dpa3_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa3_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=0.35, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Downregulated)),
#              raxlab=seq(0,max(PARN_vs_WT_dpa3_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa3")
# pyramid.plot(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated, PARN_vs_WT_dpa5_reshape$Frequency.Upregulated,labels= PARN_vs_WT_dpa5_reshape$tRNA_type,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Frequency",
#              gap=5, space=0.15, top.labels = c("Downregulated", "tRNA-AA","Upregulated"),laxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Downregulated)), 
#              raxlab=seq(0,max(PARN_vs_WT_dpa5_reshape$Frequency.Upregulated)),main="Differentially accumulated tRF in PARN_vs_WT_dpa5")
################################################################################
#-----------------------------miRNA---------------------------------------------
library(ggplot2)
miRNA_set=subset(DGE_set,DGE_set$RNA_type=="miRNA")
miRNA_set_table=as.data.frame(table(miRNA_set$set,miRNA_set$expression))
colnames(miRNA_set_table)=c("Comparison","Expression","Frequency")
miRNA_set_table_sum=aggregate(Frequency~Comparison, miRNA_set_table,sum)
colnames(miRNA_set_table_sum)=c("Comparison","sum")
miRNA_set_table=merge(miRNA_set_table,miRNA_set_table_sum,by="Comparison",all.x = TRUE)
miRNA_set_table$perc=(miRNA_set_table$Frequency/miRNA_set_table$sum)*100
bp=ggplot(data=miRNA_set_table, aes(x=" ", y=perc, group=Expression, colour=Expression, fill=Expression)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_wrap(.~ Comparison) +theme_void()+ggtitle("Differentially accumulated miRNA")
bp + scale_fill_manual(values=c("#E69F00", "#56B4E9"))+ theme(
  plot.title = element_text(size=14, face="bold",hjust = 0.5),
  axis.text=element_text(size=13),
  text = element_text(size=15)
)
?theme
Elac_vs_WT_dpa3_pie=ggplot(data=miRNA_set_table[miRNA_set_table$Comparison=="Elac_vs_WT_dpa3",], aes(x=" ", y=Frequency, group=Expression, colour=Expression, fill=Expression)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)

