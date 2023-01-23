# this script is to make manifest for downloading transcriptome files for each site primary
# made 2023/01/23

# activate package to get TCGA ID
library(stringr)
library(TCGAutils)
library(GenomicDataCommons)
library(magrittr)

# set function for searching row that have same character (this grep function can search multiple characters)
grep2 <- function(x,y,...) unlist(lapply(x,grep,y,...))

# make new direction
setwd("C:/Rdata")
dir.create("20230123_make_TCGA_manifest_for_each_site_primary")


#### make manifest of TCGA RNA-seq transcriptome ####

# convert file id to case id
setwd("c:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
t.manifest <-read.table("TCGA_transcriptome_all_files_gdc_manifest.2023-01-23.txt",sep="\t",header = T,stringsAsFactors = F)
t.uuid <-t.manifest[,1]
st.uuid <-UUIDtoUUID(t.uuid)
colnames(t.uuid)[2] <-"case_id"

# get information about site primary 
t.list <-gdc_clinical(t.uuid[,2])
t.df <-as.data.frame(t.list[["main"]])
site <-sort(unique(t.df[,5]))

# make table about number of files for each site primary
site.file.num <-as.data.frame(table(t.df[,5]),stringsAsFactors = F)
site.file.num[,3] <-NA

setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
dir.create("RNA_seq_transcriptome")
setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/RNA_seq_transcriptome")

# make manifest for each site primary
for (i in 1:length(site)) {
  # extract rows for each site primary 
  site.df <-t.df[t.df[,5]==site[i],]
  
  # extract subject files from manifest contains all files 
  site.case.id <-unique(site.df[,7])
  m <-grep2(site.case.id,t.uuid[,2])
  site.manifest <-t.manifest[m,]
  
  # edit name of site primary and create file name
  site.name <-gsub(" ","_",site[i])
  site.name <-gsub(",","",site.name)
  file.name <-paste0("TCGA_",site.name,"_transcriptome_manifest",".txt")
  
  # output manifest
  write.table(site.manifest,file.name,sep="\t",row.names = F,quote = F)
  site.file.num[i,3] <-length(m)
  }

# output table about number of files for each site primary
colnames(site.file.num) <-c("site_primary","number_of_case","number_of_file")
write.table(site.file.num,"table_about_number_of_TCGA_transcriptome_file_for_each_site_primary.txt",sep="\t",row.names = F,quote = F)

#### make manifest of TCGA miRNA-seq ####

# convert file id to case id
setwd("c:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
m.manifest <-read.table("TCGA_miRNA_all_files_gdc_manifest.2023-01-23.txt",sep="\t",header = T,stringsAsFactors = F)
m.uuid <-m.manifest[,1]
m.uuid <-UUIDtoUUID(m.uuid)
colnames(m.uuid)[2] <-"case_id"

# get information about site primary 
m.list <-gdc_clinical(m.uuid[,2])
m.df <-as.data.frame(m.list[["main"]])
site <-sort(unique(m.df[,5]))

# make table about number of files for each site primary
site.file.num <-as.data.frame(table(m.df[,5]),stringsAsFactors = F)
site.file.num[,3] <-NA

setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
dir.create("miRNA_seq")
setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/miRNA_seq")

# make manifest for each site primary
for (i in 1:length(site)) {
  # extract rows for each site primary 
  site.df <-m.df[m.df[,5]==site[i],]
  
  # extract subject files from manifest contains all files 
  site.case.id <-unique(site.df[,7])
  m <-grep2(site.case.id,m.uuid[,2])
  site.manifest <-m.manifest[m,]
  
  # edit name of site primary and create file name
  site.name <-gsub(" ","_",site[i])
  site.name <-gsub(",","",site.name)
  file.name <-paste0("TCGA_",site.name,"_miRNA_manifest",".txt")
  
  # output manifest
  write.table(site.manifest,file.name,sep="\t",row.names = F,quote = F)
  site.file.num[i,3] <-length(m)
}

# output table about number of files for each site primary
colnames(site.file.num) <-c("site_primary","number_of_case","number_of_file")
write.table(site.file.num,"table_about_number_of_TCGA_miRNAseq_file_for_each_site_primary.txt",sep="\t",row.names = F,quote = F)

#### make manifest of TCGA gene counts ####

# convert file id to case id
setwd("c:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
g.manifest <-read.table("TCGA_gene_counts_all_files_gdc_manifest.2023-01-23.txt",sep="\t",header = T,stringsAsFactors = F)
g.uuid <-g.manifest[,1]
g.uuid <-UUIDtoUUID(g.uuid)
colnames(g.uuid)[2] <-"case_id"

# get information about site primary 
g.list <-gdc_clinical(g.uuid[,2])
g.df <-as.data.frame(g.list[["main"]])
site <-sort(unique(g.df[,5]))

# make table about number of files for each site primary
site.file.num <-as.data.frame(table(g.df[,5]),stringsAsFactors = F)
site.file.num[,3] <-NA

setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
dir.create("RNA_seq_gene_counts")
setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/RNA_seq_gene_counts")

# make manifest for each site primary
for (i in 1:length(site)) {
  # extract rows for each site primary 
  site.df <-g.df[g.df[,5]==site[i],]
  
  # extract subject files from manifest contains all files 
  site.case.id <-unique(site.df[,7])
  m <-grep2(site.case.id,g.uuid[,2])
  site.manifest <-g.manifest[m,]
  
  # edit name of site primary and create file name
  site.name <-gsub(" ","_",site[i])
  site.name <-gsub(",","",site.name)
  file.name <-paste0("TCGA_",site.name,"_gene_counts_manifest",".txt")
  
  # output manifest
  write.table(site.manifest,file.name,sep="\t",row.names = F,quote = F)
  site.file.num[i,3] <-length(m)
}

# output table about number of files for each site primary
colnames(site.file.num) <-c("site_primary","number_of_case","number_of_file")
write.table(site.file.num,"table_about_number_of_TCGA_gene_counts_file_for_each_site_primary.txt",sep="\t",row.names = F,quote = F)
