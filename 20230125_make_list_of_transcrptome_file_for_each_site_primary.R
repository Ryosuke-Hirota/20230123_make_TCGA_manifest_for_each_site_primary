# this script is to make list for quantifying transcript level by salmon
# made 2023/01/25

# list manifests for each site primary
# these lists are located at "https://github.com/Ryosuke-Hirota/20230123_make_TCGA_manifest_for_each_site_primary/RNA_seq_transcriptome"
t.manifests <-list.files(path = "C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/RNA_seq_transcriptome",pattern = "manifest.txt")

# make new directory
setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary")
dir.create("list_of_transcriptome_file")


# import manifest one by one and make list
for (i in 1:length(t.manifests)) {
  setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/RNA_seq_transcriptome")
  manifest <-read.table(t.manifests[i],sep="\t",header = T,stringsAsFactors = F)
  file.list <- manifest[,c(1,2)]
  setwd("C:/Rdata/20230123_make_TCGA_manifest_for_each_site_primary/list_of_transcriptome_file")
  file.name <-gsub("manifest","process",t.manifests[i])
  write.table(file.list,file.name,sep="\t",row.names = F,col.names = F,quote = F)
}

