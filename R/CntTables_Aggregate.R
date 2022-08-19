## Example
## Stranded, H. sapiens simulation

library(data.table)
library(Rsubread)
library(tidyr)
library(stringr)
library(zeallot)




# Function to load SQuIRE count table

squire_load<-function(mypath){
  s1_EM<-fread(mypath, select=c(4, 7,9,10,11,14,16))
  print(head(s1_EM))
  s1_EM[,SQuIRE:=paste0(TE_chr,paste0("_",paste0(TE_start,paste0("_",paste0(TE_stop,paste0("_",TE_strand))))))]
  s1_EM[,c("TE_ID", "TE_chr", "TE_start", "TE_stop", "TE_strand"):=NULL]
  s1_EM[,Sample:=paste0(str_remove(Sample,"_MM_100Aligned.sortedByCoord.out.bam"),"_SQuIRE")]
  colnames(s1_EM)[2]<-as.character(s1_EM[1,"Sample"])
  s1_EM[,Sample:=NULL]
  s1_EM<-as.data.frame(s1_EM)
  return(s1_EM)
}

squire_load_all<-function(mypath,repmsk){
  temp<-list.files(path=mypath, pattern="*_TEcounts.txt")
  temp<-paste0(mypath,temp)
  myfiles<-lapply(temp,squire_load)
  tel_tab<-myfiles %>% purrr::reduce(full_join, by="SQuIRE")
  tel_tab[is.na(tel_tab)]<-0
  tel_tab<-as.data.table(tel_tab,key="SQuIRE")
  tel_tab<-repmsk[tel_tab]
  tel_tab[,c("V1", "V4", "V5","V7","fam", "class", "SQuIRE", "TrueCNT","T_ID"):=NULL]
  # subfamilies
  temp<-list.files(path=mypath, pattern="*_subFcounts.txt")
  temp<-paste0(mypath,temp)
  myfiles<-lapply(temp,squire_s_load)
  tel_tabs<-myfiles %>% purrr::reduce(full_join, by="ID")
  tel_tabs[is.na(tel_tabs)]<-0
  return(list(tel_tab,tel_tabs))
}

# Subfamily load function
# 1, 3, 7 _subFcounts.txt
squire_s_load<-function(mypath){
  s1_EM<-fread(mypath, select=c(1, 3, 7))
  s1_EM[,Sample:=paste0(str_remove(Sample,"_MM_100Aligned.sortedByCoord.out.bam"),"_SQuIRE")]
  colnames(s1_EM)<-c("Sample","ID",as.character(s1_EM[1,"Sample"]))
  s1_EM[,Sample:=NULL]
  s1_EM<-separate(s1_EM,ID,sep=":",into="ID",extra="drop")
  s1_EM<-unq(as.data.frame(s1_EM))
  return(s1_EM)
}





# Set path with alignment files
path_bam<-"hsapiens_sim/stranded/bam/"
# Path with TE annotation
gtf<-"hsapiens_sim/GRCh38_rmsk_TEinst_TEcount.gtf"


# List all alignment files
fs0<-list.files(path = "hsapiens_sim/stranded/bam/", pattern = ".bam")
fs02<-paste0("hsapiens_sim/stranded/bam/",fs0)



## Making count table for loci specific expression; featureCounts


sim1_FC_MM_F<-featureCounts(fs02, annot.ext=gtf, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="transcript_id", countMultiMappingReads=TRUE, fraction=TRUE, isPairedEnd=TRUE, nthreads=4,strandSpecific=1,requireBothEndsMapped=TRUE)
sim1_FC_MM_F_cnt<-sim1_FC_MM_F$counts

sim1_FC_MM_R<-featureCounts(fs02, annot.ext=gtf, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="transcript_id", countMultiMappingReads=TRUE, fraction=FALSE, isPairedEnd=TRUE, nthreads=4,strandSpecific=1,requireBothEndsMapped=TRUE)
sim1_FC_MM_R_cnt<-sim1_FC_MM_R$counts

sim1_FC_UM<-featureCounts(fs02, annot.ext=gtf, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="transcript_id", countMultiMappingReads=FALSE, isPairedEnd=TRUE, nthreads=4,strandSpecific=1,requireBothEndsMapped=TRUE)
sim1_FC_UM_cnt<-sim1_FC_UM$counts


# Edit colnames; generate subfamily tables, save subfamilies, sort by rowname and merge to single file
# Merge with SalmonTE, TElocal and TEcount counts

# sim1_FC_MM_F_cnt<-read.delim("simulation_hsapiens/Sim01_FeatureCnt_MMF_loci.txt")
# sim1_FC_MM_R_cnt<-read.delim("simulation_hsapiens/Sim01_FeatureCnt_MMR_loci.txt")
# sim1_FC_UM_cnt<-read.delim("simulation_hsapiens/Sim01_FeatureCnt_UM_loci.txt")


# Edit colnames

colnames(sim1_FC_MM_F_cnt)<-str_remove(colnames(sim1_FC_MM_F_cnt),".MM.100Aligned.sortedByCoord.out.bam")
colnames(sim1_FC_MM_F_cnt)<-str_replace(colnames(sim1_FC_MM_F_cnt),"\\.","_")
colnames(sim1_FC_MM_F_cnt)<-paste0(colnames(sim1_FC_MM_F_cnt),"_FC_MM_F")

colnames(sim1_FC_MM_R_cnt)<-str_replace(str_remove(colnames(sim1_FC_MM_R_cnt),".MM.100Aligned.sortedByCoord.out.bam"),"\\.","_")
colnames(sim1_FC_MM_R_cnt)<-paste0(colnames(sim1_FC_MM_R_cnt),"_FC_MM_R")
colnames(sim1_FC_UM_cnt)<-str_replace(str_remove(colnames(sim1_FC_UM_cnt),".MM.100Aligned.sortedByCoord.out.bam"),"\\.","_")
colnames(sim1_FC_UM_cnt)<-paste0(colnames(sim1_FC_UM_cnt),"_FC_UM")
# Write out the output
save(list=c("sim1_FC_MM_F", "sim1_FC_MM_R", "sim1_FC_UM","sim1_FC_UM_cnt","sim1_FC_MM_R_cnt","sim1_FC_MM_F_cnt"), file="hsapiens_sim/stranded/FeatureCnts_Sim01_st.RData")

write.table(sim1_FC_MM_F_cnt,"hsapiens_sim/stranded/Sim01_FeatureCnt_MMF_loci_st.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(sim1_FC_MM_R_cnt,"hsapiens_sim/stranded/Sim01_FeatureCnt_MMR_loci_st.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(sim1_FC_UM_cnt,"hsapiens_sim/stranded/Sim01_FeatureCnt_UM_loci_st.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")

# check if rownames are identical for the three
identical(rownames(sim1_FC_MM_F_cnt),rownames(sim1_FC_MM_R_cnt))
identical(rownames(sim1_FC_MM_F_cnt),rownames(sim1_FC_UM_cnt))

# cbind them
FCnt_loci<-do.call("cbind", list(sim1_FC_MM_F_cnt, sim1_FC_MM_R_cnt, sim1_FC_UM_cnt))
# remove singular df's to save space
remove(list=c("sim1_FC_MM_F_cnt", "sim1_FC_MM_R_cnt", "sim1_FC_UM_cnt"))
# Get subfamilies cnt
FCnt_loci<-as.data.table(FCnt_loci,keep.rownames = TRUE,key="rn")
setkey(repmsk,ID)
FCnt_sub<-repmsk[FCnt_loci]
FCnt_sub[,c("V1", "V4", "V5", "V7","fam","class","SQuIRE", "TrueCNT"):=NULL]
FCnt_sub<-FCnt_sub[,lapply(.SD,sum), by = .(T_ID),.SDcols=colnames(FCnt_sub)[-c(1,2)]]
fwrite(FCnt_sub,"hsapiens_sim/stranded/sim01_FeatureCNT_sub_st.txt",sep="\t",col.names=TRUE)





# SalmonTE, load tables
sim01_SalmonTEraw <- read.csv("/data/natalia/hsapiens_sim/stranded/sim01strand_SalmonTEraw.csv", row.names=1, stringsAsFactors=FALSE)
colnames(sim01_SalmonTEraw)<-paste0(colnames(sim01_SalmonTEraw),paste0("_","SalmonTE"))
sim01_SalmonTEraw<-sim01_SalmonTEraw[,order(colnames(sim01_SalmonTEraw))]
sim01_SalmonTEraw<-sim01_SalmonTEraw[,c(10,1:9)]
colnames(sim01_SalmonTEraw)[1]<-"sample_01_SalmonTE"

sim01_SalmonTEraw<-as.data.table(sim01_SalmonTEraw,keep.rownames = TRUE)
setkey(sim01_SalmonTEraw,rn)
sim01_SalmonTEraw<-separate(sim01_SalmonTEraw, rn, sep="_dup", into=c("subid","ID"),remove=FALSE)

# SalmonTE subfamily counts aggregate
slmonte_sub<-colnames(sim01_SalmonTEraw)[!colnames(sim01_SalmonTEraw) %in% c("ID","subid","rn")]
salmonte_sub_2<-sim01_SalmonTEraw[,lapply(.SD,sum), by = .(subid),.SDcols=slmonte_sub]
fwrite(salmonte_sub_2,"hsapiens_sim/stranded/salmonte_sub_raw.txt",sep="\t",col.names=TRUE)
sim01_SalmonTEraw[,subid:=NULL]
sim01_SalmonTEraw[,ID:=NULL]





# TElocal
# Load cnt files, merge them for samples, separate for genes and TEs.
c(sim01_TElocal_UM_loci,sim01_TElocal_UM_sub,sim01_TElocal_UM_gene)%<-%telocal_load("hsapiens_sim/stranded/TElocal_UM/")

# Change colnames to remove suffix and add suffix for procedure
colnames(sim01_TElocal_UM_loci)<-paste0(str_remove(colnames(sim01_TElocal_UM_loci),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_UM")
colnames(sim01_TElocal_UM_loci)[1]<-"gene.TE"
colnames(sim01_TElocal_UM_sub)<-paste0(str_remove(colnames(sim01_TElocal_UM_sub),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_UM")
colnames(sim01_TElocal_UM_sub)[1]<-"gene.TE"
colnames(sim01_TElocal_UM_gene)<-paste0(str_remove(colnames(sim01_TElocal_UM_gene),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_UM")
colnames(sim01_TElocal_UM_gene)[1]<-"gene.TE"

# Same for EM-generated counts
c(sim01_TElocal_EM_loci,sim01_TElocal_EM_sub,sim01_TElocal_EM_gene)%<-%telocal_load("hsapiens_sim/stranded/TElocal_EM/")
colnames(sim01_TElocal_EM_loci)<-paste0(str_remove(colnames(sim01_TElocal_EM_loci),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_MM")
colnames(sim01_TElocal_EM_loci)[1]<-"gene.TE"
colnames(sim01_TElocal_EM_sub)<-paste0(str_remove(colnames(sim01_TElocal_EM_sub),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_MM")
colnames(sim01_TElocal_EM_sub)[1]<-"gene.TE"
colnames(sim01_TElocal_EM_gene)<-paste0(str_remove(colnames(sim01_TElocal_EM_gene),"MM_100Aligned.sortedByCoord.out.bam"),"TElocal_MM")
colnames(sim01_TElocal_EM_gene)[1]<-"gene.TE"

# Write out tables for TElocal generated counts
fwrite(sim01_TElocal_UM_loci,"hsapiens_sim/stranded/sim01_TElocal_UM_loci.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TElocal_UM_sub,"hsapiens_sim/stranded/sim01_TElocal_UM_sub.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TElocal_UM_gene,"hsapiens_sim/stranded/sim01_TElocal_UM_gene.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TElocal_EM_loci,"hsapiens_sim/stranded/sim01_TElocal_EM_loci.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TElocal_EM_sub,"hsapiens_sim/stranded/sim01_TElocal_EM_sub.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TElocal_EM_gene,"hsapiens_sim/stranded/sim01_TElocal_EM_gene.txt",sep="\t",col.names=TRUE)


# TEcount
# Load, remove tables with class counts ("sub" table)
# TEcount, uniq
c(sim01_TEcount_UM_loci,sim01_TEcount_UM_sub,sim01_TEcount_UM_gene)%<-%telocal_load("hsapiens_sim/stranded/TEcount_UM/")
remove(sim01_TEcount_UM_sub)
colnames(sim01_TEcount_UM_loci)<-paste0(str_remove(colnames(sim01_TEcount_UM_loci),"MM_100Aligned.sortedByCoord.out.bam"),"TEcount_UM")
colnames(sim01_TEcount_UM_loci)[1]<-"gene.TE"
colnames(sim01_TEcount_UM_gene)<-paste0(str_remove(colnames(sim01_TEcount_UM_gene),"MM_100Aligned.sortedByCoord.out.bam"),"TEcount_UM")
colnames(sim01_TEcount_UM_gene)[1]<-"gene.TE"

# TEcount, EM
c(sim01_TEcount_EM_loci,sim01_TEcount_EM_sub,sim01_TEcount_EM_gene)%<-%telocal_load("hsapiens_sim/stranded/TEcount_EM/")
remove(sim01_TEcount_EM_sub)
colnames(sim01_TEcount_EM_loci)<-paste0(str_remove(colnames(sim01_TEcount_EM_loci),"MM_100Aligned.sortedByCoord.out.bam"),"TEcount_MM")
colnames(sim01_TEcount_EM_loci)[1]<-"rn"
colnames(sim01_TEcount_EM_gene)<-paste0(str_remove(colnames(sim01_TEcount_EM_gene),"MM_100Aligned.sortedByCoord.out.bam"),"TEcount_MM")
colnames(sim01_TEcount_EM_gene)[1]<-"gene.TE"

# Save TEcount joint files
fwrite(sim01_TEcount_UM_loci,"hsapiens_sim/stranded/sim01_TEcount_UM_sub.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TEcount_UM_gene,"hsapiens_sim/stranded/sim01_TEcount_UM_gene.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TEcount_EM_loci,"hsapiens_sim/stranded/sim01_TEcount_EM_sub.txt",sep="\t",col.names=TRUE)
fwrite(sim01_TEcount_EM_gene,"hsapiens_sim/stranded/sim01_TEcount_EM_gene.txt",sep="\t",col.names=TRUE)



setkey(repmsk,SQuIRE)
# SQuIRE
# load counts
setkey(repmsk,SQuIRE)
c(sim01_squire_EM_loci,sim01_squire_EM_sub)%<-%squire_load_all("/data/natalia/hsapiens_sim/stranded/SQuIRE/MM/",repmsk)

colnames(sim01_squire_EM_loci)[-1]<-paste0(colnames(sim01_squire_EM_loci)[-1],"_EM")
colnames(sim01_squire_EM_sub)[-1]<-paste0(colnames(sim01_squire_EM_sub)[-1],"_EM")
fwrite(sim01_squire_EM_loci,"/data/natalia/hsapiens_sim/stranded/sim01_SQuIRE_EM_loci.txt",sep="\t",col.names=TRUE)
fwrite(sim01_squire_EM_sub,"/data/natalia/hsapiens_sim/stranded/sim01_SQuIRE_EM_sub.txt",sep="\t",col.names=TRUE)

# Unique
c(sim01_squire_UM_loci,sim01_squire_UM_sub)%<-%squire_load_all("hsapiens_sim/stranded/SQuIRE/UM/",repmsk)

colnames(sim01_squire_UM_loci)[-1]<-paste0(colnames(sim01_squire_UM_loci)[-1],"_UM")
colnames(sim01_squire_UM_sub)[-1]<-paste0(colnames(sim01_squire_UM_sub)[-1],"_UM")
fwrite(sim01_squire_UM_loci,"hsapiens_sim/stranded/sim01_SQuIRE_UM_loci.txt",sep="\t",col.names=TRUE)
fwrite(sim01_squire_UM_sub,"hsapiens_sim/stranded/sim01_SQuIRE_UM_sub.txt",sep="\t",col.names=TRUE)


##########

# SQuIRE tables generated have duplicated values for some elements, that we need to get rid of

sim01_squire_UM_loci<-sim01_squire_UM_loci[complete.cases(sim01_squire_UM_loci), ]
# DT["A",mult="first"]
filter_dt<-function(sim01_squire_UM_loci){
  dp<-unique(sim01_squire_UM_loci$ID[duplicated(sim01_squire_UM_loci$ID)])
  for (i in 1:length(dp)){
    print(dp[i])
    a<-sim01_squire_UM_loci[ID == dp[i], which = TRUE]
    print(a)
    for (b in c(2:ncol(sim01_squire_UM_loci))){
      print(class(b))
      print(sum(unique(sim01_squire_UM_loci[a[1]:a[length(a)],as.integer(b)])))
      sim01_squire_UM_loci[a[1],as.integer(b)]<-sum(unique(sim01_squire_UM_loci[a[1]:a[length(a)],as.integer(b)]))
    }
    sim01_squire_UM_loci<-sim01_squire_UM_loci[-c(a[-1]),]
  }
  return(sim01_squire_UM_loci)
}
sim01_squire_UM_loci<-filter_dt(sim01_squire_UM_loci)


sim01_squire_EM_loci<-sim01_squire_EM_loci[complete.cases(sim01_squire_EM_loci), ]
sim01_squire_EM_loci<-filter_dt(sim01_squire_EM_loci)

# True counts
load("hsapiens_sim/stranded/sim_counts_matrix.rda")
# Remove genic counts, rename TE coordinates to matching IDs, mark as true, create a subfamily table
# Merge all tables, create meta-table

setkey(repmsk,TrueCNT)
c(counts_matrix,counts_matrix_s)%<-%true_load(counts_matrix,repmsk)
fwrite(counts_matrix,"hsapiens_sim/stranded/sim01_TrueCnt_TE_loci.txt",sep="\t",col.names=TRUE)
fwrite(counts_matrix_s,"hsapiens_sim/stranded/sim01_TrueCnt_TE_sub.txt",sep="\t",col.names=TRUE)
counts_matrix[,T_ID:=NULL]
# Merge tables

# FCnt_loci, FCnt_sub,sim01_SalmonTEraw,salmonte_sub_2,sim01_TElocal_UM_loci,sim01_TElocal_UM_sub,
# sim01_TElocal_EM_loci,sim01_TElocal_EM_sub,sim01_TEcount_UM_loci,sim01_TEcount_EM_loci,sim01_squire_EM_loci,
# sim01_squire_EM_sub,sim01_squire_UM_loci,sim01_squire_UM_sub,counts_matrix,counts_matrix_s

# Loci
# FCnt_loci, sim01_SalmonTEraw, sim01_TElocal_UM_loci, sim01_TElocal_EM_loci, sim01_squire_EM_loci, sim01_squire_UM_loci, counts_matrix
setkey(FCnt_loci,rn)
setkey(sim01_SalmonTEraw,rn)
colnames(sim01_TElocal_UM_loci)[1]<-"rn"
setkey(sim01_TElocal_UM_loci,rn)
colnames(sim01_TElocal_EM_loci)[1]<-"rn"
setkey(sim01_TElocal_EM_loci,rn)
colnames(sim01_squire_EM_loci)[1]<-"rn"
setkey(sim01_squire_EM_loci,"rn")
colnames(sim01_squire_UM_loci)[1]<-"rn"
setkey(sim01_squire_UM_loci,rn)
colnames(counts_matrix)[1]<-"rn"
setkey(counts_matrix,rn)


counts_matrix_c<-counts_matrix[sim01_squire_UM_loci[sim01_squire_EM_loci[sim01_TElocal_EM_loci[sim01_TElocal_UM_loci[sim01_SalmonTEraw[FCnt_loci]]]]]]
counts_matrix_c[is.na(counts_matrix_c)]<-0
counts_matrix_c[,T_ID:=NULL]
# Filter out all-zero rows
counts_matrix_c<-counts_matrix_c[rowSums(counts_matrix_c[,-1])>0,]
fwrite(counts_matrix_c,"hsapiens_sim/stranded/sim01_allCnt_loci.txt",sep="\t",col.names=TRUE)


colnames(counts_matrix_c)[12:21]<-paste0(colnames(counts_matrix_c)[12:21],"_UM")

# create meta table
allmeta_loci<-data.frame(ID=colnames(counts_matrix_c)[-1],groups=rep(c(rep(1,5),rep(2,5)),9))
allmeta_loci[allmeta_loci$ID %like% "_TC","Tool"]<-"TrueCount"
allmeta_loci[allmeta_loci$ID %like% "_SQuIRE","Tool"]<-"SQuIRE_UM"
allmeta_loci[allmeta_loci$ID %like% "_SQuIRE_EM","Tool"]<-"SQuIRE_EM"
allmeta_loci[allmeta_loci$ID %like% "_TElocal_MM","Tool"]<-"TElocal_MM"
allmeta_loci[allmeta_loci$ID %like% "_TElocal_UM","Tool"]<-"TElocal_UM"
allmeta_loci[allmeta_loci$ID %like% "_SalmonTE","Tool"]<-"SalmonTE"
allmeta_loci[allmeta_loci$ID %like% "FC_MM_F","Tool"]<-"FC_MM_F"
allmeta_loci[allmeta_loci$ID %like% "FC_MM_R","Tool"]<-"FC_MM_R"
allmeta_loci[allmeta_loci$ID %like% "FC_UM","Tool"]<-"FC_UM"
allmeta_loci[allmeta_loci$Tool %like% "EM","EM"]<-"EM"
allmeta_loci[allmeta_loci$Tool %like% "MM","EM"]<-"EM"
allmeta_loci[allmeta_loci$Tool %like% "UM","EM"]<-"UM"
allmeta_loci[allmeta_loci$Tool %like% "TrueCount","EM"]<-"TrueCnt"
allmeta_loci[allmeta_loci$Tool %like% "SalmonTE","EM"]<-"EM"
allmeta_loci[allmeta_loci$Tool %like% "FC_MM_F","EM"]<-"MM_F"
allmeta_loci[allmeta_loci$Tool %like% "FC_MM_R","EM"]<-"MM_R"
rownames(allmeta_loci)<-allmeta_loci$ID


fwrite(allmeta_loci,"hsapiens_sim/stranded/sim01_allCnt_loci_meta.txt",sep="\t",col.names=TRUE)
