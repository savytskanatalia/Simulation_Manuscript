# Simulation, human
library(polyester)
library(seqinr)
library(tidyr)
library(stringr)

# SECTION 1
# Read in Genes expressed in Substantia nigra from GTEx 8
GTEx_v8_median_tpm <- read.delim("GTEx_v8_RNASeQCv1.1.9_Substantia_nigra_gene_median_tpm.gct", comment.char="#")
# Sort by expression in decreasing order
GTEx_v8_median_tpm<-GTEx_v8_median_tpm[order(GTEx_v8_median_tpm$Brain...Substantia.nigra,decreasing=TRUE),]
# Remove everything in gene IDs after .
GTEx_v8_median_tpm$Name<- sub("*\\.[0-9]", "", GTEx_v8_median_tpm$Name)
# Subsample top 13000 expressed genes
GTEx_v8_median_tpm<-GTEx_v8_median_tpm[1:13000,]
# Write out IDs for preparing .fa with transcripts derived from these genes
write.csv(GTEx_v8_median_tpm$Name,"Substantia_nigra_HS.ID",row.names = FALSE,quote=FALSE)


# SECTION 2
# Load human repeat masker
GRCh38_rmsk_TEinst_modified_complete_sorted <- read.delim("GRCh38_rmsk_TEinst_modified_complete_sorted.bed", header=FALSE)
# Total number: 4532102
# get the length
GRCh38_rmsk_TEinst_modified_complete_sorted$V11<-GRCh38_rmsk_TEinst_modified_complete_sorted$V3-GRCh38_rmsk_TEinst_modified_complete_sorted$V2
# exclude the short ones; for all - cutoff 200 because cutoff for long non-coding RNAs - for which RNAseq is appropriate. Shorter even if expressed- not appropriately estimated from this data

GRCh38_rmsk_TEinst_modified_complete_sorted<-GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V11>=200,]

# Create region column readible for samtools faidx
GRCh38_rmsk_TEinst_modified_complete_sorted$V12<-paste0(as.character(GRCh38_rmsk_TEinst_modified_complete_sorted$V1),(paste0(":",
                                                                                                             paste0(as.character(GRCh38_rmsk_TEinst_modified_complete_sorted$V2),
                                                                                                                    paste0("-",as.character(GRCh38_rmsk_TEinst_modified_complete_sorted$V3))))))

# Remove chrY elements, not interested in those for simulation
GRCh38_rmsk_TEinst_modified_complete_sorted<-GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V1!="chrY",]


# Write out 
write.table(GRCh38_rmsk_TEinst_modified_complete_sorted,'HS_TE_forSimulation.bed', col.names=FALSE,
            row.names=FALSE,
            quote=FALSE,sep="\t")


# As for faidx extraction used regions; for simulation logFC table need to use "dictionary" total file with ID-samtools faidx region to attribute region ID to TE




# SECTION 3



# For genes; load the transcript ID file for substantia nigra; subset random 2000 genes to be overexpressed and 2000 to be downregulated; assign same logFC to multiple transcripts of individual genes
# Subset and get only longest transcripts for each gene
# 

GRCh38_all <- read.delim("GRCh38_genes_lncRNAs_RepMask.fa.fai", header=FALSE)
# Subset only for genes
GRCh38_all<-GRCh38_all[GRCh38_all$V1 %like% "ENST",]
# separate gene from transcript
GRCh38_all<-separate(GRCh38_all,V1,into=c("transcript","trangenescript"),sep="_", remove=FALSE)
# Order by length
GRCh38_all<-GRCh38_all[order(GRCh38_all$V2,decreasing=TRUE),]
# Remove duplicates - retain only first-longest transcript- instance for each gene
GRCh38_all<-GRCh38_all[!(duplicated(GRCh38_all$trangenescript)),]

# Set seed and subset to 13000 genes
set.seed(25)
Expressed_genes<-GRCh38_all[GRCh38_all$trangenescript %in% sample(GRCh38_all$trangenescript,size=13000),]
# Assign FC of 1 to two new columns- sample groups
Expressed_genes$A<-1
Expressed_genes$B<-1
Expressed_genes$seeds<-1:nrow(Expressed_genes)

# Subsample for DE
set.seed(1)
up_genes<-Expressed_genes[Expressed_genes$trangenescript %in% sample(Expressed_genes$trangenescript,size=1250),"trangenescript"]
set.seed(3)
down_genes<-Expressed_genes[Expressed_genes$trangenescript %in% sample(Expressed_genes[!(Expressed_genes$trangenescript %in% up_genes) ,"trangenescript"],size=1250),"trangenescript"]

# Assign log fold changes, random logFC, set seed from the table 
for (i in up_genes){
  set.seed(Expressed_genes[Expressed_genes$trangenescript==i,"seeds"])
  Expressed_genes[Expressed_genes$trangenescript==i,"B"]<-round(runif(1,2,8))
}
for (i in down_genes){
  set.seed(Expressed_genes[Expressed_genes$trangenescript==i,"seeds"])
  Expressed_genes[Expressed_genes$trangenescript==i,"A"]<-round(runif(1,2,8))
}

write.table(Expressed_genes,'HS_Simulation/gene_Table.txt',quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
# Subset TEs to be expressed and to be DE
# Subfamilies to be DE: L1Hs, L1PA2, L1PA6, L1PA8, HERVK, HERVH, SVA_E, SVA_F, SVA_A, AluYa, AluYb9, AluYb8, AluSg1, AluSx1, AluSx3, AluSc8, AluJr

# Subset loci for DE expression:
# First random 100 for singular loci-drivers
set.seed(0)
TE_expr<-GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %in% sample(GRCh38_rmsk_TEinst_modified_complete_sorted$V4,size=100),]


DE_TE<-c("L1HS", "L1PA2", "L1PA6", "L1PA8", "HERVK", "HERVH", "SVA_E", "SVA_F", "SVA_A", "AluYa", "AluYb9", "AluYb8", "AluSg4", "AluSx1", "AluSx3", "AluSc8", "AluJr")
for (i in c(1:length(DE_TE))){
  
  print(i)
  print(DE_TE[i])
  if (length(GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %like% DE_TE[i],"V4"])>=300){
      set.seed(i)
      TE_expr<-rbind(TE_expr,GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %in% sample(GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %like% DE_TE[i],"V4"],size=200),])
    
  } else {
      set.seed(i)
      TE_expr<-rbind(TE_expr,GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %in% sample(GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %like% DE_TE[i],"V4"],size=100),])
    
  }
}

# Subset 3900 loci as expressed
set.seed(666)
TE_expr<-rbind(TE_expr,GRCh38_rmsk_TEinst_modified_complete_sorted[GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %in% sample(GRCh38_rmsk_TEinst_modified_complete_sorted[!(GRCh38_rmsk_TEinst_modified_complete_sorted$V4 %in% TE_expr$V4),"V4"],size=3900),])
rownames(TE_expr)<-1:nrow(TE_expr)

# Create logFC table.
# Loci in first 100 rows will have individual change in expression
# Up subfamilies: L1PA2, HERVK, SVA_F,AluYa, AluSg4,AluSx3, AluJr,"AluSc8"
# Down subfamilies: L1HS, L1PA6, L1PA8, "HERVH","SVA_E","SVA_A","AluYb8","AluYb9","AluSx1"
TE_expr[,c("A","B")]<-1

for (i in 1:100){
  set.seed(i)
  a<-runif(1,0,100)
  if (a>50){
    set.seed(i)
    TE_expr[i,"B"]<-round(runif(1,2,8))
  } else {
    set.seed(i)
    TE_expr[i,"A"]<-round(runif(1,2,8))
  }
}
# Up, Down vectors; for loop
TE_up<-c("L1PA2", "HERVK", "SVA_F","AluYa", "AluSg4","AluSx3", "AluJr","AluSc8")
TE_down<-c("L1HS", "L1PA6", "L1PA8", "HERVH","SVA_E","SVA_A","AluYb8","AluYb9","AluSx1")
for (i in 1:length(TE_up)){
  set.seed(i*100)
  TE_expr[TE_expr$V4 %like% TE_up[i],"B"]<-round(runif(1,2,8))
}
for (i in 1:length(TE_down)){
  set.seed(i*50)
  TE_expr[TE_expr$V4 %like% TE_down[i],"A"]<-round(runif(1,2,8))
}
# Subset random TE loci, then assign DE for the chosen subfamilies, then form a table
TE_expr<-TE_expr[,c("V4","V10","V12","V11","A","B")]
colnames(TE_expr)<-c("TE_ID","TE_Class","FA_ID","Length","A","B")
# Write out table
write.table(TE_expr,'TE_Table.txt',quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
colnames(TE_expr)<-c("TE_ID","TE_Class","V1","Length","A","B")

# logFC table
total_logfc<-rbind(Expressed_genes[,c("V1","A","B")],TE_expr[,c("V1","A","B")])

# Get length for genes+TEs
all_length<-c(Expressed_genes$V2,TE_expr$Length)
readspertx_g = round(round(runif(1,5,20))*all_length/100)


write.table(data.frame(ID=total_logfc$V1,A=total_logfc$A,B=total_logfc$B,length=all_length,reads_per_tr=readspertx_g),'/media/savytska/natalia/Projects/Simulations/HS_Simulation01/Total_ReadPerTrns.txt',quote=FALSE,row.names=FALSE,sep="\t")
write.table(total_logfc,'Total_logFC.txt',quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.csv(total_logfc$V1,'Total_ID.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
rownames(total_logfc)<-total_logfc$V1
total_logfc_m<-as.matrix(total_logfc[,-1])

# Run simulation
library(polyester)
# unstranded
simulate_experiment(fasta='Simulation_HS_01.fa', gtf=NULL, seqpath=NULL, num_reps=c(5,5), reads_per_transcript=readspertx_g, fold_changes=total_logfc_m, paired=TRUE, seed=100, gzip=TRUE,lib_sizes=c(5,5,5,5,5,5,5,5,5,5),outdir='unstranded/')


# Run simulation with stranded data 
simulate_experiment(fasta='Simulation_HS_01.fa', gtf=NULL, seqpath=NULL, num_reps=c(5,5), reads_per_transcript=Total_ReadPerTrns$reads_per_tr, fold_changes=as.matrix(Total_ReadPerTrns[,c("A","B")]), paired=TRUE, seed=25, strand_specific=TRUE, gzip=TRUE,lib_sizes=c(5,5,5,5,5,5,5,5,5,5),outdir='stranded/')

