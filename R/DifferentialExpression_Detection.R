################################################################################################################################################

# Example for Hsapiens stranded simulation; Differential Expression Detection benchmarking



## Differential Expression and FDR for it, simmulation, non-stranded data

myPaths <- .libPaths()

myPaths <- c("/data/natalia_lib", myPaths)

.libPaths(myPaths)

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(rlang)
library(DESeq2)
library(umap)
library(ggplot2)
library(data.table)
library(ggforce)
library(zeallot)
library(tidyr)
library(stringr)




## PCA for normalized counts, outputs one dataframe with first N PCs and table with first 6 PC+percentage of Variation

pca_dt<-function(object, ntop = 500,quant_m, sgroup){
  object<-as.matrix(object)
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(object[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  nrow(pca$x)
  length(quant_m)
  d2 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3=pca$x[, 3], PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6], Tools=quant_m, Group=sgroup)
  attr(d2, "percentVar") <- percentVar[1:6]
  return(list(pca$x,d2))
}


## Function to get PCA results - first 6 PC+percentVar from normalized counts
pca_pca <- function (object, ntop = 500,quant_m, sgroup) {
  object<-as.matrix(object)
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(object[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  nrow(pca$x)
  length(quant_m)
  d2 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3=pca$x[, 3], PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6], Tools=quant_m, Group=sgroup)
  attr(d2, "percentVar") <- percentVar[1:6]
  print(conames(d2))
  return(d2)}


## Function to get PCA results - first N PC  from normalized counts that can be fed to UMAP
pca_umap <- function (object, ntop = 500,quant_m, sgroup) {
  object<-as.matrix(object)
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(object[select, ]))
  return(pca$x)}


## Plotting UMAP plot. Needs output from umapa function+directory and prefix to save the file
umaplot <- function(df,mypath,outprefix) {
  df$Group<-as.factor(as.character(df$Group))
  uplot <- ggplot(df, aes(x=UMAP1,y=UMAP2, color=Tool, shape=Group)) + geom_point(size=8) + scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic()+ theme(panel.grid.minor=element_line(size=0.5),
                           panel.grid.major=element_line(size=0.5),
                           text=element_text(size=25),
                           legend.text=element_text(size=25),
                           axis.text.x = element_text( hjust=1))
  ggsave(filename =paste0(outprefix,"_UMAP.png"), plot = uplot, device = "png", path = mypath, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  return(uplot)}


# Calculating UMAP from PCA results
umapa<- function(df,quant_m, sgroup) {
  umap.custom<-umap.defaults
  umap.custom$metric<-"pearson"
  umap_o<-umap(df,umap.custom)
  umap_df<-as.data.frame(umap_o$layout)
  umap_df[,3]<-quant_m
  umap_df[,4]<-sgroup
  colnames(umap_df)<-c("UMAP1","UMAP2","Tool","Group")
  return(umap_df)}



# Plotting PCA results
FDplotPCA<-function(df,pc1=1,pc2=2,mypath="/data/natalia_lib/",filenm="PCA12_Tool.png"){
  df$Group<-as.factor(as.character(df$Group))
  pc1_str<-paste0("PC",as.character(pc1))
  pc2_str<-paste0("PC",as.character(pc2))
  PCA_plot<-ggplot(data = df, aes_string(x = pc1_str, y = pc2_str, color = "Tools", shape="Group")) +
    xlab(paste0(pc1_str,": ", round(attr(df,"percentVar")[pc1] * 100), "% variance")) +
    ylab(paste0(pc2_str,": ", round(attr(df,"percentVar")[pc2] * 100), "% variance")) + coord_fixed() + 
    geom_point(aes(size = 8)) +
    theme_classic()+theme(panel.grid.minor=element_line(size=0.5),
                          panel.grid.major=element_line(size=0.5),
                          text=element_text(size=25),
                          legend.text=element_text(size=25),
                          legend.position="bottom")+
    guides(color=guide_legend(title.position="top",title.hjust = 0.5),size=FALSE) + scale_color_viridis(option = "inferno",discrete=TRUE)
  ggsave(filename =filenm, plot = PCA_plot, device = "png", path = mypath, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  return(PCA_plot)
}



# Function, that calls DE on the total count table, gets normalized data, calls to PCA+UMAP and plotting the results
# Example to get output
# c(sim2_rep1_PC12,sim2_rep1_PC34,sim2_rep1_PC56,sim2_rep1_umap)%<-%simul_deseq(simulation02_wmask_loci_counts_mx,simulation2_loci_meta,"/data/natalia_lib/simulation02_all/","sim2_rep1_wmask")


simul_deseq<-function(mymx,mymeta,outdir,outprefix){
  
  desobj<-DESeqDataSetFromMatrix(countData=mymx,
                                 colData=mymeta,
                                 design=~Tool)
  # Set appropriate levels for factors
  desobj$Tool <- relevel(desobj$Tool, "TC")
  desobj<-estimateSizeFactors(desobj)
  desobj_nc<-counts(desobj, normalized=TRUE)
  fwrite(as.data.frame(desobj_nc),file=paste0(outdir,paste0(outprefix,"_MedRat_norm.txt")), row.names = TRUE, sep="\t")
  desobj_nc_df<-as.data.frame(desobj_nc)
  # pca_plotobj<-pca_pca(desobj_nc_df, ntop=20000,mymeta$Tool, mymeta$Group)
  #  pca_umapa<-pca_umap(desobj_nc_df, ntop=20000)
  c(pca_umapa,pca_plotobj)%<-%pca_dt(desobj_nc_df, ntop=20000,mymeta$Tool, mymeta$Group)
  
  PC12<-FDplotPCA(pca_plotobj,pc1=1,pc2=2,mypath=outdir,filenm=paste0(outprefix,"_PC12_Tool.png"))
  PC34<-FDplotPCA(pca_plotobj,pc1=3,pc2=4,mypath=outdir,filenm=paste0(outprefix,"_PC34_Tool.png"))
  PC56<-FDplotPCA(pca_plotobj,pc1=5,pc2=6,mypath=outdir,filenm=paste0(outprefix,"_PC56_Tool.png"))
  
  umap_obj<-umapa(pca_umapa,mymeta$Tool, mymeta$Group)
  uplot<-umaplot(umap_obj,outdir,outprefix) 
  return(list(PC12,PC34,PC56,uplot))
}






# Unified function to call DESeq2 and generate needed object

# Function, that catches possible stop() at Dispersion assessment for DESeq2 and switches to estimating estimateDispersionsGeneEst instead

disp_switch<-function(obj){
  obj1<-obj
  obj1 <- try(estimateDispersions(obj1))
  if (class(obj1)=="try-error") {
    print("Gotcha!")
    obj <- estimateDispersionsGeneEst(obj)
    dispersions(obj) <- mcols(obj)$dispGeneEst
    return(obj)
  } else {
    obj<-estimateDispersions(obj)
    return(obj)
  }
  return(obj)
}


# DE function



whole_dea<-function(mymx,mymeta,myoutputdir,outprefix){
  dds_object<-DESeqDataSetFromMatrix(countData=mymx,
                                     colData=mymeta,
                                     design=~Group)
  # Set appropriate levels for factors
  dds_object$Group <- relevel(dds_object$Group, "1")
  dds_object$Group 
  
  dds_object<-estimateSizeFactors(dds_object)
  sizeFactors(dds_object)
  dds_object<-disp_switch(dds_object)
  dds_object <- nbinomWaldTest(dds_object)
  dds_object_nc<-counts(dds_object, normalized=TRUE)
  fwrite(as.data.frame(results(dds_object, pAdjustMethod = "BH")),file=paste0(myoutputdir,(paste0(outprefix,"_DEres.txt"))),sep="\t",row.names=TRUE)
  dt<-as.data.table(results(dds_object, pAdjustMethod = "BH"), keep.rownames=TRUE)
  dt[,ID:=rn]
  setkey(dt,ID)
  dt[,rn:=NULL]
  return(dt)
}




de_stats2<-function(dt,dt_abs,cutoff2, mytool){
  dt2<-rt_clean(dt)
  setkey(dt_abs,ID)
  setkey(dt2,ID)
  dt_fin<-dt2[dt_abs,]
  dt3<-empty_dt(8,8)
  print(dt3)
  colnames(dt3)<-c("Tool","TP","FP","TN","FN","F1","FDR","baseMean")
  dt3$Tool<-as.character(dt3$Tool)
  dt3[,Tool:=mytool]
  counter<-1
  for (i in c(0,5,10,20,30,50,70,100)){
    dt3$baseMean[counter]<-i
    AP<-nrow(dt_fin[(dt_fin$baseMean>i & dt_fin$padj<cutoff2) & (dt_fin$log2FoldChange>2 | dt_fin$log2FoldChange<(-2)),])
    
    print(AP)
    dt3$TP[counter]<-nrow(dt_fin[(dt_fin$baseMean>i & dt_fin$padj<cutoff2) & ((dt_fin$log2FoldChange>2 & dt_fin$log2FC_T>0) |(dt_fin$log2FoldChange<(-2) & dt_fin$log2FC_T<0) ),])
    dt3$FP[counter]<-AP-dt3$TP[counter]
    dt3$FN[counter]<-nrow(dt_fin[dt_fin$Status!=0,])-dt3$TP[counter]
    AS<-nrow(dt_fin)
    dt3$TN[counter]<-AS-dt3$FN[counter]-AP
    dt3$FDR[counter]<-dt3$FP[counter]/(dt3$FP[counter]+dt3$TP[counter])
    dt3$F1[counter]<-2*((( dt3$TP[counter]/(dt3$TP[counter]+dt3$FN[counter]))*dt3$TP[counter]/(dt3$TP[counter]+dt3$FP[counter]))/(( dt3$TP[counter]/(dt3$TP[counter]+dt3$FN[counter]))+(dt3$TP[counter]/(dt3$TP[counter]+dt3$FP[counter]))))
    counter<-counter+1
  }
  
  return(list(dt3,dt_fin))
}


# Function to collect FDR/F1 statistics for DE results based on the chosen expression cutoff for a range of TE length cutoffs
len_cut_stat<-function(dt_fin,cutoff2, mytool,exp_c,repmsk){ # input - DE table+absolute status from original table = dt_fin, obtained as dt_fin with de_stats2; cutoff2 - pval cutoff, mytool - toolname, expression cutoff-exp_c, repmasker mod table
  dt3<-empty_dt(8,7)
  print(dt3)
  colnames(dt3)<-c("Tool","TP","FP","TN","FN","F1","FDR","Length")
  dt3$Tool<-as.character(dt3$Tool)
  dt3[,Tool:=mytool]
  counter<-1
  for (i in c(0,50,100,150,200,250,300)){
    dt3$Length[counter]<-i
    cst_mx<-id_len_f(repmsk,i,dt_fin)
    dt3<-stat_collector(dt3,cst_mx,cutoff2, exp_c,counter)
    counter<-counter+1
  }
  return(dt3)
}

# Function to subset DE table to only TEs longer than cutoff i 
id_len_f<-function(repmsk,i,dt_fin){
  ID_filt<-repmsk[(repmsk$V5-repmsk$V4)>i,"ID"][[1]]
  print(length(ID_filt))
  ID_filt<-ID_filt[ID_filt %in% dt_fin$ID]
  cst_mx<-dt_fin[ID_filt,]
  return(cst_mx)
}

# Stat collector function that counts TP..FDR for DE table (dt_fin) with pval cutoff (cutoff2) and expression cutoff (exp_c), for stat table (dt3)
stat_collector<-function(dt3,dt_fin,cutoff2, exp_c,counter){
  AP<-nrow(dt_fin[(dt_fin$baseMean>exp_c & dt_fin$padj<cutoff2) & (dt_fin$log2FoldChange>2 | dt_fin$log2FoldChange<(-2)),])
  dt3$TP[counter]<-nrow(dt_fin[(dt_fin$baseMean>exp_c & dt_fin$padj<cutoff2) & ((dt_fin$log2FoldChange>2 & dt_fin$log2FC_T>0) |(dt_fin$log2FoldChange<(-2) & dt_fin$log2FC_T<0) ),])
  dt3$FP[counter]<-AP-dt3$TP[counter]
  dt3$FN[counter]<-nrow(dt_fin[dt_fin$Status!=0,])-dt3$TP[counter]
  AS<-nrow(dt_fin)
  dt3$TN[counter]<-AS-dt3$FN[counter]-AP
  dt3$FDR[counter]<-dt3$FP[counter]/(dt3$FP[counter]+dt3$TP[counter])
  dt3$F1[counter]<-2*((( dt3$TP[counter]/(dt3$TP[counter]+dt3$FN[counter]))*dt3$TP[counter]/(dt3$TP[counter]+dt3$FP[counter]))/(( dt3$TP[counter]/(dt3$TP[counter]+dt3$FN[counter]))+(dt3$TP[counter]/(dt3$TP[counter]+dt3$FP[counter]))))
return(dt3)}

# Function for getting table of all measurements
# F1 Score = 2*(Recall * Precision) / (Recall + Precision)
# Recall (TPR) = TP/(TP+FN)
# Precision = TP/(TP+FP)
# FDR = FP/(FP+TP)

empty_dt<-function(ncols, nrows){
  dt<-as.data.table(matrix(rep(1:nrows,each=ncols),ncol=ncols,byrow=TRUE))
  dt[,]<-0
  return(dt)
}
# Function to substitute NAs with specific values, in this case 0, https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
na_to_z <- function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0.5]
  return(DT)
}


# Function for removing extra unused columns 

rt_clean<-function(dt){
  na_to_z(dt)
  return(dt)
}

## Expanding the table with DE absolute status for loci or subfamily, having a table of it + full list of subf/loci in dataset
## For subfamilies-needs to additionally take care of the subfamilies that have loci expressed in both directions


tabexpand_universal <- function(fullist, etable){
  etable<-as.data.frame(etable)
  abs_dif<-fullist[!fullist %in% rownames(etable)]
  etable[abs_dif,]<-0
  etable<-as.data.table(etable, keep.rownames=TRUE)
  etable$ID<-etable$rn
  setkey(etable,ID)
  etable[,rn:=NULL]
  return(etable) }


###############################################################################################################################
###############################################################################################################################



######
######




# Loading absolute values for simulation, stranded
sim_tx_info <- read.delim("hsapiens_sim/stranded/sim_tx_info.txt", stringsAsFactors=FALSE)
rownames(sim_tx_info)<-sim_tx_info$transcriptid
sim_tx_info<-sim_tx_info[rownames(sim_tx_info) %in% rownames(sim_tx_info)[!(rownames(sim_tx_info) %in% rownames(sim_tx_info[grep("ENST",rownames(sim_tx_info),value=TRUE),]))],]
sim_tx_info$ID<-""

# Translate the IDs
# Load dictionary, just conveniently saved TE annotation .gtf column 9 into separate columns (ID, T_ID, family, class...)
TE_Table <- read.delim("/data/natalia/hsapiens_sim/TE_Table.txt", stringsAsFactors=FALSE)
identical(TE_Table$FA_ID,sim_tx_info$transcriptid)
# [1] TRUE
sim_tx_info$ID<-TE_Table$TE_ID



# Adding columns to absolute DE status for processing data further
# If DE>CTRL, Status=1, log2FC>0; If DE<CTRL, Status=-1, log2FC<0
sim_tx_info<-as.data.table(sim_tx_info,key="ID")
sim_tx_info[,c("Status","log2FC_T"):=list(0,0)]
for (i in c(1:nrow(sim_tx_info))){
  if (sim_tx_info[i,"foldchange.A"]==sim_tx_info[i,"foldchange.B"]){
    sim_tx_info[i,"Status"]=0
    sim_tx_info[i,"log2FC_T"]=0
  } else if (sim_tx_info[i,"foldchange.B"]>sim_tx_info[i,"foldchange.A"]) {
    sim_tx_info[i,"Status"]=1
    sim_tx_info[i,"log2FC_T"]=log2(sim_tx_info[i,"foldchange.B"])
  } else if (sim_tx_info[i,"foldchange.B"]<sim_tx_info[i,"foldchange.A"]) {
    sim_tx_info[i,"Status"]= -1
    sim_tx_info[i,"log2FC_T"]=-1*log2(sim_tx_info[i,"foldchange.A"])
  } else
  {print("Something wrong with this line!")
    print(sim_tx_info[i,])}
}



# Subset to have a smaller dataframe to operate with
sim_tx_info_l<-sim_tx_info[,c("ID","Status","log2FC_T")]
sim_tx_info_l<-setDF(sim_tx_info_l,rownames = sim_tx_info_l$ID)





# counts_matrix<-fread("hsapiens_sim/sim01/sim01_allCnt_loci.txt")
# # dt>df>mx
# counts_matrix<-as.data.frame(counts_matrix)
# rownames(counts_matrix)<-counts_matrix$rn


counts_matrix<-as.data.frame(counts_matrix_c)
rownames(counts_matrix)<-counts_matrix$rn
# round counts (for the EM-derived counts)
counts_matrix<-as.matrix(round(counts_matrix[,-1]))


# Now need to add non-expressed TEs to the list, indicating they are not DE
# Need all_loci<-gen_l2$ID
# Getting this list from count table
all_loci<-rownames(counts_matrix)
sim_tx_info_l<-tabexpand_universal(all_loci,sim_tx_info_l)


# checking colnames(counts_matrix)==rownames(allmeta_loci)
colnames(allmeta_loci)[2]<-"Group"

allmeta_loci$Group<-as.factor(as.character(allmeta_loci$Group))
#colnames(simulation02_wmask_loci_counts_mx)==simulation2_loci_meta$Sample
allmeta_loci[allmeta_loci$Tool=="TrueCount","Tool"]<-"TC"
c(sim01_ns_loci_PC12,sim01_ns_loci_PC34,sim01_ns_loci_PC56,sim01_ns_loci_umap)%<-%simul_deseq(counts_matrix,allmeta_loci,"hsapiens_sim/stranded/","sim01_st_loci")


PC12<-sim01_ns_loci_PC12+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
PC34<-sim01_ns_loci_PC34+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
PC56<-sim01_ns_loci_PC56+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
uplot<-sim01_ns_loci_umap+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()

p_u_batch <- plot_grid(
  PC12 + theme(legend.position="none"),
  PC34 + theme(legend.position="none"),
  PC56 + theme(legend.position="none"),
  uplot + theme(legend.position="none"),
  align = 'vh',labels = c("A", "B", "C","D"), hjust = -1,nrow = 2, ncol=2)
p_u_batch

legend_u_b <- get_legend(PC12 +guides(color = guide_legend(nrow = 3),shape = guide_legend(nrow = 2),
                                      guide = guide_legend(direction = "horizontal",  label.position="bottom",
                                                           label.hjust = 0.5, label.vjust = 0.5,label.theme = element_text(angle = 90)))+ theme(legend.position = "bottom"))
p4<-plot_grid(p_u_batch, legend_u_b,ncol = 1, rel_heights = c(1, .1))
save_plot("/data/natalia/hsapiens_sim/stranded/loci_PCA_UMAP.png", p4, ncol = 1, base_height = 10)






# colnames(simulation02_wmask_loci_counts_mx)==simulation2_loci_meta$Sample

# SalTE
DE_SalTE_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "SalmonTE"]],allmeta_loci[allmeta_loci$Tool %like% "SalmonTE",],"hsapiens_sim/stranded/","sim01_st_loci_SalTE")
c(DE_SalTE_sim2_rep1_l_stat,DE_SalTE_sim2_rep1_l2)%<-%de_stats2(DE_SalTE_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "SalTE")

DE_SalTE_len_stat_5<-len_cut_stat(DE_SalTE_sim2_rep1_l2,cutoff2=0.01, "SalTE",5, repmsk)
DE_SalTE_len_stat_50<-len_cut_stat(DE_SalTE_sim2_rep1_l2,cutoff2=0.01, "SalTE",50, repmsk)


# TElocal_MM


DE_TElocal_MM_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "TElocal_MM"]],allmeta_loci[allmeta_loci$Tool %like% "TElocal_MM",],"hsapiens_sim/stranded/","sim01_st_loci_TElocalMM")

c(DE_TElocal_MM_sim2_rep1_l_stat,DE_TElocal_MM_sim2_rep1_l2)%<-%de_stats2(DE_TElocal_MM_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "TElocal_MM")

DE_TElocal_MM_len_stat_5<-len_cut_stat(DE_TElocal_MM_sim2_rep1_l2,cutoff2=0.01, "TElocal_MM",5, repmsk)
DE_TElocal_MM_len_stat_50<-len_cut_stat(DE_TElocal_MM_sim2_rep1_l2,cutoff2=0.01, "TElocal_MM",50, repmsk)

# TElocal_UM
DE_TElocal_UM_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "TElocal_UM"]],allmeta_loci[allmeta_loci$Tool %like% "TElocal_UM",],"hsapiens_sim/stranded/","sim01_st_loci_TElocalUM")

c(DE_TElocal_UM_sim2_rep1_l_stat,DE_TElocal_UM_sim2_rep1_l2)%<-%de_stats2(DE_TElocal_UM_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "TElocal_UM")

DE_TElocal_UM_len_stat_5<-len_cut_stat(DE_TElocal_UM_sim2_rep1_l2,cutoff2=0.01, "TElocal_UM",5, repmsk)
DE_TElocal_UM_len_stat_50<-len_cut_stat(DE_TElocal_UM_sim2_rep1_l2,cutoff2=0.01, "TElocal_UM",50, repmsk)


# SQuIRE_UM
DE_SQuIRE_UM_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "SQuIRE_UM"]],allmeta_loci[allmeta_loci$Tool %like% "SQuIRE_UM",],"hsapiens_sim/stranded/","sim01_st_loci_SQuIRE_UM")
c(DE_SQuIRE_UM_sim2_rep1_l_stat,DE_SQuIRE_UM_sim2_rep1_l2)%<-%de_stats2(DE_SQuIRE_UM_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "SQuIRE_UM")

DE_SQuIRE_UM_len_stat_5<-len_cut_stat(DE_SQuIRE_UM_sim2_rep1_l2,cutoff2=0.01, "SQuIRE_UM",5, repmsk)
DE_SQuIRE_UM_len_stat_50<-len_cut_stat(DE_SQuIRE_UM_sim2_rep1_l2,cutoff2=0.01, "SQuIRE_UM",50, repmsk)

# SQuIRE_EM
DE_SQuIRE_EM_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "SQuIRE_EM"]],allmeta_loci[allmeta_loci$Tool %like% "SQuIRE_EM",],"hsapiens_sim/stranded/","sim01_st_loci_SQuIRE_EM")
c(DE_SQuIRE_EM_sim2_rep1_l_stat,DE_SQuIRE_EM_sim2_rep1_l2)%<-%de_stats2(DE_SQuIRE_EM_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "SQuIRE_EM")

DE_SQuIRE_EM_len_stat_5<-len_cut_stat(DE_SQuIRE_EM_sim2_rep1_l2,cutoff2=0.01, "SQuIRE_EM",5, repmsk)
DE_SQuIRE_EM_len_stat_50<-len_cut_stat(DE_SQuIRE_EM_sim2_rep1_l2,cutoff2=0.01, "SQuIRE_EM",50, repmsk)



# FC_MM_F
DE_FC_MM_F_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "FC_MM_F"]],allmeta_loci[allmeta_loci$Tool %like% "FC_MM_F",],"hsapiens_sim/stranded/","sim01_st_loci_FC_MM_F")
c(DE_FC_MM_F_sim2_rep1_l_stat,DE_FC_MM_F_sim2_rep1_l2)%<-%de_stats2(DE_FC_MM_F_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "FC_MM_F")

DE_FC_MM_F_len_stat_5<-len_cut_stat(DE_FC_MM_F_sim2_rep1_l2,cutoff2=0.01, "FC_MM_F",5, repmsk)
DE_FC_MM_F_len_stat_50<-len_cut_stat(DE_FC_MM_F_sim2_rep1_l2,cutoff2=0.01, "FC_MM_F",50, repmsk)


# FC_MM_R
DE_FC_MM_R_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "FC_MM_R"]],allmeta_loci[allmeta_loci$Tool %like% "FC_MM_R",],"hsapiens_sim/stranded/","sim01_st_loci1_FC_MM_R")
c(DE_FC_MM_R_sim2_rep1_l_stat,DE_FC_MM_R_sim2_rep1_l2)%<-%de_stats2(DE_FC_MM_R_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "FC_MM_R")


DE_FC_MM_R_len_stat_5<-len_cut_stat(DE_FC_MM_R_sim2_rep1_l2,cutoff2=0.01, "FC_MM_R",5, repmsk)
DE_FC_MM_R_len_stat_50<-len_cut_stat(DE_FC_MM_R_sim2_rep1_l2,cutoff2=0.01, "FC_MM_R",50, repmsk)

# FC_UM
DE_FC_UM_sim2_rep1_l<-whole_dea(counts_matrix[,colnames(counts_matrix)[colnames(counts_matrix) %like% "FC_UM"]],allmeta_loci[allmeta_loci$Tool %like% "FC_UM",],"hsapiens_sim/stranded/","sim01_st_loci_FC_UM")
c(DE_FC_UM_sim2_rep1_l_stat,DE_FC_UM_sim2_rep1_l2)%<-%de_stats2(DE_FC_UM_sim2_rep1_l,sim_tx_info_l,cutoff2=0.01, "FC_UM")

DE_FC_UM_len_stat_5<-len_cut_stat(DE_FC_UM_sim2_rep1_l2,cutoff2=0.01, "FC_UM",5, repmsk)
DE_FC_UM_len_stat_50<-len_cut_stat(DE_FC_UM_sim2_rep1_l2,cutoff2=0.01, "FC_UM",50, repmsk)

###########################################################################################################



sim2_r1_l_total_stats<-do.call("rbind", list(DE_SalTE_sim2_rep1_l_stat, DE_TElocal_MM_sim2_rep1_l_stat, DE_TElocal_UM_sim2_rep1_l_stat,
                                             DE_SQuIRE_UM_sim2_rep1_l_stat,DE_SQuIRE_EM_sim2_rep1_l_stat,DE_FC_MM_F_sim2_rep1_l_stat,DE_FC_MM_R_sim2_rep1_l_stat,DE_FC_UM_sim2_rep1_l_stat))


sim2_r1_l_total_stats$baseMean<-factor(sim2_r1_l_total_stats$baseMean, levels=c("0","5","10","20","30","50","70","100"))
fwrite(sim2_r1_l_total_stats,file="hsapiens_sim/stranded/simulation01_st_loci_DE_F1FDR.txt", sep="\t",quote="auto", col.names=TRUE)


sim2_r1_l_F1<-ggplot(sim2_r1_l_total_stats, aes(x=baseMean, y=F1,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff, baseMean")+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


sim2_r1_l_FDR<-ggplot(sim2_r1_l_total_stats, aes(x=baseMean, y=FDR,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff, baseMean")+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


F1c<-sim2_r1_l_F1+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
FDR1c<-sim2_r1_l_FDR+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()

F1FDR_batch <- plot_grid(
  F1c + theme(legend.position="none"),
  FDR1c + theme(legend.position="none"),
  align = 'vh',labels = c("F1", "FDR"), hjust = -1,ncol=2)
p_u_batch

legend_f1fdr_l <- get_legend(F1c +guides(color = guide_legend(nrow = 3),shape = guide_legend(nrow = 2),
                                         guide = guide_legend(direction = "horizontal",  label.position="bottom",
                                                              label.hjust = 0.5, label.vjust = 0.5,label.theme = element_text(angle = 90)))+ theme(legend.position = "bottom"))
plot_grid(F1FDR_batch, legend_f1fdr_l,ncol = 1, rel_heights = c(1, .1))

# 1200x700


# With length cutoffs

# Cutoff 5


sim2_r1_l_total_stats_5<-do.call("rbind", list(DE_SalTE_len_stat_5, DE_TElocal_MM_len_stat_5, DE_TElocal_UM_len_stat_5, 
                                             DE_SQuIRE_UM_len_stat_5,DE_SQuIRE_EM_len_stat_5, DE_FC_MM_F_len_stat_5, DE_FC_MM_R_len_stat_5, DE_FC_UM_len_stat_5))


sim2_r1_l_total_stats_5$Length<-factor(sim2_r1_l_total_stats_5$Length, levels=c("0", "50", "100", "150", "200", "250", "300"))
fwrite(sim2_r1_l_total_stats_5,file="hsapiens_sim/stranded/simulation01_st_loci_DE_F1FDR_Expr5_Length.txt", sep="\t",quote="auto", col.names=TRUE)


sim2_r1_l_F1<-ggplot(sim2_r1_l_total_stats_5, aes(x=Length, y=F1,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff Length, bp")+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


sim2_r1_l_FDR<-ggplot(sim2_r1_l_total_stats_5, aes(x=Length, y=FDR,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff Length, bp")+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


F1c<-sim2_r1_l_F1+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
FDR1c<-sim2_r1_l_FDR+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()

F1FDR_batch <- plot_grid(
  F1c + theme(legend.position="none"),
  FDR1c + theme(legend.position="none"),
  align = 'vh',labels = c("F1", "FDR"), hjust = -1,ncol=2)
p_u_batch

legend_f1fdr_l <- get_legend(F1c +guides(color = guide_legend(nrow = 3),shape = guide_legend(nrow = 2),
                                         guide = guide_legend(direction = "horizontal",  label.position="bottom",
                                                              label.hjust = 0.5, label.vjust = 0.5,label.theme = element_text(angle = 90)))+ theme(legend.position = "bottom"))
plot_grid(F1FDR_batch, legend_f1fdr_l,ncol = 1, rel_heights = c(1, .1))



# Expression cutoff 50

sim2_r1_l_total_stats_50<-do.call("rbind", list(DE_SalTE_len_stat_50, DE_TElocal_MM_len_stat_50, DE_TElocal_UM_len_stat_50, 
                                               DE_SQuIRE_UM_len_stat_50,DE_SQuIRE_EM_len_stat_50, DE_FC_MM_F_len_stat_50, DE_FC_MM_R_len_stat_50, DE_FC_UM_len_stat_50))


sim2_r1_l_total_stats_50$Length<-factor(sim2_r1_l_total_stats_50$Length, levels=c("0", "50", "100", "150", "200", "250", "300"))
fwrite(sim2_r1_l_total_stats_50,file="hsapiens_sim/stranded/simulation01_st_loci_DE_F1FDR_Expr50_Length.txt", sep="\t",quote="auto", col.names=TRUE)


sim2_r1_l_F1<-ggplot(sim2_r1_l_total_stats_50, aes(x=Length, y=F1,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff Length, bp")+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


sim2_r1_l_FDR<-ggplot(sim2_r1_l_total_stats_50, aes(x=Length, y=FDR,color=Tool,group=Tool)) + 
  geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+ scale_color_viridis(option = "inferno",discrete=TRUE)+
  theme_classic()+theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                        legend.text=element_text(size=15),
                        axis.text = element_text(size=15),legend.position="bottom")+
  xlab("Cutoff Length, bp")+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
  geom_hline(yintercept = 0.5, color="black", linetype="dashed")


F1c<-sim2_r1_l_F1+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
FDR1c<-sim2_r1_l_FDR+theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()

F1FDR_batch <- plot_grid(
  F1c + theme(legend.position="none"),
  FDR1c + theme(legend.position="none"),
  align = 'vh',labels = c("F1", "FDR"), hjust = -1,ncol=2)
p_u_batch

legend_f1fdr_l <- get_legend(F1c +guides(color = guide_legend(nrow = 3),shape = guide_legend(nrow = 2),
                                         guide = guide_legend(direction = "horizontal",  label.position="bottom",
                                                              label.hjust = 0.5, label.vjust = 0.5,label.theme = element_text(angle = 90)))+ theme(legend.position = "bottom"))
plot_grid(F1FDR_batch, legend_f1fdr_l,ncol = 1, rel_heights = c(1, .1))
