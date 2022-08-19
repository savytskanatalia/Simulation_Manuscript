# StatTable for quantification
# Example with Hsapiens stranded simulation





# Creating empty table for ROC metrics
# tools- vector with tool names, as represented in count table
# nsample- int, number of samples analyzed

wholempty_table <- function(tools,nsample){
  s2_roc<-as.data.frame(matrix(rep(1:(nsample*length(tools)),each=11),ncol=11,byrow=TRUE))
  s2_roc[,]<-0
  colnames(s2_roc)<-c("Samples","Tools","Group","TP","FP","TN","FN","Precision","Recall","F1","FDR")
  
  s2_roc$Samples<-rep(paste0(rep("sample_", times=nsample), str_pad(c(1:nsample),2,pad="0")), times=length(tools))
  
  s2_roc$Group<-rep(c(rep("1",times=(nsample/2)),rep("2",times=(nsample/2))),times=length(tools))
  s2_roc$Tools<-rep(tools,each=nsample)
  s2_roc<-as.data.table(s2_roc,key = c("Samples","Tools"))
  return(s2_roc)
}



# Function for getting table of all measurements
# F1 Score = 2*(Recall * Precision) / (Recall + Precision)
# Recall (TPR) = TP/(TP+FN)
# Precision = TP/(TP+FP)
# FDR = FP/(FP+TP)

wholespace_table_custom<-function(tools,all,df,nsample,cutoff=0){
  for (i in c(1:nsample)){
    for (j in tools){
      # .(2,c("A","C"))
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"]<-nrow(all[eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),"_TC")))!=0 & eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),j)))>cutoff])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FP"]<-nrow(all[eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),"_TC")))==0 & eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),j)))>cutoff])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TN"]<-nrow(all[eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),"_TC")))==0 & eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),j)))<=cutoff])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FN"]<-nrow(all[eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),"_TC")))!=0 & eval(parse(text=paste0(paste0("sample_",str_pad((i),2,pad="0")),j)))<=cutoff])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Precision"]<-df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"]/(df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"]+df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FP"])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Recall"]<-df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"]/(df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"]+df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FN"])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"F1"]<-2*(df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Recall"]*df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Precision"])/(df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Precision"]+df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"Recall"])
      df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FDR"]<-df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FP"]/(df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"FP"]+df[paste0(paste0("sample_",str_pad((i),2,pad="0")),j),"TP"])
    }}
  return(df)}
# paste0(paste0("sample_",str_pad((i),2,pad="0")),"_TC")
# paste0(paste0("sample_",str_pad((i),2,pad="0")),j)
# counts_matrix<-as.data.frame(counts_matrix_c)
# rownames(counts_matrix)<-counts_matrix$rn
# counts_matrix<-counts_matrix[,-1]




sim01_loci_quant_stat<-wholempty_table(unique(allmeta_loci$Tool)[-1],10)
sim01_loci_quant_stat<-as.data.frame(sim01_loci_quant_stat)
rownames(sim01_loci_quant_stat)<-paste0(sim01_loci_quant_stat$Samples,paste0("_",sim01_loci_quant_stat$Tools))


sim01_loci_quant_stat0<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,0)
sim01_loci_quant_stat0$Cutoff<-0
sim01_loci_quant_stat5<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,5)
sim01_loci_quant_stat5$Cutoff<-5
sim01_loci_quant_stat10<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,10)
sim01_loci_quant_stat10$Cutoff<-10
sim01_loci_quant_stat20<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,20)
sim01_loci_quant_stat20$Cutoff<-20
sim01_loci_quant_stat30<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,30)
sim01_loci_quant_stat30$Cutoff<-30
sim01_loci_quant_stat50<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,50)
sim01_loci_quant_stat50$Cutoff<-50
sim01_loci_quant_stat70<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,70)
sim01_loci_quant_stat70$Cutoff<-70
sim01_loci_quant_stat100<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),counts_matrix_c,sim01_loci_quant_stat,10,100)
sim01_loci_quant_stat100$Cutoff<-100


total_df_loci_quant<-rbind(sim01_loci_quant_stat0,rbind(sim01_loci_quant_stat5,(rbind(sim01_loci_quant_stat10,rbind(sim01_loci_quant_stat20,rbind(sim01_loci_quant_stat30,rbind(sim01_loci_quant_stat50,rbind(sim01_loci_quant_stat70,sim01_loci_quant_stat100))))))))

fwrite(total_df_loci_quant,"hsapiens_sim/stranded/sim01_quant_loci_stat.txt", sep="\t",quote="auto", col.names=TRUE,row.names=TRUE)



###########################################################################


## Plotting F1 and FDR rates for quantification at different cutoffs


## Plotting F1+FDR for Quantification step
# df - df containing TP,TN,FP,FN,FDR, F1 + Cutoff stating for which count cutoff it was calculated.
# outfir- directory to save files to , outprefix - prefix for saved final files


plot_f1fdr<-function(df,outdir,outprefix,cutoff){
  
  df<-df[order(df$Cutoff),]
  df$Cutoff<-as.character(df$Cutoff)
  if (cutoff=="count"){
    df$Cutoff<-factor(df$Cutoff, levels=c("0","5","10","20","30","50","70","100"))
  } else if (cutoff=="length"){
    df$Cutoff<-factor(df$Cutoff, levels=c("0", "50", "100", "150", "200", "250", "300"))
  } else {
    print("Something fishy with Cutoff system, input it please")
  }
  
  
  df_1 = df %>% group_by(Tools, Group, Cutoff) %>%
    summarise(mean_gr = mean(FDR), sd_gr = sd(FDR))
  df_2 = df %>% group_by(Tools,Cutoff) %>%
    summarise(mean_gr = mean(FDR), sd_gr = sd(FDR))
  df_1_group1 = subset(df_1, df_1$Group == 1)
  df_1_group2 = subset(df_1, df_1$Group == 2)
  
  
  p0<-ggplot(df_2, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
    geom_hline(yintercept = 0.5, color="black", linetype="dashed")
  ggsave(filename =paste0(outprefix,"_joint_FDR.png"), plot = p0, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  
  p1<-ggplot(df_1_group1, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
    geom_hline(yintercept = 0.5, color="black", linetype="dashed")
  ggsave(filename =paste0(outprefix,"_g1_FDR.png"), plot = p1, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  p2<-ggplot(df_1_group2, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("FDR")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)+
    geom_hline(yintercept = 0.5, color="black", linetype="dashed")
  ggsave(filename =paste0(outprefix,"_g2_FDR.png"), plot = p2, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  df_1f = df %>% group_by(Tools, Group, Cutoff) %>%
    summarise(mean_gr = mean(F1), sd_gr = sd(F1))
  df_2f = df %>% group_by(Tools,Cutoff) %>%
    summarise(mean_gr = mean(F1), sd_gr = sd(F1))
  df_1f_group1 = subset(df_1f, df_1$Group == 1)
  df_1f_group2 = subset(df_1f, df_1$Group == 2)
  
  p0f<-ggplot(df_2f, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)
  ggsave(filename =paste0(outprefix,"_joint_F1.png"), plot = p0f, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  p1f<-ggplot(df_1f_group1, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)
  ggsave(filename =paste0(outprefix,"_g1_F1.png"), plot = p1f, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  
  p2f<-ggplot(df_1f_group2, aes(x=Cutoff, y=mean_gr, group=Tools, color=Tools)) + 
    geom_line(size=2)+geom_point(shape = 21, fill = "black", size = 3, stroke = 1)+
    geom_errorbar(aes(ymin=mean_gr - sd_gr, ymax=mean_gr + sd_gr), width=.5, position=position_dodge(0.05))+ scale_color_viridis(option = "inferno",discrete=TRUE)+
    theme_classic() +theme(panel.grid.minor=element_line(size=0.5),panel.grid.major=element_line(size=0.5),text=element_text(size=25),
                           legend.text=element_text(size=15),
                           axis.text = element_text(size=15),legend.position="bottom")+
    xlab(paste0("Cutoff, ",cutoff))+ylab("F1")+guides(color=guide_legend(ncol=2,nrow=5))+ylim(0,0.95)
  ggsave(filename =paste0(outprefix,"_g2_F1.png"), plot = p2f, device = "png", path = outdir, scale = 1,width = 30 ,height = 25,dpi = 300, units="cm")
  
  p0c<-p0+ theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
  p1c<-p1+ theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
  p2c<-p2+ theme(plot.margin = margin(6, 0, 6, 0)) +  panel_border()
  p0fc<-p0f+ theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
  p1fc<-p1f+ theme(plot.margin = margin(6, 0, 6, 2)) + panel_border()
  p2fc<-p2f+ theme(plot.margin = margin(6, 0, 6, 0)) +  panel_border()
  
  p_u_batch <- plot_grid(
    p0c + theme(legend.position="none"),p1c + theme(legend.position="none"), p2c + theme(legend.position="none"),
    p0fc + theme(legend.position="none"),p1fc + theme(legend.position="none"),p2fc + theme(legend.position="none"),
    align = 'vh',labels = c("A", "B", "C","D", "E", "F"), hjust = -1,nrow = 2, ncol=3
  )
  p_u_batch
  
  legend_u_b <- get_legend(
    p0c +guides(color = guide_legend(nrow = 2),
                guide = guide_legend(direction = "horizontal",  label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,label.theme = element_text(angle = 90)))+ theme(legend.position = "bottom")
  )
  
  plot_grid(p_u_batch, legend_u_b,ncol = 1, rel_heights = c(1, .1))
  
}





plot_f1fdr(total_df_loci_quant,"hsapiens_sim/stranded/","sim01_loci_quant",cutoff="count")

########################################

# Then produce for count cutoff=5, and lengths 0,50,100,150,200,250,300
# Subset input table by length in repmsk, then 

# Set baseline table for binding and no cutoff
l_cst_stat<-sim01_loci_quant_stat5
l_cst_stat$Cutoff<-0
l_cst_stat$Sample_ID<-rownames(l_cst_stat)
rownames(l_cst_stat)<-1:nrow(l_cst_stat)
setkey(counts_matrix_c,rn)
for (i in c(50,100,150,200,250,300)){
  ID_filt<-repmsk[(repmsk$V5-repmsk$V4)>i,"ID"][[1]]
  print(length(ID_filt))
  ID_filt<-ID_filt[ID_filt %in% counts_matrix_c$rn]
  cst_mx<-counts_matrix_c[ID_filt,]
  cst_stat<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),cst_mx,sim01_loci_quant_stat,10,5)
  cst_stat$Cutoff<-i
  cst_stat$Sample_ID<-rownames(cst_stat)
  rownames(cst_stat)<-1:nrow(cst_stat)
  l_cst_stat<-rbind(l_cst_stat,cst_stat)
}

plot_f1fdr(l_cst_stat,"hsapiens_sim/stranded/","sim01_loci_quant_cutoff5_length",cutoff="length")





# Cutoff expr 50

l_cst_stat<-sim01_loci_quant_stat50
l_cst_stat$Cutoff<-0
l_cst_stat$Sample_ID<-rownames(l_cst_stat)
rownames(l_cst_stat)<-1:nrow(l_cst_stat)
setkey(counts_matrix_c,rn)
for (i in c(50,100,150,200,250,300)){
  ID_filt<-repmsk[(repmsk$V5-repmsk$V4)>i,"ID"][[1]]
  print(length(ID_filt))
  ID_filt<-ID_filt[ID_filt %in% counts_matrix_c$rn]
  cst_mx<-counts_matrix_c[ID_filt,]
  cst_stat<-wholespace_table_custom(paste0("_",unique(allmeta_loci$Tool)[-1]),cst_mx,sim01_loci_quant_stat,10,50)
  cst_stat$Cutoff<-i
  cst_stat$Sample_ID<-rownames(cst_stat)
  rownames(cst_stat)<-1:nrow(cst_stat)
  l_cst_stat<-rbind(l_cst_stat,cst_stat)
}

plot_f1fdr(l_cst_stat,"hsapiens_sim/stranded/","sim01_loci_quant_cutoff50_length",cutoff="length")







