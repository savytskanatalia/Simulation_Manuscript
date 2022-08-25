# Savytska Natalia, 2021

# Function for calculating enrichment with Fisher's test
# As input, it takes two frequency tables obtained with function table().
# TE_com_SubF_stat - tested table, hg38_subf_stat - reference or contrast table
# If TE entry is not found in reference (absent), the count for reference is 0, the total sum of all elements apart from ref - total sum of elements for ref 

enrich_te<-function(TE_com_SubF_stat,hg38_subf_stat){
  TE_com_SubF_stat$OR<-0
  TE_com_SubF_stat$pval<-1
  my_te_sum<-sum(TE_com_SubF_stat$Freq)
  hg38_te_sum<-sum(hg38_subf_stat$Freq)
  for (i in 1:nrow(TE_com_SubF_stat)){
    if (any(hg38_subf_stat[hg38_subf_stat$Var1==TE_com_SubF_stat$Var1[i],"Freq"])==TRUE){
      contab<-data.frame(TE=c(TE_com_SubF_stat[i,"Freq"],hg38_subf_stat[hg38_subf_stat$Var1==TE_com_SubF_stat[i,"Var1"],"Freq"]),
                       other=c((my_te_sum-TE_com_SubF_stat[i,"Freq"]),(hg38_te_sum-hg38_subf_stat[hg38_subf_stat$Var1==TE_com_SubF_stat[i,"Var1"],"Freq"])))
      TE_com_SubF_stat[i,"OR"]<-fisher.test(contab)$estimate
      TE_com_SubF_stat[i,"pval"]<-fisher.test(contab)$p.value
    }
    else {
      contab<-data.frame(TE=c(TE_com_SubF_stat[i,"Freq"],0),
                         other=c((my_te_sum-TE_com_SubF_stat[i,"Freq"]),hg38_te_sum))
      TE_com_SubF_stat[i,"OR"]<-fisher.test(contab)$estimate
      TE_com_SubF_stat[i,"pval"]<-fisher.test(contab)$p.value
    }
  }
  TE_com_SubF_stat$FDR<-p.adjust(TE_com_SubF_stat$pval,method="fdr")
  return(TE_com_SubF_stat)
}


# Example. I test for enrichment of specific TE subfamilies among False hits as compared to True hits.
# Load a table
sim01_str_FP5_loci <- read.delim("sim01_str_FP5_loci.txt")
# Retain only column with loci IDs
sim01_str_FP5_loci_ID<-unique(sim01_str_FP5_loci$TE)
# Get subfamilies for each element by removing "_dupN" identifier, as their base IDs are their subfamilies. Alternatively you could match those in full annotation from annotation table
sim01_str_FP5_loci_ID<-cbind(sim01_str_FP5_loci_ID,str_split_fixed(sim01_str_FP5_loci_ID,"_dup",n=2))
# save as df
sim01_str_FP5_loci_ID<-as.data.frame(sim01_str_FP5_loci_ID)
# count frequencies of the subfamilies (V2), save as a df
sim01_str_FP5_loci_ID_freq<-as.data.frame(table(sim01_str_FP5_loci_ID$V2))
# set TE IDs (Var1) as character variable
sim01_str_FP5_loci_ID_freq$Var1<-as.character(sim01_str_FP5_loci_ID_freq$Var1)


# Repeat same for True hits
sim01_str_TP5_loci <- read.delim("sim01_str_TP5_loci.txt")
sim01_str_TP5_loci_ID<-unique(sim01_str_TP5_loci$TE)
sim01_str_TP5_loci_ID<-as.data.frame(cbind(sim01_str_TP5_loci_ID,str_split_fixed(sim01_str_TP5_loci_ID,"_dup",n=2)))
sim01_str_TP5_loci_ID_freq<-as.data.frame(table(sim01_str_TP5_loci_ID$V2))
sim01_str_TP5_loci_ID_freq$Var1<-as.character(sim01_str_TP5_loci_ID_freq$Var1)


# Check enrichment with enrich_te()
TE_FP5_SubF_stat<-enrich_te(sim01_str_FP5_loci_ID_freq,sim01_str_TP5_loci_ID_freq)



# record number of simulated loci (so frequency in True hits table)
TE_FP5_SubF_stat$Simulated<-0
for (i in 1:nrow(TE_FP5_SubF_stat)){
  if (any(sim01_str_TP5_loci_ID_freq[sim01_str_TP5_loci_ID_freq$Var1==TE_FP5_SubF_stat[i,"Var1"],"Freq"])==TRUE){
    TE_FP5_SubF_stat[i,"Simulated"]<-sim01_str_TP5_loci_ID_freq[sim01_str_TP5_loci_ID_freq$Var1==TE_FP5_SubF_stat[i,"Var1"],"Freq"]
  } else {
    TE_FP5_SubF_stat[i,"Simulated"]<-0
  }
}

# Save the table with a function of choice...
# fwrite(TE_com_SubF_stat,'Enrichment.txt',sep="\t",row.names = TRUE)
