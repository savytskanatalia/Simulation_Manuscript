### POLYESTER
library(polyester)
library(seqinr)
library(tidyr)
library(stringr)
library("biomaRt")

# Load the file "long_active_TE_mm9.txt"
library(optparse)
library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
library(polyester)


option_list = list(
  make_option(c("-t", "--number_of_te"), type="numeric", default=2000, 
              help="Number of the transposable element loci to be expressed in total, [default= %default]", metavar="numeric"),
  make_option(c("-f", "--number_of_fam"), type="numeric", default=6, 
              help="Number of subfamilies per class (LINE,SINE,LTR) to be differentially expressed (up, down, in both directions), [default= %default]. Not recommended to take more than 9 (resulting in 27 DE subfamilies per calss)", metavar="numeric"),
  make_option(c("-g", "--number_of_genes"), type="numeric", default=2000, 
              help="Number of the genes to be differentially expressed, [default= %default]", metavar="numeric"),
  make_option(c("-1", "--te_file"), type="character", default="long_active_TE_mm9.txt", 
              help="Mice active transpsable elements table, includes locations and IDs; input file name", metavar="character"),
  make_option(c("-2", "--gene_file"), type="character", default="brain_gene_transcript_id.txt", 
              help="List of longest transcripts (one per gene) of the genes active in mice brain; file name [default= %default]", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


active_te<-read.table(file=opt$te_file,stringsAsFactors=FALSE)

# Hard-coded. List of active TE elements in mice based on published study Molaro, et al (2014)

te_ltr<- c("MTD-int" , "ORR1B1-int" , "ORR1A2-int" , "RMER19B" , "IAP-d-int" , "IAPEz-int" , "MTA_Mm-int" , "RLTR10-int" , "ETnERV3-int" , "ORR1A0-int" , "ORR1A4-int" , "ORR1A3-int" , "IAPEY3-int" , "RMER17C-int" , 
          "MTC-int" ,  "ORR1A1-int" , "MMETn-int" , "IAPEy-int" ,  "ORR1B2-int" , "ETnERV2-int" , "IAPLTR3-int" , "MTEa-int" , "RMER15-int" , "RMER17A2" , "RLTR13B3" , "ORR1D1-int" , "RMER17B" , "ORR1E-int" , 
          "MTB_Mm-int" , "ORR1C1-int" , "MTE-int" , "RMER6-int" , "MTB-int" , "RLTR1B-int" , "RMER17D2" , "ORR1A2" , "ORR1C2-int" , "ORR1D2-int" , "MTE2b-int" , "ERVL-B4-int" , "MTEb-int" , "ETnERV-int" , 
          "RMER6C" , "IAPLTR4_I" , "ORR1D-int" , "MLT1H-int" , "MLT1A1-int" , "MTE2a-int" , "MMTV-int" , "RMER17A-int" , "RMER17D" , "RMER19A" , "RMER6D" , "MLT1D-int" , "MLT1A0-int" , "MT-int" , 
          "MLT1F-int" , "RMER6A" , "MLT1A-int" , "MLT1-int" , "MLT1J-int" ,  "MLT1G3-int" , "MLT1H2-int" , "MLT1E1-int" , "MLT1F2-int" , "MLT1F1-int" ,  "RMER17A" , "MLT1C-int" , "MLT1B-int" , "RLTR13B4" , 
          "MLT1G1-int" , "MLT1G-int" , "MTD"  , "MLT1I-int" , "RMER19C" , "MLT1E1A-int")

te_sine<- c("ID_B1" , "B4A" , "MIRb" , "B4"  ,  "B3" , "B3A"  ,  "MIR" ,  "MIRc" ,  "RSINE1" , "B2_Mm1t" , "B2_Mm2" , "B2_Mm1a" , "MIRm"  , "B1_Mur4" , "MIR3"  , "B1_Mus1" , "B1_Mus2" , "B1_Mur3" , "B1_Mm" , "B1F1" , 
             "PB1D10" , "B1F" , "B1F2" , "PB1D9" , "B1_Mur1" , "B1_Mur2")

te_line<- c("L1Md_A" , "Lx9" , "L1Md_T" , "L1_Mus2" , "Lx7" , "Lx6" , "Lx4A" , "L1Md_F2" , "L1M2" , "L1_Mus3" , "L1_Mur3" , "L1Md_F3" , "L1VL1" , "L1_Mm" , "L1Md_Gf" , "Lx3B" , "Lx2B" , "L1VL4" , "L1_Mus1" , "Lx2" , 
            "Lx3C" , "L1_Rod" , "Lx8" , "L1_Mur2" , "Lx3_Mus" , "Lx" , "L1_Mur1" , "L1Md_F" , "Lx4B" , "Lx5" , "L1VL2" , "L1_Mus4" , "MusHAL1" , "L1MA6" , "Lx3A" , "L1MC" , "L1MD" , "Lx2A" , "L1MA4" , "L1MB1" , 
            "L1ME3B" , "Lx2A1" , "L1MC3" , "L1MA9" , "L1MD1" , "L1MB4" , "L1MC1" , "L1M4" ,  "L1MB7" , "L2"  ,  "L1MDa" ,  "L2a"  ,  "L1MEf" , "L1MA5A" , "L1M4c" ,  "L1M3"  ,  "L1M3e" , "L1MA4A" , "L1MCa" , "L1MA5" , 
            "HAL1" ,  "L1MA10" , "L1MEc" , "L1M5" , "L1MA8" ,  "L1MA7" ,  "L1MC4" , "L1M3a"  ,  "L1ME5" , "L1MC2"  , "L1ME1" , "L1MB5"  , "L1MB2"  , "L1MB8" ,  "L1ME3A" , "L1MEg"  , "L1M4b" ,  "L1MB3" , "L1MEe" ,  "L1MCb" , 
            "L1MD2" ,  "L1M3c" , "L1PB4" ,  "L1ME3" ,  "L1MC4a" , "L1M3f" ,  "L1M3d" ,  "L1MEd" , "L1ME2" , "L1M2c"  ,  "L2c" ,  "L1ME4a" , "L1MD3"  , "L1M2a" , "L1MC5" , "L1ME2z" , "L1MCc"  ,  "L1M6" , "L1M3b" , "HAL1b" , 
            "L1MEb" , "L1MEg1" , "L1MDb" , "L1MEg2" , "L1M3de" ,  "L2b" , "L1MEa")

# Sampling by number of families given by input or default

# Sampling LINEs
te_line_up<-sample(te_line,opt$number_of_fam)
te_line_down<-sample(te_line[!(te_line %in% te_line_up)],opt$number_of_fam)
te_line_both<-sample(te_line[!(te_line %in% c(te_line_up,te_line_down))],opt$number_of_fam)

# Sampling SINEs
te_sine_up<-sample(te_sine,opt$number_of_fam)
te_sine_down<-sample(te_sine[!(te_sine %in% te_sine_up)],opt$number_of_fam)
te_sine_both<-sample(te_sine[!(te_sine %in% c(te_sine_up,te_sine_down))],opt$number_of_fam)

# Sampling LTRs
te_ltr_up<-sample(te_ltr,opt$number_of_fam)
te_ltr_down<-sample(te_ltr[!(te_line %in% te_ltr_up)],opt$number_of_fam)
te_ltr_both<-sample(te_ltr[!(te_ltr %in% c(te_ltr_up,te_ltr_down))],opt$number_of_fam)


# Joining UP, DOWN, MIXED into common lists

te_up<-c(te_line_up,te_sine_up,te_ltr_up)
te_down<-c(te_line_down,te_sine_down,te_ltr_down)
te_both<-c(te_line_both,te_sine_both,te_ltr_both)




# Sampling 6000 loci
active_te_sampled<-active_te[sample(nrow(active_te), opt$number_of_te),]


# Creating DF with fold change values

count_TE<-data.frame(active_te_sampled$V7, active_te_sampled$V5, active_te_sampled$V2, active_te_sampled$V2)
colnames(count_TE)<-c("ID","Subfamily","CTRL","DIS")
count_TE$CTRL<-0
count_TE$DIS<-0





# Populating FC-like DF with values



n=0    #Setting a variable number, that will tell the for loop how to populate the TE var subfamilies
for (i in 1:nrow(count_TE)) {
  if (count_TE$Subfamily[i] %in% te_up) {
    count_TE$CTRL[i]=1
    count_TE$DIS[i]=round(runif(1,5,16))
  } else if (count_TE$Subfamily[i] %in% te_down) {
    count_TE$CTRL[i]=round(runif(1,5,16))
    count_TE$DIS[i]=1
  } else if (count_TE$Subfamily[i] %in% te_both) {
    if (n==0) {
      count_TE$CTRL[i]=round(runif(1,5,16))
      count_TE$DIS[i]=1
      n=1
    } else {
      count_TE$CTRL[i]=1
      count_TE$DIS[i]=round(runif(1,5,16)) 
      n=0
    }
  } else {
    count_TE$CTRL[i]=1
    count_TE$DIS[i]=1
  }
}



gene_file<-opt$gene_file
active_genes<-read.table(file=gene_file,stringsAsFactors=FALSE, header=TRUE)


de_genes<-active_genes$ensembl_transcript_id[sample(nrow(active_genes),opt$number_of_genes)]
de_genes_up<-de_genes[1:round(opt$number_of_genes/2)]
de_genes_down<-de_genes[(round(opt$number_of_genes/2)+1):opt$number_of_genes]

## Creating Fold Change table for the genes (transcripts)
gene_fc<-active_genes[,c(3,2,2)]
gene_fc[,c(2,3)]<-0
colnames(gene_fc)<-c("ID","CTRL","DIS")

# Populating the FC table

for (i in 1:nrow(gene_fc)) {
  if (gene_fc$ID[i] %in% de_genes_up) {
    gene_fc$CTRL[i]=1
    gene_fc$DIS[i]=round(runif(1,5,16))
  } else if (gene_fc$ID[i] %in% de_genes_down) {
    gene_fc$CTRL[i]=round(runif(1,5,16))
    gene_fc$DIS[i]=1
  }  else {
    gene_fc$CTRL[i]=1
    gene_fc$DIS[i]=1
  }
}


# Creating common FC table for both genes and TE
total_fc<-rbind(gene_fc,count_TE[,c(1,3,4)])
rownames(total_fc)<-total_fc$ID












## Merging back ID and class for TE
#df$x <- paste0(df$n, "-", df$s)

active_te_sampled$TE <-paste0(active_te_sampled[,5], ";",active_te_sampled[,7], ";", active_te_sampled[,9], ";", active_te_sampled[,10])

te_bed_f<-active_te_sampled[,c(1,2,3,13,4,12)]
te_bed_f[,5]<-0
te_bed_f[,1]<-str_remove(te_bed_f[,1],"chr")


### Writing all the output files I need for simulation
#unnecessary, can be obtained with loaded .fa later+will be more accurate

write.table(total_fc,file="foldchange_for_sim.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep = "\t")

write.table(active_genes$ensembl_transcript_id,file="ENST_id_for_sim.txt", row.names = FALSE, col.names=FALSE, quote=FALSE, sep = "\t")

write.table(te_bed_f,file="TE_bed_for_sim.bed", row.names = FALSE, col.names=FALSE, quote=FALSE, sep = "\t")







fasta="total_simulation.fasta"
sim_fa<-read.fasta(file=fasta)
length(sim_fa)
readspertx<-round(getLength(sim_fa)/2)

fc_count="foldchange_for_sim.txt"
fold <- read.delim(fc_count, header=TRUE, stringsAsFactors=FALSE)

fold_m<-as.matrix(fold[,c(2,3)])

simulate_experiment(fasta=fasta, gtf=NULL, seqpath=NULL, outdir='/media/savytska/28789cbc-3aa5-41b9-80e9-4b6eed772cb7/savytska/Projects/Banchmarking_TE/GTF/01_necessary_files/out', num_reps=c(5,5), reads_per_transcript=readspertx, fold_changes=fold_m, gzip=TRUE)


# Mmusculus, stranded
# Get .fasta with TE+genic transcripts to simulate from 

fasta='Mmusculus_stranded/total_simulation.fasta'
sim_fa<-read.fasta(file=fasta)
length(sim_fa)
readspertx<-round(getLength(sim_fa)/2)

fc_count='Mmusculus_stranded/foldchange_for_sim.txt'
fold <- read.delim(fc_count, header=TRUE, stringsAsFactors=FALSE)

fold_m<-as.matrix(fold[,c(2,3)])
# gains equivalent results to the previous simulation in regards to library size and read sim - 2*readspertx and 0.5 library == readspertx and 1. library
simulate_experiment(fasta=fasta, gtf=NULL, seqpath=NULL, num_reps=c(5,5), reads_per_transcript=2*readspertx, fold_changes=fold_m, gzip=TRUE, paired=TRUE, seed=25, strand_specific=TRUE, gzip=TRUE,lib_sizes=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),outdir='Mmusculus_stranded/raw')
