# Load the list of active elements
mm9_active_TE <- read.table("mm9_active_TE.txt", quote="\"", comment.char="")
mm9_active_TE$V1<-as.character(mm9_active_TE$V1)



# Load GTF of "active elements" to clean it up for truncated LTR and LINEs


mm9_activeTE_s <- read.delim("mm9_activeTE_s.bed", header=FALSE)
mm9_activeTE_s$V5<-mm9_activeTE_s$V3-mm9_activeTE_s$V2

#gene_id "MLT1E2"; transcript_id "MLT1E2_dup79"; family_id "MaLR"; class_id "LTR";


mm9_activeTE_s<-separate(mm9_activeTE_s, V4, into=c("gene","transcript","family","class"), sep="; ", extra="merge")
mm9_activeTE_s$class<-str_remove(mm9_activeTE_s$class,"class_id")
mm9_activeTE_s$class<-str_remove(mm9_activeTE_s$class,";")
mm9_activeTE_s$class<-str_remove(mm9_activeTE_s$class, "\"")
mm9_activeTE_s$class<-str_remove(mm9_activeTE_s$class, "\"")


## Extract either SINEs with over 200 bp length or non-sines over 900)
n_mm9_TE<-mm9_activeTE_s[ mm9_activeTE_s$class == "SINE" & mm9_activeTE_s$V5>200 | mm9_activeTE_s$class != "SINE" & mm9_activeTE_s$V5>900, ]
n_mm9_TE$V1<-as.character(n_mm9_TE$V1)
n_mm9_TE<-n_mm9_TE[!(n_mm9_TE$V1 %in% c("chr13_random","chr16_random","chr17_random","chr1_random","chr3_random","chr4_random","chr5_random","chr7_random","chr8_random","chr9_random","chrUn_random","chrX_random","chrY_random")), ]
## Save, prepare bed, prepare the sampler of 2000 elements
write.table(n_mm9_TE,file="long_active_TE_mm9.txt", row.names = FALSE, col.names=FALSE, quote=FALSE, sep = "\t")

n_mm9_TE_bed<-n_mm9_TE
n_mm9_TE_bed$gene<-str_remove(n_mm9_TE_bed$gene,"gene_id ")
n_mm9_TE_bed$gene<-str_remove(n_mm9_TE_bed$gene,"\"")
n_mm9_TE_bed$gene<-str_remove(n_mm9_TE_bed$gene,"\"")



####################################
# Obtaining gene list and respective tables


# List of mm9 transcripts
genes_mm9 <- read.delim("Mus_musculus.NCBIM37.67.cdna.all.fa.fai", header=FALSE, stringsAsFactors=FALSE)

genes_mm9$V1<-as.character(genes_mm9$V1)

# List of brain expressed genes
list_brain_genes <- read.table("list_of_expressed_inmouse_brain_genes.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
list_brain_genes<-separate(list_brain_genes, V1, into=c("gene","rest"), sep="\\.",extra="merge")
#ntab<- genes_mm9[genes_mm9$V1 %in% list_brain_genes$gene,]
list_brain_genes<-list_brain_genes[,1]






# Getting only longest transcript for each possible gene

# Get Transcript to gene mapping for the cDNA .fa

ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",mirror = "useast")
mart <- useMart(biomart = "ENSEMBL_MART_MOUSE", dataset = "mmusculus_gene_ensembl")
res1 <- getBM(attributes = c('ensembl_transcript_id',
                            'ensembl_gene_id',
                            'transcript_length'),
             filters = 'ensembl_transcript_id', 
             values = genes_mm9$V1, 
             mart = ensembl)

# Retaining only single longest isoform per gene
res1.agg<-aggregate(transcript_length ~ ensembl_gene_id, res1, max)
res1<-merge(res1.agg,res1)

# Retaining only "genes for brain dataset"
res1<-res1[res1$ensembl_gene_id %in% list_brain_genes,]

write.table(res1,file="brain_gene_transcript_id.txt", row.names = FALSE, col.names=TRUE, quote=FALSE, sep = "\t")
