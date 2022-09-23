# Get gene ID's and expression values from GTEx
cat '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct' | awk 'BEGIN {FS="\t";OFS="\t"};{print $1,$22}' > '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GTEx_v8_RNASeQCv1.1.9_Substantia_nigra_gene_median_tpm.gct'
# Delete first two rows in the file - those are comments. Proceed with Section 1 of Human.R script




# Merge transcripts .fa for lncRNAs and protein-coding genes.
cat '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/gencode.v38.pc_transcripts.fa' '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/gencode.v38.lncRNA_transcripts.fa' > '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome.fa'
# Index the total .fa to extract sequences easier later
cat /media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome.fa | awk -F _ '{FS="|";OFS="_"}; /^>/ {print $1,$2; next} 1' > /media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome_fin.fa

samtools faidx '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome_fin.fa'

# Subset only overlapping genes from GTEx and hg38 coding-non-coding transcriptome
while read i; do grep $i '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome_fin.fa' | sed 's/>//g' >> "/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/Substantia_nigra_HS.txt" ; done < "/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/Substantia_nigra_HS.ID"
samtools faidx '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/GRCh38_Coding_Transcriptome_fin.fa' -r "/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/Substantia_nigra_HS.txt" > "/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/Substantia_nigra_HS.fa"
 
 
 
 
 
 ## TE section
 # Get lists of regions of TEs of interest from SECTION 2 in Human.R script
 # For TEs on forward strand it`s ok to just exctract the sequences
samtools faidx '/media/savytska/natalia/Annotation/hg38/fasta/genome.fa' -r '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_ForwardStr.txt' > '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_ForwardStr.fa'
 # For TEs on reverse strand, one need to get reverse complement sequences
samtools faidx '/media/savytska/natalia/Annotation/hg38/fasta/genome.fa' -r '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_RevStr.txt' -i --mark-strand no > '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_RevStr.fa'
# Merge both files
cat '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_ForwardStr.fa' '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_RevStr.fa' > '/media/savytska/natalia/Projects/Simulations/GRCh38_HOMO_SAPIENS/TE/HS_TE_forSimulation_total.fa'
