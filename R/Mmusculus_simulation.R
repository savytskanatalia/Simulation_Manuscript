
### POLYESTER
library(polyester)
library(seqinr)


# Mmusculus, stranded
# Get .fasta with TE+genic transcripts to simulate from 

fasta='Mmusculus_stranded/total_simulation.fasta'
sim_fa<-read.fasta(file=fasta)
length(sim_fa)
readspertx<-round(getLength(sim_fa)/2)

fc_count='Mmusculus_stranded/foldchange_for_sim.txt'
fold <- read.delim(fc_count, header=TRUE, stringsAsFactors=FALSE)

fold_m<-as.matrix(fold[,c(2,3)])

simulate_experiment(fasta=fasta, gtf=NULL, seqpath=NULL, num_reps=c(5,5), reads_per_transcript=2*readspertx, fold_changes=fold_m, gzip=TRUE, paired=TRUE, seed=25, strand_specific=TRUE, gzip=TRUE,lib_sizes=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),outdir='Mmusculus_stranded/raw')
