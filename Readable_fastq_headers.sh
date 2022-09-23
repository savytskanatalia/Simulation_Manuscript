
# Edit header to resemble SRA headers (needed for SalmonTE)
awk '{FS="";OFS="_"};/^>/{print ">SSRR02", substr($0,2);next}{print}' | sed 's/\//\t/g' > fastq/sample_02_1.fasta < sample_02_1.fasta
# .fa to .fq conversion
perl fasta_to_fastq.pl sample_02_1.fasta > sample_02_1.fastq
