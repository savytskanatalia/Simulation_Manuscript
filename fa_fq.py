# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 12:41:43 2019

@author: savytska
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-


import subprocess

for i in range(1,9,1):
	print('started copying')
#	subprocess.call("cp sample_0"+str(i)+"_1.fasta.gz fastq/", shell=True)
	print('finished copying')
	subprocess.call("gzip -d sample_0"+str(i)+"_1.fasta.gz", shell=True)
	print('successfully unzipped')
	com_list = ["awk '",'{FS="";OFS="_"};/^>/{print ', '">SSRR0'+str(i), '",substr($0,2); next}', "{print}'> fastq/sample_0"+str(i)+"_1.fasta", "< sample_0"+str(i)+"_1.fasta"]
	com_list2 = ["awk '",'{FS="";OFS="_"};/^>/{print ', '">SSRR0'+str(i), '",substr($0,2); next}', "{print}'> fastq/sample_0"+str(i)+"_2.fasta", "< sample_0"+str(i)+"_2.fasta"]
	com_com = " ".join(com_list)
	print(com_com)
	com_com2 = " ".join(com_list2)
	print('starting fasta modification')
	subprocess.call(com_com, shell=True)
	print('starting fasta to fastq conversion')
	subprocess.call("perl /data/scripts/fasta_to_fastq.pl fastq/sample_0"+str(i)+"_1.fasta > fastq/sample_0"+str(i)+"_1.fastq", shell=True)
	subprocess.call("gzip fastq/sample_0"+str(i)+"_1.fastq", shell=True)
#	subprocess.call("rm fastq/sample_0"+str(i)+"_1.fasta", shell=True)
#	subprocess.call("cp sample_0"+str(i)+"_2.fasta.gz fastq/", shell=True)
	subprocess.call("gzip -d sample_0"+str(i)+"_2.fasta.gz", shell=True)
	subprocess.call(com_com2, shell=True)
	subprocess.call("perl /data/scripts/fasta_to_fastq.pl fastq/sample_0"+str(i)+"_2.fasta > fastq/sample_0"+str(i)+"_2.fastq", shell=True)
	subprocess.call("gzip fastq/sample_0"+str(i)+"_2.fastq", shell=True)
#	subprocess.call("rm fastq/sample_0"+str(i)+"_2.fasta", shell=True)
