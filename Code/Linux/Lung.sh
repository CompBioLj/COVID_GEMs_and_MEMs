#!/bin/bash

mkdir LUNG
mkdir LUNG/LUNG_mock
mkdir LUNG/LUNG_infected


for i in `seq 25  40`;
do
	prefetch -v SRR115177${i};
	fastq-dump --outdir LUNG/  SRR115177${i}/SRR115177${i}.sra;
	rm -fr SRR115177${i};
	
done 

declare -a avg_length
avg_length=(133.0 133.0 133.0 133.0 134.0 134.0 134.0 134.0)

declare -a std_deviation
std_deviation=(23.5 23.5 23.5 23.5 21.9 21.8 21.8 21.8)

for i in `seq 25  32`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o LUNG/LUNG_mock/SRR115177${i} --single -l ${avg_length[${i}-25]} -s  ${std_deviation[${i}-25]} LUNG/SRR115177${i}.fastq
	rm LUNG/SRR115177${i}.fastq;
done


declare -a avg_length_p
avg_length_p=(141.0 141.0 141.0 141.0)
declare -a std_deviation_p
std_deviation_p=(22.0 22.0 22.1 22.1)

for i in `seq 37  40`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o LUNG/LUNG_infected/SRR115177${i} --single -l ${avg_length_p[${i}-37]} -s ${std_deviation_p[${i}-37]} LUNG/SRR115177${i}.fastq
	rm LUNG/SRR115177${i}.fastq;
done
