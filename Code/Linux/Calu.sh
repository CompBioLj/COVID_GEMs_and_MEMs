#!/bin/bash

mkdir CALU
mkdir CALU/CALU_mock
mkdir CALU/CALU_infected


for i in `seq 44 49`;
do
	prefetch -v SRR115177${i};
	fastq-dump --outdir CALU/  SRR115177${i}/SRR115177${i}.sra;
	rm -fr SRR115177${i};
	
done 

declare -a avg_length
avg_length=(135.0 137.0 137.0)

declare -a std_deviation
std_deviation=(22.1 20.8 21.1)

for i in `seq 44  46`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o CALU/CALU_mock/SRR115177${i} --single -l ${avg_length[${i}-44]} -s  ${std_deviation[${i}-44]} CALU/SRR115177${i}.fastq
	rm CALU/SRR115177${i}.fastq;
done


declare -a avg_length_p
avg_length_p=(137.0 140.0 135.0)
declare -a std_deviation_p
std_deviation_p=(21.6 17.7 23.4)

for i in `seq 47  49`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o CALU/CALU_infected/SRR115177${i} --single -l ${avg_length_p[${i}-47]} -s ${std_deviation_p[${i}-47]} CALU/SRR115177${i}.fastq
	rm CALU/SRR115177${i}.fastq;
done

