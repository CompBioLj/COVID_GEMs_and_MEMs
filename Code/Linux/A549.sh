#!/bin/bash

mkdir A549
mkdir A549/A549_mock
mkdir A549/A549_infected


for i in `seq 39 62`;
do
	prefetch -v SRR114122${i};
	fastq-dump --outdir A549/  SRR114122${i}/SRR114122${i}.sra;
	rm -fr SRR114122${i};
	
done 

declare -a avg_length
avg_length=(145.0 145.0 145.0 145.0 145.0 145.0 145.0 145.0 146.0 146.0 146.0 146.0)

declare -a std_deviation
std_deviation=(13.2 13.1 13.1 13.0 12.4 12.4 12.3 12.3 12.4 12.4 12.6 12.5)

for i in `seq 39  50`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o A549/A549_mock/SRR114122${i} --single -l ${avg_length[${i}-39]} -s  ${std_deviation[${i}-39]} A549/SRR114122${i}.fastq
	rm A549/SRR114122${i}.fastq;
done


declare -a avg_length_p
avg_length_p=(146.0 146.0 146.0 146.0 144.0 144.0 144.0 144.0 146.0 146.0 146.0 146.0)
declare -a std_deviation_p
std_deviation_p=(12.6 12.6 12.6 12.5 15.3 15.2 15.3 15.2 12.3 12.3 12.2 12.2)

for i in `seq 51  62`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o A549/A549_infected/SRR114122${i} --single -l ${avg_length_p[${i}-51]} -s ${std_deviation_p[${i}-51]} A549/SRR114122${i}.fastq
	rm A549/SRR114122${i}.fastq;
done
