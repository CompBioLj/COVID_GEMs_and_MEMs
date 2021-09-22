#!/bin/bash

mkdir NHBE
mkdir NHBE/NHBE_mock
mkdir NHBE/NHBE_infected


for i in `seq 15 38`;
do
	prefetch -v SRR114122${i};
	fastq-dump --outdir NHBE/  SRR114122${i}/SRR114122${i}.sra;
	rm -fr SRR114122${i};
	
done 

declare -a avg_length
avg_length=(131.0 131.0 131.0 131.0 128.0 128.0 128.0 128.0 120.0 120.0 120.0 120.0)

declare -a std_deviation
std_deviation=(21.1 21.0 21.0 21.0 22.2 22.2 22.2 22.2 25.0 25.0 24.9 24.9)

for i in `seq 15  26`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o NHBE/NHBE_mock/SRR114122${i} --single -l ${avg_length[${i}-15]} -s  ${std_deviation[${i}-15]} NHBE/SRR114122${i}.fastq
done


declare -a avg_length_p
avg_length_p=(135.0 135.0 135.0 135.0 134.0 134.0 134.0 133.0 132.0 132.0 132.0 132.0)
declare -a std_deviation_p
std_deviation_p=(20.2 20.2 20.2 20.2 20.7 20.6 20.6 20.6 20.6 20.6 20.5 20.5)

for i in `seq 27  38`
do 
	
	kallisto quant -i Index/Homo_sapien.idx -o NHBE/NHBE_infected/SRR114122${i} --single -l ${avg_length_p[${i}-27]} -s ${std_deviation_p[${i}-27]} NHBE/SRR114122${i}.fastq

done
