#!/bin/bash

mkdir 293T
mkdir 293T/293T_mock
mkdir 293T/293T_infected

#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789547/293_24hr_control_1_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789547/293_24hr_control_1_fq2.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789548/293_24hr_control_2_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789548/293_24hr_control_2_fq2.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789549/293_24hr_control_3_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789549/293_24hr_control_3_fq2.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789556/293_24hr_COVID_1_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789556/293_24hr_COVID_1_fq2.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789557/293_24hr_COVID_2_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789557/293_24hr_COVID_2_fq2.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789558/293_24hr_COVID_3_fq1.fq.gz.1
#wget https://storage.googleapis.com/nih-sequence-read-archive/sra-src/SRR12789558/293_24hr_COVID_3_fq2.fq.gz.1


#for i in `seq 47  49`
#do
#	for j in `seq 1 2`
#	do
#		wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/0${i}/SRR127895${i}/SRR127895${i}_${j}.fastq.gz
#	done
#done
	 
#for i in `seq 56  58`
#do
#	for j in `seq 1 2`
#	do
#		wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/0${i}/SRR127895${i}/SRR127895${i}_${j}.fastq.gz
#	done
#done	


for i in `seq 47  49`
do
	
	kallisto quant -i Index/Homo_sapien.idx -o 293T/293T_mock/SRR127895${i} SRR127895${i}_1.fastq.gz SRR127895${i}_2.fastq.gz
	#rm 293T/SRR127895${i}.fastq;
done

for i in `seq 56  58`
do
	

	
	kallisto quant -i Index/Homo_sapien.idx -o 293T/293T_mock/SRR127895${i} SRR127895${i}_1.fastq.gz SRR127895${i}_2.fastq.gz
	#rm 293T/SRR127895${i}.fastq;
	
done



#for i in `seq 51  62`
#do 
	
#	kallisto quant -i Index/Homo_sapien.idx -o 293T/293T_infected/SRR127895${i} --single -l ${avg_length_p[${i}-51]} -s ${std_deviation_p[${i}-51]} 293T/SRR127895${i}.fastq
#	rm 293T/SRR127895${i}.fastq;
#done
