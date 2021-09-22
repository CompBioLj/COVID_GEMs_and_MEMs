#!/bin/bash

mkdir Index

wget /home/gemfinal/Index http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

kallisto index -i Index/Homo_sapien.idx Homo_sapiens.GRCh38.cdna.all.fa.gz



