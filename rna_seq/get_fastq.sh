#!/bin/bash

sras=`tail -n +2 info/sra_run.tsv | cut -f1 | tr -d '"'`

for sra in $sras
do
root=${sra:0:6}
id=${sra:9:2}
ftp="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$root/0$id/$sra/$sra.fastq.gz"
wget -P fastq $ftp
done


