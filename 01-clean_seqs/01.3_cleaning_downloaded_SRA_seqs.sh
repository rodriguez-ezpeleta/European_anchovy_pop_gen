#!/bin/bash

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
prj=/share/projects/ANE # project directory
log_folder=$prj/log_files # directory to save the output of each Stacks command/program
clean=$prj/clean # output directory

folder=$prj/SRA_download # folder were the downloaded 'fastq.gz' files are located
sample=(`ls ${folder}/*.fastq.gz | rev| cut -d "/" -f 1 | rev | sed 's/.fastq.gz//'`)
nsample=(`ls ${folder}/*.fastq.gz | rev| cut -d "/" -f 1 | rev | sed 's/.fastq.gz//' | wc -l`)
n=(`expr $nsample - 1`)

for i in `seq 0 $n`; do
	echo ${sample[i]} >> $prj/clean/samples_SRA.txt
	process_radtags -i gzfastq -f $folder/{sample[i]}.fastq.gz -o $prj/clean -e sbfI -c -q -r --score-limit 28 -t 95
	echo $(zcat $prj/clean/${sample[i]}.fq.gz | wc -l) / 4 | bc >> $prj/clean/Samples_reads_SRA.txt
done 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_process_radtags_SRA.log

paste $prj/clean/samples_SRA.txt $prj/clean/Samples_reads_SRA.txt > $prj/clean/Samples_reads_before_ustacks_SRA_"$TIMESTAMP".txt #info file
rm -rf $prj/clean/samples_SRA.txt $prj/clean/Samples_reads_SRA.txt
