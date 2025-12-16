#!/bin/bash

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
prj=/share/projects/ANE # project directory
log_folder=$prj/log_files # directory to save the output of each Stacks command/program
clean=$prj/clean 

mv $clean/ANE*/*.fq.gz $clean/
names=(`cat $prj/scripts/anchovy_samples_pools_ANE01_to_ANE19.txt`)
nsamples=(`wc -l $prj/scripts/anchovy_samples_pools_ANE01_to_ANE19.txt`)

n=(`expr $nsamples - 1`)
mkdir $clean/rem_reads

# CLONE FILTER. (for paired-end reads). 
for i in `seq 0 $n`; do
	echo ${names[i]} >> $clean/samples.txt
	echo $(zcat $clean/${names[i]}.1.fq.gz | wc -l) / 4 | bc >> $clean/Reads_pre_clone_filter.txt
	clone_filter -1 $clean/${names[i]}.1.fq.gz -2 $clean/${names[i]}.2.fq.gz -o $clean -i gzfastq
	echo $(zcat $clean/${names[i]}.1.1.fq.gz | wc -l) / 4 | bc >> $clean/Reads_after_clone_filter.txt
	echo $(zcat $clean/${names[i]}.rem.1.fq.gz | wc -l) / 4 | bc >> $clean/Reads_rem1.txt
	echo $(zcat $clean/${names[i]}.rem.2.fq.gz | wc -l) / 4 | bc >> $clean/Reads_rem2.txt
	mv $clean/${names[i]}.2.2.fq.gz $clean/${names[i]}_notclean.2.2.fq.gz
	#Remove the first base of read2 (T). #HEADCROP: Cut the specified number of bases from the start of the read.
	java -jar trimmomatic-0.38.jar SE -threads 16 $clean/${names[i]}_notclean.2.2.fq.gz $clean/${names[i]}.2.2.fq.gz HEADCROP:1 
	mv $clean/${names[i]}.1.fq.gz $clean/rem_reads/${names[i]}_preclonefilter.1.fq.gz
	mv $clean/${names[i]}.2.fq.gz $clean/rem_reads/${names[i]}_preclonefilter.2.fq.gz
	mv $clean/${names[i]}.rem.1.fq.gz $clean/rem_reads/
	mv $clean/${names[i]}.rem.2.fq.gz $clean/rem_reads/
	mv $clean/${names[i]}.1.1.fq.gz $clean/${names[i]}.fq.gz
	mv $clean/${names[i]}.2.2.fq.gz $clean/${names[i]}.2.fq.gz
	rm -rf $clean/${names[i]}_notclean.*
	rm -rf $clean/rem_reads/${names[i]}_preclonefilter_*
done 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_clonefilter_pools_ANE01_to_ANE19.log

#Create a table with cleaning info. 
paste $clean/samples.txt $clean/Reads_pre_clone_filter.txt  $clean/Reads_after_clone_filter.txt $clean/Reads_rem1.txt $clean/Reads_rem2.txt  | awk '{print $1, $2, $3, $3*100 / $2, $4, $5}' > $clean/Samples_reads_table.txt
echo "Names   Reads_pre_clonefilter   Reads_after_clonefilter   Percent_no_clones   Rem1   Rem2" > $clean/Samples_reads_header.txt
cat $clean/Samples_reads_header.txt $clean/Samples_reads_table.txt > $clean/Samples_reads_before_ustacks_"$TIMESTAMP".txt #info file
rm -rf $clean/Samples_reads_table.txt $clean/Samples_reads_header.txt $clean/Reads_pre_clone_filter.txt $clean/Reads_after_clone_filter.txt $clean/Reads_rem1.txt $clean/Reads_rem2.txt $clean/samples.txt
