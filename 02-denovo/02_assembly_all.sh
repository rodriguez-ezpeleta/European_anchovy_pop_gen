
#!/bin/bash

TIMESTAMP=$(date +%Y-%m-%d)
prj=/share/projects/ANE # project directory
log_folder=$prj/log_files
clean=$prj/clean # input directory

ls ${clean}/*.2.fq.gz | rev| cut -d "/" -f 1 | rev |  sed 's/\.2\.fq\.gz//'> ${clean}/AZTI_anchovy_names.txt  
ls ${clean}/SRR*.fq.gz | rev| cut -d "/" -f 1 | rev | sed 's/\.fq\.gz//'> ${clean}/SRA_anchovy_names.txt 
cat ${clean}/AZTI_anchovy_names.txt ${clean}/SRA_anchovy_names.txt > ${clean}/all_anchovy_names.txt
sample=(`cat ${clean}/all_anchovy_names.txt`)
nsample=(`cat ${clean}/all_anchovy_names.txt | wc -l`)

### 1. ustacks ###

#Build loci de novo in each sample.
m=5
M=2
n=3
mkdir $prj/denovo_m${m}M${M}n${n}
out=$prj/denovo_m${m}M${M}n${n}

id=1
s=(`expr $nsample - 1`)
for i in `seq 0 $s`; do
	echo -e "\n\n##### USTACKS: Treating individual $id: ${sample[i]}"
	ustacks -f $clean/${sample[i]}.fq.gz -o $out -i $id -m $m -M $M --disable-gapped -p 24
	let "id+=1"
done 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_ustacks.log

#Write a table with ustacks output info. 
cd $out
echo "Names   RetainedReads   Primary   Secondary   Tags   Alleles SNPs" > $out/ustacks_stats_"$TIMESTAMP".txt

for i in `seq 0 $s`; do
	echo -n ${sample[i]}
	echo ${sample[i]} >> $out/names_"$TIMESTAMP".txt
        echo -n ${sample[i]} >> $out/ustacks_stats_"$TIMESTAMP".txt
        echo -n "       "  >> $out/ustacks_stats_"$TIMESTAMP".txt
        use="$(zcat $out/${sample[i]}.tags.tsv.gz | grep -c 'primary\|secondary')"
        echo -n $use >> $out/ustacks_stats_"$TIMESTAMP".txt
        echo -n "       "  >> $out/ustacks_stats_"$TIMESTAMP".txt
	primary="$(zcat $out/${sample[i]}.tags.tsv.gz | grep -c primary)"
	echo -n $primary >> $out/ustacks_stats_"$TIMESTAMP".txt
	echo -n "       "  >> $out/ustacks_stats_"$TIMESTAMP".txt
	secondary="$(zcat $out/${sample[i]}.tags.tsv.gz | grep -c secondary)"
        echo -n $secondary >> $out/ustacks_stats_"$TIMESTAMP".txt
        echo -n "       "  >> $out/ustacks_stats_"$TIMESTAMP".txt
        tag="$(zcat $out/${sample[i]}.tags.tsv.gz | grep -c consensus)"
        echo -n $tag >> $out/ustacks_stats_"$TIMESTAMP".txt
        echo -n "       "  >> $out/ustacks_stats_"$TIMESTAMP".txt
        all="$(zcat $out/${sample[i]}.alleles.tsv.gz | wc -l)"
        echo -n $all >> $out/ustacks_stats_"$TIMESTAMP".txt
        echo -n "       " >> $out/ustacks_stats_"$TIMESTAMP".txt
        snp="$(zcat $out/${sample[i]}.snps.tsv.gz | grep -c 'E')"
        echo $snp  >> $out/ustacks_stats_"$TIMESTAMP".txt
done

#Coverage info
echo "Names   Coverage" > $out/coverage_titles.txt
grep "Final coverage" $log_folder/"$TIMESTAMP"_stacks_ustacks.log | grep -oE "mean=[0-9.]+" | sed -E 's/mean=//' > $out/coverage.txt
paste $out/names_"$TIMESTAMP".txt $out/coverage.txt > $out/ustacks_final_coverage.txt
cat $out/coverage_titles.txt $out/ustacks_final_coverage.txt > $out/ustacks_coverage_results_"$TIMESTAMP".txt
rm $out/coverage.txt $out/coverage_titles.txt $out/ustacks_final_coverage.txt

### 2. cstacks ###

# Build the catalog of loci from the samples contained in the population map.
#ls *alleles* | sed 's/\.alleles\.tsv\.gz//' > $prj/scripts/popmap_all_inds.txt
#cut -c 1-3 $prj/scripts/popmap_all_inds.txt > $prj/scripts/popmap_all_pops.txt
#paste $prj/scripts/popmap_all_inds.txt $prj/scripts/popmap_all_pops.txt > $prj/scripts/popmap_all.txt
#rm -rf $prj/scripts/popmap_all_inds.txt $prj/scripts/popmap_all_pops.txt

#To build the catalog from a subset of individuals, supply a separate population map only containing those samples.
awk '$5>25000' $out/ustacks_stats_"$TIMESTAMP".txt | grep -v "^Names" | cut -d " " -f 1  | cut -c1-3 > $prj/scripts/pops.txt
awk '$5>25000' $out/ustacks_stats_"$TIMESTAMP".txt | grep -v "^Names" | cut -d " " -f 1  > $prj/scripts/inds.txt
paste $prj/scripts/inds.txt $prj/scripts/pops.txt > $prj/scripts/popmap_m5M2n3_25000Tags_ALL_502M.txt
rm $prj/scripts/inds.txt $prj/scripts/pops.txt

popmap=$prj/scripts/popmap_m5M2n3_25000Tags_ALL_502M.txt
cstacks -P $out -M $popmap -n $n -p 24 --disable-gapped 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_cstacks.log

### 3. sstacks ###

#Match all samples supplied in the population map against the catalog.
sstacks -P $out -M $popmap -p 24 --disable-gapped 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_sstacks.log

#info
cd $out
for f in *matches.tsv.gz ;  do 
	echo -n ${f%.matches.tsv.gz} ; echo -n " " ; zcat $f | cut -f 3 |sort | uniq |wc -l 
done > cstacks_stats_"$TIMESTAMP".txt
cd ..

### 4. tsv2bam ### 

#Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include paired-end reads using tsv2bam. 
#Samples with R1 and R2
popmap=$prj/scripts/popmap_m5M2n3_25000Tags_AZTI_375M.txt
tsv2bam -P $out -M $popmap -t 24 -R $clean 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_tsv2bam_AZTI.log
#Samples with only R1
popmap=$prj/scripts/popmap_m5M2n3_25000Tags_SRA_127M.txt
tsv2bam -P $out -M $popmap -t 24 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_tsv2bam_SRA.log


### 5. gstacks ###

#Build a paired-end contig from the metapopulation data (if paired-reads provided),align reads per sample, call variant sites in the population, genotypes in each individual.
popmap=$prj/scripts/popmap_m5M2n3_25000Tags.txt
gstacks -P $out -M $popmap -t 8 --ignore-pe-reads 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_gstacks_onlyR1.log

### 6. populations ###

popmap=$prj/scripts/popmap_m5M2n3_25000Tags_ALL_502M.txt
outpops=/share/projects/ANE/denovo_m${m}M${M}n${n}/populations_r75
mkdir $outpops
populations -P $out -M $popmap -O $outpops -t 24 -R 0.75 --plink --vcf 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_populations.log

#set different popmaps for analysis of subsets. 
