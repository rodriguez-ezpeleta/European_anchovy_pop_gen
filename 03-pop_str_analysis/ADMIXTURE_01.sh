
#ADMIXTURE SCRIPT 1
file=outlier_snps_AZTI_SRA_L1L4_502M_anchovy_oneSNPperTAG_LD02
# output dir
outfold=/share/projects/ANE/denovo_m5M2n3/populations_r75_AZTI_SRA_L1L4_502M/ADMIXTURE/
# SNPs number
snp=(`wc -l ${file}.map`)
#give the SNPs a distance of 100000. 
snps=$(($snp*100000))

# FILE MODIFICATION
#Modify the ".map" file (to add the distance between SNPs).
seq 1 100000 $snps > maps_column.txt
cut -f 1-3 ${file}.map > new_map.txt
paste new_map.txt maps_column.txt > ${outfold}/${file}.map
rm -rf ./new_map.txt
rm -rf maps_column.txt
#Modify the ".ped" file.  To change the alleles from "ATGC" to the alleles as 1 and 2.
plink --noweb --file ${file} --recode12 --out ${file}


# CROSS VALIDATION (to set the correct value of "K" and "-c"). "-c" value is the number of the steps for admixture to run over each iteration. First run one iteration and explore the number of steps it needs to reach convergence, and then, use a value for c which ensure convergence for each step.
#mkdir cross_validation
folder="$outfold/cross_validation"
cd $folder
for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv ../$file.ped $K | tee $folder/log${K}_outlier.out; done
grep -h CV $folder/log*_outlier.out > $folder/cross_validation_outlier.txt
cd ..


#ADMIXTURE WITH BOOTSTRAPING
# Check the cross-validation results to set the "-c".
K=2
admixture -j16 ${file}.ped ${K} -B1000 -c14
K=3
admixture -j16 ${file}.ped ${K} -B1000 -c17
K=4
admixture -j16 ${file}.ped ${K} -B1000 -c21
K=5
admixture -j16 ${file}.ped ${K} -B1000 -c22
K=6
admixture -j16 ${file}.ped ${K} -B1000 -c24
K=7
admixture -j16 ${file}.ped ${K} -B1000 -c23

# Now we run the R scritp "ADMIXTURE_02_plot.R" to get the plot.





