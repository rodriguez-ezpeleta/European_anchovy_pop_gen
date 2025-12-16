# SNP filtering using plink #

#Check the missingness.
plink --file populations.plink --noweb --missing --out populations.plink --recode

#-Missingness per marker. Step to exclude SNPs on the basis of missing genotype rate. 
#cat populations.plink.lmiss | awk '$5 >0.1' | wc -l
plink --file populations.plink --noweb --geno 0.1 --out r75_geno10 --recode
plink --file r75_geno10 --noweb --missing --out r75_geno10_missing --recode


#-Missingness per individual. Step to exclude individuals with too much missing genotype data.
#cat r75_geno10_missing.imiss | awk '$6 >0.1' | wc -l
plink --file r75_geno10 --noweb --mind 0.1 --out r75_geno10_mind10 --recode
plink --file r75_geno10_mind10 --noweb --missing --out r75_geno05_mind10_missing --recode

#-Minor Allele frequency
plink --file r75_geno10_mind10 --noweb --maf 0.05 --out r75_geno10_mind10_maf05 --recode

#-Repeat to change some thresholds. 
plink --file r75_geno10_mind10_maf05 --noweb --geno 0.05 --out r75_geno10_mind10_maf05_geno05 --recode
plink --file r75_geno10_mind10_maf05_geno05 --noweb --mind 0.1 --out r75_geno10_mind10_maf05_geno05_mind10 --recode

#-Relatedness using vcf
plink19 --file r75_geno10_mind10_maf05_geno05_mind10 --noweb --recode vcf --out r75_geno10_mind10_maf05_geno05_mind10 # change from plink to vcf
vcftools --relatedness --vcf r75_geno10_mind10_maf05_geno05_mind10.vcf --out r75_geno10_mind10_maf05_geno05_mind10
plink --file r75_geno10_mind10_maf05_geno05_mind10 --noweb --remove indiv_to_remove_relatedness.txt --out r75_geno10_mind10_maf05_geno05_mind10_rel --recode

##  ONE SNP PER TAG. 
sed 's/ \+/\t/g' r75_geno10_mind10_maf05_geno05_mind10_rel.map | sed 's/_/\t/g' | sort -n -k 2,2 -u | sed 's/\t/_/2' | cut -f 2 > oneSNPperTAG_r75_geno10_mind10_maf05_geno05_mind10_rel.txt
plink --noweb --file r75_geno10_mind10_maf05_geno05_mind10_rel --extract oneSNPperTAG_r75_geno10_mind10_maf05_geno05_mind10_rel.txt --out r75_geno10_mind10_maf05_geno05_mind10_oneSNPperTAG_rel --recode
# another option: stacks "populations" ("--write-single-snp").

#-Linkage disequilibrium
#sed 's/0/un/1' r75_geno10_mind10_maf05_geno05_mind10_oneSNPperTAG_rel.map > all_9209SNP_LD.map
#cp r75_geno10_mind10_maf05_geno05_mind10_oneSNPperTAG_rel.ped all_9209SNP_LD.ped
#"r75_geno10_mind10_maf05_geno05_mind10_oneSNPperTAG_rel" = "all_9209SNP_LD".
plink --file all_9209SNP_LD --allow-extra-chr --indep-pairwise 9209 1 0.2 --out all_9209SNP_LD_results02 --recode
plink --file all_9209SNP_LD_results02 --allow-extra-chr --exclude all_9209SNP_LD_results02.prune.out --out plink_filtered_LD02_allSNPs --recode

File from "map/ped" to "structure".
./ped2struc.sh plink_filtered_LD02_allSNPs
	
	
	
