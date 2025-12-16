#!/bin/bash


TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
prj=/share/projects/ANE # project directory
log_folder=$prj/log_files # directory to save the output of each Stacks command/program
clean=$prj/clean # output directory
file_rutas=$prj/scripts/ruta_pools_ANE01_to_ANE19.txt
info_pools=$prj/scripts/info_pools_ANE01_to_ANE19.txt

npools=(`grep -v "^#" $file_rutas | wc -l`)


# PROCESS_RADTAGS.
for f in `seq 1 $npools`; do 
	pool=(`grep -v "^#" $file_rutas | head -n $f | tail -n 1 | cut -f 2`)
	R1pool=(`grep -v "^#" $file_rutas | head -n $f | tail -n 1 | cut -f 3`)
	R2pool=(`grep -v "^#" $file_rutas | head -n $f | tail -n 1 | cut -f 4`)
	grep -w $pool $info_pools | cut -f 2 > $clean/${pool}_bars.txt
	grep -w $pool $info_pools | cut -f 1 > $clean/${pool}_names.txt
	paste $clean/${pool}_bars.txt $clean/${pool}_names.txt > $clean/${pool}_barcodes.txt
	rm -rf $clean/${pool}_bars.txt $clean/${pool}_names.txt
	mkdir $clean/${pool}
	process_radtags -i gzfastq -1 $R1pool -2 $R2pool -o $clean/$pool/ -e sbfI -b $clean/${pool}_barcodes.txt -c -q -r --score-limit 28 -t 95 --adapter-1 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG --adapter-2 GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
done 2>&1 | tee $log_folder/"$TIMESTAMP"_stacks_process_radtags_t95bp_remadap_AZTI_pools_01to19.log

grep "ANE" $info_pools | cut -f 1 > $prj/scripts/anchovy_samples_pools_01to19.txt
