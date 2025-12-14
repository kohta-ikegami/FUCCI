#!/bin/bash

# runSTAR_batch.sh

if [ $# -lt 10 ]; then # '-lt' stands for less than

	echo -e "\tThis script runs runSTAR.sh given a lookup table in CCHMC bsub.\n"

	echo -e "\tUsage: `basename $0` [lookup] [genome] [STAR index dir] [gtf] [exon|gene] [cpu] [mem] [time] [parent dir] [subdir name]\n"

	echo -e "\trunSTAR.sh runs STAR on single-end stranded RNA-seq data, with following key features:"
	echo -e "\t1) Softclips Truseq Read1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA."
	echo -e "\t2) Compute Per gene read count based on overlapping with exons or genes.\n"
	
	echo -e "\t[lookup]          ID[tab]desc[tab]fastq(fullpath)"
	echo -e "\t[genome]          Genome name to be used in the output filename (such as hg38)"
	echo -e "\t[STAR index dir]  STAR index directory. Use index without annotation."
	echo -e "\t[gtf]             GTF file"		
	echo -e "\t[exon|gene]       Reads overlapping exons [exon] or the entire gene [gene] should be computed for gene count file."
	echo -e "\t[cpu]             Number of CPU to be requested"
	echo -e "\t[mem]             Total memory (in GB) to request."
	echo -e "\t[time]            Walltime in hour."
	echo -e "\t[parent dir]      Results will be saved in [parent dir]/[ID]/[subdir name]"
	echo -e "\t[subdir name]     See above\n"
	
	exit 
fi


# 1) Global variables 

	lookup=${1};
	genome=${2};
	index=${3};
	gtf=${4};
	gtf_feature=${5};
	cpu=${6};
	mem=${7};
	runtime=${8};
	patdir=${9};
	subdir=${10};	


# 2) Process lookup file
	
	while IFS=$' \t\n' read -r idname desc fastq
	do

		# 1) Local variables	
		id_dir=$patdir"/"$idname;
		outdir=$patdir"/"$idname"/"$subdir;
			
		bsub_file=$outdir"/"$idname"_STAR.bat"
		prefix=$idname"_"$desc"_"$genome

		# 2) Make directories
		if [ ! -d $id_dir ]; 
		then mkdir $id_dir; 
		else echo -e "Directory $id_dir exits. Skip making it."
		fi
		
		if [ ! -d $outdir ]; 
		then mkdir $outdir;
#		else echo -e "Directory $outdir exits. Skip making it."		
		else rm -r $outdir; mkdir $outdir;			
		fi

		# 2) Make and run bsub

		if [ -f $bsub_file ]; then 
			echo "$bsub_file already exists. Skip to next data."
		else
			echo "Will run $bsub_file"

			# 2.1) Create pbs file 
			command="runSTAR.sh $fastq $outdir $cpu $index $gtf $gtf_feature $prefix"

			makeBat_nosub.sh $idname"_STAR" $cpu $mem $runtime $outdir "$command"

			# 2.2) Qsub .pbs 
			bsub < $bsub_file

		fi
		
	done < $lookup

exit;

