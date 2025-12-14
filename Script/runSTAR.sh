#!/bin/bash

# runSTAR.sh

if [ $# -lt 7 ]; then # '-lt' stands for less than

	echo -e "\n\tUsage: `basename $0` [fastq] [outdir] [cpu] [STAR index dir] [gtf] [exon|gene] [prefix]\n"

	echo -e "\tThis script runs STAR on single-end stranded RNA-seq data, with following key features:"
	echo -e "\t1) Softclips Truseq Read1 adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA."
	echo -e "\t2) Compute Per gene read count based on overlapping with exons or genes."
	
	echo -e "\t[fastq]           If more than one, separate by commas."	
	echo -e "\t[outdir]          Output directory in full path."
	echo -e "\t[cpu]             Num of CPU cores requested."
	echo -e "\t[STAR index dir]  STAR index directory. Use index without annotation."
	echo -e "\t[gtf]             GTF file."		
	echo -e "\t[exon|gene]       Reads overlapping exons [exon] or the entire gene [gene] should be computed for gene count file."
	echo -e "\t[prefix]          Prefix for begraph and bigwig files. Recommended to use ID_desc_genome.\n"
	
	exit 
fi

# 0) Modules
	module load STAR/2.7.9 samtools/1.18.0 ucsctools/v466

# 1) Variables
	fastq=$1;
	outdir=$2"/";
	cpu=$3;
	index=$4;
	gtf=$5;
	gtf_feature=$6;
	seed=$7
	myrandom=$RANDOM.$RANDOM;	
	chrom=chrom.$myrandom.temp
	
	echo -e "\n1) Finished setting variables.\n";
	
# 2) Outdir
	if [ ! -d $outdir ]; then mkdir $outdir; fi

	echo -e "2) Finished checking outdir.\n";
			
# 3) Run STAR
	STAR --runMode alignReads --runThreadN $cpu --genomeDir $index --readFilesIn $fastq \
	--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $cpu \
	--quantMode GeneCounts --sjdbGTFfile $gtf --sjdbGTFfeatureExon $gtf_feature \
	--outWigType bedGraph --outWigNorm RPM \
	--clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	--outFileNamePrefix $outdir

	echo -e "3) Completed STAR alignment\n";

# 4) Make bam index
	samtools index Aligned.sortedByCoord.out.bam

	echo -e "4) Made bam index\n";

# 5) Make compressed bedGraphs and bigwig

	## 5.1) Make chrome size file from bam header
	bamToChromfile.sh Aligned.sortedByCoord.out.bam | sort -k1,1 -k2,2n > $chrom		# bedGraphToBigWig requires -k1,1 -k2,2n sorted bg

	## 5.2) Sort bedgraph files
	sort -k1,1 -k2,2n Signal.Unique.str1.out.bg > $seed"_k1sort_STAR_str1.bg" 
	sort -k1,1 -k2,2n Signal.Unique.str2.out.bg > $seed"_k1sort_STAR_str2.bg" 

	## 5.3) Make bigwig files
	bedGraphToBigWig $seed"_k1sort_STAR_str1.bg" $chrom $seed"_k1sort_STAR_str1.bw"
	bedGraphToBigWig $seed"_k1sort_STAR_str2.bg" $chrom $seed"_k1sort_STAR_str2.bw"
	
	## 5.4) Gzip bedGraphs
	gzip $seed"_k1sort_STAR_str1.bg"
	gzip $seed"_k1sort_STAR_str2.bg"

	## 5.5) Remove unnecessary bedgraphs
	rm Signal.Unique.str1.out.bg Signal.Unique.str2.out.bg Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str2.out.bg $chrom
	
	echo -e "5) Made bedGraphs and bigwigs\n";

exit;


