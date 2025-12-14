#!/bin/bash

## process_GROseq_v1.sh

if [ $# -lt 3 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2025"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tTaking paired-end GRO-seq bam file as input and return fragment count per transcript.\n"
	echo -e "Usage:  \t`basename $0` [bam] [seed] [transcript]\n"

	echo -e "        \t[bam]              \tSorted paired-end bam. MAPQ >20 will be used."
	echo -e "        \t[seed]             \tseed.rpkm and seed.count will be generated. Will also be used in bg/bw file names."	
	echo -e "        \t[transcript]       \tbed6+ for which read count should be computed. Col4 => ID. Col6 => strand. No header.\n"
	
	echo -e "        \tTemp files will be generated in current working directory.\n"
		
	exit 
fi

# History
# v1.0: 


# 1) Modules
	module unload samtools
	module unload ucsctools
	module load bedtools/2.30.0 samtools/1.18.0 ucsctools/v466 #v380 stopped working.

# 2) Variables

	bam=$1
	seed=$2
	transcript=$3
	mapq=20
	RAN=$RANDOM.$RANDOM
	
# 3) File names

	logfile=logfile.$RAN

	chrom=chrom.$RAN.temp

	plustxpt=plustxpt.$RAN.temp
	minustxpt=minustxpt.$RAN.temp

	plustxptcov=plustxptcov.$RAN.temp
	minustxptcov=minustxptcov.$RAN.temp
	
	plusbg="$seed"_5end_p.bg
	minusbg="$seed"_5end_m.bg
	plusbw="$seed"_5end_p.bw
	minusbw="$seed"_5end_m.bw
	norm_plusbg="$seed"_5end10ext_np.bg
	norm_minusbg="$seed"_5end10ext_nm.bg
	norm_plusbw="$seed"_5end10ext_np.bw
	norm_minusbw="$seed"_5end10ext_nm.bw
	
	transcript_count=per_transcript_count.txt
	
#########################################
# Make 5' end bedgraph
#########################################

# 4) Make genome file from bam
	bamToChromfile.sh $bam > ./$chrom
	echo -e "Made chromosome length file:" $chrom >> $logfile

# 5) Make single-base bedgraph for 5' position of READ2 in paired-end fragments
	# In NEBNext directional RNA-seq library prep:
	# READ2 reads from the 5'end of fragmented RNAs
	# READ1 reads from the 3' end of the reverse complement of the fragmented RNAs.
	# So, we will extract READ2 5' end because this is the position where Pol II was located before run-on.
	# SAM flag 131: proper pair, second in pair 

	samtools view -u -q $mapq -f 131 $bam | 
	bedtools genomecov -ibam stdin -bg -5 -strand + > ./$plusbg
	
	samtools view -u -q $mapq -f 131 $bam | 
	bedtools genomecov -ibam stdin -bg -5 -strand - > ./$minusbg
	
	echo -e "Made 5'-end single-base bedgraph for plus and minus strand:" ./$plusbg ./$minusbg >> $logfile

#########################################
# Compute 5' end coverage per transcript
#########################################

# 5) Split transcripts by strand
	cat $transcript |
	awk '{ if($6 == "+") {OFS="\t"; print $1, $2, $3, $4 > "./'$plustxpt'"} else {OFS="\t"; print $1, $2, $3, $4 > "./'$minustxpt'"} }'
	
	echo -e "Completed splitting transcripts by strand" >> $logfile
	
# 6) Count 5' position per transcript for each strand
	bedtools sort -faidx $chrom -i ./$plustxpt |
	bedtools map -a stdin -b ./$plusbg -c 4 -o sum -null 0 -g $chrom > ./$plustxptcov
	
	bedtools sort -faidx $chrom -i ./$minustxpt |
	bedtools map -a stdin -b ./$minusbg -c 4 -o sum -null 0 -g $chrom > ./$minustxptcov
	
	echo -e "Computed 5'end counts for plus and minus transcripts" >> $logfile

	cat ./$plustxptcov ./$minustxptcov | cut -f 4,5 | sort -k1,1 > ./$transcript_count 
	
	echo -e "Printed 5'end counts for all transcripts:" $transcript_count >> $logfile
		

#####################################################################################################
# Make plus/minus 10 bp-extended, depth-normalized 5' position bedgraph (5end10ext) for visualization
#####################################################################################################

# 8) Get total fragment count to compute scale factor
	totcount=$(samtools view -c -f 131 -q $mapq $bam)
	echo -e "Total number of mapped fragments in "$bam" at MAPQ "$mapq" with FLAG "$flag":" $totcount >> $logfile

# 9) Scale factor (normalized by total fragment counts in millions)
	scale=$(awk 'BEGIN{printf "%.3f", 1000000/'$totcount'}')
	echo -e "Here is the scale factor:" $scale >> $logfile

# 10) Single-base bedgraph => +/- 10-bp-extend bed => scaled bedgraph
	awk '{if($2-10 < 0){s=0} else {s=$2-10} } {OFS="\t"; k=$4; for(i=0;i<k;i++) print $1, s, $2+10}' ./$plusbg |
	bedtools genomecov -i stdin -bg -g $chrom -scale $scale > ./$norm_plusbg

	awk '{if($2-10 < 0){s=0} else {s=$2-10} } {OFS="\t"; k=$4; for(i=0;i<k;i++) print $1, s, $2+10}' ./$minusbg |
	bedtools genomecov -i stdin -bg -g $chrom  -scale $scale > ./$norm_minusbg

	echo -e "Made +/-10-bp extended, depth-normalized 5' position bedgraph" >> $logfile
	
# 11) Convert to bw
	bedGraphToBigWig ./$plusbg $chrom ./$plusbw	
	bedGraphToBigWig ./$minusbg $chrom ./$minusbw	

	bedGraphToBigWig ./$norm_plusbg $chrom ./$norm_plusbw	
	bedGraphToBigWig ./$norm_minusbg $chrom ./$norm_minusbw	

	echo -e "Finished writing bigwig file." >> $logfile
		
	
#############################################
# Gzip count bedgraph and remove temp files
###############################################

# 12) Gzip
	gzip ./$plusbg
	gzip ./$minusbg
	

# 9) Remove temp files
	rm ./$plustxpt ./$minustxpt ./$plustxptcov ./$minustxptcov ./$norm_plusbg ./$norm_minusbg


exit;
	

#########################################
# Make bedgraph for 5'-end coverage
#########################################

	awk '{OFS="\t"; k=$4; for(i=0;i<k;i++) print $1, $2-10, $2+10}' ./$plusbg |
	awk '{OFS="\t"; k=$4; for(i=0;i<k;i++) print $1, $2-10, $2+10}' ./$minusbg |


# 7) Get total fragment count

	totcount=$(samtools view -c -f 131 -q $mapq $bam)
	echo -e "Total number of mapped fragments in "$bam" at MAPQ "$mapq" with FLAG "$flag":" $totcount >> $logfile

# 8) Make normalized 5'-end count bigwig	
	norm_plusbg="$seed"_readcount_np.bg
	norm_minusbg="$seed"_readcount_nm.bg
	norm_plusbw="$seed"_readcount_np.bw
	norm_minusbw="$seed"_readcount_nm.bw
	
# 9) Scale factor (normalized by total fragment counts in millions)
	scale=$(awk 'BEGIN{printf "%.3f", 1000000/'$totcount'}')
	echo -e "Here is the scale factor:" $scale >> $logfile

# 10) Write scaled bg
	awk '{printf "%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4*'$scale'}' ./$plusbg > ./$norm_plusbg
	echo -e "Finished writing plus-strand bg file." >> $logfile
	
	awk '{printf "%s\t%s\t%s\t%.3f\n", $1, $2, $3, $4*'$scale'}' ./$minusbg > ./$norm_minusbg
	echo -e "Finished writing minus-strand bg file." >> $logfile
	
# 11) Convert to bw
	bedGraphToBigWig ./$norm_plusbg $chrom ./$norm_plusbw	
	echo -e "Finished writing plus-strand bw file." >> $logfile

	bedGraphToBigWig ./$norm_minusbg $chrom ./$norm_minusbw	
	echo -e "Finished writing minus-strand bw file." >> $logfile

	
###########################################
# Author: Kohta Ikegami
# Contact: kohta.ikegami@cchmc.org
# Copyright: CC-BY-NC-SA 
###########################################	

