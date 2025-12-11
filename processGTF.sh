#!/bin/bash

if [ $# -lt 1 ]  
then
	echo -e "\n"
	echo -e "Tool:   \t`basename $0`"
	echo -e "Author: \tKohta Ikegami, 2021"
	echo -e "Version:\tv1.0"
	echo -e "Summary:\tFrom GTF, make a bed file (*_pergene.bed) with a single longest transcript by exon sum.\n"
	echo -e "Usage:  \t`basename $0` [GTF]"
	echo -e "\n"
	exit 
fi


# Modules
	module load ucsctools/v380 bedtools/2.30.0 

# Variables
	gtf=$1

# Filenames
	myrandom=$RANDOM.$RANDOM
	predtemp=pred.$myrandom.temp
	txtemp=tx.$myrandom.temp
	genetemp=gene.$myrandom.temp
	geneBed=${gtf/.gtf/_pergene.bed}


# 1) Make genePred from gtf 
	gtfToGenePred -genePredExt $gtf $predtemp

# 2) Get pergene transcript from genePred. A single transcript with a largest sum of exons will be retained. 
	awk '{OFS="\t"; print $2,$4,$5,$1,$8,$3,$12,$9,$10}' $predtemp |
	bedtools expand -c 8,9 |
	awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$9-$8}' |
	bedtools groupby -g 1,2,3,4,5,6,7 -c 8 -o sum |
	sort -k7,7 -k8,8nr | 
	bedtools groupby -g 7 -c 8 -o first -full |
	awk '{OFS="\t"; print $7,$1,$2,$3,$4,$8,$6}' |
	sort -k1,1 > $txtemp

# 3) Get gene list from gtf
	cat $gtf |
	sed 's/"//g' | sed 's/;//g' |
	awk -F "\t| " '($3 == "gene"){OFS="\t"; print $1,$4,$5,".",$7,$10,$14,$12}' |
	sort -k6,6 > $genetemp

# 4) Check gene number  
	txtempLen=$(wc -l $txtemp | awk '{print $1}')
	genetempLen=$(wc -l $genetemp | awk '{print $1}')

	if [ $txtempLen != $genetempLen ]; then
			echo -e "Pergene transcript number does not match gene number. Exit."
			exit				    
		else 
			echo -e "Pergene transcript number matches gene number. Proceed."
	fi 
	
# 5) Join gene file with pergene transcript file 
	join $genetemp $txtemp -1 6 -2 1 |
	awk 'BEGIN{print "#chr\tstart\tstop\tscore\tgene_id\tstrand\tgene_name\tgene_type\ttxchr\ttxstart\ttxstop\ttranscript_id\ttxlength\ttxstrand"} {OFS="\t"; print $2,$3,$4,$5,$1,$6,$7,$8,$9,$10,$11,$12,$13,$14}' |
	sort -k1,1 -k2,2n > $geneBed

# 6) Remove temp files
	rm $predtemp $txtemp $genetemp

exit