if [ $# -lt 1 ]; then # '-lt' stands for less than
  echo "Usage: `basename $0` <file.bam>"
  exit 
fi

samtools view -H $1 | awk '($1=="@SQ"){OFS="\t"; print $2, $3}' | sed 's/SN://' | sed 's/LN://'

exit

