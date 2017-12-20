echo "Extracting barcodes..."

for i in rLP5*.fastq; do awk 'NR%4==2' $i > ${i/.fastq/}.txt; done

echo "Counting Barcodes..."

for i in rLP5*.txt; do cut -c 1-20 $i | rev | tr ACGT TGCA | sort -T ./ --parallel=20 | uniq -c > counts_$i; done

echo "Finalizing..."

mkdir -p ./Final_BC; mv counts_* $_

cd Final_BC

#removes counts_ from the name

rename 's/^counts_//' counts*
