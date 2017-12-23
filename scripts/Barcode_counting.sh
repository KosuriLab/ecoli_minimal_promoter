echo "Extracting barcodes..."

dir='../processed_data/intermediate'

for i in ../rawdata/rLP5*.fastq
do 
	awk 'NR%4==2' $i > ${dir}/$(basename ${i/.fastq/}).txt
done

echo "Counting Barcodes..."

for i in ../processed_data/intermediate/rLP5*.txt
do 
	cut -c 1-20 $i | rev | tr ACGT TGCA | sort -T ./ --parallel=20 | uniq -c > \
	../processed_data/counts_$(basename $i)
done
# remove leftover sort files
rm -f sort*

# remove intermediate text files with barcodes
rm -f ../processed_data/intermediate/rLP5*.txt

# remove counts from name
cd ../processed_data/
rename 's/^counts_//' counts*
