for i in 82 83
do
	echo "Aligning sample $i"
	pear -j 30 -f ../rawdata/pMin_$i_R2.fastq* -r ../rawdata/pMin_$i_R1.fastq* -o ../processed_data/$i.aligned.txt
	echo "Combining Replicates"
	cat ../processed_data/*.aligned.txt.assembled.fastq* | awk 'NR%4==2' > ../processed_data/min_assembled_reads.txt
done
