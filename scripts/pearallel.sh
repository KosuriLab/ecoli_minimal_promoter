for i in 82 83
	do
		echo "Aligning sample $i"
pear -j 30 -f *$i*R2* -r *$i*R1* -o $i.aligned.txt

echo "Combining Replicates"

cat *.aligned.txt.assembled.fastq | awk 'NR%4==2' > min_assembled_reads.txt

done;
