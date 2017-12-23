all: merge map rna_barcode_count

merge: ../processed_data/min_assembled_reads.txt

map: ../processed_data/min_mappedv2.txt

# merge paired end DNA reads with PEAR
../processed_data/min_assembled_reads.txt: ./pearallel.sh
	@ ./pearallel.sh
	@ mkdir ../processed_data/intermediate
	@ mv ../processed_data/*assembled*.fastq ../processed_data/intermediate
	@ mv ../processed_data/*discarded.fastq ../processed_data/intermediate


# run barcode mapping
# for parameter descriptions: python barcode_mapping_ecoli.py -h
../processed_data/min_mappedv2.txt: ../processed_data/min_assembled_reads.txt
	@ python barcode_mapping_ecoli.py $< \
	../ref/min_promoter_lib_up_20161118_controls_clean.txt \
	tab 150 24 end 20 $@
	@ mv barcode_statistics.txt ../processed_data
	@ mv variant_statistics.txt ../processed_data


# count RNA barcode reads
# this produces several text files, so not assigning a single target
rna_barcode_count: ./Barcode_counting.sh





