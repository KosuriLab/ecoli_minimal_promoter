all: merge map rna_barcode_count

merge: ../processed_data/min_assembled_reads.txt

map: ../processed_data/min_mappedv2.txt

rna_barcode_count: $(patsubst %,../processed_data/rLP5_min_%.txt,$(BASE))

BASE := DNA_rep1.1 DNA_rep1.2 DNA_rep2 DNA_rep3 RNA_rep1.1 RNA_rep1.2 RNA_rep2 RNA_rep3


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
$(patsubst %,../processed_data/rLP5_min_%.txt,$(BASE)): ./Barcode_counting.sh





