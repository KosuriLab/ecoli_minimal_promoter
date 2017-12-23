all: process_seq analysis

process_seq: merge map rna_barcode_count

analysis: expression fig2 fig3 fig4 fig5 fig6 supp_figs

merge: ../processed_data/min_assembled_reads.txt

map: ../processed_data/min_mappedv2.txt
	
rna_barcode_count: $(patsubst %,../processed_data/rLP5_min_%.txt,$(BASE))

expression: ../processed_data/revLP5_Min_MOPS_glu_expression.txt

fig2: ../figs/Figure2C.pdf ../figs/Figure2E.pdf

fig3: ../figs/Figure3A.pdf ../figs/Figure3B.pdf

fig4: ../figs/Figure4A.pdf ../figs/Figure4B.pdf ../figs/Figure4C.pdf ../figs/Figure4D.pdf

fig5: ../figs/Figure5A.pdf ../figs/Figure5B.pdf ../figs/Figure5C.pdf ../figs/Figure5D.pdf

fig6: ../figs/Figure6A.pdf ../figs/Figure6B.pdf

# list one supp fig for brevity
supp_figs: ../figs/SF10_ECK726.pdf

BASE := DNA_rep1.1 DNA_rep1.2 DNA_rep2 DNA_rep3 RNA_rep1.1 RNA_rep1.2 RNA_rep2 RNA_rep3


# merge paired end DNA reads with PEAR
../processed_data/min_assembled_reads.txt: ./pearallel.sh
	@ echo "Merging reads..."
	@ ./pearallel.sh
	@ mkdir ../processed_data/intermediate
	@ mv ../processed_data/*assembled*.fastq ../processed_data/intermediate
	@ mv ../processed_data/*discarded.fastq ../processed_data/intermediate


# run barcode mapping
# for parameter descriptions: python barcode_mapping_ecoli.py -h
../processed_data/min_mappedv2.txt: ../processed_data/min_assembled_reads.txt
	@ echo "Mapping barcodes..."
	@ python barcode_mapping_ecoli.py $< \
	../ref/min_promoter_lib_up_20161118_controls_clean.txt \
	tab 150 24 end 20 $@
	@ mv barcode_statistics.txt ../processed_data
	@ mv variant_statistics.txt ../processed_data


# count RNA barcode reads
$(patsubst %,../processed_data/rLP5_min_%.txt,$(BASE)): $(patsubst %,../rawdata/rLP5_min_%.fastq,$(BASE))
	@ echo "Counting RNA barcode reads..."
	@ ./Barcode_counting.sh


../processed_data/revLP5_Min_MOPS_glu_expression.txt: $(patsubst %,../processed_data/rLP5_min_%.txt,$(BASE))
	@ echo "Calculating expression..."
	@ Rscript ./Barcode\ Counts\ to\ Expression.R


../figs/Figure2C.pdf ../figs/Figure2E.pdf: ../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Figure 2..."
	@ Rscript ./Figure2.R


../figs/Figure3A.pdf ../figs/Figure3B.pdf: ../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Figure 3..."
	@ Rscript ./Figure3.R


../figs/Figure4A.pdf ../figs/Figure4B.pdf ../figs/Figure4C.pdf ../figs/Figure4D.pdf: \
../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Figure 4..."
	@ Rscript ./Figure4.R


../figs/Figure5A.pdf ../figs/Figure5B.pdf ../figs/Figure5C.pdf ../figs/Figure5D.pdf: \
../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Figure 5..."
	@ Rscript ./Figure5.R


../figs/Figure6A.pdf ../figs/Figure6B.pdf: ../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Figure 6..."
	@ Rscript ./Figure6.R


../figs/SF10_ECK726.pdf: ../processed_data/revLP5_Min_MOPS_glu_expression.txt
	@ echo "Making Supplementary Figures..."
	@ Rscript ./Supplementary\ figures.R


print-%  : ; @echo $* = $($*)


