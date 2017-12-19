This directory contains all data files to be used in conjuction with our R scripts available on Github.

Brief Descriptions:

barcode_statistics.txt: This file contains mapping information for all barcodes. Most notably, this includes the barcode sequence, 
what variant sequence it mapped to,
the name of that variant sequence,
and the number of reads mapping to that variant sequence

variant_statistics.txt: This file is essentially a summary of barcode_statistics.txt and describes mapping information for each variant.

SpacerSequences.txt:
This file contains each spacer sequence and the name. In the paper, spacers are generally referred to by their last 3 digits.
It is used in the R script to calculate GC content of each spacer.


revLP5_Min_MOPS_glu_expression.txt:

This file contains the DNA, RNA, and expression measurements for all promoters tested with more than 3 barcodes integrated.
Most importantly is RNA_exp_12, which is an average of our two expression measurements.

rlp5Min_SplitVariants.txt:

This file is the same as rLP5....txt except the names are split by minus10, minus35, e.t.c, to allow for facile analysis of each variant/element.


Thanks for your interest in this study! 

