# This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter

# Install required packages
# install.packages("dplyr")
# install.packages("tidyr")


library("dplyr")
library("tidyr")

options(stringsAsFactors = F)

# SET TO WORKING DIRECTORY CONTAINING BARCODE RNA AND DNA COUNTS 


# Read in all barcode counts and normalize by Reads per Million

filelist = list.files(pattern = 'rLP5_*')
for(i in filelist) {
    x <- read.table(i, col.names=c(i, 'barcode'), header = F)
    x[[i]] <- 1000000*x[[i]] / sum(x[[i]])  #Normalizes by RPM
    assign(i,x)  
}

# combine reads for all barcodes 
rLP5_min <- full_join(rLP5_min_DNA1.txt, rLP5_min_DNA2.txt, by='barcode') %>%
    full_join(., rLP5_min_RNA1.txt, by='barcode') %>%
    full_join(., rLP5_min_RNA2.txt, by='barcode')

# rename all colummns that were named after text file
names(rLP5_min) = sub(".txt","", names(rLP5_min)) 
rm(list = c(filelist))
rm(x)

# Combine barcode counts with their promoter identity
barcode_stats_min <- read.table("./barcode_statistics.txt", header = T)
Compare_barcode_Reps <-  barcode_stats_min[!is.na(barcode_stats_min$name),] %>%
    # Filter out unmapped barcodes and assign read counts to mapped barcodes
    left_join(., rLP5_min , by ='barcode') 
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0
  

# Evaluate how many barcodes appear in integrated sequencing reads
nrow(filter(Compare_barcode_Reps, 
            rLP5_min_DNA1 > 0 | 
            rLP5_min_DNA2 > 0 | 
            rLP5_min_RNA1 > 0 | 
            rLP5_min_RNA2 > 0)) # 304,073/342,746 (89%) Barcodes show up

# Determine Expression as a function of the sum of all RNA counts for all barcodes 
# of a promoter divided by the sum of all DNA counts for all barcodes of a promoter
# Furthermore, filter out all promoters with fewer than 3 barcodes integrated
min_MOPS <- Compare_barcode_Reps %>% group_by(name) %>% 
    mutate(num_barcodes = n()) %>%
    filter(., rLP5_min_RNA1 > 0 | 
             rLP5_min_RNA2 > 0 | 
             rLP5_min_DNA1 > 0 | 
             rLP5_min_DNA2 > 0) %>%
    mutate(num_barcodes_integrated = n()) %>%
    # Filter out promoters with fewer than 3 barcodes integrated
    filter(num_barcodes_integrated >= 4) %>% 
    mutate(RNA_exp_1 = sum(rLP5_min_RNA1) / (sum(rLP5_min_DNA1) + sum(rLP5_min_DNA2)),
           RNA_exp_2 = sum(rLP5_min_RNA2) / (sum(rLP5_min_DNA2) + sum(rLP5_min_DNA1)),
           DNA_sum = (sum(rLP5_min_DNA2) + sum(rLP5_min_DNA1)),
           # Calculate expression averaged across replicates
           RNA_exp_12 = ((sum(rLP5_min_RNA1) + sum(rLP5_min_RNA2))/2) / (sum(rLP5_min_DNA1) + sum(rLP5_min_DNA2))) %>% 
    ungroup() %>% 
  select(name, RNA_exp_1, RNA_exp_2, RNA_exp_12, DNA_sum, num_barcodes, num_barcodes_integrated) %>% 
  distinct() 

write.table(min_MOPS, "./revLP5_Min_MOPS_glu_expression.txt", quote = F, row.names = F)
