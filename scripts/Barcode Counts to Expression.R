#This script takes Barcode RNA and DNA counts, 
# merges them to their respective promoters, 
# and calculates the strength of each promoter

#Install required packages
# install.packages("dplyr")
# install.packages("tidyr")


library("dplyr")
library("tidyr")

options(stringsAsFactors = F)

#Read in all barcode counts and normalize by Reads per Million

filelist = list.files('../processed_data', full.names = T)
filelist = filelist[grepl('rLP5', filelist)]
for(i in filelist) {
    name <- basename(i)
    x <- read.table(i, col.names=c(name, 'barcode'), header = F)
    x[[name]] <- 1000000*x[[name]] / sum(x[[name]])  #Normalizes by RPM
    assign(name,x)  
}

#combine reads for all barcodes 
rLP5_min <- full_join(rLP5_min_DNA_rep1.1.txt, rLP5_min_DNA_rep1.2.txt, by='barcode') %>%
  full_join(., rLP5_min_RNA_rep1.1.txt, by='barcode') %>%
  full_join(., rLP5_min_RNA_rep1.2.txt, by='barcode') %>%
  full_join(., rLP5_min_DNA_rep2.txt, by='barcode') %>%
  full_join(., rLP5_min_RNA_rep2.txt, by='barcode') %>%
  full_join(., rLP5_min_DNA_rep3.txt, by='barcode') %>%
  full_join(., rLP5_min_RNA_rep3.txt, by='barcode')
  
  
#rename all colummns that were named after text file
names(rLP5_min) = sub(".txt","", names(rLP5_min)) 
rm(list = c(filelist))
rm(x)


#Combine barcode counts with their promoter identity

barcode_stats_min <- read.table("../processed_data/barcode_statistics.txt", header = T)
Compare_barcode_Reps <-  barcode_stats_min[!is.na(barcode_stats_min$name),] %>%
  left_join(., rLP5_min , by ='barcode') #Filter out unmapped barcodes and assign read counts to mapped barcodes
Compare_barcode_Reps[is.na(Compare_barcode_Reps)] <- 0
  

#Evaluate how many barcodes appear in integrated sequencing reads
nrow(filter(Compare_barcode_Reps, rLP5_min_DNA_rep1.1 > 0 |
              rLP5_min_DNA_rep1.2 > 0 |
              rLP5_min_RNA_rep1.1 > 0 |
              rLP5_min_RNA_rep1.2 > 0 |
              rLP5_min_DNA_rep2 > 0 |
              rLP5_min_RNA_rep2 > 0 |
              rLP5_min_DNA_rep3 > 0 |
              rLP5_min_RNA_rep3 > 0))  #318,825/351,275 (90.5%) Barcodes show up

# Determine Expression as a function of the sum of all RNA counts for all barcodes 
# of a promoter divided by the sum of all DNA counts for all barcodes of a promoter
# Furthermore, filter out all promoters with fewer than 3 barcodes integrated

min_MOPS <- Compare_barcode_Reps %>% group_by(name) %>% 
  mutate(num_barcodes = n()) %>%
  filter(rLP5_min_DNA_rep1.1 > 0 |
           rLP5_min_DNA_rep1.2 > 0 |
           rLP5_min_RNA_rep1.1 > 0 |
           rLP5_min_RNA_rep1.2 > 0 |
           rLP5_min_DNA_rep2 > 0 |
           rLP5_min_RNA_rep2 > 0 |
           rLP5_min_DNA_rep3 > 0 |
           rLP5_min_RNA_rep3 > 0) %>%
  mutate(num_barcodes_integrated = n()) %>%
  filter(num_barcodes_integrated >= 4) %>% #Filter out promoters with fewer than 3 barcodes integrated
  mutate(RNA_exp_rep1.1 = sum(rLP5_min_RNA_rep1.1)/sum(rLP5_min_DNA_rep1.1),
         RNA_exp_rep1.2 = sum(rLP5_min_RNA_rep1.2)/sum(rLP5_min_DNA_rep1.2),
         RNA_exp_rep2 = sum(rLP5_min_RNA_rep2)/sum(rLP5_min_DNA_rep2),
         RNA_exp_rep3 = sum(rLP5_min_RNA_rep3)/sum(rLP5_min_DNA_rep3),
         DNA_sum = (sum(rLP5_min_DNA_rep1.1)+sum(rLP5_min_DNA_rep1.2)+sum(rLP5_min_DNA_rep2)+sum(rLP5_min_DNA_rep3)),
         RNA_exp_average = ((((RNA_exp_rep1.1+RNA_exp_rep1.2)/2)+RNA_exp_rep2+RNA_exp_rep3)/3)) %>% #Calculate expression averaged across replicates
  ungroup() %>% 
  select(name, RNA_exp_rep1.1, RNA_exp_rep1.2, RNA_exp_rep2,RNA_exp_rep3, RNA_exp_average, DNA_sum, num_barcodes, num_barcodes_integrated) %>% 
  distinct() 

write.table(min_MOPS, "../processed_data/revLP5_Min_MOPS_glu_expression.txt", quote = F, row.names = F)










