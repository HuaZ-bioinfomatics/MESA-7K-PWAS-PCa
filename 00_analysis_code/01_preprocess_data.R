###################################
### 1. generate annotation file ###
###################################

library(biomaRt)
library(SomaDataIO)
library(data.table)
library(dplyr)
library(stringr)
my_adat <- read_adat('00_rawdata/UCSF-Ganz-MESA-Clinical_and_SomaScan_Data/Final_combined_ADAT/SS-2212385_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat')
Feature_Data <- data.frame(getAnalyteInfo(my_adat))
Feature_Data <- Feature_Data[1:12]
nrow(Feature_Data) #7596
write.table(Feature_Data, 'Feature_Data.txt', row.names = F, sep = '\t', quote = F)
Feature_Data <- filter(Feature_Data, Type == 'Protein')
nrow(Feature_Data) #7524
Feature_Data <- filter(Feature_Data,Organism == 'Human')
nrow(Feature_Data) #7289


sum(str_detect(Feature_Data$UniProt, "[|]")) #80 -- mapped to multiple targets
# mapped to single targets
Feature_Data_single <- Feature_Data[!str_detect(Feature_Data$UniProt, "[|]"),]
nrow(Feature_Data_single) # 7209 -- mapped to single targets

annotation_single <- dplyr::select(Feature_Data_single, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol)
nrow(annotation_single) # 7209
#######################################################################################
### first search by uniprot_gn_id in the ensembl database and matched entrezgene_id ###
#######################################################################################
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110)
allgene <- unique(annotation_single$UniProt)
length(allgene) # 6334 -- unique UniProt

# extract from ensembl database based on uniprot_gn_id
annotation_try1 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "uniprot_gn_id",
    values = allgene,
    mart = ensembl)
nrow(annotation_try1) # 7349
annotation_try1 <- unique(annotation_try1)
nrow(annotation_try1) # 7349
# Update TSS column based on the condition
annotation_try1$TSS <- ifelse(annotation_try1$strand == 1, annotation_try1$start_position, annotation_try1$end_position)

# ann_Feature_Data_single <- dplyr::select(ann_Feature_Data_single, uniprot_gn_id, hgnc_symbol, chromosome_name, TSS)
annotation_try1 <- dplyr::inner_join(annotation_single, annotation_try1, by = c('UniProt' = 'uniprot_gn_id'))
annotation_try1 <- unique(annotation_try1)
annotation_try1 <- annotation_try1[annotation_try1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
nrow(annotation_try1) # 7175
length(unique(annotation_try1$AptName)) # found 7045 AptName in autosome and sex chromosomes

correct0 <- annotation_try1[!annotation_try1$AptName %in% annotation_try1[duplicated(annotation_try1$AptName),]$AptName,]
nrow(correct0) # 6969
duplicated_rows_df <- annotation_try1[annotation_try1$AptName %in% annotation_try1[duplicated(annotation_try1$AptName),]$AptName,]


correct1 <- duplicated_rows_df[duplicated_rows_df$EntrezGeneID == duplicated_rows_df$entrezgene_id & duplicated_rows_df$EntrezGeneSymbol == duplicated_rows_df$hgnc_symbol,]
correct1 <- correct1[!duplicated(correct1$AptName),] # X, Y duplicate, force remain first row
nrow(correct1) # 67

unique(duplicated_rows_df[!duplicated_rows_df$AptName %in% correct1$AptName,]$AptName)
# [1] "seq.13944.3"  "seq.15542.19" "seq.21535.5"  "seq.22468.54" "seq.23586.32"
# [6] "seq.24465.28" "seq.2781.63"  "seq.7097.8"   "seq.7842.52" 
# cannot matched by both EntrezGeneID and EntrezGeneSymbol
# manually select
correct2 <- duplicated_rows_df[!duplicated_rows_df$AptName %in% correct1$AptName,][c(1, 4, 5, 11, 12, 14, 17, 19, 20),]
nrow(correct2) # 9

annotation_res1 <- rbind(correct0, correct1, correct2)

##############################################################
### second search by entrezgene_id in the ensembl database ###
##############################################################

not_found_in_annotation_try1 <- annotation_single[!annotation_single$AptName %in% annotation_res1$AptName,]
nrow(not_found_in_annotation_try1) # 164 AptName not found

not_found_in_annotation_try1 <- dplyr::filter(not_found_in_annotation_try1, EntrezGeneID != '')
nrow(not_found_in_annotation_try1) # 162 AptName not found

allgene <- unique(not_found_in_annotation_try1$EntrezGeneID)
annotation_try2 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "entrezgene_id",
    values = allgene,
    mart = ensembl)

nrow(annotation_try2) #2326
annotation_try2 <- unique(annotation_try2)
nrow(annotation_try2) #2326
annotation_try2 <- annotation_try2[annotation_try2$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
nrow(annotation_try2) #547
# Update TSS column based on the condition
annotation_try2$TSS <- ifelse(annotation_try2$strand == 1, annotation_try2$start_position, annotation_try2$end_position)
annotation_try2$entrezgene_id <- as.character(annotation_try2$entrezgene_id)
annotation_try2 <- inner_join(not_found_in_annotation_try1, annotation_try2, by = c('EntrezGeneID' = 'entrezgene_id'))
nrow(annotation_try2) # 644

correct3 <- annotation_try2[!annotation_try2$AptName %in% annotation_try2[duplicated(annotation_try2$AptName),]$AptName,]
nrow(correct3) # 41
duplicated_rows_df <- annotation_try2[annotation_try2$AptName %in% annotation_try2[duplicated(annotation_try2$AptName),]$AptName,]
duplicated_rows_df$uniprot_gn_id <- NA
duplicated_rows_df <- unique(duplicated_rows_df)
correct4 <- duplicated_rows_df
length(unique(correct4$AptName)) # 102
annotation_res2 <- rbind(correct3, correct4)
nrow(annotation_res2) # 143
###########################################################
### third search by hgnc_symbol in the ensembl database ###
###########################################################
not_found_in_annotation_try2 <- annotation_single[!annotation_single$AptName %in% c(annotation_res1$AptName, annotation_res2$AptName),]
nrow(not_found_in_annotation_try2) # 21

annotation_try3 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "hgnc_symbol",
    values = unique(not_found_in_annotation_try2$EntrezGeneSymbol),
    mart = ensembl)

nrow(annotation_try3) #138
annotation_try3 <- unique(annotation_try3)
nrow(annotation_try3) #138
annotation_try3 <- annotation_try3[annotation_try3$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
nrow(annotation_try3) #13
# Update TSS column based on the condition
annotation_try3$TSS <- ifelse(annotation_try3$strand == 1, annotation_try3$start_position, annotation_try3$end_position)
annotation_try3$entrezgene_id <- as.character(annotation_try3$entrezgene_id)
annotation_try3 <- inner_join(not_found_in_annotation_try2, annotation_try3, by = c('EntrezGeneSymbol' = 'hgnc_symbol'))
annotation_try3$uniprot_gn_id <- NA
annotation_res3 <- unique(annotation_try3)


#############################
### forth search manually ###
#############################
not_found_in_annotation_try3 <- annotation_single[!annotation_single$AptName %in% c(annotation_res1$AptName, annotation_res2$AptName, annotation_res3$AptName),]
nrow(not_found_in_annotation_try3) # 14
manually_df <- not_found_in_annotation_try3[not_found_in_annotation_try3$EntrezGeneSymbol %in% c('CARD17', 'TXNRD3NB', 'GGT2', 'C5orf38'),]
manually_df$chromosome_name <- c(11, 3, 22, 5)
manually_df$start_position <- c(105092486, 126571779, 21207973, 2752131)
manually_df$end_position <- c(105101414, 126655124, 21225554, 2755397)
manually_df$strand <- c(-1, -1, -1, 1)
manually_df$TSS <- ifelse(manually_df$strand == 1, manually_df$start_position, manually_df$end_position)
# not_found_in_annotation_try3$EntrezGeneSymbol
#  [1] "KIR2DS2"  "CARD17"   "GSTT1"    "TXNRD3NB" "FLJ44635" "KIR3DS1" 
#  [7] "GGT2"     "C5orf38"  "LILRA3"   "BAGE3"    "HLA-DRB3" "KIR2DL2" 
# [13] "KIR2DL5A" "KIR2DL5A"
# KIR2DS2 Scaffold HSCHR19KIR_LUCE_BDEL_HAP_CTG3_1: 130,021-144,329 reverse strand.
# CARD17 Chromosome 11: 105,092,486-105,101,414 reverse strand.
# GSTT1 Scaffold HSCHR22_1_CTG7: 270,314-278,855 reverse strand.
# TXNRD3NB Chromosome 3: 126,571,779-126,655,124 reverse strand.
# FLJ44635 not found
# KIR3DS1 Scaffold HSCHR19KIR_T7526_BDEL_HAP_CTG3_1: 65,476-80,089 reverse strand.
# GGT2 Chromosome 22: 21,207,973-21,225,554 reverse strand.
# C5orf38 Chromosome 5: 2,752,131-2,755,397 forward strand.
# LILRA3 Scaffold HSCHR19_4_CTG3_1: 268,758-275,852 reverse strand.
# BAGE3 not found
# HLA-DRB3 Scaffold HSCHR6_MHC_COX_CTG1: 3,934,009-3,947,126 reverse strand.
# KIR2DL2 Scaffold HSCHR19KIR_G248_BA2_HAP_CTG3_1: 133,429-147,923 reverse strand.
# KIR2DL5A Scaffold HSCHR19KIR_FH15_B_HAP_CTG3_1: 60,105-69,601 reverse strand.
# KIR2DL5A Scaffold HSCHR19KIR_FH15_B_HAP_CTG3_1: 60,105-69,601 reverse strand.

annotation_res1 <- select(annotation_res1, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS)
annotation_res2 <- select(annotation_res2, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS)
annotation_res3 <- select(annotation_res3, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS)

annotation_single_all <- rbind(annotation_res1, annotation_res2, annotation_res3, manually_df )
nrow(annotation_single_all) # 7199
# 7209 - 7199 = 10. ten are in Scaffold are could not be found 

# mapped to multiple targets
library(tidyverse)
Feature_Data_multiple <- Feature_Data[str_detect(Feature_Data$UniProt, "[|]"),]
nrow(Feature_Data_multiple) # 80 -- mapped to multiple targets

annotation_multiple <- dplyr::select(Feature_Data_multiple, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol)
nrow(annotation_multiple) # 80

annotation_multiple <- data.frame(separate_rows(annotation_multiple, UniProt, sep = "[|]"))
allgene <- unique(annotation_multiple$UniProt)
length(allgene) # 117 -- unique UniProt

# extract from ensembl database based on uniprot_gn_id
annotation_multiple_try1 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "uniprot_gn_id",
    values = allgene,
    mart = ensembl)

annotation_multiple_try1 <- annotation_multiple_try1[annotation_multiple_try1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ] 
# manually remove duplicate
duplicated(annotation_multiple_try1$uniprot_gn_id)
annotation_multiple_try1 <- annotation_multiple_try1[c(-47, -52, -53, -92),]
nrow(annotation_multiple_try1) # 114
# write.table(annotation_multiple_try1, 'annotation_multiple_try1.txt', row.names = F, sep  = '\t', quote = F)


annotation_multiple[!annotation_multiple$UniProt %in% annotation_multiple_try1$uniprot_gn_id,]

annotation_multiple_try2 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "hgnc_symbol",
    values = c('ITGB2', 'C8B', 'C1QB'),
    mart = ensembl)
annotation_multiple_try2$uniprot_gn_id <- NA
annotation_multiple_try2 <- unique(annotation_multiple_try2)
annotation_multiple_try2$uniprot_gn_id <- c('P02746', 'P07358', 'P05107')

annotation_multiple_res <- rbind(annotation_multiple_try1 , annotation_multiple_try2)

annotation_multiple_all <- left_join(annotation_multiple, annotation_multiple_res, by = c('UniProt' = 'uniprot_gn_id'))
annotation_multiple_all$TSS <- ifelse(annotation_multiple_all$strand == 1, annotation_multiple_all$start_position, annotation_multiple_all$end_position)
annotation_multiple_all <- dplyr::select(annotation_multiple_all, AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS )

# rbind single and multiple annotation result
annotation_all <- rbind(annotation_single_all, annotation_multiple_all)
length(unique(annotation_all$AptName)) # 7279

if (!file.exists("01_preprocess_data/annotation")) {
  dir.create("01_preprocess_data/annotation", recursive = TRUE)
}

write.table(annotation_all, '01_preprocess_data/annotation/SomaScan_annotation.txt', sep = '\t', quote = F, row.names =F)

annotation_TSS <- select(annotation_all, AptName, chromosome_name, TSS, start_position, end_position, strand)
colnames(annotation_TSS) <- c('AptName', 'Chr', 'TSS', 'Start', 'End', 'Strand')
annotation_TSS <- annotation_TSS[annotation_TSS$Chr %in% c(as.character(1:22)), ]
write.table(annotation_TSS, '01_preprocess_data/annotation/SomaScan_TSS.txt', sep = '\t', quote = F, row.names =F)


#############################################################
### 2. get overlapped subjects of covariate and phenotype ###
#############################################################

library(SomaDataIO)
library(data.table)
library(dplyr)


covariate <- fread('/mnt/lvm_vol_2/sliu/project/MESA_Olink_3K_PWAS/00_rawdata/MESA_Olink3Kdata/MESA_PhenoLangWu_20230811/MESA_PhenoLangWu_20230811.txt', data.table = F)
subject_with_covariate <- covariate$sidno

my_adat <- read_adat('00_rawdata/UCSF-Ganz-MESA-Clinical_and_SomaScan_Data/Final_combined_ADAT/SS-2212385_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat')
phenotype <- data.frame(my_adat) #nrow(phenotype) 17895
phenotype <- filter(phenotype, RowCheck == 'PASS' & SampleType == 'Sample') #nrow(phenotype) 15728

phenotype <- phenotype[, c(12, 34:ncol(phenotype))]
mapping <- read.table('/mnt/lvm_vol_2/sliu/project/MESA_Olink_3K_PWAS/00_rawdata/MESA_Olink3Kdata/MESA_A370_BridgeID_SIDNO_20231025.txt')

# exam 1, total 6212, D means duplicate, only select Barcode_2D1_E1 for further analysis
length(intersect(phenotype$SampleId, mapping$Barcode_2D1_E1[!is.na(mapping$Barcode_2D1_E1)])) #5923
length(intersect(phenotype$SampleId, mapping$Barcode_2D2_E1[!is.na(mapping$Barcode_2D2_E1)])) #286
length(intersect(phenotype$SampleId, mapping$Barcode_2D3_E1[!is.na(mapping$Barcode_2D3_E1)])) #1
length(intersect(phenotype$SampleId, mapping$Barcode_2D4_E1[!is.na(mapping$Barcode_2D4_E1)])) #2

# exam 4, total 5256
length(intersect(phenotype$SampleId, mapping$Barcode_2D1_E4[!is.na(mapping$Barcode_2D1_E4)])) #4983
length(intersect(phenotype$SampleId, mapping$Barcode_2D2_E4[!is.na(mapping$Barcode_2D2_E4)])) #271
length(intersect(phenotype$SampleId, mapping$Barcode_2D3_E4[!is.na(mapping$Barcode_2D3_E4)])) #2

# exam 5, total 4362
length(intersect(phenotype$SampleId, mapping$Barcode_2D1_E5[!is.na(mapping$Barcode_2D1_E5)])) #4183
length(intersect(phenotype$SampleId, mapping$Barcode_2D2_E5[!is.na(mapping$Barcode_2D2_E5)])) #179

# select only exam 1 for further analysis
mapping1 <- dplyr::select(mapping, sidno, Barcode_2D1_E1)
colnames(mapping1)[2] <- 'Barcode_2D'
mapping2 <- dplyr::select(mapping, sidno, Barcode_2D2_E1)
colnames(mapping2)[2] <- 'Barcode_2D'
mapping3 <- dplyr::select(mapping, sidno, Barcode_2D3_E1)
colnames(mapping3)[2] <- 'Barcode_2D'
mapping4 <- dplyr::select(mapping, sidno, Barcode_2D4_E1)
colnames(mapping4)[2] <- 'Barcode_2D'

# mapping_exam1 <- rbind(mapping1, mapping2, mapping3, mapping4)
mapping_exam1 <- mapping1
mapping_exam1 <- filter(mapping_exam1, sidno != 'NA' & Barcode_2D != 'NA')
# nrow(mapping_exam1) 5774
subjects_with_SomaScan <- mapping_exam1$sidno

subjects_with_covariate_and_SomaScan <- intersect(subject_with_covariate, subjects_with_SomaScan)

write.table(subjects_with_covariate_and_SomaScan, '01_preprocess_data/subjects_with_covariate_and_SomaScan.txt', row.names = F, sep = '\t', quote = F, col.names =F)


#################################
### 3. preprare genotype data ###
#################################
library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)


races <- c('AFA', 'CHN', 'CAU', 'HIS')

for (race in races){
    ## Create the output folder if it doesn't exist
    if (!dir.exists(file.path('01_preprocess_data', race))) {
        dir.create(file.path('01_preprocess_data', race))
    }
    if (!dir.exists(file.path('01_preprocess_data', race, 'Genotyped'))) {
        dir.create(file.path('01_preprocess_data', race, 'Genotyped'))
    }
    if (!dir.exists(file.path('01_preprocess_data', race, 'Imputed'))) {
        dir.create(file.path('01_preprocess_data', race, 'Imputed'))
    }

    files <- Sys.glob(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/00_rawdata/TOPMed_Imputation_V2/', race, '/chr*.dose.vcf.gz'))
    # Initialize parallel processing
    num_cores <- 30
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # for (file in files){
    foreach(file = files, .packages = c("data.table", "dplyr", "readr", "stringr"), .combine = rbind) %dopar% {
        
        ## extract Imputed SNP
        ## all SNPs R2 >= 0.8, SNP only, rename SNP ID as Chr:pos, remove duplicated SNPs ID
        ## plink2 --vcf chr10.dose.vcf.gz --exclude-if-info "R2 < 0.8" --snps-only --set-all-var-ids "@:#" --rm-dup exclude-all --make-bed --out chr10_SNP_unique
        outfile_name <- paste0('01_preprocess_data/', race, '/Imputed/',gsub('.dose.vcf.gz', '', last(str_split(file, '/')[[1]])), '_SNP')
        system2(
            command = "plink2",
            args = c(
            "--vcf", file,
            "--exclude-if-info", "\"R2 < 0.8\"",
            "--snps-only",
            "--set-all-var-ids", "\"@:#\"",
            "--rm-dup", "exclude-all",
            "--make-bed",
            "--out", outfile_name
            )
        )
        outfile_name2 <- paste0('01_preprocess_data/', race, '/Imputed/',gsub('.dose.vcf.gz', '', last(str_split(file, '/')[[1]])), '_SNP_unique')
        system2(
            command = "plink2",
            args = c(
            "--bfile", outfile_name,
            "--rm-dup", "exclude-all",
            "--keep", '01_preprocess_data/subjects_with_covariate_and_SomaScan.txt',
            "--make-bed",
            "--out", outfile_name2
            )
        )
        
        ## extract Genotyped SNP
        ## extract all SNPs R2 >= 0.99, SNP only, rename SNP ID as Chr:pos, remove duplicated SNPs ID
        outfile_name <- paste0('01_preprocess_data/', race, '/Genotyped/',gsub('.dose.vcf.gz', '', last(str_split(file, '/')[[1]])), '_SNP')
        system2(
            command = "plink2",
            args = c(
            "--vcf", file,
            "--require-info", "\"TYPED\"",
            "--exclude-if-info", "\"R2 < 0.99\"",
            "--snps-only",
            "--set-all-var-ids", "\"@:#\"",
            "--make-bed",
            "--out", outfile_name
            )
        )
        outfile_name2 <- paste0('01_preprocess_data/', race, '/Genotyped/',gsub('.dose.vcf.gz', '', last(str_split(file, '/')[[1]])), '_SNP_unique')
        system2(
            command = "plink2",
            args = c(
            "--bfile", outfile_name,
            "--rm-dup", "exclude-all",
            "--keep", '01_preprocess_data/subjects_with_covariate_and_SomaScan.txt',
            "--make-bed",
            "--out", outfile_name2
            )
        )
    }
    stopCluster(cl)
}



#####################################################
### 4. prepare phenotype data and covariate files ###
#####################################################

if (!file.exists("01_preprocess_data/covariate")) {
  dir.create("01_preprocess_data/covariate")
}


library(SomaDataIO)
library(data.table)
my_adat <- read_adat('00_rawdata/UCSF-Ganz-MESA-Clinical_and_SomaScan_Data/Final_combined_ADAT/SS-2212385_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat')
phenotype <- data.frame(my_adat) 
#nrow(phenotype) 17895
phenotype <- filter(phenotype, RowCheck == 'PASS' & SampleType == 'Sample') 
#nrow(phenotype) 15728

phenotype <- phenotype[, c(12, 34:ncol(phenotype))]
mapping <- read.table('/mnt/lvm_vol_2/sliu/project/MESA_Olink_3K_PWAS/00_rawdata/MESA_Olink3Kdata/MESA_A370_BridgeID_SIDNO_20231025.txt')

# exam 1, total 6212, D means duplicate, only select Barcode_2D1_E1 for further analysis
length(intersect(phenotype$Barcode2d, mapping$Barcode_2D1_E1[!is.na(mapping$Barcode_2D1_E1)])) #5896
length(intersect(phenotype$Barcode2d, mapping$Barcode_2D2_E1[!is.na(mapping$Barcode_2D2_E1)])) #285
length(intersect(phenotype$Barcode2d, mapping$Barcode_2D3_E1[!is.na(mapping$Barcode_2D3_E1)])) #1
length(intersect(phenotype$Barcode2d, mapping$Barcode_2D4_E1[!is.na(mapping$Barcode_2D4_E1)])) #2

# select only exam 1 for further analysis
mapping1 <- dplyr::select(mapping, sidno, Barcode_2D1_E1)

mapping_exam1 <- mapping1
mapping_exam1 <- filter(mapping_exam1, sidno != 'NA' & Barcode_2D1_E1 != 'NA')
mapping_exam1$Barcode_2D1_E1 <- as.character(mapping_exam1$Barcode_2D1_E1)
phenotype_matrix <- inner_join(mapping_exam1, phenotype, by = c('Barcode_2D1_E1' = 'Barcode2d'))
phenotype_matrix <- phenotype_matrix[, -2]

# get SomaScan annotation
Feature_Data <- data.frame(getAnalyteInfo(my_adat))
Feature_Data <- filter(Feature_Data, Type == 'Protein')
Feature_Data <- dplyr::select(Feature_Data, AptName, SeqId, SeqIdVersion, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol)
write.table(Feature_Data, 'SomaScan_7K_annotation.txt', row.names = F, sep = '\t', quote = F)

# extract protein level matrix
phenotype_matrix <- phenotype_matrix[, colnames(phenotype_matrix) %in% c('sidno', Feature_Data$AptName)]

####################################
### prepare phenotype 
####################################
annotation_TSS <- fread('01_preprocess_data/annotation/SomaScan_TSS.txt', data.table = F)
if (!file.exists("01_preprocess_data/phenotype")) {
  dir.create("01_preprocess_data/phenotype")
}
phenotype_matrix_extract <- phenotype_matrix[, colnames(phenotype_matrix) %in% annotation_TSS$AptName]
phenotype_matrix_extract <- cbind(sidno = phenotype_matrix$sidno, phenotype_matrix_extract)
write.table(phenotype_matrix_extract, '01_preprocess_data/phenotype/Soma_SMP.txt', row.names = F, sep = '\t', quote = F)


# covariate
covariate <- fread('/mnt/lvm_vol_2/sliu/project/MESA_Olink_3K_PWAS/00_rawdata/MESA_Olink3Kdata/MESA_PhenoLangWu_20230811/MESA_PhenoLangWu_20230811.txt', data.table = F)
covariate$kg <- covariate$wtlb1 * 0.453592
covariate$BMI <- covariate$kg / (covariate$htcm1/100)^2

# covariate without diabetes
covariate_no_DM <- filter(covariate, dm031c == 0)
covariate_no_DM <- filter(covariate_no_DM, covariate_no_DM$sidno %in% phenotype_matrix$sidno )

races <- c('AFA', 'CAU', 'CHN', 'HIS')

# male
for (race in races){
    # read sample with genotype
    sample <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/01_preprocess_data/', race, '/Genotyped/chr1_SNP_unique.fam'), data.table = F)
    covariate_in_race <- filter(covariate_no_cancer, covariate_no_cancer$sidno %in% sample$V2 & covariate_no_cancer$gender1 == 1)
    covariate_in_race <- dplyr::select(covariate_in_race, sidno, age1c, site1c, BMI, cig1c, pkyrs1c) # discussed on Oct 31, 2023, remove curalc1 because only 80.07% data
    covariate_in_race <- covariate_in_race[complete.cases(covariate_in_race), ]
    message ('race information from covariate is ',unique(covariate_in_race$race1c))
    message ('number of male samples without any cancers in ',race, ' is ', nrow(covariate_in_race))
    write.table(covariate_in_race, paste0('01_preprocess_data/covariate/covariate_for_', race, '_male.txt'), quote = F, row.names = F, sep = '\t')
}



