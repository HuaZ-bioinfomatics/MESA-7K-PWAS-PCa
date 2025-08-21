# validate non-Hispanic White models using INTERVAL data

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(stringr)

if (!dir.exists('07_validation_using_INTERVAL')){
    dir.create('07_validation_using_INTERVAL')
}

# prepare INTERVAL genotype data
if (!dir.exists('07_validation_using_INTERVAL/01_INTERVAL_genotype')){
    dir.create('07_validation_using_INTERVAL/01_INTERVAL_genotype')
}

# extract only SNP
files <- Sys.glob('/data/sliu/database/INTERVAL/EGAD00010001544/INTERVAL_SOMALOGIC_POSTQC_chrom*final_v1.bgen.bim')
for (file in files){
    chr_number <- str_split(last(str_split(file, '/')[[1]]), '_')[[1]][[5]]
    infile <- gsub('.bgen.bim', '.bgen', file)
    system2(
      command = "00_analysis_code/bin/plink2",
      args = c(
        "--bfile", infile,
        "--snps-only",
        "--make-bed",
        "--out", paste0('07_validation_using_INTERVAL/01_INTERVAL_genotype/chr',chr_number)
      )
    )
}

# generate rsID convert file
rsID <- fread('/data/sliu/database/INTERVAL/EGAD00010001544/Variant_info.tsv', data.table = F)
rsID <- select(rsID, VARIANT_ID, rsID)
write.table(rsID, '07_validation_using_INTERVAL/01_INTERVAL_genotype/rsID_convert.txt', col.names = F, row.names = F, sep = '\t', quote = F)

# convert to rsID
files <- Sys.glob('07_validation_using_INTERVAL/01_INTERVAL_genotype/chr*.bim')
for (file in files){
    chr_number <- last(str_split(file, '/')[[1]]) %>% gsub('.bim', '', .)
    infile <- gsub('.bim', '', file)
    system2(
      command = "00_analysis_code/bin/plink2",
      args = c(
        "--bfile", infile,
        "--snps-only",
        "--update-name", "07_validation_using_INTERVAL/01_INTERVAL_genotype/rsID_convert.txt",
        "--rm-dup", "exclude-all",
        "--make-bed",
        "--out", paste0('07_validation_using_INTERVAL/01_INTERVAL_genotype/',chr_number, '_rsID')
      )
    )
}

# merge all chromosomes
# Combine chromosomes 1-22 into a single file
all_chromosome <- file.path("07_validation_using_INTERVAL/01_INTERVAL_genotype/all_chromosome.txt")
writeLines(paste0("07_validation_using_INTERVAL/01_INTERVAL_genotype/chr", 2:22, "_rsID"), con = all_chromosome)

system2(
    command = "00_analysis_code/bin/plink",
    args = c(
        "--bfile", '07_validation_using_INTERVAL/01_INTERVAL_genotype/chr1_rsID',
        "--merge-list", all_chromosome,
        "--make-bed",
        "--allow-no-sex",
        "--out", paste0("07_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID")
    )
)




annotation_MESA <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
# select only single target
annotation_MESA <- annotation_MESA[!duplicated(annotation_MESA$AptName),]
annotation_MESA <- select(annotation_MESA, UniProt, AptName)


annotation_INTERVAL <- fread('/data/sliu/database/INTERVAL/EGAD00010001544/001_SOMALOGIC_GWAS_protein_info.csv', data.table = F)
annotation_INTERVAL <- select(annotation_INTERVAL, UniProt, SOMAMER_ID)

intersect <- inner_join(annotation_MESA, annotation_INTERVAL, by = c('UniProt'))

measured_pheno_all <- fread('/data/sliu/database/INTERVAL/EGAD00010001544/INTERVAL_SOMALOGIC_POSTQC_GWASIN_PROTEINDATA_v1.tsv', data.table = F)
male_sample <- filter(measured_pheno_all, gender == 'M')$Sample_Name
measured_pheno_all <- select(measured_pheno_all, c('Sample_Name',intersect$SOMAMER_ID))


# prepare INTERVAL measured protein level
if (!dir.exists('07_validation_using_INTERVAL/02_INTERVAL_protein')){
    dir.create('07_validation_using_INTERVAL/02_INTERVAL_protein')
}


for(i in 2:ncol(measured_pheno_all)){
    tmp <- measured_pheno_all[,c(1,i)]
    id <- colnames(measured_pheno_all)[i]
    write.table(tmp, paste0('07_validation_using_INTERVAL/02_INTERVAL_protein/', id, '.pheno'), sep = '\t', quote = F, row.names = F)
}



#------------- models in CAU male with Rsq >= 0.01 
models <- fread('03_establish_models/CAU_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
models$ID <- gsub('.wgt.RDat', '', models$model)

intersect_in_MESA <- filter(intersect, AptName %in% models$ID)
# select only unique pairs
intersect_in_MESA <- filter(intersect_in_MESA, !AptName %in% unique(intersect_in_MESA[duplicated(intersect_in_MESA$AptName),]$AptName))


# Initialize parallel processing
num_cores <- 50
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# generate predicted protein levels and calcuate Rsq
if (!dir.exists('07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male')){
    dir.create('07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male')
}

Rsq_res <- foreach(i = 1:nrow(intersect_in_MESA),.packages = c("data.table", "dplyr"), .combine = rbind)%dopar%{
    ID_MESA <- intersect_in_MESA[i,]$AptName
    ID_INTERVAL <- intersect_in_MESA[i,]$SOMAMER_ID
    system(paste0('Rscript 00_analysis_code/bin/make_score.R 03_establish_models/CAU_male_5e-9/08_establish_prediction_models/models/', ID_MESA, '.wgt.RDat > 07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male/', ID_MESA ))
    system(paste0('00_analysis_code/bin/plink --bfile 07_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID --score 07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male/', ID_MESA, ' 1 2 4 --out 07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male/', ID_MESA, '_predicted_pheno'))
    if (file.exists(paste0('07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male/', ID_MESA, '_predicted_pheno.profile'))){
        predicted_pheno <- fread(paste0('07_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level_CAU_male/', ID_MESA, '_predicted_pheno.profile'), data.table = F)
        predicted_pheno <- select(predicted_pheno, IID, SCORE)
        measured_pheno <- fread(paste0('07_validation_using_INTERVAL/02_INTERVAL_protein/', ID_INTERVAL, '.pheno'), data.table = F)
        tmp <- left_join(predicted_pheno, measured_pheno, by = c('IID' = 'Sample_Name'))
        # only extract male
        tmp <- tmp[tmp$IID %in% male_sample,]
        colnames(tmp)[3] <- 'protein_abundance'
        fit <- lm(SCORE~protein_abundance, tmp)
        Rsq <- summary(fit)$adj.r.squared
        res_tmp <- c(ID_MESA, ID_INTERVAL, Rsq)
    }
}

Rsq_res <- data.frame(Rsq_res)
colnames(Rsq_res) <- c('AptName', 'INTERVAL_ID', 'External_Rsq')
Rsq_res$External_Rsq <- as.numeric(Rsq_res$External_Rsq)
paste0((round(nrow(filter(Rsq_res, External_Rsq >= 0.01))/nrow(Rsq_res), 4) * 100 ), '%') # "69.19%"
stopCluster(cl)
# get internal Rsq
Internal_R2 <- select(models, ID, Rsq)
Rsq_res <- left_join(Rsq_res, Internal_R2, by = c('AptName' = 'ID'))

colnames(Rsq_res)[4] <- 'Internal_Rsq'
write.table(Rsq_res, '07_validation_using_INTERVAL/validation_Rsq_External_Internal_CAU_male.txt', row.names = F, sep = '\t', quote = F)




# CAU male models classification
library(data.table)
library(dplyr)
models <- fread('03_establish_models/CAU_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
models$ID <- gsub('.wgt.RDat', '', models$model)


models_cis <- filter(models,  N_trans == 0)
models_trans <- filter(models, N_cis ==0)
models_cis_trans <- filter(models, N_cis != 0 & N_trans != 0)

df <- fread('07_validation_using_INTERVAL/validation_Rsq_External_Internal_CAU_male.txt', data.table = F)
df$External_Rsq <- as.numeric(df$External_Rsq)

df_all <- filter(df, AptName %in% models$ID)
paste0(nrow(filter(df_all, External_Rsq >= 0.01)), '/', nrow(df_all), '=',(round(nrow(filter(df_all, External_Rsq >= 0.01))/nrow(df_all), 4) * 100 ), '%')

df_cis <- filter(df, AptName %in% models_cis$ID)
paste0(nrow(filter(df_cis, External_Rsq >= 0.01)), '/', nrow(df_cis), '=',(round(nrow(filter(df_cis, External_Rsq >= 0.01))/nrow(df_cis), 4) * 100 ), '%')


df_cis_trans <- filter(df, AptName %in% models_cis_trans$ID)
paste0(nrow(filter(df_cis_trans, External_Rsq >= 0.01)), '/', nrow(df_cis_trans), '=',(round(nrow(filter(df_cis_trans, External_Rsq >= 0.01))/nrow(df_cis_trans), 4) * 100 ), '%')

df_trans <- filter(df, AptName %in% models_trans$ID)
paste0(nrow(filter(df_trans, External_Rsq >= 0.01)), '/', nrow(df_trans), '=',(round(nrow(filter(df_trans, External_Rsq >= 0.01))/nrow(df_trans), 4) * 100 ), '%')



# [1] "69.19%"
# [1] "584/844=69.19%"
# [1] "331/471=70.28%"
# [1] "95/115=82.61%"
# [1] "158/258=61.24%"
