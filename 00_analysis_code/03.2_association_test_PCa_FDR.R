library(data.table)
library(dplyr)
library(stringr)
res_summary <- data.frame(
    race_treshold = as.character(),
    N_models = as.numeric(),
    N_test = as.numeric(),
    N_sig_FDR = as.numeric(),
    N_sig_Bonf = as.numeric()
)


# Set the path to the directory you want to explore
directory_path <- "04_association"

# Use list.dirs to get a list of all directories in the specified path
all_directories <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)
matched_dirs <- all_directories[grepl("^association_prostate_cancer_2023_NG_", basename(all_directories))]



for (dir in matched_dirs){
    res <- data.frame()
    files <- Sys.glob(paste0(dir, '/*.txt'))
    for (file in files){
        df <- fread(file, data.table = F)
        res <- rbind(res, df)
    }
    if (grepl('prostate', dir)){
        if (grepl('AFA', dir)){
            posfile <- paste0('03_establish_models/AFR_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
        }else{
            posfile <- paste0('03_establish_models/', last(str_split(dir, '_')[[1]]), '_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
        }
    }   
    # get models with Rsq >= 0.01
    models <- fread(posfile, data.table = F)
    models$ID <- gsub('.wgt.RDat', '', models$model)

    # remove PWAS with NA
    res <- res[!is.na(res$PWAS.P),]

    # filter models with Rsq >= 0.01
    res <- filter(res, ID %in% models$ID)
    res <- res[order(res$PWAS.P),]
    res$FDR <- p.adjust(res$PWAS.P, method = 'fdr')
    res$Bonferroni_p <- p.adjust(res$PWAS.P, method = 'bonferroni')
    res_summary[nrow(res_summary) + 1,] <- c(last(str_split(dir, '/')[[1]]), nrow(models), nrow(res), sum(res$FDR < 0.05), sum(res$Bonferroni_p < 0.05))
    write.table(res, paste0(dir, '/all_association.out'), row.names = F, quote = F, sep = '\t')
   
}

print(res_summary)
