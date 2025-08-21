library(GenABEL)
library(data.table)
library(bacon)

races <- c('AFA', 'CHN', 'CAU', 'HIS')
diseases <- c('prostate_cancer_2023_NG')

for (disease in diseases){
    for (race in races){

         df <- fread(paste0('04_association/association_', disease, '_', race, '/all_association.out'), data.table = F)


        
        bc <- bacon(df$PWAS.Z)

        # control inflation using bacon
        df$P.bacon <- pval(bc)
        df$FDR.bacon <- p.adjust(df$P.bacon, method = 'fdr')

        # calculate  lambdas for raw data and bacon corrected 

        df = df[complete.cases(df), ]
        df$CHISQ <- qchisq(df$PWAS.P,1,lower.tail=FALSE)
        lambda_raw <- median(df$CHISQ) / qchisq(0.5,1)

        df$CHISQ_bacon_controled <- qchisq(df$P.bacon,1,lower.tail=FALSE)
        lambda_bacon_controled <- median(df$CHISQ_bacon_controled) / qchisq(0.5,1)


        # lambda_raw <- round(estlambda(df$PWAS.P, method="median")$estimate, 2)
        # lambda_bacon_controled <- round (estlambda(df$P.bacon, method="median")$estimate, 2)
        message (paste(disease, race, lambda_raw, sum(df$FDR < 0.05), lambda_bacon_controled, sum(df$FDR.bacon < 0.05)))
    }
}


# transfer
races <- c('AFA', 'CHN', 'CAU', 'HIS')
diseases <- c('T2D')

for (disease in diseases){
    for (race in races){
        df <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS_Transfer/04_association/', disease, '_', race, '/all_T2D_PWAS.txt'), data.table = F)
        bc <- bacon(df$PWAS.Z)

        # control inflation using bacon
        df$P.bacon <- pval(bc)
        df$FDR.bacon <- p.adjust(df$P.bacon, method = 'fdr')

        # calculate lambdas for raw data and bacon corrected 
        lambda_raw <- round(estlambda(df$PWAS.P, method="median")$estimate, 2)
        lambda_bacon_controled <- round (estlambda(df$P.bacon, method="median")$estimate, 2)
        message (paste(disease, race, lambda_raw, sum(df$FDR < 0.05), lambda_bacon_controled, sum(df$FDR.bacon < 0.05)))
    }
}



# ------- check T2D GWAS summary lambdas
library(GenABEL)
library(data.table)
races <- c('AFR', 'EAS', 'EUR', 'AMR')
for (race in races){
    df <- fread(paste0('/mnt/lvm_vol_2/sliu/database/GWAS_summary/T2D_2023/', race, '_MetalFixed_LDSC-CORR_Results1TBL.gz'), data.table = F)
    lambda <- round(estlambda(df$'P-value', method="median")$estimate, 2)
    print(paste(race, lambda))
}

