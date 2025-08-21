

# 2ScML analysis
library(data.table)
library(dplyr)
library(BEDMatrix)
library(stringr)
library(TScML)
library(stats)

if(!dir.exists('06_robustness_using_2ScML')){
    dir.create('06_robustness_using_2ScML')
}

convert_37 <- fread('/data/sliu/database/dbSNP/GRCh37_to_rsID.txt', data.table = F, header = F)

# #----------- PCa 2ScML analysis
# # generate input file for AFA
# PCa_AFA <- fread('/data/sliu/database/GWAS_summary/prostate_multiethnic_2023_NG_Anqi/GCST90274715.tsv', data.table = F)
# PCa_AFA <- select(PCa_AFA, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value)
# PCa_AFA <- filter(PCa_AFA, PCa_AFA$effect_allele %in% c('A', 'T', 'C', 'G') & PCa_AFA$other_allele %in% c('A', 'T', 'C', 'G'))
# PCa_AFA$Chr_Pos <- paste0(PCa_AFA$chromosome, ':', PCa_AFA$base_pair_location)
# PCa_AFA <- inner_join(PCa_AFA, convert_37, by = c('Chr_Pos' = 'V1'))
# PCa_AFA <- select(PCa_AFA, V2, effect_allele, other_allele, beta, standard_error, p_value)
# colnames(PCa_AFA) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
# PCa_AFA$N <- 19391 + 61608
# fwrite(PCa_AFA, '06_robustness_using_2ScML/PCa_AFA.txt' , sep = '\t')


# # generate input file for CHN
# PCa_CHN <- fread('/data/sliu/database/GWAS_summary/prostate_multiethnic_2023_NG_Anqi/GCST90274716.tsv', data.table = F)
# PCa_CHN <- select(PCa_CHN, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value)
# PCa_CHN <- filter(PCa_CHN, PCa_CHN$effect_allele %in% c('A', 'T', 'C', 'G') & PCa_CHN$other_allele %in% c('A', 'T', 'C', 'G'))
# PCa_CHN$Chr_Pos <- paste0(PCa_CHN$chromosome, ':', PCa_CHN$base_pair_location)
# PCa_CHN <- inner_join(PCa_CHN, convert_37, by = c('Chr_Pos' = 'V1'))
# PCa_CHN <- select(PCa_CHN, V2, effect_allele, other_allele, beta, standard_error, p_value)
# colnames(PCa_CHN) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
# PCa_CHN$N <- 10809 + 95790
# fwrite(PCa_CHN, '06_robustness_using_2ScML/PCa_CHN.txt' , sep = '\t')

# # generate input file for CAU
# PCa_CAU <- fread('/data/sliu/database/GWAS_summary/prostate_multiethnic_2023_NG_Anqi/GCST90274714.tsv', data.table = F)
# PCa_CAU <- select(PCa_CAU, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value)
# PCa_CAU <- filter(PCa_CAU, PCa_CAU$effect_allele %in% c('A', 'T', 'C', 'G') & PCa_CAU$other_allele %in% c('A', 'T', 'C', 'G'))
# PCa_CAU$Chr_Pos <- paste0(PCa_CAU$chromosome, ':', PCa_CAU$base_pair_location)
# PCa_CAU <- inner_join(PCa_CAU, convert_37, by = c('Chr_Pos' = 'V1'))
# PCa_CAU <- select(PCa_CAU, V2, effect_allele, other_allele, beta, standard_error, p_value)
# colnames(PCa_CAU) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
# PCa_CAU$N <- 122188 + 604640
# fwrite(PCa_CAU, '06_robustness_using_2ScML/PCa_CAU.txt' , sep = '\t')


# # generate input file for HIS
# PCa_HIS <- fread('/data/sliu/database/GWAS_summary/prostate_multiethnic_2023_NG_Anqi/GCST90274717.tsv', data.table = F)
# PCa_HIS <- select(PCa_HIS, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value)
# PCa_HIS <- filter(PCa_HIS, PCa_HIS$effect_allele %in% c('A', 'T', 'C', 'G') & PCa_HIS$other_allele %in% c('A', 'T', 'C', 'G'))
# PCa_HIS$Chr_Pos <- paste0(PCa_HIS$chromosome, ':', PCa_HIS$base_pair_location)
# PCa_HIS <- inner_join(PCa_HIS, convert_37, by = c('Chr_Pos' = 'V1'))
# PCa_HIS <- select(PCa_HIS, V2, effect_allele, other_allele, beta, standard_error, p_value)
# colnames(PCa_HIS) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
# PCa_HIS$N <- 3931 + 26405 
# fwrite(PCa_HIS, '06_robustness_using_2ScML/PCa_HIS.txt' , sep = '\t')


# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

#allele qc
allele.qc = function(a1,a2,ref1,ref2) {
    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"

	flip1 = flip
	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
	return(snp)
}



perform_2ScML <- function(sig_protein, n1, n2, n.ref, sumstat, race){
    out <- data.frame(
        Protein = character(),
        beta = numeric(),
        se = numeric(),
        Z = numeric(),
        pval = numeric(),
        stringsAsFactors = FALSE
    )

    for (protein in sig_protein){
        res <- rep(NA, 5)
        res[1] <- protein
        if (race == 'AFA'){
            race <- 'AFR'
        }
        weight <- load(paste0('03_establish_models/', race, '_male_5e-9/08_establish_prediction_models/models/', protein, '.wgt.RDat'))
        # #count blup number, if blup is the best model, while number of SNP used in model is bigger than 500, change the p value to 1, we will not use this model
        # count_blup_SNP <- sum(!is.na(wgt.matrix[,'blup']))
        # if (count_blup_SNP > 500){
        #     cv.performance['pval', 'blup'] <- 1
        # }
        #change the p-value of top1 to 1, we will not use this model
        cv.performance['pval', 'top1'] <- 1
        #change the p-value of blup to 1, we will not use this model
        cv.performance['pval', 'blup'] <- 1
        #change the p-value of lasso to 1, we will not use this model
        cv.performance['pval', 'lasso'] <- 1

        # which rows have rsq
        row.rsq = grep( "rsq" , rownames(cv.performance) )
        # which rows have p-values
        row.pval = grep( "pval" , rownames(cv.performance) )	
        # Identify the best model
        mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))

        best = which.min(cv.performance[2,])

        if ( names(best) == "lasso" || names(best) == "enet" ) {
        keep = wgt.matrix[,best] != 0
        } else if ( names(best) == "top1" ) {
        keep = which.max(wgt.matrix[,best]^2)
        } else { 
        keep = 1:nrow(wgt.matrix)
        }

        #get SNP rsID used in models
        SNPS_used_in_model <- (snps[,c(2,5,6)])[keep,]$V2

        #get SNP in reference panel
        Z.ref.bim <- fread(paste0('03_establish_models/', race, '_male_5e-9/08_establish_prediction_models/ref_ld/', protein, '_1000G.bim'))
        SNPS_in_ref <- Z.ref.bim$V2

        #get SNP GWAS summary
        SNPS_in_GWAS_summary <- sumstat$SNP

        #common SNPs in model, reference panel, and GWAS summary
        common_SNPs <- intersect(intersect(SNPS_used_in_model, SNPS_in_ref), SNPS_in_GWAS_summary)

        #filter 

        snps <- filter(snps, snps$V2 %in% common_SNPs)
        Z.ref.bim <- filter(Z.ref.bim, Z.ref.bim$V2 %in% common_SNPs)


        pheno <- fread(paste0('03_establish_models/', race, '_male_5e-9/05_pheno_adjustment/', protein, '.pheno'), data.table = F)


        #reference panel cor, non-scale same as scale, here use non-scale
        Z.ref <- BEDMatrix(paste0('03_establish_models/', race, '_male_5e-9/08_establish_prediction_models/ref_ld/', protein, '_1000G'))
        Z.ref  <- PatchUp(Z.ref)
        Z.ref<- as.matrix(Z.ref)
        Z.ref <- as.data.frame(Z.ref)
        ref_SNP_ID <- str_split(colnames(Z.ref),'_') %>% unlist() %>% matrix(., ncol = 2, byrow = TRUE)
        colnames(Z.ref) <- ref_SNP_ID[,1]
        Z.ref <- Z.ref [,common_SNPs]
        Z.ref.original <- Z.ref
        cor.Z.ref.original = cor(Z.ref.original) + diag(0.00001,ncol(Z.ref.original))



        #calculate cor.D1Z1.original
        
        Z.Stage1 = BEDMatrix(paste0('03_establish_models/', race, '_male_5e-9/07_extract_potential_SNP_predictors/', protein, '_filtered'))
        Z.Stage1   <- PatchUp(Z.Stage1)
        Z.Stage1 <- as.matrix(Z.Stage1)
        Z.Stage1 <- as.data.frame(Z.Stage1)
        tmp_SNP_ID <- str_split(colnames(Z.Stage1),'_') %>% unlist() %>% matrix(., ncol = 2, byrow = TRUE)
        colnames(Z.Stage1) <- tmp_SNP_ID[,1]
        Z.Stage1 <- Z.Stage1[,common_SNPs]
        Z <- Z.Stage1 + diag(0.00001,ncol(Z.Stage1))
        Z1 = scale(Z,scale = F)
        D <- pheno[,protein]
        D1 = scale(D,scale = F)
        cor.D1Z1.original = as.numeric(cor(D1,Z1))



        ### Stage1 with individual-level data
        Stage1FittedModel = 
        TScMLStage1(cor.D1Z1 = cor.D1Z1.original,
                    Cap.Sigma.stage1 = cor(Z1),
                    n1 = n1,
                    p = nrow(Z.Stage1),
                    ind.stage1 = 1:ncol(Z1))
        #flip 

        qc1 = allele.qc( snps$V5 , snps$V6 , Z.ref.bim$V5 , Z.ref.bim$V6 )
        Stage1FittedModel[ qc1$flip ] = -1 * Stage1FittedModel[qc1$flip]

        # Match summary data to input, drop NA where summary data is missing
        m = match( Z.ref.bim$V2 , sumstat$SNP )
        sumstat_select = sumstat[m,]

        # QC / allele-flip the input and output
        qc = allele.qc( sumstat_select$A1 , sumstat_select$A2 , Z.ref.bim$V5 , Z.ref.bim$V6 )

        # Flip BETA-scores for mismatching alleles
        sumstat_select$BETA[ qc$flip ] = -1 * sumstat_select$BETA[ qc$flip ]
        sumstat_select$A1[ qc$flip ] = Z.ref.bim[ qc$flip]$V5
        sumstat_select$A2[ qc$flip ] = Z.ref.bim[ qc$flip]$V6
        cor.Y2Z2.original <- sumstat_select$BETA/sqrt(sumstat_select$BETA^2 + ( sumstat_select$N-2)*sumstat_select$SE^2)

        ### TScML Stage2 with summary data and reference panel
        lm1 = summary(lm(D1~Z1))
        Est.Sigma1Square = (lm1$sigma)^2
        if (!is.na(Est.Sigma1Square)){
            Est.Sigma2Square = as.numeric(1 - cor.Y2Z2.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.Y2Z2.original)
            if(Est.Sigma1Square<=0)
            {
            Est.Sigma1Square = 1
            }
            if(Est.Sigma2Square<=0)
            {
            Est.Sigma2Square = 1
            }
            number <- floor(ncol(Z.Stage1)/2)
            TScMLStage2.Ref =
                TScMLStage2(gamma.hat.stage1 = Stage1FittedModel,
                            cor.Y2Z2 = cor.Y2Z2.original,
                            Estimated.Sigma = cor.Z.ref.original,
                            n1 = n1,
                            n2 = n2,
                            p = ncol(Z.Stage1),   # may be different
                            K.vec.stage2 = 0:number,
                            Est.Sigma1Square = Est.Sigma1Square,
                            Est.Sigma2Square = Est.Sigma2Square)

            TScML.Summary.Var = 
                TScMLVar(Z.ref.original = Z.ref.original,
                        Stage1FittedModel = Stage1FittedModel,
                        betaalpha.hat.stage2 = TScMLStage2.Ref$betaalpha.hat.stage2,
                        Est.Sigma1Square = Est.Sigma1Square,
                        Est.Sigma2Square = Est.Sigma2Square,
                        n1 = n1,
                        n2 = n2,
                        n.ref = n.ref)
            res.beta <- TScMLStage2.Ref$betaalpha.hat.stage2[1] 
            res.se <- sqrt(TScML.Summary.Var )
            res.Z <- res.beta/res.se
            res.pval <- pnorm(-abs(res.beta/res.se))*2
            res[2] <- res.beta
            res[3] <- res.se
            res[4] <- res.Z
            res[5] <- res.pval
            print (res)
        }
        out[nrow(out) + 1, ] <- res
	
    }
    return(out)
}


# AFA 
AFA_sig <- fread('04_association/association_prostate_cancer_2023_NG_AFA/all_association.out', data.table = F)
AFA_sig <- filter(AFA_sig, FDR < 0.05)
n1 <- 450   	#sample size used to build models
n2 <- 19391 + 61608		#sample size of GWAS summary
n.ref <- 319
sumstat <- fread('06_robustness_using_2ScML/PCa_AFA.txt', data.table = F)
race <- 'AFA'
res_AFA <- perform_2ScML(AFA_sig$ID, n1, n2, n.ref, sumstat, race)
res_AFA$race <- 'AFA'
res_AFA <- inner_join(AFA_sig, res_AFA, by = c('ID' = 'Protein'))
res_AFA <- select(res_AFA, ID, race, PWAS.Z, PWAS.P, FDR, beta, se, Z, pval)
colnames(res_AFA) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')

# CHN 
CHN_sig <- fread('04_association/association_prostate_cancer_2023_NG_CHN/all_association.out', data.table = F)
CHN_sig <- filter(CHN_sig, FDR < 0.05)
n1 <- 289   	#sample size used to build models
n2 <- 10809 + 95790	#sample size of GWAS summary
n.ref <- 244
sumstat <- fread('06_robustness_using_2ScML/PCa_CHN.txt', data.table = F)
race <- 'CHN'
res_CHN <- perform_2ScML(CHN_sig$ID, n1, n2, n.ref, sumstat, race)
res_CHN$race <- 'CHN'
res_CHN <- inner_join(CHN_sig, res_CHN, by = c('ID' = 'Protein'))
res_CHN <- select(res_CHN, ID, race, PWAS.Z, PWAS.P, FDR, beta, se, Z, pval)
colnames(res_CHN) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')


# CAU 
CAU_sig <- fread('04_association/association_prostate_cancer_2023_NG_CAU/all_association.out', data.table = F)
CAU_sig <- filter(CAU_sig, FDR < 0.05)
n1 <- 758   	#sample size used to build models
n2 <- 122188 + 604640	#sample size of GWAS summary
n.ref <- 239
sumstat <- fread('06_robustness_using_2ScML/PCa_CAU.txt', data.table = F)
race <- 'CAU'
res_CAU <- perform_2ScML(CAU_sig$ID, n1, n2, n.ref, sumstat, race)
res_CAU$race <- 'CAU'
res_CAU <- inner_join(CAU_sig, res_CAU, by = c('ID' = 'Protein'))
res_CAU <- select(res_CAU, ID, race, PWAS.Z, PWAS.P, FDR, beta, se, Z, pval)
colnames(res_CAU) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')


# HIS 
HIS_sig <- fread('04_association/association_prostate_cancer_2023_NG_HIS/all_association.out', data.table = F)
HIS_sig <- filter(HIS_sig, FDR < 0.05)
n1 <- 474   	#sample size used to build models
n2 <- 3931 + 26405	#sample size of GWAS summary
n.ref <- 170
sumstat <- fread('06_robustness_using_2ScML/PCa_HIS.txt', data.table = F)
race <- 'HIS'
res_HIS <- perform_2ScML(HIS_sig$ID, n1, n2, n.ref, sumstat, race)
res_HIS$race <- 'HIS'
res_HIS <- inner_join(HIS_sig, res_HIS, by = c('ID' = 'Protein'))
res_HIS <- select(res_HIS, ID, race, PWAS.Z, PWAS.P, FDR, beta, se, Z, pval)
colnames(res_HIS) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')




res_all <- rbind(res_AFA, res_CHN, res_CAU, res_HIS)
res_all$p <- as.numeric(res_all$p)

res_all <- res_all[order(res_all$p), ]
res_all$fdr <- p.adjust(res_all$p, method = 'fdr')
res_all$Z <- as.numeric(res_all$Z)
res_all$z <- as.numeric(res_all$z)
res_all$Consistent_direction <- ifelse(sign(res_all$Z) == sign(res_all$z), 'TRUE', 'FALSE')
write.table(res_all, '06_robustness_using_2ScML/robustness_analysis_using_2ScML_PCa.txt', row.names = F, sep = '\t', quote = F)


library(data.table)
library(dplyr)
df <- fread('06_robustness_using_2ScML/robustness_analysis_using_2ScML_PCa.txt')

ann <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt')
df <- left_join(df, ann[, c(1, 3)], by = c("ID" = "AptName"))
df <- select(df, ID, TargetFullName, race, Z, P, FDR, beta, se, z, p, fdr, Consistent_direction)
colnames(df) <- c('ID', 'Target Full Name', 'Race', 'PWAS.Z', 'PWAS.P', 'PWAS.FDR', 'Beta', 'Se', 'Z', 'P', 'FDR', 'Consistent direction')
library(writexl)
write_xlsx(df, "06_robustness_using_2ScML/TableS2.xlsx")