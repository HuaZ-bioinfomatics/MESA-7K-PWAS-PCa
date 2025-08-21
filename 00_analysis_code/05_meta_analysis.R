library(data.table)
library(dplyr)

if (!dir.exists('05_meta_analysis')){
    dir.create('05_meta_analysis')
}


# ------ Prostate Cancer
races <- c('AFA', 'CHN', 'CAU', 'HIS')

sample_size <- data.frame(
    AFA = 80999,
    CHN = 106599,
    CAU = 726828,
    HIS = 30336
)
for (race in races){
    association_df <- fread(paste0('04_association/association_prostate_cancer_2023_NG_', race, '/all_association.out'), data.table = F)
    association_df <- select(association_df, ID, PWAS.Z, PWAS.P)
    association_df$WEIGHT <- sample_size[race][[1]]
    write.table(association_df, paste0('05_meta_analysis/prostate_cancer_', race, '.txt'), sep = '\t', quote = F, row.names = F)
}
# generate metal config files manually
# meta for prostate cancer 
system('/data/sliu/software/generic-metal/metal 05_meta_analysis/metal_prostate_cancer.txt')


# prostate cancer combine meta and ethinic specific results
meta <- fread('05_meta_analysis/prostate_cancer_meta1.tbl', data.table = F)
meta <- select(meta, MarkerName, Zscore, 'P-value', 'Direction', 'HetISq', 'HetPVal')
colnames(meta) <- c('ID', 'Z_meta', 'P_meta', 'Direction', 'HetISq', 'HetPVal')
meta <- meta[order(meta$P_meta),]
meta$FDR_meta <- p.adjust(meta$P_meta, method = 'fdr')

annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, AptName, TargetFullName, Target, EntrezGeneSymbol)
annotation <- annotation[!duplicated(annotation$AptName),]
meta <- left_join(meta, annotation, by = c('ID' = 'AptName'))
write.table(meta, '05_meta_analysis/meta_PCa_association.txt', row.names = F, sep = '\t', quote = F)

AFA <- fread('04_association/association_prostate_cancer_2023_NG_AFA/all_association.out', data.table = F)
AFA <- select(AFA, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, FDR)
CHN <- fread('04_association/association_prostate_cancer_2023_NG_CHN/all_association.out', data.table = F)
CHN <- select(CHN, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, FDR)
CAU <- fread('04_association/association_prostate_cancer_2023_NG_CAU/all_association.out', data.table = F)
CAU <- select(CAU, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, FDR)
HIS <- fread('04_association/association_prostate_cancer_2023_NG_HIS/all_association.out', data.table = F)
HIS <- select(HIS, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, FDR)

res <- full_join(meta, AFA, by = 'ID')
res <- full_join(res, CHN, by = 'ID')
res <- full_join(res, CAU, by = 'ID')
res <- full_join(res, HIS, by = 'ID')
colnames(res) <- c('ID', 'Z_meta', 'P_meta', 'Direction', 'HetISq', 'HetPVal', 'FDR_meta','TargetFullName', 'Target', 'EntrezGeneSymbol', 'MODEL_AFA', 'Rsq_AFA', 'Z_AFA', 'P_AFA', 'FDR_AFA', 'MODEL_CHN', 'Rsq_CHN', 'Z_CHN', 'P_CHN', 'FDR_CHN', 'MODEL_CAU', 'Rsq_CAU', 'Z_CAU', 'P_CAU', 'FDR_CAU', 'MODEL_HIS', 'Rsq_HIS', 'Z_HIS', 'P_HIS', 'FDR_HIS')
write.table (res, '05_meta_analysis/combine_meta_ethinc_specific.txt', row.names = F, sep = '\t', quote = F)
res <- filter(res, FDR_meta < 0.05 | FDR_AFA < 0.05 | FDR_CHN < 0.05 | FDR_CAU < 0.05 | FDR_HIS < 0.05 )

write.table (res, '05_meta_analysis/combine_meta_ethinc_specific_Sig.txt', row.names = F, sep = '\t', quote = F)






# ------ T2D
races <- c('AFA', 'CHN', 'CAU', 'HIS')

sample_size <- data.frame(
    AFA = 156738,
    CHN = 427504,
    CAU = 1812017,
    HIS = 88743
)
for (race in races){
    association_df <- fread(paste0('04_association/association_T2D_2023_', race, '_1e-12/all_association.out'), data.table = F)
    association_df <- select(association_df, ID, PWAS.Z, PWAS.P)
    association_df$WEIGHT <- sample_size[race][[1]]
    write.table(association_df, paste0('05_meta_analysis/T2D_', race, '.txt'), sep = '\t', quote = F, row.names = F)
}
# generate metal config files manually
# meta for T2D 
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 05_meta_analysis/metal_T2D.txt')


# T2D combine meta and ethinic specific results
meta <- fread('05_meta_analysis/T2D_meta1.tbl', data.table = F)
meta <- select(meta, MarkerName, Zscore, 'P-value', 'Direction', 'HetISq', 'HetPVal')
colnames(meta) <- c('ID', 'Z_meta', 'P_meta', 'Direction', 'HetISq', 'HetPVal')
meta <- meta[order(meta$P_meta),]
meta$Bonf_meta <- p.adjust(meta$P_meta, method = 'bonferroni')

annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, AptName, TargetFullName, Target, EntrezGeneSymbol)
annotation <- annotation[!duplicated(annotation$AptName),]
meta <- left_join(meta, annotation, by = c('ID' = 'AptName'))
write.table(meta, '05_meta_analysis/meta_T2D_association.txt', row.names = F, sep = '\t', quote = F)

AFA <- fread('04_association/association_T2D_2023_AFA_1e-12/all_association.out', data.table = F)
AFA <- select(AFA, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, Bonferroni_p)
CHN <- fread('04_association/association_T2D_2023_CHN_1e-12/all_association.out', data.table = F)
CHN <- select(CHN, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, Bonferroni_p)
CAU <- fread('04_association/association_T2D_2023_CAU_1e-12/all_association.out', data.table = F)
CAU <- select(CAU, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, Bonferroni_p)
HIS <- fread('04_association/association_T2D_2023_HIS_1e-12/all_association.out', data.table = F)
HIS <- select(HIS, ID, MODEL, MODELCV.R2, PWAS.Z, PWAS.Z, PWAS.P, Bonferroni_p)

res <- full_join(meta, AFA, by = 'ID')
res <- full_join(res, CHN, by = 'ID')
res <- full_join(res, CAU, by = 'ID')
res <- full_join(res, HIS, by = 'ID')
colnames(res) <- c('ID', 'Z_meta', 'P_meta', 'Direction', 'HetISq', 'HetPVal', 'Bonf_meta','TargetFullName', 'Target', 'EntrezGeneSymbol', 'MODEL_AFA', 'Rsq_AFA', 'Z_AFA', 'P_AFA', 'Bonf_AFA', 'MODEL_CHN', 'Rsq_CHN', 'Z_CHN', 'P_CHN', 'Bonf_CHN', 'MODEL_CAU', 'Rsq_CAU', 'Z_CAU', 'P_CAU', 'Bonf_CAU', 'MODEL_HIS', 'Rsq_HIS', 'Z_HIS', 'P_HIS', 'Bonf_HIS')
write.table (res, '05_meta_analysis/T2D_combine_meta_ethinc_specific.txt', row.names = F, sep = '\t', quote = F)
res <- filter(res, Bonf_meta < 0.05 | Bonf_AFA < 0.05 | Bonf_CHN < 0.05 | Bonf_CAU < 0.05 | Bonf_HIS < 0.05 )

write.table (res, '05_meta_analysis/T2D_combine_meta_ethinc_specific_Sig.txt', row.names = F, sep = '\t', quote = F)


