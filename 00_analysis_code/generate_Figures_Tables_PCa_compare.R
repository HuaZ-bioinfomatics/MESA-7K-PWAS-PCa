library(data.table)
library(dplyr)
library(writexl)
library(readxl)
ann <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt')
ann <- ann[!duplicated(ann$AptName),]
ann <- select(ann, AptName, SomaId, TargetFullName, Target, EntrezGeneSymbol, UniProt)
AFA <- fread('04_association/association_prostate_cancer_2023_NG_AFA/all_association.out')
AFA <- filter(AFA, AFA$FDR < 0.05 )
AFA <- select(AFA, ID, PWAS.Z, PWAS.P)
AFA$Race <- 'AFA'
CHN <- fread('04_association/association_prostate_cancer_2023_NG_CHN/all_association.out')
CHN <- filter(CHN, CHN$FDR < 0.05 )
CHN <- select(CHN, ID, PWAS.Z, PWAS.P)
CHN$Race <- 'CHN'
CAU <- fread('04_association/association_prostate_cancer_2023_NG_CAU/all_association.out')
CAU <- filter(CAU, CAU$FDR < 0.05 )
CAU <- select(CAU, ID, PWAS.Z, PWAS.P)
CAU$Race <- 'CAU'
HIS <- fread('04_association/association_prostate_cancer_2023_NG_HIS/all_association.out')
HIS <- filter(HIS, HIS$FDR < 0.05 )
HIS <- select(HIS, ID, PWAS.Z, PWAS.P)
HIS$Race <- 'HIS'


all <- rbind(AFA, CHN, CAU, HIS)
all <- left_join(all, ann, by = c('ID' = 'AptName'))


# Wu's Excel file
wu_df <- as.data.table(read_excel("Figures_Tables/reported_proteins/Wu.CancerRes_2019_protein.xlsx"))
wu_df <- wu_df[, .(`Protein full name`, `OR (95% CI)b`, P)]
setnames(wu_df, old = c("Protein full name", "OR (95% CI)b", "P"),
         new = c("TargetFullName", "OR", "P"))
wu_df[, (names(wu_df)) := lapply(.SD, function(col) {
  if (is.character(col)) gsub("^\\p{Z}+|\\p{Z}+$", "", col, perl = TRUE) else col
})]

all_joined_1 <- inner_join(all, wu_df, by = "TargetFullName")
all_joined_1$Study <- 'Wu et.al 2019'


# zhong
zhong_df <- as.data.table(read_excel("Figures_Tables/reported_proteins/zhong2023HMG.xlsx"))

zhong_df <- zhong_df[, .(`Protein full name`, `Z score`, `P-value`)]
setnames(zhong_df, old = c("Protein full name", "Z score", "P-value"),
         new = c("TargetFullName", "Z", "P"))
zhong_df[, (names(zhong_df)) := lapply(.SD, function(col) {
  if (is.character(col)) gsub("^\\p{Z}+|\\p{Z}+$", "", col, perl = TRUE) else col
})]
all_joined_2 <- inner_join(all, zhong_df, by = "TargetFullName")
all_joined_2$Study <- 'Zhong et.al 2023'

# Ren
Ren_df <- as.data.table(read_excel("Figures_Tables/reported_proteins/Ren2023JournalofTranslationalMedicine.xlsx"))

Ren_df <- Ren_df[, .(`UniProt`, `OR (95%)`, `P`)]
setnames(Ren_df, old = c("UniProt", "OR (95%)", "P"),
         new = c("UniProt", "OR", "P"))
Ren_df[, (names(Ren_df)) := lapply(.SD, function(col) {
  if (is.character(col)) gsub("^\\p{Z}+|\\p{Z}+$", "", col, perl = TRUE) else col
})]
all_joined_3 <- inner_join(all, Ren_df, by = "UniProt")
all_joined_3$Study <- 'Ren et.al 2023'


# Desai
Desai_df <- as.data.table(read_excel("Figures_Tables/reported_proteins/Desai2023eBioMedicine.xlsx"))
Desai_df$Study<- 'NA'
Desai_df$Study[c(1:17)] <- 'Deasi et.al 2023 Overall'
Desai_df$Study[c(18:28)] <- 'Deasi et.al 2023 Early onset'
Desai_df$Study[c(29:35)] <- 'Deasi et.al 2023 Aggressive'


Desai_df <- Desai_df[, .(`Uniprot ID`, `Odds ratio (95% CI) practical`, `p-value (unadjusted)`, `Study`)]
setnames(Desai_df, old = c("Uniprot ID", "Odds ratio (95% CI) practical", "p-value (unadjusted)", 'Study'),
         new = c("UniProt", "OR", "P", 'Study'))
Desai_df[, (names(Desai_df)) := lapply(.SD, function(col) {
  if (is.character(col)) gsub("^\\p{Z}+|\\p{Z}+$", "", col, perl = TRUE) else col
})]
all_joined_4 <- inner_join(all, Desai_df, by = "UniProt")



# JianWu
JianWu_df <- as.data.table(read_excel("Figures_Tables/reported_proteins/JianWu2025HumanGenomic.xlsx", 
                                      sheet = "TableS1PRACTICAL"))
JianWu_df$Z <- JianWu_df$beta/JianWu_df$se
JianWu_df <- JianWu_df[, .(`exposure`, `Z`, `pval`)]
setnames(JianWu_df, old = c("exposure", "Z", "pval"),
         new = c("EntrezGeneSymbol", "Z", "P"))
JianWu_df <- filter(JianWu_df, JianWu_df$P < 3.3e-5)      
JianWu_df[, (names(JianWu_df)) := lapply(.SD, function(col) {
  if (is.character(col)) gsub("^\\p{Z}+|\\p{Z}+$", "", col, perl = TRUE) else col
})]
all_joined_5 <- inner_join(all, JianWu_df, by = "EntrezGeneSymbol")
all_joined_5$Study <- "Wu et.al 2025"




final_col_order <- c("ID", "PWAS.Z", "PWAS.P", "SomaId",
                     "TargetFullName", "Target", "EntrezGeneSymbol", "UniProt","Race",
                     "OR", "Z", "P", "Study")


standardize_df <- function(df) {

  if (!"OR" %in% names(df)) df[, OR := NA_character_]

  if (!"Z" %in% names(df)) df[, Z := NA_real_]

  if (!"P" %in% names(df)) df[, P := NA_real_]
  if (!"Study" %in% names(df)) df[, Study := NA_character_]
  
  df <- df[, ..final_col_order]
  return(df)
}


df_list <- lapply(list(all_joined_1, all_joined_2, all_joined_3, all_joined_4, all_joined_5), standardize_df)


final_df <- rbindlist(df_list, use.names = TRUE, fill = TRUE)

final_df <- final_df[order(final_df$PWAS.P, ID),]
fwrite(final_df, 'Compared.txt', sep = '\t')

