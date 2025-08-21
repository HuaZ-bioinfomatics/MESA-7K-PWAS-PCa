library(data.table)
library(dplyr)
library(ggvenn)
library(ggpubr)

dir_name <- 'Figures_Tables'
if(!exists(dir_name)){
    dir.create(dir_name)
}

#############################################
################## Figure 1 #################
#############################################
#---------- venn plot for models with Rsq >= 0.01
# Rsq 0.01 models venn plot
# AFA
AFA_models <- fread('03_establish_models/AFR_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
AFA_models <- AFA_models$model
# CHN
CHN_models <- fread('03_establish_models/CHN_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
CHN_models <- CHN_models$model
# CAU
CAU_models <- fread('03_establish_models/CAU_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
CAU_models <- CAU_models$model
# HIS
HIS_models <- fread('03_establish_models/HIS_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
HIS_models <- HIS_models$model


x <- list(
  'African\n(1,578)' = AFA_models,
  'Asian (1,218)' = CHN_models,
  'European (1,993)' = CAU_models,
  'Hispanic/\nLatino\n(1,390)' = HIS_models
)

# pdf(file="Figures_Tables/prostate_cancer_venn_models_Rsq_0.01.pdf", width = 8, height = 8)
venn_plot <- ggvenn(
  x, 
  fill_color = c('#BC3C29', '#0072B5', '#E18727', '#20854E'),
  stroke_size = 0.8, stroke_color = 'white', set_name_size = 5, digits = 2, fill_alpha = 0.6
)
# dev.off()



#------ Rsq plot
library(data.table)
library(ggplot2)
library(dplyr)
races <- c('AFA', 'CHN', 'CAU', 'HIS')
merge <- data.frame()
for (race in races){
    if (race == 'AFA'){
      df <- fread(paste0('03_establish_models/AFR_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt'), data.table = F)
    }else{
      df <- fread(paste0('03_establish_models/', race, '_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt'), data.table = F)

    }
    df <- select(df, Rsq)
    df$race <- race
    merge <- rbind(merge, df)
}
# replace race name
merge <- merge %>%
  mutate(race = case_when(
    race == 'AFA' ~ 'African',
    race == 'CHN' ~ 'Asian',
    race == 'CAU' ~ 'European',
    race == 'HIS' ~ 'Hispanic/Latino',
    TRUE ~ race  # If none of the above conditions are met, keep the original value
  ))
mean_values <- merge %>%
  group_by(race) %>%
  summarise(mean_Rsq = mean(Rsq))




merge$race <- factor(merge$race , levels=c("African", "Asian", "European", "Hispanic/Latino"))

Rsq_plot <- ggplot(data = merge, aes(race, Rsq)) + 
  ggdist::stat_halfeye(aes(color=race,fill=race),adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA, show.legend = FALSE) + 
  geom_boxplot(aes(color=race),width = .2, outlier.shape = NA, show.legend = FALSE) + 
  coord_cartesian(ylim=c(0,1))+
#   geom_jitter(aes(color=race),width = .05, alpha = .3) +
  scale_fill_manual(values=c('#BC3C29', '#0072B5', '#E18727', '#20854E'))+
  scale_color_manual(values=c('#BC3C29', '#0072B5', '#E18727', '#20854E'))+
  ylab(model ~ performance ~ (italic(R)^2))+ xlab ('')+
  geom_text(data = mean_values, aes(label=sprintf("%.2f", mean_Rsq), y =mean_Rsq), color = 'white', vjust = 4.2, hjust = -0.7, show.legend = FALSE) +
  theme_classic()



# ------------- pie chart
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggpmisc)
mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
# AFA
AFA <- fread('03_establish_models/AFR_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
count.data <- data.frame(
  class = c("Cis", "Trans", "Cis+Trans"),
  n = c(nrow(filter(AFA, AFA$N_trans == 0)), nrow(filter(AFA, AFA$N_cis == 0)), nrow(filter(AFA, AFA$N_trans != 0 & AFA$N_cis != 0 )))
)
count.data$prop <- count.data$n/sum(count.data$n) * 100
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
pie_AFA <- ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white", show.legend = FALSE) +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(class, '\n', n)), size = 3, color = "black")+
  scale_fill_manual(values = c('#bc3c29', '#d07669', '#e4b1a9')) +
  theme_void()

# CHN
CHN <- fread('03_establish_models/CHN_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
count.data <- data.frame(
  class = c("Cis", "Trans", "Cis+Trans"),
  n = c(nrow(filter(CHN, CHN$N_trans == 0)), nrow(filter(CHN, CHN$N_cis == 0)), nrow(filter(CHN, CHN$N_trans != 0 & CHN$N_cis != 0 )))
)
count.data$prop <- count.data$n/sum(count.data$n) * 100
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
pie_CHN <- ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white", show.legend = FALSE) +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(class, '\n', n)), size = 3, color = "black")+
  scale_fill_manual(values = c('#0072B5', '#4c9ccb', '#99c6e1')) +
  theme_void()

# CAU
CAU <- fread('03_establish_models/CAU_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
count.data <- data.frame(
  class = c("Cis", "Trans", "Cis+Trans"),
  n = c(nrow(filter(CAU, CAU$N_trans == 0)), nrow(filter(CAU, CAU$N_cis == 0)), nrow(filter(CAU, CAU$N_trans != 0 & CAU$N_cis != 0 )))
)
count.data$prop <- count.data$n/sum(count.data$n) * 100
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
pie_CAU <- ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white", show.legend = FALSE) +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(class, '\n', n)), size = 3, color = "black")+
  scale_fill_manual(values = c('#E18727', '#eaab67', '#f3cfa8')) +
  theme_void()

# HIS
HIS <- fread('03_establish_models/HIS_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
count.data <- data.frame(
  class = c("Cis", "Trans", "Cis+Trans"),
  n = c(nrow(filter(HIS, HIS$N_trans == 0)), nrow(filter(HIS, HIS$N_cis == 0)), nrow(filter(HIS, HIS$N_trans != 0 & HIS$N_cis != 0 )))
)
count.data$prop <- count.data$n/sum(count.data$n) * 100
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
pie_HIS <- ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white", show.legend = FALSE) +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = paste0(class, '\n', n)), size = 3, color = "black")+
  scale_fill_manual(values = c('#20854E', '#62a983', '#a5ceb8')) +
  theme_void()

pie_plot <- ggarrange(pie_AFA, pie_CHN, pie_CAU, pie_HIS, ncol = 4, nrow = 1)

# ggsave('test.pdf', pie_plot, width = 8, height = 2)

plot_b <- Rsq_plot+
  geom_plot(data=tibble(x=0.65,y=1.4,plot=list(pie_plot)),
            aes(x=x,y=y,label=plot),
            vp.width=0.89,vp.height=0.89)
        

p <- ggarrange(venn_plot, plot_b, ncol = 2, nrow = 1, labels = c('a', 'b'))


ggsave('Figures_Tables/Figure1.pdf', p, height = 8, width = 16)





#############################################
################## Figure 2 #################
#############################################
rm(list=ls())

# ----------- cis+trans
library(OmicCircos)
data(UCSC.hg19)
library(data.table)
library(dplyr)
library(stringr)
library(draw)
chr_info <- UCSC.hg19[c(1:22),c(1,6,7)]
chr_info$seg.name <- gsub('chr', '', chr_info$seg.name)
chr_info$the.v <- 'NA'
chr_info$NO <- 'NA'
chr_info <- segAnglePo(chr_info, seg=chr_info[[1]])
chr_info <- data.frame(chr_info)
chr_info$seg.name <- paste0('chr', chr_info$seg.name)
chr_info$seg.name <- as.factor(chr_info$seg.name)
UCSC.hg19 <- chr_info



# p value

AFA <- fread('04_association/association_prostate_cancer_2023_NG_AFA/all_association.out', data.table = F)
sig1 <- -log10(last(AFA[AFA$FDR < 0.05,])$PWAS.P)
AFA <- select(AFA, ID, PWAS.P, FDR)
AFA$Race <- 'AFA'


CAU <- fread('04_association/association_prostate_cancer_2023_NG_CAU/all_association.out', data.table = F)
sig2 <- -log10(last(CAU[CAU$FDR < 0.05,])$PWAS.P)
CAU <- select(CAU, ID, PWAS.P, FDR)
CAU$Race <- 'CAU'


CHN <- fread('04_association/association_prostate_cancer_2023_NG_CHN/all_association.out', data.table = F)
sig3 <- -log10(last(CHN[CHN$FDR < 0.05,])$PWAS.P)
CHN <- select(CHN, ID, PWAS.P, FDR)
CHN$Race <- 'CHN'


HIS <- fread('04_association/association_prostate_cancer_2023_NG_HIS/all_association.out', data.table = F)
sig4 <- -log10(last(HIS[HIS$FDR < 0.05,])$PWAS.P)
HIS <- select(HIS, ID, PWAS.P, FDR)
HIS$Race <- 'HIS'



annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, AptName, EntrezGeneSymbol, chromosome_name, TSS)
annotation <- annotation[!duplicated(annotation$AptName),]
all <- rbind(AFA, CAU, CHN, HIS)
all <- left_join(all, annotation, by = c('ID' = 'AptName'))
all$chromosome_name <- paste0('chr', all$chromosome_name)
all$color <- ifelse(all$Race == 'AFA', '#BC3C29',
                    ifelse(all$Race == 'CAU', '#E18727',
                            ifelse(all$Race == 'CHN', '#0072B5', '#20854E')))

all$P <- -log10(all$PWAS.P)

pvalue <- select(all, chromosome_name, TSS, EntrezGeneSymbol, PWAS.P, Race)
pvalue$PWAS.P <- -log10(pvalue$PWAS.P)
pvalue_AFA <- filter(pvalue, Race == 'AFA')
pvalue_AFA <- pvalue_AFA[,1:4]
colnames(pvalue_AFA) <- c('chr', 'po', 'gene', 'pvalue')
pvalue_CAU <- filter(pvalue, Race == 'CAU')
pvalue_CAU <- pvalue_CAU[,1:4]
colnames(pvalue_CAU) <- c('chr', 'po', 'gene', 'pvalue')
pvalue_CHN <- filter(pvalue, Race == 'CHN')
pvalue_CHN <- pvalue_CHN[,1:4]
colnames(pvalue_CHN) <- c('chr', 'po', 'gene', 'pvalue')
pvalue_HIS <- filter(pvalue, Race == 'HIS')
pvalue_HIS <- pvalue_HIS[,1:4]
colnames(pvalue_HIS) <- c('chr', 'po', 'gene', 'pvalue')



AFA_sig <- all[all$FDR < 0.05 & all$Race == 'AFA',]
CAU_sig <- all[all$FDR < 0.05 & all$Race == 'CAU',]
CHN_sig <- all[all$FDR < 0.05 & all$Race == 'CHN',]
HIS_sig <- all[all$FDR < 0.05 & all$Race == 'HIS',]
all_sig <- rbind(AFA_sig, CAU_sig, CHN_sig, HIS_sig)

# # only select duplicated in 2 or more ethnic

# all_sig <- filter(all_sig, ID %in% all_sig[duplicated(all_sig$ID),]$ID)

all_sig <- select(all_sig, chromosome_name, TSS, EntrezGeneSymbol, P, color)
colnames(all_sig)[1:3] <- c('chr', 'po', 'Gene')
all_sig$chromosome <- as.numeric(gsub('chr', '', all_sig$chr))
all_sig <- all_sig[order(all_sig$chromosome, all_sig$po),]
# all_sig <- all_sig[!duplicated(all_sig[,c(1,2,3,5)]),]

pdffile  <- "Figures_Tables/Figure2.pdf";
pdf(pdffile, 8, 8);
par(mar=c(1, 1, 1, 1));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", xpd=1);
legend(350,430, box.col = "white",
       bg ="white",  cex = 0.5,
       legend=c("African", "Asian", 'European', 'Hispanic/Latino'),  
       fill = c("#BC3C29", '#0072B5', "#E18727",  '#20854E'),border = c("#BC3C29", '#0072B5', "#E18727",  '#20854E')) 
circos(R=290, cir=UCSC.hg19, type="chr", W=10, col = 'black')
circos(R=320, cir=UCSC.hg19, W=10, mapping=all_sig, col.v=4, type="label", side="out", cex=0.4, col = all_sig$color);
circos(R=50, cir=UCSC.hg19, W=60, mapping=pvalue_AFA, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#BC3C29')
circos(R=110, cir=UCSC.hg19, W=60, mapping=pvalue_CHN, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#0072B5')
circos(R=170, cir=UCSC.hg19, W=60, mapping=pvalue_CAU, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#E18727')
circos(R=230, cir=UCSC.hg19, W=60, mapping=pvalue_HIS, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#20854E')



start1 <- 0.545
end1 <- 0.86
start4 <- 2.12
end4 <- 2.443
start2 <- (start4 - start1)/3 + start1
start3 <- (start4 - start1)/3 + start2
end2 <- (end4 - end1)/3 + end1
end3 <- (end4 - end1)/3 + end2

length1 <- 13.3 # need change accordingly
length2 <- 20.1 # need change accordingly
length3 <- 273.7 # need change accordingly
length4 <- 16.6 # need change accordingly

r1 <- (end1-start1)/length1 * sig1 + start1
r2 <- (end2-start2)/length2 * sig2 + start2
r3 <- (end3-start3)/length3 * sig3 + start3
r4 <- (end4-start4)/length4 * sig4 + start4
r3 <- 1.606
drawCircle(x = 4, y = 4, radius = r1, lineType = 'dashed', lineColor = 'red', lineWidth = 0.5)
drawCircle(x = 4, y = 4, radius = r2, lineType = 'dashed', lineColor = 'red', lineWidth = 0.5)
drawCircle(x = 4, y = 4, radius = r3, lineType = 'dashed', lineColor = 'red', lineWidth = 0.5)
drawCircle(x = 4, y = 4, radius = r4, lineType = 'dashed', lineColor = 'red', lineWidth = 0.5)


dev.off()




#############################################
################## FigureS1 #################
#############################################
# ------
library(data.table)
library(dplyr)
library(ggpmisc)
library(ggplot2)
data <- read.table('07_validation_using_INTERVAL/validation_Rsq_External_Internal_CAU_male.txt',header = T,row.names  = 1)
data <- filter(data, External_Rsq >= 0.01)
my.formula <- y ~ x
Figure1B <- ggplot(data = data, aes(x = Internal_Rsq, y = External_Rsq)) +      
  geom_point(size=3,
             alpha=1,
             color="#E18727")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab(External~italic(R)^2~(EUR))+
  xlab(Cross~validation~italic(R)^2~(EUR))+
  # stat_poly_eq(use_label(c("eq"))) 
  annotate("text", x = 0.13, y = 0.75, label = 'y = 0.809 x - 0.000545',color="black")

ggsave('Figures_Tables/FigureS1.pdf', Figure1B, width = 6, height = 6)




#############################################
################## Table S1 #################
#############################################
library(data.table)
library(dplyr)
library(writexl)
res_list <- list()
races <- c('AFR', 'CAU', 'CHN', 'HIS')

for (race in races) {
  filepath <- paste0('03_establish_models/', race, '_male_5e-9/08_establish_prediction_models/model_Rsq_0.01.txt')
  tmp <- fread(filepath)
  
  race_label <- ifelse(race == 'AFR', 'AFA', race)
  tmp$race <- race_label
  
  res_list[[race]] <- tmp
}

res <- rbindlist(res_list)
res$ID <- gsub('.wgt.RDat$', '', res$model)
ann <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt')
ann <- ann[!duplicated(ann$AptName),]
ann <- select(ann, AptName, SomaId, TargetFullName, Target, EntrezGeneSymbol)
res <- left_join(res, ann, by = c('ID' = 'AptName'))
res <- select(res, model, ID, SomaId, TargetFullName, Target, EntrezGeneSymbol, best_model, Rsq, N_SNP, N_cis, N_trans, race)
write_xlsx(as.data.frame(res), "Figures_Tables/TableS1.xlsx")

#############################################
################## Table S2 #################
#############################################
library(data.table)
library(dplyr)
library(writexl)

# Step 1: Read in data
df <- fread('05_meta_analysis/combine_meta_ethinc_specific_Sig.txt')

# Step 2: Reorder columns
df <- df %>%
  relocate(TargetFullName, Target, EntrezGeneSymbol, .after = ID)

# Step 3: Define population codes and their labels
pop_codes <- c("AFA", "CHN", "CAU", "HIS")
pop_labels <- c("African", "Asian", "European", "Hispanic/Latino")  

# Step 4: Add population_specific column — FDR < 0.05 in any population
df <- df %>%
  rowwise() %>%
  mutate(population_specific = {
    fdr_vals <- c_across(all_of(paste0("FDR_", pop_codes)))

    # Include all populations where FDR < 0.05
    sig_flags <- pop_labels[!is.na(fdr_vals) & fdr_vals < 0.05]

    paste(sig_flags, collapse = ", ")
  }) %>%
  ungroup()

# Step 5: Create output directory if needed
if (!dir.exists("Figures_Tables")) dir.create("Figures_Tables")

# Step 6: Write to Excel
write_xlsx(as.data.frame(df), "Figures_Tables/TableS2.xlsx")







####### Table S3
library(data.table)
library(dplyr)
library(writexl)

df <- fread('Figures_Tables/TableS4_from_AnqiWang.txt')

  