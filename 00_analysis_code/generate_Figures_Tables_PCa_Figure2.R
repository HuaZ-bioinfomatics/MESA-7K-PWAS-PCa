# First, make sure you have the plotrix package installed. If not, run this line once:
# install.packages("plotrix")

# ----------- cis+trans
library(OmicCircos)
data(UCSC.hg19) 
library(data.table)
library(dplyr)
library(circlize)
library(stringr)
library(plotrix) # <-- 1. LOAD THE PLOTRIX LIBRARY

# --- [Your data loading and processing code remains the same] ---
# (Code from fread('04_association/...') down to ordering all_sig)

AFA <- fread('04_association/association_prostate_cancer_2023_NG_AFA/all_association.out', data.table = F)
sig1_pvalue <- max(AFA[AFA$FDR < 0.05, "PWAS.P"])
sig1 <- -log10(sig1_pvalue)
AFA <- select(AFA, ID, PWAS.P, FDR)
AFA$Race <- 'AFA'

CAU <- fread('04_association/association_prostate_cancer_2023_NG_CAU/all_association.out', data.table = F)
sig2_pvalue <- max(CAU[CAU$FDR < 0.05, "PWAS.P"])
sig2 <- -log10(sig2_pvalue)
CAU <- select(CAU, ID, PWAS.P, FDR)
CAU$Race <- 'CAU'

CHN <- fread('04_association/association_prostate_cancer_2023_NG_CHN/all_association.out', data.table = F)
sig3_pvalue <- max(CHN[CHN$FDR < 0.05, "PWAS.P"])
sig3 <- -log10(sig3_pvalue)
CHN <- select(CHN, ID, PWAS.P, FDR)
CHN$Race <- 'CHN'

HIS <- fread('04_association/association_prostate_cancer_2023_NG_HIS/all_association.out', data.table = F)
sig4_pvalue <- max(HIS[HIS$FDR < 0.05, "PWAS.P"])
sig4 <- -log10(sig4_pvalue)
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

all_sig <- select(all_sig, chromosome_name, TSS, EntrezGeneSymbol, P, color)
colnames(all_sig)[1:3] <- c('chr', 'po', 'Gene')
all_sig$chromosome <- as.numeric(gsub('chr', '', all_sig$chr))
all_sig <- all_sig[order(all_sig$chromosome, all_sig$po),]

# --- [Plotting Code] ---

pdffile  <- "Figures_Tables/Figure2.pdf";
pdf(pdffile, 8, 8);
par(mar=c(1, 1, 1, 1));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", xpd=1);
legend(350,430, box.col = "white",
       bg ="white",  cex = 0.5,
       legend=c("African", "Asian", 'European', 'Hispanic/Latino'),
       fill = c("#BC3C29", '#0072B5', "#E18727",  '#20854E'),border = c("#BC3C29", '#0072B5', "#E18727",  '#20854E'))
circos(R=290, cir=UCSC.hg19, type="chr", W=10, col = 'black')

# Track definitions
track_R_AFA <- 50
track_W_AFA <- 60
track_R_CHN <- 110
track_W_CHN <- 60
track_R_CAU <- 170
track_W_CAU <- 60
track_R_HIS <- 230
track_W_HIS <- 60

# First, draw the label track
circos(R=320, cir=UCSC.hg19, W=10, mapping=all_sig, col.v=4, type="label", side="out", cex=0.4, col = all_sig$color);

# Now, draw the p-value tracks
circos(R=track_R_AFA, cir=UCSC.hg19, W=track_W_AFA, mapping=pvalue_AFA, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#BC3C29')
circos(R=track_R_CHN, cir=UCSC.hg19, W=track_W_CHN, mapping=pvalue_CHN, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#0072B5')
circos(R=track_R_CAU, cir=UCSC.hg19, W=track_W_CAU, mapping=pvalue_CAU, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#E18727')
circos(R=track_R_HIS, cir=UCSC.hg19, W=track_W_HIS, mapping=pvalue_HIS, col.v=4, type="b", B=F, scale=TRUE, cex=0.5, col = '#20854E')


# Calculate p-value ranges for each track
min_p_AFA <- min(pvalue_AFA$pvalue, na.rm = TRUE)
max_p_AFA <- max(pvalue_AFA$pvalue, na.rm = TRUE)

min_p_CHN <- min(pvalue_CHN$pvalue, na.rm = TRUE)
max_p_CHN <- max(pvalue_CHN$pvalue, na.rm = TRUE)

min_p_CAU <- min(pvalue_CAU$pvalue, na.rm = TRUE)
max_p_CAU <- max(pvalue_CAU$pvalue, na.rm = TRUE)

min_p_HIS <- min(pvalue_HIS$pvalue, na.rm = TRUE)
max_p_HIS <- max(pvalue_HIS$pvalue, na.rm = TRUE)


# Function to map a data value to a radial position within a track
map_to_radius <- function(value, min_data, max_data, inner_R, width_W) {
  if (is.infinite(value) || is.na(value)) { # Handle cases where significance threshold is outside data range
      return(NA) 
  }
  if (max_data == min_data) {
    return(inner_R + width_W / 2)
  }
  radius <- inner_R + ((value - min_data) * width_W) / (max_data - min_data)
  return(radius)
}

# Calculate the radii for the significance lines
r1 <- map_to_radius(sig1, min_p_AFA, max_p_AFA, track_R_AFA, track_W_AFA)
r2 <- map_to_radius(sig3, min_p_CHN, max_p_CHN, track_R_CHN, track_W_CHN) # <-- Mapped sig3 to CHN track
r3 <- map_to_radius(sig2, min_p_CAU, max_p_CAU, track_R_CAU, track_W_CAU) # <-- Mapped sig2 to CAU track
r4 <- map_to_radius(sig4, min_p_HIS, max_p_HIS, track_R_HIS, track_W_HIS)


# --- CORRECTED DRAWING CALLS ---
# Define plot center
plot_center_x <- 400.5
plot_center_y <- 400.5

# 2. Use draw.circle, correct center coordinates, and updated arguments
if(!is.na(r1)) draw.circle(x = plot_center_x, y = plot_center_y, radius = r1, lty = 'dashed', border = 'red', lwd = 0.5)
if(!is.na(r2)) draw.circle(x = plot_center_x, y = plot_center_y, radius = r2, lty = 'dashed', border = 'red', lwd = 0.5)
if(!is.na(r3)) draw.circle(x = plot_center_x, y = plot_center_y, radius = r3, lty = 'dashed', border = 'red', lwd = 0.5)
if(!is.na(r4)) draw.circle(x = plot_center_x, y = plot_center_y, radius = r4, lty = 'dashed', border = 'red', lwd = 0.5)

dev.off()