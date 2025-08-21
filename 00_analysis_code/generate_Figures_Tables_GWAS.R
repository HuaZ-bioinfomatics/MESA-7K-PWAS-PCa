library(data.table)
library(dplyr)
library(writexl)

# Step 1: Load meta-analysis result
meta_df <- fread('05_meta_analysis/combine_meta_ethinc_specific_Sig.txt')

# Step 2: Load and deduplicate annotation by AptName
annot <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt') %>%
  distinct(AptName, .keep_all = TRUE)

# Step 3: Select and rename needed columns
annot_clean <- annot %>%
  select(AptName, chromosome_name, start_position, end_position) %>%
  rename(chr = chromosome_name,
         gene_start = start_position,
         gene_end = end_position) %>%
  mutate(chr = as.character(chr))

# Step 4: Load and clean GWAS loci (GRCh38)
gwas_df <- fread('Figures_Tables/TableS4_from_AnqiWang.txt') %>%
  filter(!is.na(`Position (GRCh38)`)) %>%
  rename(rsID = `rsID`,
         chr = Chromosome,
         pos = `Position (GRCh38)`) %>%
  mutate(chr = as.character(chr))

# Step 5: Join gene position info to meta_df by ID
meta_df <- meta_df %>%
  left_join(annot_clean, by = c("ID" = "AptName"))

# Step 6: Define function using gene start and end
get_gwas_overlap_info <- function(chr, start_pos, end_pos) {
  if (is.na(chr) | is.na(start_pos) | is.na(end_pos)) return(c("No", NA, NA))
  
  subset <- gwas_df %>% filter(chr == chr)
  if (nrow(subset) == 0) return(c("No", NA, NA))
  
  # Check overlap with gene ±500kb
  overlap <- subset %>%
    filter(pos >= (start_pos - 500000) & pos <= (end_pos + 500000))
  
  # Distance to gene region
  subset <- subset %>%
    mutate(distance = pmin(abs(pos - start_pos), abs(pos - end_pos)))
  
  nearest <- subset %>%
    filter(distance == min(distance)) %>%
    slice(1)
  
  if (nrow(overlap) > 0) {
    return(c("Yes", nearest$rsID, round(nearest$distance / 1000, 2)))
  } else {
    return(c("No", nearest$rsID, round(nearest$distance / 1000, 2)))
  }
}

# Step 7: Apply the function
gwas_info <- t(mapply(get_gwas_overlap_info, meta_df$chr, meta_df$gene_start, meta_df$gene_end))
meta_df$GWAS_overlap <- gwas_info[, 1]
meta_df$Nearest_GWAS_SNP <- gwas_info[, 2]
meta_df$Distance_to_SNP_kb <- as.numeric(gwas_info[, 3])

# Step 8: Export updated table
write_xlsx(as.data.frame(meta_df), "Figures_Tables/TableS2_with_GWAS_overlap_minimal.xlsx")
