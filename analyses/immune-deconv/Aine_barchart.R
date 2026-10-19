x <- readRDS("analyses/immune-deconv/results/quantiseq_output.rds")
xcell <- readRDS("analyses/immune-deconv/results/xcell_output.rds")

#print(names(x))

### MB rows
mb <- x[grepl("Medulloblastoma", x$cancer_group, ignore.case = TRUE), ]
# print(head(mb, n=20))
#write.table(mb, file = "analyses/immune-deconv/aine_scratch/MB_xcell.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#mini_MB <- read.table("analyses/immune-deconv/aine_scratch/MB_xcell.tsv", sep = "\t", header = TRUE)


### unique samples per subtype
# unique_MB_types <- unique(mb[c("molecular_subtype")])
# unique_MB_samples2 <- unique(mb[c("Kids_First_Biospecimen_ID")])
# unique_MB_samples <- unique(mb[c("Kids_First_Biospecimen_ID","molecular_subtype")])
#n_1 <- nrow(unique_MB_samples)
#n_2 <- nrow(unique_MB_samples2)
#print(paste0(n_1, ", ", n_2))   
# they match! yippee QC!

#print(table(unique_MB_samples$molecular_subtype, useNA = "ifany"))

### Result: same across both
# MB, Group3 MB, Group4    MB, SHH    MB, WNT 
#         68        102         84         24 

available_cell_types <- (unique(mb$cell_type))
#print(available_cell_types)

### Result for xcell:
#  [1] "Myeloid dendritic cell activated" "B cell"                          
#  [3] "T cell CD4+ memory"               "T cell CD4+ naive"               
#  [5] "T cell CD4+ (non-regulatory)"     "T cell CD4+ central memory"      
#  [7] "T cell CD4+ effector memory"      "T cell CD8+ naive"               
#  [9] "T cell CD8+"                      "T cell CD8+ central memory"      
# [11] "T cell CD8+ effector memory"      "Class-switched memory B cell"    
# [13] "Common lymphoid progenitor"       "Common myeloid progenitor"       
# [15] "Myeloid dendritic cell"           "Endothelial cell"                
# [17] "Eosinophil"                       "Cancer associated fibroblast"    
# [19] "Granulocyte-monocyte progenitor"  "Hematopoietic stem cell"         
# [21] "Macrophage"                       "Macrophage M1"                   
# [23] "Macrophage M2"                    "Mast cell"                       
# [25] "B cell memory"                    "Monocyte"                        
# [27] "B cell naive"                     "Neutrophil"                      
# [29] "NK cell"                          "T cell NK"                       
# [31] "Plasmacytoid dendritic cell"      "B cell plasma"                   
# [33] "T cell gamma delta"               "T cell CD4+ Th1"                 
# [35] "T cell CD4+ Th2"                  "T cell regulatory (Tregs)"       
# [37] "immune score"                     "stroma score"                    
# [39] "microenvironment score"          

### Quaniseq
#  [1] "B cell"                       "Macrophage M1"               
#  [3] "Macrophage M2"                "Monocyte"                    
#  [5] "Neutrophil"                   "NK cell"                     
#  [7] "T cell CD4+ (non-regulatory)" "T cell CD8+"                 
#  [9] "T cell regulatory (Tregs)"    "Myeloid dendritic cell"      
# [11] "uncharacterized cell"        

mini_metadata <- unique(mb[c("Kids_First_Biospecimen_ID", "molecular_subtype", "cohort")])
#print(with(mini_metadata, table(molecular_subtype, cohort, useNA = "ifany")))

# Results: same across both
#                  cohort
# molecular_subtype PBTA PPTC
#        MB, Group3   60    8
#        MB, Group4  100    2
#        MB, SHH      78    6
#        MB, WNT      22    2

# there's like nothing in PPTC
mb_2 <- mb[mb$cohort == "PBTA", ]
mb_unique_samples <- unique(mb_2[c("Kids_First_Biospecimen_ID", "molecular_subtype")])
#print(table(mb_unique_samples$molecular_subtype)) 

# --> 260 samples
# MB, Group3 MB, Group4    MB, SHH    MB, WNT 
#         60        100         78         22 

# ------------------------------ PLOTTING -------------------------------- #
# Plan: make one of those stacked bar charts showing immune composition!
# except using quantiseq instead of xcell so i can go straight to absolute fractions...
# 4 bars; one per subtype
# each bar is the same height; but will also list number of contributing samples on top

library(ggplot2)
library(dplyr)

subtypes <- c("MB, WNT", "MB, SHH", "MB, Group3", "MB, Group4")

immune_samples <- mb_2 %>%
  filter(molecular_subtype %in% subtypes, sample_type == "Tumor", cell_type != "uncharacterized cell")

stopifnot(!anyDuplicated(immune_samples[c("Kids_First_Biospecimen_ID", "cell_type")]), all(is.finite(immune_samples$fraction)), all(immune_samples$fraction >= 0))

immune_samples <- immune_samples %>%
  group_by(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  mutate(total_immune = sum(fraction),n_cell_types = n_distinct(cell_type)) %>%
  ungroup()

sample_checks <- immune_samples %>%
  distinct(Kids_First_Biospecimen_ID, molecular_subtype, total_immune, n_cell_types)

# print(table(sample_checks$n_cell_types))
# print(table(sample_checks$total_immune == 0))

stopifnot(all(sample_checks$n_cell_types == 10))

immune_samples <- immune_samples %>%
  filter(total_immune > 0) %>%
  mutate(proportion = fraction / total_immune)

immune <- immune_samples %>%
  group_by(molecular_subtype, cell_type) %>%
  summarise(mean_proportion = mean(proportion), .groups = "drop")

immune %>%
  group_by(molecular_subtype) %>%
  summarise(total = sum(mean_proportion), .groups = "drop") %>%
  print()

# Result:
#   molecular_subtype total
#   <chr>             <dbl>
# 1 MB, Group3            1
# 2 MB, Group4            1
# 3 MB, SHH               1
# 4 MB, WNT               1

# ============================ ACTUAL PLOTTING ============================ #
immune$molecular_subtype <- factor(immune$molecular_subtype, levels = subtypes)

sample_counts <- immune_samples %>%
  distinct(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  count(molecular_subtype, name = "n")

sample_counts$molecular_subtype <- factor(sample_counts$molecular_subtype, levels = subtypes)

plot <- ggplot(immune, aes(x = molecular_subtype, y = mean_proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = sample_counts, aes(x = molecular_subtype, y = 1.02, label = paste0("n = ", n)), size = 3, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_discrete(labels = c(
    "MB, WNT" = "WNT",
    "MB, SHH" = "SHH",
    "MB, Group3" = "Group 3",
    "MB, Group4" = "Group 4")) +
  labs(x = "MB Subtype", y = "Mean Proportion of Immune Cells", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))

print(plot)
#dir.create("analyses/immune-deconv/aine_scratch/plots")
#ggsave("analyses/immune-deconv/aine_scratch/plots/MB_immune_composition_quantiseq.png", plot = plot, width = 8, height = 6, dpi = 300)
ggsave("analyses/immune-deconv/results/Aine_MB_immune_composition.png", plot = plot, width = 8, height = 6, dpi = 300)