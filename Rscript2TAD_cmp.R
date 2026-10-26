

library(dplyr)
library(SpectralTAD)
library(TADCompare)
library(ggplot2)

# Read in data (E14-E18)
cool_CEBP <- read.table("/home/orsalia/Desktop/HiC_PHiC/Final/TAD_Compare/chr1/CEBP_chr1.txt")
cool_HNF4A <- read.table("/home/orsalia/Desktop/HiC_PHiC/Final/TAD_Compare/chr1/HNF4A_chr1.txt")
cool_P60 <- read.table("/home/orsalia/Desktop/HiC_PHiC/Final/TAD_Compare/chr1/P60_chr1.txt")


# Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
sparse_CEBP <- HiCcompare::cooler2sparse(cool_CEBP)
sparse_HNF4A <- HiCcompare::cooler2sparse(cool_HNF4A)
sparse_P60 <- HiCcompare::cooler2sparse(cool_P60)


plot_CEBP_P60 = TADCompare(cont_mat1 = sparse_CEBP, cont_mat2 = sparse_P60, resolution = 20000, z_thresh = 1.5) 
plot_CEBP_P60
write.csv(plot_CEBP_P60$TAD_Frame,"CEBP_P60intertablechr1.csv")

plot_HNF4A_P60 = TADCompare(cont_mat1 = sparse_HNF4A, cont_mat2 = sparse_P60, resolution = 20000, z_thresh = 1.5) 
plot_HNF4A_P60
write.csv(plot_HNF4A_P60$TAD_Frame,"HNF4A_P60intertablechr1.csv")

# Visualizing the results
DiffPlot(tad_diff    = TD_Compare, 
         cont_mat1   = sparse_E14,
         cont_mat2   = sparse_E18,
         resolution  = 20000,
         start_coord = 115000000,
         end_coord   = 125000000,
         show_types  = TRUE, 
         point_size  = 3, palette = "RdYlBu",
         rel_heights = c(1, 2))
