rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(tidyverse) # for data cleaning
library(pheatmap) # for making heatmap
library(viridis) # for heatmap color scheme
library(table1) # for making cytokine summary table
library(ggfortify) # for PCA
library(factoextra) # for PCA 
library(ggplot2) # for plotting
library(readxl) # for data import
library(ggpubr) # for making figure panel

# Import data
MDMCytokineData <- read.csv("2021_05_20 MSD MDM Table.csv")
MDMCytokineData_df <- data.frame(MDMCytokineData, row.names = MDMCytokineData$Clean.Sample.ID)
MDMCytokineData_df <- MDMCytokineData_df[, 3:30]

# Create summary data frame with means by group
MDM_mean_df = MDMCytokineData_df %>% group_by(Polarization) %>% summarize(across(everything(), mean))
MDM_mean_df <- data.frame(MDM_mean_df[,-1], row.names = MDM_mean_df$Polarization)

###############################################
#### 1. Summary table
###############################################

# Function so that mean (SEM) can be shown in table
my.render.cont.mean.sem <- function(x) {
  s <- stats.default(x)
  s$SEM <- with(s, SD/sqrt(N))
  with(stats.apply.rounding(s), c("",
                                  "Mean (SEM)"=sprintf("%s (%s)", MEAN, SEM)))
}

# Example code that will create a list that can be copied and pasted into the table1 function
paste0(" ", names(MDMCytokineData_df[2:ncol(MDMCytokineData_df)]), " ", collapse="+")

# Make table. Row labels, including Greek symbols, as well as excess spacing, were fixed in Word. 
# Significance markers were also added in Word using output from analysis in GraphPad Prism 9. 
table1(~ Eotaxin + Eotaxin.3 + GM.CSF + IL.1a + IL.1B + IL.2 + IL.5 + IL.6 + IL.7 + IL.8 + IL.10 
       + IL.12p40 + IL.12p70 + IL.13 + IL.15 + IL.16 + IL.17 + IP.10 + MCP.1 + MCP.4 + MDC
       + MIP.1a + MIP.1B + TARC + TNF.a + TNF.B + VEGF | Polarization, 
       data = MDMCytokineData_df, 
       render.continuous = my.render.cont.mean.sem,
       render.missing = NULL,
       overall = NULL)

###############################################
#### 2. Heat map
###############################################

# Delete periods from all of the columns
names(MDM_mean_df) <- gsub("\\.", "", names(MDM_mean_df))

# Turn data frame into matrix for pheatmap input
MDM_mean_dm <- data.matrix(MDM_mean_df)

# Make heat map 
pheatmap(t(MDM_mean_dm), 
         color = viridis(100), # sets color scheme
         angle_col = c("0"), # makes column labels horizontal
         cellwidth = 75, # sets dimensions of cells so that they don't change with viewer pane size
         cellheight = 10, # sets dimensions of cells so that they don't change with viewer pane size
         border_color = "black", # adds black border around cells
         treeheight_col = 10, # sets dims of trees so that they don't change with viewer pane size
         fontsize_row = 9, # sets dims of trees so that they don't change with viewer pane size
         scale = 'row', # scales data by row
         fontsize_col = 12, # sets font size for column labels
         cutree_rows = 5) # indicates how many clusters to show separated by spaces

###############################################
#### 3. PCA
###############################################

# Remove periods from column names
names(MDMCytokineData_df) <- gsub("\\.", "", names(MDMCytokineData_df))

# Create a data frame with only mediator columns
MDMDataPCA <- MDMCytokineData_df[, 2:28]

# Log transform mediator data so that PCA plot isn't as bunched up.
# Pseudocount of 1 added because some values are zero and you can't take the log of 0.
MDMDataPCA_log <- log(MDMDataPCA+1)

# Run PCA
# Centering ensures components are only looking at variance within the dataset.
# Scaling brings all variables to the same magnitude so that high abundance variables 
# do not influence results more than low abundance variables. 
pca.res <- prcomp(MDMDataPCA_log, center = TRUE, scale = TRUE)

# Set theme
theme_set(theme_bw())

# Cluster plot. Export as 3 x 4 PDF.
fviz_pca_ind(pca.res,
             label = "none",
             habillage = MDMCytokineData_df$Polarization, 
             palette = c("#0A0A0A", "#0253F5", "#02D40D"),
             addEllipses = TRUE) + 
  ggtitle("PCA - Clustering") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.border = element_rect(fill = NA, color = "black", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", 
                                         size = 0.3, 
                                         linetype = "solid",
                                         color = "black"),
        legend.text = element_text(size = 10)) +
  xlim(-10, 10)

# To decide which contribution variables to label in M1 direction, 
# took top 5 from axes 1 and 2 contributions.
fviz_contrib(pca.res, choice = "var", axes = 1, top = 10)
fviz_contrib(pca.res, choice = "var", axes = 2, top = 10)

# Choose which variables to show on label by creating vector. 
# Did this because all of the M1-contributing variables overlapped a lot.
labels <- c("TARC", "MCP4", "MDC", "IL7", "MIP1a", "TNFa", "IL8", "IL6", "IL15", "IL10", "Eotaxin3", "IL5")

# Contribution plot. Export as 3 x 5 PDF.
fviz_pca_var(pca.res, col.var = "contrib",
             repel = TRUE, col.circle = NA, labelsize = 4,
             select.var = list(name = labels)) +
  ggtitle("PCA - Contributions") +
  scale_color_viridis(begin = 0, end = 0.6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.border = element_rect(fill = NA, color = "black", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10), 
        aspect.ratio=1/1.75) +
  labs(color = "Contribution") +
  ylim(-1.25, 1.25) +
  xlim(-1.25, 1.25)

###############################################
#### 4. Bioenergetics PCA
###############################################

# Import data
bioenergetics_data <- read_excel("2022_05_27 hMDM Bioenergetics Data.xlsx")
bioenergetics_data <- data.frame(bioenergetics_data, row.names = bioenergetics_data$RowID)
bioenergetics_data <- bioenergetics_data[, -1]

# Set Factors
bioenergetics_data$Sample.Type <- factor(bioenergetics_data$Sample.Type, 
                                         levels = c("M0 hMDM", "M1 hMDM", "M2 hMDM"))

# Create a data frame with only bioenergetic columns
bioenergetics_data_pca <- bioenergetics_data[, -1]

# Run PCA on each set of data
# Centering ensures components are only looking at variance within the dataset.
# Scaling brings all variables to the same magnitude so that high abundance variables 
# do not influence results more than low abundance variables. 
pca.res <- prcomp(bioenergetics_data_pca, center = TRUE, scale = TRUE)

# Set theme
theme_set(theme_bw())

# Cluster Plots
bioenergetics.pca <- fviz_pca_ind(pca.res,
             label = "none",
             habillage = bioenergetics_data$Sample.Type, 
             palette = c("#0A0A0A", "#0253F5", "#02D40D"),
             addEllipses = TRUE) + 
  ggtitle("PCA - Bioenergetics Parameters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.border = element_rect(fill = NA, color = "black", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", 
                                         size = 0.3, 
                                         linetype = "solid",
                                         color = "black"),
        legend.text = element_text(size = 7)) +
  xlim(-10, 10)

# Contribution Plot
bioenergetics.contrib <- fviz_pca_var(pca.res, col.var = "contrib",
             repel = TRUE, col.circle = NA, labelsize = 4) +
  ggtitle("PCA - Contributions of Bioenergetic Parameters") +
  scale_color_viridis(begin = 0, end = 0.6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.border = element_rect(fill = NA, color = "black", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7), 
        aspect.ratio=1/1.75) +
  labs(color = "Contribution") +
  ylim(-1.25, 1.25) +
  xlim(-1.25, 1.25)

# Make Panel
dev.off()

pdf(file = "hMDM Bioenergetics PCA.pdf",
    width = 10,
    height = 5)

ggarrange(bioenergetics.pca, bioenergetics.contrib,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()
