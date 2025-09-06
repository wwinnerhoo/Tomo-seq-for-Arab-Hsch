
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../"))

rm(list = ls())

getwd()
source("dependency/0__Functions_publication.R")
dev.off()

library(viridis)  # Provides colorblind-friendly palettes
library(tibble)

#################################  Load the data 
##//load count tables for Arabidopsis
s78_Atha <- read.csv("data/78_Atha_processed__deduplicated.csv", row.names = 1, check.names = FALSE) # colnames are numbers
s78_Atha

s69_Atha <- read.csv("data/69_Atha_processed__deduplicated.csv", row.names = 1, check.names = FALSE) # colnames are numbers
s69_Atha

#load count tables for H. schachtii
s78_Hsch <- read.csv("data/78_Hsch_processed__deduplicated.csv", row.names = 1, check.names = FALSE) # colnames are numbers
s78_Hsch

s69_Hsch <- read.csv("data/69_Hsch_processed__deduplicated.csv", row.names = 1, check.names = FALSE) # colnames are numbers
s69_Hsch

################################## pre-process
#################### S69 Sample

s69_Atha_noERCC <- prepocessMTX(s69_Atha, geneIDcol = 'Gene_id', chromosomeCol = 'Chromosome')
s69_Hsch_noERCC <- prepocessMTX(s69_Hsch, geneIDcol = 'Gene_id', chromosomeCol = 'Chromosome')

##// SCT normalize data
sobj.s69 <- runSCT(s69_Atha_noERCC, s69_Hsch_noERCC, combineMTX = TRUE, min.cells = c(1), min.features = c(300), regressVars = FALSE)

#// List of oflength 1; unlist
sobj.s69 <- sobj.s69[[1]]

#Otherwise, i.e. select the highly variable genes for the NEM and A.thaliana combined matrix, run below code;
#// Select highly variable features
sojb.s69_At_var <- FindVariableFeatures(sobj.s69, nfeatures = 10000)
s69_At_var.features <- c(sojb.s69_At_var@assays$SCT@var.features)

#// Extract scaled and normalized matrix
s69.full.mtx <- as.matrix(sobj.s69@assays$SCT@scale.data)
s69.var.mtx <- as.matrix(sobj.s69@assays$SCT@scale.data)[s69_At_var.features,]
s69.hs.mtx <- as.matrix(sobj.s69@assays$SCT@scale.data)[grepl('^transcript:Hs', rownames(sobj.s69@assays$SCT@scale.data)),]
s69.at.mtx <- as.matrix(sojb.s69_At_var@assays$SCT@scale.data)[grepl('^AT', rownames(sojb.s69_At_var@assays$SCT@scale.data)),] #26165

s69_hvgs <- s69.at.mtx[s69_At_var.features[grep('^AT', s69_At_var.features)],]
s69_hvgs


#// Extract non-scaled matrices species matrices
s69.hs.mtx_unscaled <- as.matrix(sobj.s69@assays$SCT@counts)[grepl('^transcript:Hs', rownames(sobj.s69@assays$SCT@counts)),]
s69.at.mtx_unscaled <- as.matrix(sobj.s69@assays$SCT@counts)[grepl('^AT', rownames(sobj.s69@assays$SCT@counts)),]


####################### S78 sample
s78_Atha_noERCC <- prepocessMTX(s78_Atha, geneIDcol = 'Gene_id', chromosomeCol = 'Chromosome')
s78_Hsch_noERCC <- prepocessMTX(s78_Hsch, geneIDcol = 'Gene_id', chromosomeCol = 'Chromosome')

##// SCT normalize data
sobj.s78 <- runSCT(s78_Atha_noERCC, s78_Hsch_noERCC, combineMTX = TRUE, min.cells = c(1), min.features = c(300), regressVars = FALSE)

#// List of oflength 1; unlist
sobj.s78 <- sobj.s78[[1]]

#Otherwise, i.e. select the highly variable genes for the NEM and A.thaliana combined matrix, run below code;
#// Select highly variable features
sojb.s78_At_var <- FindVariableFeatures(sobj.s78, nfeatures = 10000)
s78_At_var.features <- c(sojb.s78_At_var@assays$SCT@var.features)

#// Extract scaled and normalized matrix
s78.full.mtx <- as.matrix(sobj.s78@assays$SCT@scale.data)
s78.var.mtx <- as.matrix(sobj.s78@assays$SCT@scale.data)[s78_At_var.features,]
s78.hs.mtx <- as.matrix(sobj.s78@assays$SCT@scale.data)[grepl('^transcript:Hs', rownames(sobj.s78@assays$SCT@scale.data)),]
s78.at.mtx <- as.matrix(sojb.s78_At_var@assays$SCT@scale.data)[grepl('^AT', rownames(sojb.s78_At_var@assays$SCT@scale.data)),] #26165

s78_hvgs <- s78.at.mtx[s78_At_var.features[grep('^AT', s78_At_var.features)],]
s78_hvgs

dim(s78.hs.mtx)

#// Extract non-scaled matrices species matrices
s78.hs.mtx_unscaled <- as.matrix(sobj.s78@assays$SCT@counts)[grepl('^transcript:Hs', rownames(sobj.s78@assays$SCT@counts)),]
s78.at.mtx_unscaled <- as.matrix(sobj.s78@assays$SCT@counts)[grepl('^AT', rownames(sobj.s78@assays$SCT@counts)),]


#################################  Fig2 barpolots
# Figure 2A
#getGPS(list('At' = s69_Atha_noERCC, 'Hs' = s69_Hsch_noERCC), sid = 'S69') 
getGPS(list('At' = s69_Atha_noERCC, 'Hs' = s69_Hsch_noERCC), plot_title = expression(bold("1 dpi")), sid = 'S69') 
# Figure 2C
#getGPS(list('At' = s69.at.mtx_unscaled, 'Hs' = s69.hs.mtx_unscaled), sid = 'S69') 
getGPS(list('At' = s69.at.mtx_unscaled, 'Hs' = s69.hs.mtx_unscaled), plot_title = expression(bold("1 dpi")), sid = 'S69') 
# Figure 2B
#getGPS(list('At' = s78_Atha_noERCC, 'Hs' = s78_Hsch_noERCC), sid = 'S78')
getGPS(list('At' = s78_Atha_noERCC, 'Hs' = s78_Hsch_noERCC), plot_title = expression(bold("2 dpi")), sid = 'S78')
# Figure 2D
#getGPS(list('At' = s78.at.mtx_unscaled, 'Hs' = s78.hs.mtx_unscaled), sid = 'S78')
getGPS(list('At' = s78.at.mtx_unscaled, 'Hs' = s78.hs.mtx_unscaled), plot_title = expression(bold("2 dpi")), sid = 'S78')

### To be deleted
# Figure 2A
Figure_2A <- getGPS(list('At' = s69_Atha_noERCC, 'Hs' = s69_Hsch_noERCC), plot_title = expression(bold("1 dpi")), sid = 'S69') 
ggsave(file.path(getwd(), 'plots/Fig2A_barplot.png'), dpi = 600, width = 6, height = 4)
# Figure 2C
Figure_2C <- getGPS(list('At' = s69.at.mtx_unscaled, 'Hs' = s69.hs.mtx_unscaled), plot_title = expression(bold("1 dpi")), sid = 'S69') 
ggsave(file.path(getwd(), 'plots/Fig2C_barplot.png'), dpi = 600, width = 6, height = 4)
# Figure 2B
Figure_2B <- getGPS(list('At' = s78_Atha_noERCC, 'Hs' = s78_Hsch_noERCC), plot_title = expression(bold("2 dpi")), sid = 'S78')
ggsave(file.path(getwd(), 'plots/Fig2B_barplot.png'), dpi = 600, width = 6, height = 4)
# Figure 2D
Figure_2D <- getGPS(list('At' = s78.at.mtx_unscaled, 'Hs' = s78.hs.mtx_unscaled), plot_title = expression(bold("2 dpi")), sid = 'S78')
ggsave(file.path(getwd(), 'plots/Fig2D_barplot.png'), dpi = 600, width = 6, height = 4)


## The number of detected Hs genes with log10 values using selected raw Hs matrix
s69_df_logHs <- s69.hs.mtx_unscaled %>%
  as.data.frame() %>%  # Convert matrix to data frame
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%  # Convert nonzero values to 1
  summarise(across(everything(), sum)) %>%  # Column sums
  pivot_longer(cols = everything(), names_to = "Section", values_to = "genecount") %>%  # Reshape to long format
  mutate(
    log10GeneCounts = log10(genecount),  # Compute log10 of gene counts
    Section = factor(Section, levels = unique(Section))  # Order factor levels
  ) %>%
  mutate(Section = as.numeric(as.character(Section))) %>%
  mutate(HighGeneCount = ifelse(genecount > 75, "Yes", "No")) ## Gene threshold: 75 

# Create ggplot: Figure 2E
ggplot(s69_df_logHs, aes(x = Section, y = genecount, fill = HighGeneCount)) +
  geom_bar(stat = "identity") +  # Bar plot with actual counts
  theme_bw() +
  #ggtitle(paste('Genes per section', "S69", sep = ' - ')) +
  xlab('Section') +
  ylab(expression(paste(italic('H. schachtii'), ' gene counts'))) +
  # scale_fill_brewer(palette = 'Greys') + ## what we used previously
  scale_fill_manual(values = c("#E1BE6A","#D35FB7")) +  # Custom colors
  labs(x = "Section", fill = "Meeting threshold") +
  ggtitle(expression(bold("1 dpi"))) +
  presentation +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 75, color = "red") 

ggsave(file.path(getwd(), 'plots/Fig2E_1dpi_S69_barplot.png'), dpi = 600, width = 6, height = 4)

##############  Figure 2 E/F
## The number of detected Hs genes with log10 values using selected raw Hs matrix
s78_df_logHs <- s78.hs.mtx_unscaled %>%
  as.data.frame() %>%  # Convert matrix to data frame
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%  # Convert nonzero values to 1
  summarise(across(everything(), sum)) %>%  # Column sums
  pivot_longer(cols = everything(), names_to = "Section", values_to = "genecount") %>%  # Reshape to long format
  mutate(
    log10GeneCounts = log10(genecount),  # Compute log10 of gene counts
    Section = factor(Section, levels = unique(Section))  # Order factor levels
  ) %>%
  mutate(Section = as.numeric(as.character(Section))) %>%
  mutate(HighGeneCount = ifelse(genecount > 75, "Yes", "No")) ## Gene threshold: 75 

# Create ggplot: Figure 2F
ggplot(s78_df_logHs, aes(x = Section, y = genecount, fill = HighGeneCount)) +
  geom_bar(stat = "identity") +  # Bar plot with actual counts
  theme_bw() +
  #ggtitle(paste('Genes per section', "S69", sep = ' - ')) +
  xlab('Section') +
  ylab(expression(paste(italic('H. schachtii'), ' gene counts'))) +
  # scale_fill_brewer(palette = 'Greys') + ## what we used previously
  scale_fill_manual(values = c("#E1BE6A","#D35FB7" )) +  # Custom colors
  #scale_fill_viridis_d(option = "D") +
  labs(x = "Section", fill = "Meeting threshold") +
  ggtitle(expression(bold("2 dpi"))) +
  presentation +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 75, color = "red") 

ggsave(file.path(getwd(), 'plots/Fig2F_2dpi_S78_barplot.png'), dpi = 600, width = 6, height = 4)

################################ Sup Fig 1 A-D

################  Sup .Fig1, number of detected Ara and Hs genes
##// plot gene counts in ascending order
s78_Atha_Hsch_comb_noERCC <- rbind(s78_Atha_noERCC, s78_Hsch_noERCC)

sup_fig1_b_s78 <- getGPS(list('At + Hs' = s78_Atha_Hsch_comb_noERCC), sid = 'S78', plotAscending = TRUE) 

## If we prefer to use darker grey for the bars
p <- ggplot(sup_fig1_b_s78[[2]][[1]], aes(x = Section, y = GPS, fill = dataset)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  #ggtitle(paste('Genes per section', 'S78', sep = ' - ')) +
  ggtitle(expression(bold("2 dpi"))) +
  xlab('Section') +
  ylab('Total Gene count') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("At + Hs" = "#0072B2")) + 
  geom_hline(yintercept = 300, color = "red") +
  theme(axis.text.x = element_text(size=5, angle = 90, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_text(size=8, color="black"),
        axis.title.y = element_text(size=8, color="black"),
        strip.text.x = element_text(size=8, color="black"),
        strip.text.y = element_text(size=8, color="black"),
        plot.title = element_text(size=12),
        panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "grey80"), #80
        panel.grid.minor = element_blank(),
        legend.position = "right")
p  # Sup Fig 1B

ggsave(file.path(getwd(), 'plots/Sup_fig1B_1dpi_barplot.png'), dpi = 600, width = 6, height = 4)

## The number of detected Hs genes with log10 values using selected raw Hs matrix
s78_df_logHs <- s78_Hsch_noERCC %>%
  as.data.frame() %>%  # Convert matrix to data frame
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%  # Convert nonzero values to 1
  summarise(across(everything(), sum)) %>%  # Column sums
  pivot_longer(cols = everything(), names_to = "Section", values_to = "genecount") %>%  # Reshape to long format
  #mutate(Section = forcats::fct_reorder(Section, genecount, .fun = min))  # reorder factor by genecount
  arrange(genecount)   %>%
  mutate(Section = factor(Section, levels = Section)) # %>%
print(s78_df_logHs, n = 100)

# s78_Hs_barplot <- ggplot(s78_df_logHs, aes(x = Section, y = genecount)) + geom_bar(stat = "identity", fill = "grey70") + theme_bw() +
#   ggtitle(paste('H. schachtii genes per section - S78')) + xlab('Section') + ylab('Total Gene count') + 
#   theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size = 6)
#   ) +presentation + geom_hline(yintercept = 75, color = "red")
# print(s78_Hs_barplot)   

s78_Hs_barplot <- ggplot(s78_df_logHs, aes(x = Section, y = genecount)) + geom_bar(stat = "identity", fill = "#0072B2") + theme_bw() +
  #ggtitle(paste('H. schachtii genes per section - S78')) + 
  ggtitle(expression(bold("2 dpi"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Section') + ylab('Total Gene count') + 
  geom_hline(yintercept = 75, color = "red") + 
  theme(axis.text.x = element_text(size=5, angle = 90, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_text(size=8, color="black"),
        axis.title.y = element_text(size=8, color="black"),
        strip.text.x = element_text(size=8, color="black"),
        strip.text.y = element_text(size=8, color="black"),
        plot.title = element_text(size=12),
        panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "grey80"), #80
        panel.grid.minor = element_blank(),
        legend.position = "right")
print(s78_Hs_barplot) ### Sup Fig 1D

ggsave(file.path(getwd(), 'plots/Sup_fig1D_2dpi_Hsch_barplot.png'), dpi = 600, width = 6, height = 4)

##################################  Sup fig 1 A and C
s69_Atha_Hsch_comb_noERCC <- rbind(s69_Atha_noERCC, s69_Hsch_noERCC)

sup_fig1_a_s69 <- getGPS(list('At + Hs' = s69_Atha_Hsch_comb_noERCC), sid = 'S69', plotAscending = TRUE) 

## If we prefer to use darker grey for the bars
p <- ggplot(sup_fig1_a_s69[[2]][[1]], aes(x = Section, y = GPS, fill = dataset)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  #ggtitle(paste('Genes per section', 'S69', sep = ' - ')) +
  ggtitle(expression(bold("1 dpi"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Section') +
  ylab('Total Gene count') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("At + Hs" = "#0072B2")) + 
  geom_hline(yintercept = 300, color = "red") +
  theme(axis.text.x = element_text(size=5, angle = 90, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_text(size=8, color="black"),
        axis.title.y = element_text(size=8, color="black"),
        strip.text.x = element_text(size=8, color="black"),
        strip.text.y = element_text(size=8, color="black"),
        plot.title = element_text(size=12),
        panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "grey80"), #80
        panel.grid.minor = element_blank(),
        legend.position = "right")
p # Sup Fig 1A 

ggsave(file.path(getwd(), 'plots/Sup_fig1A_1dpi_barplot.png'), dpi = 600, width = 6, height = 4)

s69_df_logHs <- s69_Hsch_noERCC %>%
  as.data.frame() %>%  # Convert matrix to data frame
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) %>%  # Convert nonzero values to 1
  summarise(across(everything(), sum)) %>%  # Column sums
  pivot_longer(cols = everything(), names_to = "Section", values_to = "genecount") %>%  # Reshape to long format
  #mutate(Section = forcats::fct_reorder(Section, genecount, .fun = min))  # reorder factor by genecount
  arrange(genecount)   %>%
  mutate(Section = factor(Section, levels = Section)) # %>%
print(s69_df_logHs, n = 100)

s69_Hs_barplot <- ggplot(s69_df_logHs, aes(x = Section, y = genecount)) + geom_bar(stat = "identity", fill = "#0072B2") + theme_bw() +
  #ggtitle(paste('H. schachtii genes per section - S69')) + xlab('Section') + ylab('Total Gene count') + 
  ggtitle(expression(bold("1 dpi"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 75, color = "red") + 
  theme(axis.text.x = element_text(size=5, angle = 90, color="black"),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_text(size=8, color="black"),
        axis.title.y = element_text(size=8, color="black"),
        strip.text.x = element_text(size=8, color="black"),
        strip.text.y = element_text(size=8, color="black"),
        plot.title = element_text(size=12),
        panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "grey80"), #80
        panel.grid.minor = element_blank(),
        legend.position = "right")
print(s69_Hs_barplot) ### Sup Fig 1C

ggsave(file.path(getwd(), 'plots/Sup_fig1C_1dpi_Hsch_barplot.png'), dpi = 600, width = 6, height = 4)
#ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/modified/Sup_Fig1C_Hsch_barplot_darkgrey.png'), dpi = 600, width = 6, height = 4)



###################################### Sup fig 1 A-D

###### Elbod method 
# S78 sample: s78.hs.mtx
dim(s78.hs.mtx) #10655, 92

s78.wcss_mt <- s78.hs.mtx

wcss <- function(k) {
  set.seed(42)
  cluster_var <- kmeans(s78.wcss_mt, centers = k, nstart = 100) #$tot.withinss
  return(list(
    tot_withiness = cluster_var$tot.withinss,
    clusters = cluster_var$cluster
  ))
}

# Set range for k (e.g., 1 to 10)
k_values <- 1:15

# Use lapply to get a list of results for different k values
s78_wcss_results <- lapply(k_values, wcss)

# Extract tot.withinss values for plotting
s78_wcss_values <- sapply(s78_wcss_results, function(x) x$tot_withiness)
# Extract cluster assignments for each k (optional)
s78_cluster_assignments <- lapply(s78_wcss_results, function(x) x$clusters)

### Section 2: The second derivative method

wcss_diff <- diff(s78_wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element

# plot- Elbow curve  
elbow_plot <- data.frame(k = k_values, WCSS = s78_wcss_values)
p_elbow_red <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") 
p_elbow_red    ### Sup Fig 2B

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(s78_wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;  
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(4, 8, 12), limits = c(min(acc_plot$k), max(acc_plot$k)))  # <-- this line sets custom x-axis ticks
p_acc       ### Sup Fig 2D

print(paste("Optimal number of clusters based on acceleration method:", optimal_k))   ## 9 for s78, NEM genes

dev.off()
# plot
# Arabidopsis heatmaps combined matrix_plots_clusters
        print(p_elbow_red)
        ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3B_9clusters_elbow.png'), dpi = 600, width = 6, height = 4)
        print(p_acc)
        ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3B_9clusters_elbow_second_der.png'), dpi = 600, width = 6, height = 4)
        
        ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/modified/Fig3B_9clusters_elbow_second_der.png'), dpi = 600, width = 6, height = 4)


######################### S69
s69.hs.mtx
dim(s69.hs.mtx) #

s69.wcss_mt <- s69.hs.mtx

wcss <- function(k) {
  set.seed(42)
  cluster_var <- kmeans(s69.wcss_mt, centers = k, nstart = 100) #$tot.withinss
  return(list(
    tot_withiness = cluster_var$tot.withinss,
    clusters = cluster_var$cluster
  ))
}

# Set range for k (e.g., 1 to 10)
k_values <- 1:15

# Use lapply to get a list of results for different k values
s69_wcss_results <- lapply(k_values, wcss)

# Extract tot.withinss values for plotting
s69_wcss_values <- sapply(s69_wcss_results, function(x) x$tot_withiness)
# Extract cluster assignments for each k (optional)
s69_cluster_assignments <- lapply(s69_wcss_results, function(x) x$clusters)

### Section 2: The second derivative method

wcss_diff <- diff(s69_wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element

# plot- Elbow curve
elbow_plot <- data.frame(k = k_values, WCSS = s69_wcss_values)
p_elbow_red <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_elbow_red ### Sup Fig 2A

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(s69_wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(4, 8, 12), limits = c(min(acc_plot$k), max(acc_plot$k)))  # <-- this line sets custom x-axis ticks
p_acc       ### Sup Fig 2C

print(paste("Optimal number of clusters based on acceleration method:", optimal_k))   

dev.off()
# plot
##################################################################################   Nematode heatmap
#################################  Fig3. NEM spatial heatmap
#Fig 3A

######### preparation 
actual_col_order <- c(1:36, 38:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:56)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 57:68)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 69:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 57 & actual_col_order <= 68, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#### cluster and genes
vals <- c(as.matrix(s69_cluster_assignments[[6]])) #only values
table(vals)

opt_clusters <- s69_cluster_assignments[[6]]
opt_clusters

############### S69 NEM spatial map,  Fig 3A
##

#png("plots/Fig3A_advanced_heatmap__improved_dendro.png", width = 12, height = 10, units = "in", res = 600)

png("plots/Fig3A_advanced_heatmap__improved_dendro___wid3_8-8.png", width = 8, height = 8, units = "in", res = 600)


ht_final <- ComplexHeatmap::Heatmap(s69.hs.mtx, cluster_columns = FALSE, border = TRUE, 
                                    row_split = vals,
                                    row_dend_width = unit(3,"cm"),
                                    clustering_method_rows = "average",
                                    show_row_names = FALSE, name = 'Pearson residuals', 
                                    column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                                    column_names_gp = gpar(fontsize = 7),
                                    column_split = factor(column_groups, levels = unique (column_groups)), top_annotation = col_ha
)

ht_final <- draw(ht_final, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE) 
ht_final

dev.off()

      
      ######################### Save the gene in 6 clusters for S69
      dim(s69.hs.mtx) #2609  95
      
      rownames(s69.hs.mtx) <- gsub("transcript:", "",rownames(s69.hs.mtx))
      
      # Extract row cluster assignments
      row_clusters <- row_order(ht_final)  # A list of row indices for each cluster
      
      # Extract the names of the genes in each cluster
      S69_gene_names_by_cluster <- lapply(row_clusters, function(idx) rownames(s69.hs.mtx)[idx])
      S69_gene_names_by_cluster <- S69_gene_names_by_cluster[order(as.numeric(names(S69_gene_names_by_cluster)))]
      table(vals)
      
      file_path <- "Arabidopsis heatmaps combined matrix_plots_clusters/Final/S69_Hs_2609genes_6clusters_spatial_map.xlsx"
      Save_clusters(S69_gene_names_by_cluster,file_path)
      
      dev.off()

################################# S78 NEM spatial sample， Figure 3B
s78.hs.mtx 

## preparation 
actual_col_order <- c(1:14, 16:36, 38:69, 71, 72, 74:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:19)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 20:39)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 40:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 20 & actual_col_order <= 39, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#### cluster and genes
vals <- c(as.matrix(s78_cluster_assignments[[9]])) #only values
table(vals)

opt_clusters <- s78_cluster_assignments[[9]]
opt_clusters


### Spatial heatmap, for S78. NEM part:  _10-6

png("plots/Fig3B_advanced_heatmap__improved_dendro___wid3_8-8.png", width = 8, height = 8, units = "in", res = 600)


ht_final <- ComplexHeatmap::Heatmap(s78.hs.mtx, cluster_columns = FALSE, border = TRUE, 
                                    row_split = vals,
                                    row_dend_width = unit(3,"cm"),
                                    clustering_method_rows = "average",
                                    show_row_names = FALSE, name = 'Pearson residuals', 
                                    column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                                    column_names_gp = gpar(fontsize = 7),
                                    column_split = factor(column_groups, levels = unique (column_groups)), top_annotation = col_ha)

ht_final <- draw(ht_final, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE) 
ht_final

dev.off()

      
      ######################### Save the gene in 6 clusters for S78
      rownames(s78.hs.mtx) <- gsub("transcript:", "",rownames(s78.hs.mtx))
      
      dim(s78.hs.mtx)
      # Extract row cluster assignments
      row_clusters <- row_order(ht_final)  # A list of row indices for each cluster
      
      # Extract the names of the genes in each cluster
      S78_gene_names_by_cluster <- lapply(row_clusters, function(idx) rownames(s78.hs.mtx)[idx])
      S78_gene_names_by_cluster <- S78_gene_names_by_cluster[order(as.numeric(names(S78_gene_names_by_cluster)))]
      table(vals)
      
      file_path <- "Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3B_9clusters_heatmap.xlsx"
      Save_clusters(S78_gene_names_by_cluster,file_path)
      
######################## Sug fig  2E-H

### 1DPI sample, s69

############ Construct the gland matrix
#// Read genes
effectorome <- read_excel("data/Effertorome_paper__Supplementary_Table_S1_latest.xlsx",skip=4)
#//filter out known genes
effectorome <- subset(effectorome, `Known effector or putative effector`== 'known effector' )

#// Get Hs genes present in Sample
detected.HsGenes <- rownames(s69.hs.mtx)  #s78.hs.mtx

#// Remove the "transcript:" part
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, ':')[[1]][2]})) 

#// Replace "-" with "_"
detected.HsGenes <- gsub('-', '_', detected.HsGenes)

#// Remove ".t1"
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, '\\.')[[1]][1]})) 
rownames(s69.hs.mtx) <- detected.HsGenes #s78.hs.mtx

#// Check which effector genes in our dataset
inDataset <- subset(effectorome, `Gene ID` %in% detected.HsGenes)

#// Subset to those with annotation
inDataset <- subset(inDataset, !is.na(`Expected gland (updated)`))

#// Write to file
#write.table(inDataset, file = file.path(getwd(), 'S78KnownGlandGenes_v2.tsv'), sep = '\t', quote = FALSE, row.names = TRUE)

#// Drop useless columns
inDataset <- inDataset[,c('Gene ID', 'Expected gland (updated)')]

#// Cleanup annotation
inDataset$`Expected gland (updated)` <- unlist(lapply(inDataset$`Expected gland (updated)`, function(x) strsplit(x, '-')[[1]][1]))


########### Apply the elbow method on detected gland genes's matrix
s69_gland_matrix <- s69.hs.mtx[inDataset$`Gene ID`,]
s69_gland_matrix

k_values <- 1:15

wcss_output <- elbow_analysis(s69_gland_matrix)
wcss_values <- wcss_output$wcss_values

# Plot the Elbow Curve
wcss_diff <- diff(wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element
optimal_k

# plot- Elbow curve  
elbow_plot <- data.frame(k = k_values, WCSS = wcss_values)
p_elbow <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_elbow  ### Sup Fig 2E

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;  
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(4, 8, 12), limits = c(min(acc_plot$k), max(acc_plot$k)))  # <-- this line sets custom x-axis ticks
p_acc       ### Sup Fig 2G

dev.off()
# Output the optimal K. 
print(paste("Optimal number of clusters based on acceleration method:", optimal_k))  #11

# # plot
# print(p_elbow)

####################### Figure 3C: S69 gland heatmap
opt_clusters <- wcss_output$cluster_assignments[[11]]
vals <- c(as.matrix(wcss_output$cluster_assignments[[11]])) #only values

actual_col_order <- c(1:36, 38:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:56)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 57:68)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 69:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 57 & actual_col_order <= 68, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#opt_clusters <- s69_cluster_assignments[[6]]

hb <- rowAnnotation(Gland = inDataset$`Expected gland (updated)`,
                    col = list('Gland' = c("DG " = "#0072B2", "SvG " = "#009E73", "amphid " = "#56B4E9")))

#// Create heatmap S69
#png(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3C_11clusters_heatmap.png'), width=6,height=6,units="in",res=600)

ht <- ComplexHeatmap::Heatmap(s69_gland_matrix, cluster_columns = FALSE, border = TRUE,
                              row_split = vals,
                              show_row_names = FALSE, name = 'Pearson residuals', 
                              column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                              left_annotation = hb,column_names_gp = gpar(fontsize = 7), 
                              column_split = factor(column_groups, levels = unique(column_groups)), top_annotation = col_ha)

ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht
dev.off()



#### sup fig 2 F and H
#### 2PDI, s78

############ Construct the gland matrix

#// Read genes
effectorome <- read_excel("data/Effertorome_paper__Supplementary_Table_S1_latest.xlsx",skip=4)
#//filter out known genes
effectorome <- subset(effectorome, `Known effector or putative effector`== 'known effector' )

#// Get Hs genes present in Sample
detected.HsGenes <- rownames(s78.hs.mtx)  #s78.hs.mtx

#// Remove the "transcript:" part
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, ':')[[1]][2]})) 

#// Replace "-" with "_"
detected.HsGenes <- gsub('-', '_', detected.HsGenes)

#// Remove ".t1"
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, '\\.')[[1]][1]})) 
rownames(s78.hs.mtx) <- detected.HsGenes #s78.hs.mtx

#// Check which effector genes in our dataset
inDataset <- subset(effectorome, `Gene ID` %in% detected.HsGenes)

#// Subset to those with annotation
inDataset <- subset(inDataset, !is.na(`Expected gland (updated)`))

#// Write to file
#write.table(inDataset, file = file.path(getwd(), 'S78KnownGlandGenes_v2.tsv'), sep = '\t', quote = FALSE, row.names = TRUE)

#// Drop useless columns
inDataset <- inDataset[,c('Gene ID', 'Expected gland (updated)')]

#// Cleanup annotation
inDataset$`Expected gland (updated)` <- unlist(lapply(inDataset$`Expected gland (updated)`, function(x) strsplit(x, '-')[[1]][1]))


################# Apply the elbow method on detected gland genes's matrix
s78_gland_matrix <- s78.hs.mtx[inDataset$`Gene ID`,]
s78_gland_matrix

k_values <- 1:15
wcss_output <- elbow_analysis(s78_gland_matrix)
wcss_values <- wcss_output$wcss_values

# Plot the Elbow Curve

wcss_diff <- diff(wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element
optimal_k

# plot- Elbow curve  
elbow_plot <- data.frame(k = k_values, WCSS = wcss_values)
p_elbow <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_elbow  ### Sup Fig 2F

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;  
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")+
  scale_x_continuous(breaks = c(4, 8, 12), limits = c(min(acc_plot$k), max(acc_plot$k)))  # <-- this line sets custom x-axis ticks
p_acc       ### Sup Fig 2H

dev.off()
      # plot
      print(p_elbow)
      ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3D_13clusters_elbow.png'), dpi = 600, width = 6, height = 4)
      
      print(p_acc)
      ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3D_13clusters_elbow_second_der.png'), dpi = 600, width = 6, height = 4)
      
      ggsave(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/modified/Fig3H_13clusters_elbow_second_der.png'), dpi = 600, width = 6, height = 4)


# Output the optimal K. 
print(paste("Optimal number of clusters based on acceleration method:", optimal_k))  #13




################ Figure 3C and 3D, Gland spatial heatmap
### Figure 3D, the order of clusters is inconsistant. 


####### ###### Figure 3D,gland spatial heatmap, s78

opt_clusters <- wcss_output$cluster_assignments[[13]]
vals <- c(as.matrix(wcss_output$cluster_assignments[[13]])) #only values

## Prepare for the legend bar
actual_col_order <- c(1:14, 16:36, 38:69, 71, 72, 74:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:19)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 20:39)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 40:96)))

test_column_groups <- ifelse(actual_col_order >= 20 & actual_col_order <= 39, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

hb <- rowAnnotation(Gland = inDataset$`Expected gland (updated)`,
                    col = list('Gland' = c("DG " = "#0072B2", "SvG " = "#009E73", "amphid" = "#56B4E9")))

#// Create heatmap S78
#png(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3D_13clusters_heatmap.png'), width=6,height=6,units="in",res=600)


ht <- ComplexHeatmap::Heatmap(s78_gland_matrix, cluster_columns = FALSE, border = TRUE,
                              row_split = vals,
                              show_row_names = FALSE, name = 'Pearson residuals', 
                              column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                              left_annotation = hb,column_names_gp = gpar(fontsize = 7), 
                              column_split = factor(column_groups, levels = unique(column_groups)), top_annotation = col_ha)

ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht

dev.off()



        ################# Save the detected gland genes and cluster info into a excel file
        row_clusters <- row_order(ht)  # A list of row indices for each cluster
        
        # Extract the names of the genes in each cluster
        S78_gene_names_by_cluster <- lapply(row_clusters, function(idx) rownames(s78_gland_matrix)[idx])
        S78_gene_names_by_cluster <- S78_gene_names_by_cluster[order(as.numeric(names(S78_gene_names_by_cluster)))]
        table(opt_clusters)
        
        file_path <- "Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3D_13clusters.xlsx"
        Save_clusters(S78_gene_names_by_cluster,file_path)
        
        dev.off()
        
################# Save the detected gland genes and cluster info into a excel file


getwd()
##// Find similarly expressed genes
effectorome <- read_excel("data/Effertorome_paper__Supplementary_Table_S1_latest.xlsx",skip=4)
#//filter out known genes

## 1 DPI sample, s69

#// Get Hs genes present in Sample
detected.HsGenes <- rownames(s69.hs.mtx)

#// Replace "-" with "_"
detected.HsGenes <- gsub('-', '_', detected.HsGenes)

#// Remove ".t1"
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, '\\.')[[1]][1]})) 
detected.HsGenes <- gsub("transcript:", "",detected.HsGenes)
rownames(s69.hs.mtx) <- detected.HsGenes


rm(simGenes)
simGenes <- findSimilar('Hsc_gene_21727', s69.hs.mtx, method = 'spearman', n = 200, forceN = TRUE, threshold = 0.90)
rownames(simGenes) <- simGenes$Gene

#// Add reference gene back in
simGenes['Hsc_gene_21727',] <- c(NA, 'Hsc_gene_21727')

#// Annotate simGenes if present
simGenes$Gland <- NA
simGene.annot <- subset(effectorome, `Gene ID` %in% simGenes$Gene)

# If the overlapped genes is known (DG or SvG), replace its simGene$Gland  by the "Expected gland (updated)"
for (anno in simGene.annot$`Gene ID`) { # anno, the effecter name
  anno_row <- simGene.annot %>%
    filter(`Gene ID` == anno) #%>%
  #pull(`Expected gland (updated)`)
  
  if (!is.na(anno_row$`Expected gland (updated)`[[1]])) { # if known from which gland
    # replace it by the "Expected gland (updated)" value --> #### Annotate the putative genes.
    simGenes[anno,]$Gland <- anno_row$`Expected gland (updated)`[[1]]
    
  } else { #if not but exist 
    simGenes[anno,]$Gland <- anno_row$`Known effector or putative effector`[[1]]
  }
}


simGenes$Gland <- sapply(strsplit(simGenes$Gland, " "), `[`, 1)
simGenes$Gland <- gsub("putative", "Putative", simGenes$Gland)

#// Reformat annotation

#// Select gene for highlighting
which.gene <- which(grepl("Hsc_gene_21727", rownames(s69.hs.mtx[simGenes$Gene,])))
split <-  data.frame(x = c(rep("A", which.gene - 1), "B",
                           rep("C", nrow(s69.hs.mtx[simGenes$Gene,]) - which.gene)))

ha <- rowAnnotation(Gland = simGenes$Gland,
                    col = list('Gland' = c("DG" = "#0072B2", "SvG" = "#009E73", "amphid" = "#56B4E9","Putative" = "#F0E442")))

actual_col_order <- c(1:36, 38:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:56)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 57:68)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 69:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 57 & actual_col_order <= 68, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#// Create heatmap S69
#png(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/S69_Hs_gene9937__SIM200genes_heatmap.png'), width=6,height=6,units="in",res=600)
set.seed(42)
#png(file.path(getwd(), 'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3F_Hs_SvG_g21727_SIM200_heatmap.png'), width=6,height=6,units="in",res=600)

ht <- ComplexHeatmap::Heatmap(s69.hs.mtx[simGenes$Gene,], cluster_columns = FALSE, border = TRUE, 
                              clustering_distance_rows = 'euclidean', 
                              clustering_method_rows = 'average', show_row_names = FALSE, 
                              left_annotation = ha, row_split = split, column_names_gp = gpar(fontsize = 7),
                              name = 'Pearson residuals', column_split = factor(column_groups, levels = unique(column_groups)),
                              column_title = "Section", column_title_side = "bottom", top_annotation = col_ha, col = viridis(100) )
ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht  ### Fig 3F
dev.off()


#// write similarly expressed genes to file
write.table(simGenes, file = file.path(getwd(), 
                                       'Arabidopsis heatmaps combined matrix_plots_clusters/Final/paper/Fig3F_Hs_SvG_g2749_SIM200_heatmap.tsv'), sep = '\t', quote = FALSE, row.names = TRUE)



######################## 2DPI sample, Figure 3G

## 2 DPI sample, S78
##// Find similarly expressed genes
effectorome <- read_excel("data/Effertorome_paper__Supplementary_Table_S1_latest.xlsx",skip=4)

#// Get Hs genes present in Sample
detected.HsGenes <- rownames(s78.hs.mtx)

#// Replace "-" with "_"
detected.HsGenes <- gsub('-', '_', detected.HsGenes)

#// Remove ".t1"
detected.HsGenes <- unlist(lapply(detected.HsGenes, function(x){strsplit(x, '\\.')[[1]][1]})) 
detected.HsGenes <- gsub("transcript:", "",detected.HsGenes)
rownames(s78.hs.mtx) <- detected.HsGenes



s78.hs.mtx['Hsc_gene_21727',]
rm(simGenes)
simGenes <- findSimilar('Hsc_gene_9937', s78.hs.mtx, method = 'spearman', n = 200, forceN = TRUE, threshold = 0.90) #Hsc_gene_9937, first time
rownames(simGenes) <- simGenes$Gene

#// Add reference gene back in
simGenes['Hsc_gene_9937',] <- c(NA, 'Hsc_gene_9937')

#// Annotate simGenes if present
simGenes$Gland <- NA
simGene.annot <- subset(effectorome, `Gene ID` %in% simGenes$Gene)

# If the overlapped genes is known (DG or SvG), replace its simGene$Gland  by the "Expected gland (updated)"
for (anno in simGene.annot$`Gene ID`) { # anno, the effecter name
  anno_row <- simGene.annot %>%
    filter(`Gene ID` == anno) #%>%
  #pull(`Expected gland (updated)`)
  
  if (!is.na(anno_row$`Expected gland (updated)`[[1]])) { # if known from which gland
    # replace it by the "Expected gland (updated)" value --> #### Annotate the putative genes.
    simGenes[anno,]$Gland <- anno_row$`Expected gland (updated)`[[1]]
    
  } else { #if not but exist 
    simGenes[anno,]$Gland <- anno_row$`Known effector or putative effector`[[1]]
  }
}


simGenes$Gland <- sapply(strsplit(simGenes$Gland, " "), `[`, 1)
simGenes$Gland <- gsub("putative", "Putative", simGenes$Gland)

#// Select gene for highlighting
which.gene <- which(grepl("Hsc_gene_9937", rownames(s78.hs.mtx[simGenes$Gene,])))
split <-  data.frame(x = c(rep("A", which.gene - 1), "B",
                           rep("C", nrow(s78.hs.mtx[simGenes$Gene,]) - which.gene)))

ha <- rowAnnotation(Gland = simGenes$Gland,
                    col = list("Gland" = c("DG" = "#0072B2", "SvG" = "#009E73", "amphid" = "#56B4E9","Putative" = "#F0E442"))
)

#mat <- matrix(rnorm(92*10), nrow=10, ncol = 92)
actual_col_order <- c(1:14, 16:36, 38:69, 71, 72, 74:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:19)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 20:39)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 40:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 20 & actual_col_order <= 39, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#// Create heatmap S78
set.seed(42)

ht <- ComplexHeatmap::Heatmap(s78.hs.mtx[simGenes$Gene,], cluster_columns = FALSE, border = TRUE, 
                              clustering_distance_rows = 'euclidean', 
                              clustering_method_rows = 'average', show_row_names = FALSE, 
                              left_annotation = ha, row_split = split, column_names_gp = gpar(fontsize = 7),
                              name = 'Pearson residuals', column_split = factor(column_groups, levels = unique(column_groups)),
                              column_title = "Section", column_title_side = "bottom", top_annotation = col_ha, col = viridis(100) )
ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht    ### Fig 3G
dev.off()


##################################################################################   Plant heatmap
######################################### Elbow method, S69， sup fig 3A and C

s69.wcss_mt <- s69.at.mtx[s69_At_var.features[grep('^AT', s69_At_var.features)],]

wcss <- function(k) {
  set.seed(42)
  cluster_var <- kmeans(s69.wcss_mt, centers = k, nstart = 100) #$tot.withinss
  return(list(
    tot_withiness = cluster_var$tot.withinss,
    clusters = cluster_var$cluster
  ))
}

# Set range for k (e.g., 1 to 10)
k_values <- 1:15

# Use lapply to get a list of results for different k values
s69_wcss_results <- lapply(k_values, wcss)

# Extract tot.withinss values for plotting
s69_wcss_values <- sapply(s69_wcss_results, function(x) x$tot_withiness)
# Extract cluster assignments for each k (optional)
s69_cluster_assignments <- lapply(s69_wcss_results, function(x) x$clusters)

# Plot the Elbow Curve
### Section 2: The second derivative method

wcss_diff <- diff(s69_wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element

# plot- Elbow curve  
elbow_plot <- data.frame(k = k_values, WCSS = s69_wcss_values)
p_elbow_red <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_elbow_red  ### Sup Fig 3A

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(s69_wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;  
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(4, 8, 12), limits = c(min(acc_plot$k), max(acc_plot$k)))  # <-- this line sets custom x-axis ticks
p_acc      ### Sup Fig 3C

print(paste("Optimal number of clusters based on acceleration method:", optimal_k))   ## also 14, same as S78

dev.off()
# plot

######################################### Elbow method
s78.wcss_mt <- s78.at.mtx[s78_At_var.features[grep('^AT', s78_At_var.features)],]

wcss <- function(k) {
  set.seed(42)
  cluster_var <- kmeans(s78.wcss_mt, centers = k, nstart = 100) #$tot.withinss
  return(list(
    tot_withiness = cluster_var$tot.withinss,
    clusters = cluster_var$cluster
  ))
}

# Set range for k (e.g., 1 to 10)
k_values <- 1:15

# Use lapply to get a list of results for different k values
s78_wcss_results <- lapply(k_values, wcss)

# Extract tot.withinss values for plotting
s78_wcss_values <- sapply(s78_wcss_results, function(x) x$tot_withiness)
# Extract cluster assignments for each k (optional)
s78_cluster_assignments <- lapply(s78_wcss_results, function(x) x$clusters)

# Plot the Elbow Curve

### Section 2: The second derivative method

wcss_diff <- diff(s78_wcss_values)
# Acceleration of change
wcss_acc <- diff(wcss_diff)
optimal_k <- which.min(wcss_acc) + 1  # plus 1 because diff() adds one element

# plot- Elbow curve  
elbow_plot <- data.frame(k = k_values, WCSS = s78_wcss_values)
p_elbow_red <- ggplot(elbow_plot, aes(x = k, y = WCSS)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Elbow Method for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_elbow_red ### Sup Fig 3B

# plot  - Second Derivative curve
acc_plot <- data.frame(k = 2:(length(s78_wcss_values)-1), Acceleration = wcss_acc)  # 2：19， 1 --> 1-2; 2-3;  
p_acc <- ggplot(acc_plot, aes(x = k, y = Acceleration)) +
  geom_point(size = 3) +
  geom_line() +
  theme_bw() +
  labs(title = "Second Derivative of WCSS (Acceleration Method)",
       x = "Number of Clusters (k)",
       y = "Acceleration (Second Derivative)") +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red")
p_acc    ### Sup Fig 3D

print(paste("Optimal number of clusters based on acceleration method:", optimal_k))    ## 14




########################################## Figure 4, spatial plant heatmap 

############################################ spatial maps
#################### basic heatmap  4A

######### preparation 
actual_col_order <- c(1:36, 38:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:56)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 57:68)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 69:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 57 & actual_col_order <= 68, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#### cluster and genes
vals <- c(as.matrix(s69_cluster_assignments[[14]])) #only values
table(vals)

opt_clusters <- s69_cluster_assignments[[14]]
opt_clusters



## REAL code
panel_personal = function(index, nm) {
  # nm: the cluster name
  each_cluster <- tibble::as_tibble(opt_clusters, rownames = "Gene") %>%
    dplyr::filter(value == nm) %>%
    dplyr::pull(Gene)
  
  num_genes <- length(each_cluster)
  
  cluster_plot <- tibble::as_tibble(s69_hvgs) %>%    ### Don't forget to change the s69_hvgs if you use other samples
    mutate(Gene = rownames(s69_hvgs), .before = 1) %>%
    dplyr::filter(Gene %in% each_cluster) %>%
    gather(Section, value, -Gene) %>%
    mutate(Section = as.numeric(as.character(Section)))
  
  # Calculate mean expression explicitly
  mean_expr <- cluster_plot %>%
    group_by(Section) %>%
    summarise(mean_val = mean(value))
  
  # ggplot adjusted to minimize distortion and enhance readability
  p <- ggplot(mean_expr, aes(x = Section, y = mean_val)) +
    geom_line(color = "red", linewidth = 0.7) +
    labs(title = paste0("Gene Set:  ", num_genes)) +
    theme_void() +  # remove axes, grids, etc.
    theme(
      plot.title = element_text(size = 9, hjust = 0.5) ,
      plot.margin = margin(1, 1, 1, 1, "pt"),  # small margin 
      panel.border = element_rect(color = "black",fill = NA, linewidth = 0.5) 
      
    ) +
    coord_cartesian(expand = FALSE)
  
  # Convert to grob and draw
  g <- ggplotGrob(p)
  grid.draw(g)
}



## text anntoations
S69_annotations <- read_excel("data/GO_Enrichment_ShinyGO 0.82_S69.xlsx",sheet = "TOP3_Figure_S69")
text_final <- S69_annotations %>% 
  filter(!if_all(everything(), ~is.na(.) | . == ""))  %>% # remove the empty rows
  mutate(`Pathways (click for details)` = sub("GO:\\d+ ", "", `Pathways (click for details)`))   %>% #remove the GO:: terms
  group_by(Clusters) %>%
  summarise(`Pathways (click for details)` = list(`Pathways (click for details)`)) %>%
  deframe()

## left annotation color
color_palette <- viridis(14, option = "D")  # "D" = colorblind-friendly

# anno_link with fixed size to avoid vertical stretching/distortion
anno = anno_link(
  align_to = vals,
  which = 'row',
  #side = ifelse( "left"),
  panel_fun = panel_personal,
  size = unit(1.5, "cm"),   # absolute height (adjustable)
  gap = unit(4, "cm"),
  width = unit(3, "cm") ,  # width ensures readability
  link_gp = gpar(fill = color_palette, col=NA) # fill = "grey80"
)

png("plots/Fig4A_advanced_heatmap__9732hvgs_updated_color_last.png", width = 12, height = 10, units = "in", res = 600)

## 
ht_final <- ComplexHeatmap::Heatmap(s69_hvgs, cluster_columns = FALSE, border = TRUE, 
                                    show_row_dend = FALSE,
                                    row_split = vals,
                                    right_annotation = rowAnnotation(textbox = anno_textbox(
                                      #vals, text_final[c("1","4","6", "7")]     ### display sepcific annotations
                                      vals, text_final,                       ### show all annotations
                                      add_new_line = TRUE,
                                      gp = gpar(col =  c("black") )    #     "#E69F00", "black","#009E73"
                                    )),
                                    left_annotation = rowAnnotation(foo = anno),
                                    
                                    show_row_names = FALSE, name = 'Pearson residuals', 
                                    column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                                    column_names_gp = gpar(fontsize = 7),
                                    column_split = factor(column_groups, levels = unique (column_groups)), top_annotation = col_ha
                                    
)

ht_final <- draw(ht_final, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht_final

dev.off()


#################### basic heatmap  4B

## preparation : deleted sections are 15,37,70,73  
actual_col_order <- c(1:14, 16:36, 38:69, 71, 72, 74:96)
column_groups <- c(rep("No nematode mRNA", sum(actual_col_order %in% 1:19)), 
                   rep("Nematode mRNA", sum(actual_col_order %in% 20:39)), 
                   rep("No nematode mRNA ", sum(actual_col_order %in% 40:96)))

## Prepare for the legend bar
test_column_groups <- ifelse(actual_col_order >= 20 & actual_col_order <= 39, "Nematode-containing", "Others")
test_group_colors <- c("Others" = "#CC79A7", "Nematode-containing" = "#E69F00", "Others" = "#CC79A7")
col_ha <- HeatmapAnnotation(Group = test_column_groups, col = list(Group = test_group_colors))

#### cluster and genes
vals <- c(as.matrix(s78_cluster_assignments[[14]])) #only values
table(vals)

opt_clusters <- s78_cluster_assignments[[14]]
opt_clusters

## REAL code
panel_personal = function(index, nm) {
  # nm: the cluster name
  each_cluster <- tibble::as_tibble(opt_clusters, rownames = "Gene") %>%
    dplyr::filter(value == nm) %>%
    dplyr::pull(Gene)
  
  num_genes <- length(each_cluster)
  
  cluster_plot <- tibble::as_tibble(s78_hvgs) %>%
    mutate(Gene = rownames(s78_hvgs), .before = 1) %>%
    dplyr::filter(Gene %in% each_cluster) %>%
    gather(Section, value, -Gene) %>%
    mutate(Section = as.numeric(as.character(Section)))
  
  # Calculate mean expression explicitly
  mean_expr <- cluster_plot %>%
    group_by(Section) %>%
    summarise(mean_val = mean(value))
  
  # ggplot adjusted to minimize distortion and enhance readability
  p <- ggplot(mean_expr, aes(x = Section, y = mean_val)) +
    geom_line(color = "red", linewidth = 0.7) +
    labs(title = paste0("Gene Set:  ", num_genes)) +
    theme_void() +  # remove axes, grids, etc.
    theme(
      plot.title = element_text(size = 9, hjust = 0.5) ,
      plot.margin = margin(1, 1, 1, 1, "pt"),  # small margin 
      panel.border = element_rect(color = "black",fill = NA, linewidth = 0.5) 
      
    ) +
    coord_cartesian(expand = FALSE)
  
  # Convert to grob and draw
  g <- ggplotGrob(p)
  grid.draw(g)
}



# Generate distinct colors from a colorblind-friendly palette
color_palette <- viridis(14, option = "D")  # "D" = colorblind-friendly

# anno_link with fixed size to avoid vertical stretching/distortion
anno = anno_link(
  align_to = vals,
  which = 'row',
  #side = ifelse( "left"),
  panel_fun = panel_personal,
  size = unit(1.5, "cm"),   # absolute height (adjustable)
  gap = unit(4, "cm"),
  width = unit(3, "cm") ,  # width ensures readability
  link_gp = gpar(fill = color_palette , col=NA) # fill = "grey80"
)


getwd()



S78_annotations <- read_excel("data/GO_Enrichment_ShinyGO 0.82_S78.xlsx",sheet = "TOP3_Figure")
text_final <- S78_annotations %>% 
  filter(!if_all(everything(), ~is.na(.) | . == ""))  %>% # remove the empty rows
  mutate(`Pathways (click for details)` = sub("GO:\\d+ ", "", `Pathways (click for details)`))   %>% #remove the GO:: terms
  group_by(Clusters) %>%
  summarise(`Pathways (click for details)` = list(`Pathways (click for details)`)) %>%
  deframe()

png("plots/Fig4B_advanced_heatmap__9053hvgs_updated_color_last.png", width = 12, height = 10, units = "in", res = 600)

## 
ht_final <- ComplexHeatmap::Heatmap(s78_hvgs, cluster_columns = FALSE, border = TRUE, 
                                    show_row_dend = FALSE,
                                    row_split = vals,
                                    right_annotation = rowAnnotation(textbox = anno_textbox(
                                      #vals, text_final[c("1","4","6", "7")]     ### display sepcific annotations
                                      vals, text_final,                         ### show all annotations
                                      add_new_line = TRUE,
                                      gp = gpar(col = c("black"))    #     "#E69F00", "black","#009E73"
                                    )),
                                    #left_annotation = rowAnnotation(Cluster = vals, Expression = left_exp_plots),
                                    
                                    left_annotation = rowAnnotation(foo = anno),
                                    
                                    show_row_names = FALSE, name = 'Pearson residuals', 
                                    column_title = "Section",column_title_side = "bottom",  col = viridis(100),
                                    column_names_gp = gpar(fontsize = 7),
                                    column_split = factor(column_groups, levels = unique (column_groups)), top_annotation = col_ha
                                    
)

ht_final <- draw(ht_final, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
ht_final

dev.off()

#sessionInfo()
getwd()
rmarkdown::render("Plots.Rmd", output_format = "html_document")
