## This is a script for visualizing exported alpha and beta diversity from qiime2.
# Created 25-03-2024

### Content
## 1 Dependencies and arguments from snakefile
## 2 Alpha diversity heatmap
## 3.1 Beta diversity preparation
## 3.2 Beta diversity dendrogram
## 3.3 Beta diversity pcoa

#-------------------------------------------------------------------------------------------------------------
  
### Diversity metrics visualizations
## 1 Dependencies and arguments from snakefile
library(ggplot2)
library(ggrepel)
library(ggdendro)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)


## 2 Alpha diversity heatmap
# Load alpha diversity data
alpha_tsv <- args[1]
alpha_div <- read.delim(alpha_tsv)

# Visualize faith pd in ggplot heatmap
ggplot(alpha_div, aes(x = X.SampleID, y=1, fill = faith_pd)) +
  geom_tile(color = 'black') + 
  scale_fill_gradient(low = "palegreen", high = "palegreen4") +
  coord_fixed() + 
  labs(x = "Sample ID", y = '', fill = "Faith PD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Save plot
heatmap <- args[3]
ggsave(filename = heatmap, plot = last_plot(), path = getwd())


## 3.1 Beta diversity preparation
# Load beta diversity
beta_tsv <- args[2]
beta_div <- read.table(beta_tsv, header = TRUE, row.names = 1, sep = "\t")

# Create distance matrix
dm <- as.matrix(beta_div)
dist_matrix <- as.dist(dm)


## 3.2 Beta diversity dendrogram
# Cluster with hclust and create dendrogram
hc <- hclust(dist_matrix, method = "average")
dend_data <- as.dendrogram(hc)

# Create dendrogram with rotated axis and extra padding around the plot
dendroplot <- ggdendrogram(dend_data, rotate = TRUE, size = 0.5) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  labs(x = "", y = "")


# Save dendrogram
dendrogram <- args[4]
ggsave(filename = dendrogram, plot = dendroplot, width = 8, height = 6, dpi = 300, bg = "white")

## 3.3 Beta diversity PCOA
# Perform pcoa
pcoa_data <- cmdscale(dist_matrix)
colnames(pcoa_data) <- c("PC1", "PC2")

# Visualize pcoa with ggplot2 and ggrepel against overlap
pcoaplot <- pcoa_data %>% 
  as_tibble(rownames = "samples") %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point() + 
  geom_text_repel(aes(x = PC1, y = PC2), label = rownames(pcoa_data))
  
print(pcoaplot)
  
# Save pcoa plot
pcoa <- args[5]
ggsave(filename = pcoa, plot = last_plot(), path = getwd())

#-----------------------------------------------------------------------------------------------------------------------------------

