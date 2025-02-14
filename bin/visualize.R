#### This is a script for visualizing exported alpha and beta diversity from qiime2.
# Created 25-03-2024
# Needs integration into workflow, implementation of ITS dataset and metadata,
# and better/adaptable dimensions for visualization.

### Content
## 1 Dependencies
## 2 Load metadata
## 3 Alpha diversity
      # 3.1 Get data and prepare metadata
      # 3.2 Visualize faith pd in heatmap (ggplot2)
## 4 Beta diversity
      # 4.1 Get data and prepare metadata (ggdendro)
      # 4.2 Visualize distance matrix in dendrogram (ggplot2)
      # 4.3 Visualize pcoa (ggplot2)


#-------------------------------------------------------------------------------------------------------------
  
### Diversity metrics visualizations
## 1 Dependencies & working directory
#install.packages("ggplot2", "ggdendro", "dplyr")
library(ggplot2)
library(ggdendro)
library(dplyr)

## 3 Alpha diversity
## 3.1 Get data and prepare metadata
# Get alpha diversity faith pd
res_dir <- dirname(snakemake@input[[1]])
setwd(res_dir)

alpha_tsv <- basename(snakemake@input[[1]])
alpha_div <- read.delim(alpha_tsv)


# 3.2 Visualize faith pd in ggplot heatmap
# You need to specify the size of everything better
ggplot(alpha_div, aes(x = X.SampleID, y=1, fill = faith_pd)) +
  geom_tile(color = 'black') + 
  scale_fill_gradient(low = "palegreen", high = "palegreen4") +
  coord_fixed() + 
  labs(x = "Sample ID", y = '', fill = "Faith PD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

#heatmap_dir <- dirname(snakemake@output[[1]])
heatmap <- basename(snakemake@output[[1]])
ggsave(heatmap, plot = last_plot(), path = res_dir)

## 4 Beta diversity
# 4.1 Get data and prepare metadata
# Get beta diversity distance matrix from somewhere and cluster with hclust
beta_tsv <- basename(snakemake@input[[2]])
beta_raw <- read.delim(beta_tsv)
beta_div <- beta_raw[,-c(1)]
dist_matrix <- as.dist(beta_div)
hc <- hclust(dist_matrix, method = "average")
dend <- as.dendrogram(hc)

# Extract data from tree with ggdendro
dend_data <- dendro_data(dend, type = "rectangle")


## 4.2 Visualize distance matrix in dendrogram (ggplot2)
# You can make default plot with ggdendro, but mainly useful for extraction of data
#ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)
# Integrate with ggplot2, add coloration, mind the fixed ylim
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3) +
  ylim(-0.5, 1)

dendrogram <- snakemake@output[[2]]
ggsave('beta_dendro.png', plot = last_plot(), path = getwd())

## 4.3 PCoA and visualization
# Perform pcoa
pcoa <- cmdscale(dist_matrix)
colnames(pcoa) <- c("pcoa1", "pcoa2")

# Visualize pcoa with ggplot2
# You can also insert percentages with glue package
# Make ellipses for interquartile range?
pcoaplot <- pcoa %>% 
  as_tibble(rownames = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2)) + 
  geom_point() + 
  labs(x = "PCo 1", y = "PCo 2") +
  geom_text(hjust = 0, vjust = 0, label = rownames(pcoa))
print(pcoaplot)
  
#pcoa_dir <- dirname(snakemake@output[[3]])
pcoa <- basename(snakemake@output[[3]])
ggsave(pcoa, plot = last_plot(), path = res_dir)

#-----------------------------------------------------------------------------------------------------------------------------------

