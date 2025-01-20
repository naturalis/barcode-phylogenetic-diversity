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


## 2 Load metadata
meta_dir <- snakemake@input[[1]]
meta_raw <- read.delim(meta_dir, header = TRUE)
meta <- meta_raw[-c(1),]


## 3 Alpha diversity
## 3.1 Get data and prepare metadata
# Get alpha diversity faith pd
alpha_dir <- snakemake@input[[2]]
ad <- read.delim(alpha_dir)

# Adjust metadata to alpha diversity dataset
# Not all metadata is needed because of truncation and it needs ordering
alpha_ind <- which(meta$sample.id %in% ad$X.SampleID)
meta_afilt <- meta[alpha_ind,]

orda <- match(ad$X.SampleID, meta_afilt$sample.id)
meta_orda <- meta_afilt[orda,]

# 3.2 Visualize faith pd in ggplot heatmap
# Get vertical labeling for x-axis
# You need to specify the size of everything better
ggplot(ad, aes(x = X.SampleID, y = meta_orda$body.site, fill = faith_pd)) +
  geom_tile(color = 'black') + 
  scale_fill_gradient(low = "palegreen", high = "palegreen4") +
  coord_fixed() + 
  labs(y = "Body site", x = "Sample ID", fill = "Faith PD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

heatmap_dir <- dirname(snakemake@output[[1]])
ggsave('alpha_heatmap.png', plot = last_plot(), path = heatmap_dir)

## 4 Beta diversity
# 4.1 Get data and prepare metadata
# Get beta diversity distance matrix from somewhere and cluster with hclust
beta_dir <- snakemake@input[[3]]
beta_df <- read.delim(beta_dir)
dm_df <- beta_df[,-c(1)]
dm <- as.dist(dm_df)
hc <- hclust(dm, method = "average")
dend <- as.dendrogram(hc)

# Extract data from tree with ggdendro
dend_data <- dendro_data(dend, type = "rectangle")

# Adjust metadata to beta diversity dataset
# Not all metadata is needed because of truncation and it needs ordering
# Extraction of (meta)data for customization (ggdendro)
beta_ind <- which(meta$sample.id %in% dend_data$labels$label)
meta_filtb <- meta[beta_ind,]

ordb <- match(dend_data$labels$label, meta_filtb$sample.id)
meta_ordb <- meta_filtb[ordb,]

## 4.2 Visualize distance matrix in dendrogram (ggplot2)
# You can make default plot with ggdendro, but mainly useful for extraction of data
#ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)
# Integrate with ggplot2, add coloration, mind the fixed ylim
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label, col = meta_ordb$body.site),
            hjust = 1, angle = 90, size = 3) +
  labs(col = "Body site") +
  ylim(-0.5, 1)

dendro_dir <- dirname(snakemake@output[[2]])
ggsave('beta_dendro.png', plot = last_plot(), path = dendro_dir)

## 4.3 PCoA and visualization
# Perform pcoa
pcoa <- cmdscale(dm)
colnames(pcoa) <- c("pcoa1", "pcoa2")

# Visualize pcoa with ggplot2
# You can also insert percentages with glue package
# Make ellipses for interquartile range?
pcoaplot <- pcoa %>% 
  as_tibble(rownames = "samples") %>%
  ggplot(aes(x = pcoa1, y = pcoa2, color = meta_ordb$body.site)) + 
  geom_point() + labs(col = "Body site", x = "PCo 1", y = "PCo 2")
print(pcoaplot)
  
pcoa_dir <- dirname(snakemake@output[[3]])
ggsave('beta_pcoa.png', plot = last_plot(), path = pcoa_dir)

#-----------------------------------------------------------------------------------------------------------------------------------

