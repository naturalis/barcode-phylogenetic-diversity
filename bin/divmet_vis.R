#### This is a script for visualizing exported alpha and beta diversity from qiime2.
# Made 25-03-2024

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

### Unused extras for now
## 5 Base R dendrogram
## 6 ape package
## 7 dendextend and piping
## 8 Bootstrap and p-values
#-------------------------------------------------------------------------------------------------------------
  
### Diversity metrics visualizations
## 1 Dependencies 
#install.packages("ggplot2", "ggdendro", "dplyr")
library(ggplot2)
library(ggdendro)
library(dplyr)

## 2 Load metadata
meta_raw <- read.delim("C:/Users/chuis/Documents/Biology_BioSus/Masterstage_2/MP divmet/sample-metadata.tsv", header = TRUE)
meta <- meta_raw[-c(1),]


## 3 Alpha diversity
## 3.1 Get data and prepare metadata
# Get alpha diversity faith pd
ad <- read.delim("C:/Users/chuis/Documents/Biology_BioSus/Masterstage_2/MP divmet/alpha-diversity.tsv")

# Adjust metadata to alpha diversity dataset
# Not all metadata is needed because of truncation and it needs ordering
alpha_ind <- which(meta$sample.id %in% ad$X.SampleID)
meta_afilt <- meta[alpha_ind,]

orda <- match(ad$X.SampleID, meta_afilt$sample.id)
meta_orda <- meta_afilt[orda,]

# 3.2 Visualize faith pd in ggplot heatmap
# Get vertical labeling for x-axis
ggplot(ad, aes(x = X.SampleID, y = meta_orda$body.site, fill = faith_pd)) +
  geom_tile(color = 'black') + 
  scale_fill_gradient(low = "palegreen", high = "palegreen4") +
  coord_fixed() + 
  labs(y = "Body site", x = "Sample ID", fill = "Faith PD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# 3.3 Maybe another way? See DMB Notes

## 4 Beta diversity
# 4.1 Get data and prepare metadata
# Get beta diversity distance matrix from somewhere and cluster with hclust
beta_df <- read.delim("C:/Users/chuis/Documents/Biology_BioSus/Masterstage_2/MP divmet/distance-matrix.tsv")
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
# Integrate with ggplot2, add coloration
# Mind the fixed ylim
ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label, col = meta_ordb$body.site),
            hjust = 1, angle = 90, size = 3) +
  labs(col = "Body site") +
  ylim(-0.5, 1)



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
  
#-----------------------------------------------------------------------------------------------------------------------------------
##Extra/leftovers
# Dependencies
# install.packages("ape", "dendextend", "pvclust")
library(ape)
library(dendextend)
library(pvclust)

# 5 Base R plot
# Extra customization, type can be "cladogram", "rectangle"/ommit, "unrooted", "fan", "radial" probably
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
edgePar <- list(col = 2:3, lwd = 2:1)

plot(dend, type = "rectangle", ylab = "Height", nodePar = nodePar, leaflab = "none", horiz = TRUE, edgePar = edgePar)

## 6 With ape package
# I'd rather just do ggplot at this point qq
# Colored cluster customization
# Type can be "cladogram", "rectangle"/ommit, "unrooted", "fan", "radial"
colors = c("red", "blue", "green")
clus3 = cutree(hc, 3)

phylo <- as.phylo(hc)
plot(phylo, type = "cladogram", cex = 0.5, label.offset = 0.5, tip.color = colors[clus3], edge.color = "steelblue", edge.width = 2, edge.lty = 2) 

## 7 Use dendextend
# Dendextend customization went awry, so no
#bodsit_fac <- as.factor(meta_ord$body.site) 
#bodsit_num <- as.numeric(bodsit_fac)

# Node, edges and leaves points/shapes
dendgg <- dend %>%
  #set("labels", c("a, "b")) %>%
  set("labels_colors", bodsit_fac, order_value = TRUE)
#set("nodes_pch", 19) %>%  
#set("nodes_cex", 1) %>%  
#set("nodes_col", "red") %>% 
#set("leaves_cex", 2) %>%  
#set("leaves_col", "blue") %>% 
#set("leaves_pch", c(17, 18, 19)) %>%
#set("branches_k_color", value = c("purple", "black", "cyan3"), k = 3)

# Pipe into ggplot2
# Theme can also be null
ggd1 <- as.ggdend(dendgg)

## 8 In case you need bootstrap shit
rartab_raw <- read.delim("C:/Users/chuis/Documents/Biology_BioSus/Masterstage_2/MP divmet/feature-table.tsv", skip = 1)
rartab <- rartab_raw[,-1]
rownames(rartab) <- rartab_raw[,1]

head(rartab)

set.seed(1234)
result <- pvclust(rartab, method.dist="cor", 
                  method.hclust="average", nboot=10)
plot(result)
#pvrect(result)

# Throw in some dendextend, can customize further
result %>% as.dendrogram %>% 
  set("branches_k_color", k = 3, value = c("purple", "orange", "cyan3")) %>%
  plot
result %>% text
#result %>% pvrect
  

