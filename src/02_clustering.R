library(glue)
library(NbClust)
library(tidyverse)
library(stringr)
library(hrbrthemes)
library(cowplot)
age_group = "adult"
partitions = lapply(c("pediatric","adult"), function(age_group){
  load(glue("useful_data/{age_group}_info.RData"))
  sle = rownames(anncol %>% filter(Disease == "SLE"))
  tfmatrix = get(glue("TFActivities{str_to_title(age_group)}"))[,sle]
  bn = NbClust(t(tfmatrix),method = "average",index = "ch")
  partition = bn$Best.partition
  bn = bn$Best.nc
  print(bn)
  
  library(pheatmap)
  
  dhc <- as.dendrogram(hclust(dist(t(tfmatrix)),method = "average"))
  ddata <- dendro_data(dhc, type = "rectangle")
  
  labels_order <- ddata$labels$label
  partition <- partition[labels_order]
  
  anncolsle = anncol[labels_order,] %>% mutate(cluster = paste0("cluster",partition))
  library(ggplotify)
  plt_heat <- as.ggplot(pheatmap::pheatmap(tfmatrix, clustering_method = "average", treeheight_row = 0,
                     show_rownames = F, show_colnames = F,
                     annotation_col = anncolsle[,c("cluster","Neu","Lym")],
                     annotation_colors = list("cluster" = c("cluster1" = "gold", "cluster2" = "skyblue4"))))
  
  return(list)
})
