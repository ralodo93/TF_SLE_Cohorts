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
  
  anncolsle = anncol[names(partition),] %>% mutate(cluster = paste0("cluster",partition)) %>% 
    rownames_to_column("sample") %>% arrange(cluster) %>% mutate(sample = factor(sample, levels = unique(sample)))
  
  gene_order = as.dendrogram(hclust(dist(tfmatrix)))
  gene_order = rownames(tfmatrix)[order.dendrogram(gene_order)]
  
  dhc <- as.dendrogram(hclust(dist(t(tfmatrix[,anncolsle$sample])),method = "average"))
  dend_data <- dendro_data(dhc, type = "rectangle")
  
  dend_data$labels
  
  library(ggdendro)
  
  p = ggdendrogram(dhc)
  
  
  tf_tidy = tfmatrix %>% as.data.frame() %>% rownames_to_column("tf") %>% pivot_longer(cols = anncolsle$sample,names_to = "sample",values_to = "activity") %>%
    mutate(tf = factor(tf, levels = gene_order)) %>%
    left_join(anncolsle) %>% arrange(cluster) %>% mutate(sample = factor(sample, levels = unique(sample)))
  
  
  heat = ggplot(tf_tidy, aes(x = sample, y = tf))+
    geom_tile(aes(fill = activity))+
    scale_fill_gradient2()+
    theme_ipsum(plot_margin = margin(2,5,0,5))+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 3, angle = 90))
  
  plot_grid(plotlist = list(p,heat), align = "hv",nrow = 2)


  neu = ggplot(tf_tidy %>% dplyr::select(sample,Neu) %>% unique(), aes(x = sample, y = 1))+
    geom_tile(aes(fill = Neu))+
    scale_fill_gradient(low = "#ebf9fa", high = "#2ca162")+
    theme_ipsum(plot_margin = margin(0,5,0,5))+
    theme(axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank())
  lym = ggplot(tf_tidy %>% dplyr::select(sample,Lym) %>% unique(), aes(x = sample, y = 1))+
    geom_tile(aes(fill = Lym))+
    scale_fill_gradient(low = "#e9f9ff", high = "#8a57aa")+
    theme_ipsum(plot_margin = margin(0,5,0,5))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  clust = ggplot(tf_tidy %>% dplyr::select(sample,cluster) %>% unique(), aes(x = sample, y = 1))+
    geom_tile(aes(fill = cluster))+
    scale_fill_manual(values = c("cluster1" = "gold", "cluster2" = "skyblue4"))+
    theme_ipsum(plot_margin = margin(0,5,0,5))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  

  ful = plot_grid(plotlist = list(p,clust,neu,lym,heat), align = "v", nrow = 5, rel_heights = c(0.3,0.1,0.1,0.1,0.8))
  
  ggsave(ful, filename = "a.tiff", width = 16.5, height = 10, units = "cm")
  
  
  return(partition)
})
