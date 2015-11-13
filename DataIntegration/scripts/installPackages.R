list.of.packages <- c("ggplot2", "pheatmap", "dplyr", "reshape2", "igraph", "linkcomm", "httr", "cluster.datasets")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)
