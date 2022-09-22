install.packages("devtools")
install.packages("tibble")
install.packages("igraph")
install.packages("limma")
install.packages('Cairo')


library(Cairo)
library(nichenetr)
library(Seurat) 
library(tidyverse)
library(circlize)
library(openxlsx)
library(dplyr)
library(dittoSeq)

#ligand-target matrix 
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

#Transform human lr into mouse 
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()


#Read in expression data from RDS
seuratObj= readRDS('Seurat_objects/SeuratObject_astrocyte_0.75clust_Dim30_finalcelltype.rds')
Idents(object = seuratObj) <- seuratObj@meta.data$seurat_clusters
seuratObj_AS <- subset(seuratObj, subset = P2ry12 < 0.1)
seuratObj_AS <- subset(seuratObj_AS, subset = Plp1 < 0.1)


#######---------Performing NichenetR----------#######
#Select the Reciever cell 
Idents(object = seuratObj_AS) <- seuratObj_AS@meta.data$seurat_clusters
receiver = c("3","5") # Astrocyte Clusters
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj_AS, pct = 0.1)  #setting threshold for the number of genes 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#Define sender cells 
expressed_genes_sender = DifExpGenes_APP_APOE4_age4_pval05 %>% pull(gene) #DEGs from APP/PS1 mice (4months old) APOE4cKO vs APOE4KI 

######-------Defining potentially differentially expressed cells------########
seurat_obj_receiver= subset(seuratObj_AS, idents = receiver)
seurat_obj_receiver =SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Condition"]])

condition_ApoeKO = "KO"
condition_ApoeKI = "KI" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_ApoeKO, ident.2 = condition_ApoeKI, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


######------Define a set of potential ligands---------#########
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

potential_ligands
expressed_ligands
expressed_receptors

#########--------- Perform NicheNet ligand activity analysis ---------###########
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = base::rank(dplyr::desc(pearson)))
ligand_activities


best_upstream_ligands = ligand_activities %>% top_n(25, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
best_upstream_ligands

########----------Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis------###########
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 30) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.50)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_ligands
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
order_targets
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.001,0.0090))
p_ligand_target_network

write.xlsx(active_ligand_target_links_df,'ligand-target_mg-all_as_0.1plp1,0.1p2ry12.xlsx')

#Receptors of top ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

write.xlsx(lr_network_top_df_large,'ligand-receptor_mg_as_0.1plp1,0.1p2ry12.xlsx')






