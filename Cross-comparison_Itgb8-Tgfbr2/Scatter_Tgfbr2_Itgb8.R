####### Making Venn Diagram ######
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

df <- Tgfbr2_Itgb8_KO_comparison
ggVennDiagram(df)

######## Scatterplot Tgfbr2 and Itgb8 - KO in microglia ######
install.packages("ggpubr")
library(ggplot2)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggrepel)

df <- Scatter_comparison 
df <- df %>% remove_missing() %>% 
  filter(ITGB8_LOG2FC > 0.5 | ITGB8_LOG2FC < -0.5) %>%
  filter(TGFBR2_LOG2FC >0.5 | TGFBR2_LOG2FC < -0.5)

class(df$ITGB8_LOG2FC)
class(df$TGFBR2_LOG2FC)

theme <- theme(aspect.ratio = 1, 
             panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_line(linetype = "dashed"),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"),
             plot.title = )
?element_line

MGnD_label <- df2$`Overlap MGnD-DEGs`
M0_label <- df3$`Overlap M0-DEGs`
Total_overlap <- df4$Gene

MGnD_data <- df %>% filter(GENE %in% c(MGnD_label))
M0_data <- df %>% filter(GENE %in% c(M0_label))
Total_overlap_data <- df %>% filter(GENE %in% c(Total_overlap))

MGnD_label_selected <- df2$`Selected MGnD-DEGs` 
M0_label_selected <- df3$`Selected M0-DEGs`
MGnD_data_selected <- df %>% filter(GENE %in% c(MGnD_label_selected))
M0_data_selected <- df %>% filter(GENE %in% c(M0_label_selected))

legend <- tibble(x = c(-10,-5),
                 y = c(5,10))

scatter <- ggplot(df, aes(x= ITGB8_LOG2FC, y= TGFBR2_LOG2FC, label = GENE)) +
  geom_smooth(data = MGnD_data, method = "lm", se=TRUE, color="red", formula = y ~ x) +
  stat_cor(data = MGnD_data, color = 'red', label.x.npc = "left", label.y.npc = "top") +
  geom_smooth(data = M0_data, method = "lm", se=TRUE, color="dodgerblue", formula = y ~ x) +
  stat_cor(data = M0_data, color = 'dodgerblue', label.x.npc = "center", label.y.npc = "bottom") +
  geom_point(size=0.5, color = 'grey') +
  geom_point(data = Total_overlap_data, color = 'black', size = 0.5)+
  geom_point(data = MGnD_data, color = 'red', size = 1) +
  geom_point(data = M0_data, color = 'dodgerblue', size = 1)+
  theme +
  labs(title = "Comparison of DEGs in Microglia from Itgb8-KO and Tgfbr2-KO",
       x = "Log2FC (Itgb8-KO in Astrocytes)",
       y = "Log2FC (Tgfbr2-KO in Microglia)") +
  #geom_label(data = MGnD_data_selected, color = 'red', size = 2, label.padding = unit(0.1,"lines")) +
  #geom_label(data = M0_data_selected, color = 'blue', size = 2, label.padding = unit(0.1, "lines")) +
  geom_label_repel(data = MGnD_data_selected, 
                   color = 'red', size = 3, 
                   min.segment.length = 0,
                   max.overlaps = Inf,
                   label.padding = unit(0.1,"lines")) +
  geom_label_repel(data = M0_data_selected, 
                   color = 'dodgerblue',size = 3, 
                   min.segment.length = 0,
                   max.overlaps = Inf,
                   label.padding = unit(0.1,"lines"))
scatter 
rm(scatter)


#### Creating Circlos plot #####
######-------Load required packages and data-----#####
library(tidyverse)
library(circlize)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(ComplexHeatmap)

#########---------Making Circlos Plots---------########

#open data
df <- Circlos

#Select Data
df <- df   %>% 
  select(c(`Ligands`, `Receptor`, `Combined/10`))
df <- df %>% arrange(desc(`Combined/10`))

#Setting Colour scheme
col_fun = colorRamp2(c(-2,0,2), c("blue","whitesmoke","red"))
col_fun(seq(-2,2, by = 0.01))

#Customize graphic parameters before initializing the circlos: 
circos.clear()
circos.par(canvas.xlim = c(-1, 1), 
           canvas.ylim = c(-1, 1),
           track.margin= c(0.01, 0.01),
           start.degree = -90,     #rotating circlos plot # of °
           #gap.degree = 0.3,
           "track.height" = 0.1
           #gap.after = c("Wfdc17"=113)
           )

#circos.initialize(): allocates sectors on the circle.
#circos.track(): creates plotting regions for cells in one single track.
grid.col = c(Tgfbr2_KO_DOWN = "purple", Tgfbr2_KO_UP = "red", Itgb8_KO_UP = "orange", `Itgb8_KO_DOWN` = "dodgerblue")
             
#Assigning grid and annotation regions / size
chordDiagram(df, big.gap = 10, 
             annotationTrack = c('grid','names'), 
             annotationTrackHeight = mm_h(2), 
             preAllocateTracks = list(track.height = mm_h(4)),h.ratio=0.7,
             transparency = 0,
             grid.col = grid.col,
             directional = -1)
             # direction.type = c("diffHeight", "arrows"),
             # link.arr.type = "big.arrow"
             # scale = TRUE)
#abline(v = 0, lty = 2, col = "#00000080")
#circos.par

#Assign Annotations 
circos.track(track.index = 1, panel.fun = function(x,y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+ mm_y(1), CELL_META$sector.index,
              cex = 0.3,
              facing = 'clockwise',
              niceFacing = T, 
              adj = c(0,0.9))
}, bg.border = NA)


#adding legend 
# lgd_links = Legend(at = c( -2,0,2), col_fun = col_fun, 
#                    title_position = "topleft", title = "")
#draw(lgd_links, x = unit(6, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))

#highlighting sectors 
highlight.sector(df$Gene, track.index = 2, col = 'lightblue', text = 'Differentially expressed genes', cex = 0.4, text.col = 'black', facing = 'bending.inside', niceFacing = T)
highlight.sector(df$Receptor, track.index = 1, col = 'chartreuse3', text = 'ASTROCYTE TARGET GENES', cex = 0.7, text.col = 'black', facing = 'bending.inside', niceFacing = T)










