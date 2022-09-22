######-------Load required packages and data-----#####
library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(ComplexHeatmap)
sessionInfo()

#########---------Making Circlos Plots---------########

#open data
library(readxl)
df <- read_excel("ligand-receptor_mg_as.xlsx",sheet = "final")

#Select Data
df <- df   %>% 
  dplyr::select(c("from","to", "Combined/10+weight"))
df <- df %>% arrange(desc("Combined/10+weight"))

#Setting Colour scheme
col_fun = colorRamp2(c(-2,0,2), c("blue","whitesmoke","red"))
col_fun(seq(-1.5,1, by = 0.01))

#Customize graphic parameters before initializing the circlos: 
circos.clear()
circos.par(canvas.xlim = c(-2, 2), 
           canvas.ylim = c(-2, 2),
           track.margin= c(0.01, 0.01),
           start.degree = -90,     #rotating circlos plot # of °
           #gap.degree = 0.8,
           "track.height" = 0.1)


#rotating circlos plot by 90°
circos.par(start.degree = -90)

#Assigning grid and annotation regions / size
chordDiagram(df, big.gap = 20, small.gap = 3,
             annotationTrack = c('grid','names'), 
             col=col_fun, 
             annotationTrackHeight = mm_h(1), 
             preAllocateTracks = list(track.height = mm_h(2)),h.ratio=0.6,
             transparency = 0,  
             directional = 1, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.zindex = rank(df[[3]]))

#Assign Annotations 
circos.track(track.index = 1, panel.fun = function(x,y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+ mm_y(5), CELL_META$sector.index,
              cex = 1,
              facing = 'clockwise',
              niceFacing = T, 
              adj = c(0,0.9))
}, bg.border = NA)


#adding legend 
lgd_links = Legend(at = c( -2,0, 2), col_fun = col_fun, 
                   title_position = "topleft", title = "LR_weight")
draw(lgd_links, x = unit(6, "mm"), y = unit(6, "mm"), just = c("left", "bottom"))

#highlighting sectors 
highlight.sector(df$from, track.index = 1, col = 'orange', text = 'Microglia Ligands', cex = 0.6, text.col = 'black', facing = 'bending.inside', niceFacing = T)
highlight.sector(df$to, track.index = 1, col = 'lightblue', text = 'Astrocyte Receptors', cex = 0.6, text.col = 'black', facing = 'bending.inside', niceFacing = T)



