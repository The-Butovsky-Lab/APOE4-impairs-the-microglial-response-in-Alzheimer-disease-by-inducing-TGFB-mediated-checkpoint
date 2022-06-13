

rm(list = ls())
setwd(system("pwd", intern = T) )
#### _____________________________________________________Manual Input _____________________________________________________ ####

# Type of Species 
Species <- "mouse"            #<- Either "mouse" or "human" <- make sure to put within quotation marks


# Location of Files for analysis
Location_of_quant_Data_Files <-""  #  <--- Find the location of the folder (sometimes called transcripts_quant) that contains folders for all samples with quant.sf file in them
Name_of_Metadata_File <- ".xlsx"  # name of your metadata xlsx file. 


# Extracting data from Metadata excel file 
Position_Column_Name <- "Well"        #  <--- Well ID/Position on Plate
Number_of_variables_to_compare <- 1   #  <--- Can be 1, 2 or 3 different variables to compare. e.g., looking at effects of disease and cre is a 2 comparison 
Metadata1_Column_Name <- ""       #  <--- The more significant comparison (example: APP vs WT) 

# ____________Variable 2 ONLY for multiple comparison (if you have 2 or 3 different variables to compare) ____________
Metadata2_Column_Name <- ""   #  <--- The sub-comparison/secondary comparison (example:  APOE3 vs APOE4)

#  ____________Variable 3 ONLY for multiple comparison (if you have 3 different variables to compare)____________
Metadata3_Column_Name <- ""         #  <--- The last-comparison/most "sub" comparison (example:  KI vs KO)


# Q:  Do you have column for Order in your Data Set?
Is_Order_in_excelfile <- 'No' # Is there a column name that groups each group (from 1-x), number as you like to see on heatmap, with 1 to the far left 
Order_column <- 'Order' # name of order column in excel sheet

# Q:  Do you have Outliers in your Data Set?
Outliers_in_Data <- "No"              #  <--- Are there outliers in your data? "Yes" or "No"
Outliers_Column_Name <- "Removed"     # <--- For any outliers, indicate in metadata Excel Table by creating a new column Titled "Removed". Input "Yes" for any samples removed


# Q:  Do you have different Genders in your Data Set and do not want gender effect to influence your DEGs? 
Genders_in_Data <- "No"              #  <--- Are there different genders in your data? "Yes" or "No"
Gender_Column_Name <- "Gender"       #  <--- Name of your gender column. Recommend to name "Gender"

# Q:  Do you have Genes of Interest in your Data Set?
Select_Genes_in_Data <- "No"         #  <--- Are there different genders in your data? "Yes" or "No"
Select_Genes_Column_Name <- "Genes"  #  <--- Name of your "genes of interest" column. Recommend to name "Genes"


# Q:  Do you have 2 different plates for this run that requires batch effect removal?
Batch_effect <- "No"         #  <--- Do your samples come from seperate plates? "Yes" or "No"
Plate_Column_Name <- "Plate"  #  <--- Name of your "plate" column. Recommend to name "Plate"
#Location of Files for analysis (Plate 2 if applicable)
# PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# 
Location_of_quant_Data_Files_PLATE2 <-"/Users/Dania/Desktop/SK-4W7B Rafa, Wesley, Zhuoran, Neta, Xiaoming 5_03_2022/DM Analysis/transcripts_quant"  #  <--- Find the location of the folder (sometimes called transcripts_quant) that contains folders for all samples with quant.sf file in them
Name_of_Metadata_File_PLATE2 <- "metadata_APP_Lgals3 inj_Astro_plate2.xlsx"  # name of your metadata xlsx file. 
# PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# # PLATE 2 ONLY!!!# 



# Statistics for Data Analysis (genes have to pass these requirements to be considered "DEGs")
Statistics <- "padj"               # <---- either "pvalue" or "padj"
alpha <- 0.05                        # <---- alpha cutoff for statistical signifigance 
logFC_cutoff <- 0                    # <--- logFC_cutoff (eg. 0.5)


# Plots & Figures
Do_you_want_heatmaps <- "Yes"        # <----  "Yes" or "No"
Do_you_want_results_files <- "Yes"    # <----  "Yes" or "No"
Do_you_want_barplots_for_select_genes <- "No"      # <----  "Yes" or "No"

#Hearmaps Specific Inputs
heatmap_number_of_groups <- 2
number_of_top_genes <- 100 # number of top genes in heatmap

#Barplots Specific Inputs
ErrorBars<-'SE'  # <- 'SE' : standard error, 'SD': standard deviation (for barplots of specific genes)
# add as many colors as you want. Colors go in order. If you have 2 samples, color 1 and color 2 will be the colors for your bars. Some color options in link below
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
color1<- 'lightsteelblue'
color2<-'salmon1'
color3<-'blue'
color4<-'darkolivegreen3'
color5<- 'goldenrod3'
color6<- 'aquamarine1'
color7<-'brown4'
color8<-'lightcoral'
color9<-'palegreen2'
color10<- 'palevioletred'


#Order of metadata


Z1 <- TRUE # <- If order of metadata 1 in your figures is flipped, change (from TRUE to FALSE) or  vice versa
Z2 <- TRUE # <- If order of metadata 2 in your figures is flipped, change (from TRUE to FALSE) or  vice versa
Z3 <- TRUE # <- If order of metadata 3 in your figures is flipped, change (from TRUE to FALSE) or  vice versa
G1 <- TRUE # <- If order of gender in your figures is flipped, change (from TRUE to FALSE) or  vice versa

# _____________________________________________________End of manual Input ___________________________________________________ 







############################## Step 1:  INSTALL/LOAD REQUIRED PACKAGES ###########################

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.13")
}

list.of.packages <- c('tidyverse', 'DESeq2', 'readxl', 'readr', 'tximport', 'tximeta', 'apeglm', 'pheatmap', 'RColorBrewer',
                      'biomaRt', 'EnhancedVolcano', 'sva', 'writexl', 'openxlsx','ggpubr', 'qpcR')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

#setwd("~/Dropbox (Partners HealthCare)/  1-CAF Project/3-CAF RESULTS/2-RAW DATA/2-RNAseq/Lgals_Injection/Astrocytes") <- find code for file directory 



#### _____________________________________________________Loading in Data _____________________________________________________ ####
dir.create("plots")
dir.create("results")

if (Batch_effect == "No") {
  
  b <- function(name){
    paste("metadata$",name, sep="")
  }
  
  Variable1 <- b(Metadata1_Column_Name)
  if (Number_of_variables_to_compare >= 2) {
    Variable2 <- b(Metadata2_Column_Name)
  }
  
  if (Number_of_variables_to_compare == 3) {
    Variable3 <- b(Metadata3_Column_Name)
  }
  
  if (Genders_in_Data == "Yes") {
    Gender <- b(Gender_Column_Name)
  }
  
  if (Is_Order_in_excelfile == "Yes") {
    order <- b(Order_column)
  }
  
  
  Position <- b(Position_Column_Name)
  Outliers<- b(Outliers_Column_Name)
  Genes<- b(Select_Genes_Column_Name)
  metadata <- read_xlsx(Name_of_Metadata_File)
  if (Outliers_in_Data == "Yes") {
    idxOutliers <- grep("Yes", eval(parse(text = Outliers)))
    metadata <-metadata[-(c(idxOutliers)), ]
  }
  
  if (Is_Order_in_excelfile == "No") {
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)), variableA = eval(parse(text = Variable1)) )
      }
    }
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)))
      }
    }
    
    
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)) )
      }
    }
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)))
      }
    }
    
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)),variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)), variableC = eval(parse(text = Variable3)) )
      }
    }
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)), variableC = eval(parse(text = Variable3)) )
      }
    }
  }
  
  if (Is_Order_in_excelfile == "Yes") {
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)), variableA = eval(parse(text = Variable1)),Order=eval(parse(text = order)) )
      }
    }
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)),Order=eval(parse(text = order)))
      }
    }
    
    
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)) ,Order=eval(parse(text = order)))
      }
    }
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)),Order=eval(parse(text = order)))
      }
    }
    
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "No") {
        sample_data <- data.frame(position=eval(parse(text = Position)),variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)), variableC = eval(parse(text = Variable3)) ,Order=eval(parse(text = order)))
      }
    }
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "Yes") {
        sample_data <- data.frame(position=eval(parse(text = Position)),gender=eval(parse(text = Gender)), variableA = eval(parse(text = Variable1)), variableB = eval(parse(text = Variable2)), variableC = eval(parse(text = Variable3)),Order=eval(parse(text = order)) )
      }
    }
  }
  all_files <- list.files(Location_of_quant_Data_Files, full.names= T, pattern=NULL, all.files= FALSE)
  quant_files <- file.path(all_files, "quant.sf")
  position_list <- paste0(sample_data$position, collapse = '|')
  sample_files <- grep(position_list, quant_files, value = TRUE)
  sample_data <- na.omit(sample_data) 
}




if (Batch_effect == "Yes") {
  b1 <- function(name){
    paste("metadata1$",name, sep="")
  }
  b2 <- function(name){
    paste("metadata2$",name, sep="")
  }
  
  Variable1_1 <- b1(Metadata1_Column_Name)
  Variable1_2 <- b2(Metadata1_Column_Name)
  
  if (Number_of_variables_to_compare >= 2) {
    Variable2_1 <- b1(Metadata2_Column_Name)
    Variable2_2 <- b2(Metadata2_Column_Name)
  }
  
  if (Number_of_variables_to_compare == 3) {
    Variable3_1 <- b1(Metadata3_Column_Name)
    Variable3_2 <- b2(Metadata3_Column_Name)
  }
  
  if (Genders_in_Data == "Yes") {
    Gender1 <- b1(Gender_Column_Name)
    Gender2 <- b2(Gender_Column_Name)
  }
  
  if (Is_Order_in_excelfile == "Yes") {
    order1 <- b1(Order_column)
    order2 <- b2(Order_column)
  }
  
  
  Position1 <- b1(Position_Column_Name)
  Outliers1<- b1(Outliers_Column_Name)
  Genes<- b1(Select_Genes_Column_Name)
  Position2 <- b2(Position_Column_Name)
  Outliers2<- b2(Outliers_Column_Name)
  Batch1<- b1(Plate_Column_Name)
  Batch2<- b2(Plate_Column_Name)
  
  
  metadata1 <- read_xlsx(Name_of_Metadata_File)
  metadata2 <- read_xlsx(Name_of_Metadata_File_PLATE2)
  if (Outliers_in_Data == "Yes") {
    idxOutliers1 <- grep("Yes", eval(parse(text = Outliers1)))
    idxOutliers2 <- grep("Yes", eval(parse(text = Outliers2)))
    metadata1 <-metadata1[-(c(idxOutliers1)), ]
    metadata2 <-metadata2[-(c(idxOutliers2)), ]
  }
  
  if (Is_Order_in_excelfile == "No") {
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)), variableA = eval(parse(text = Variable1_1)) ,batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)), variableA = eval(parse(text = Variable1_2)) ,batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2) 
        
      }
    }
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)),batch=eval(parse(text = Batch2)))
        sample_data2<- na.omit(sample_data2)
        
      }
    }
    
    
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)) ,batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)) ,batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)), variableC = eval(parse(text = Variable3_1)) ,batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)), variableC = eval(parse(text = Variable3_2)) ,batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)), variableC = eval(parse(text = Variable3_1)) ,batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)), variableC = eval(parse(text = Variable3_2)) ,batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
  }
  
  if (Is_Order_in_excelfile == "Yes") {
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)), variableA = eval(parse(text = Variable1_1))  ,Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)), variableA = eval(parse(text = Variable1_2))  ,Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)),Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)),Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    
    
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)) ,Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)) ,Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)),Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)),Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "No") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)), variableC = eval(parse(text = Variable3_1)) ,Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)), variableC = eval(parse(text = Variable3_2)) ,Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "Yes") {
        sample_data1 <- data.frame(position=eval(parse(text = Position1)),gender=eval(parse(text = Gender1)), variableA = eval(parse(text = Variable1_1)), variableB = eval(parse(text = Variable2_1)), variableC = eval(parse(text = Variable3_1)) ,Order=eval(parse(text = order1)),batch=eval(parse(text = Batch1)))
        sample_data1 <- na.omit(sample_data1)
        sample_data2 <- data.frame(position=eval(parse(text = Position2)),gender=eval(parse(text = Gender2)), variableA = eval(parse(text = Variable1_2)), variableB = eval(parse(text = Variable2_2)), variableC = eval(parse(text = Variable3_2)) ,Order=eval(parse(text = order2)),batch=eval(parse(text = Batch2)))
        sample_data2 <- na.omit(sample_data2)
        
      }
    }
  }
  
  sample_data <- rbind(sample_data1,sample_data2)
  all_files1 <- list.files(Location_of_quant_Data_Files, full.names= T, pattern=NULL, all.files= FALSE)
  quant_files1 <- file.path(all_files1, "quant.sf")
  position_list1 <- paste0(sample_data1$position, collapse = '|')
  sample_files1 <- grep(position_list1, quant_files1, value = TRUE)
  all_files2 <- list.files(Location_of_quant_Data_Files_PLATE2, full.names= T, pattern=NULL, all.files= FALSE)
  quant_files2 <- file.path(all_files2, "quant.sf")
  position_list2 <- paste0(sample_data2$position, collapse = '|')
  sample_files2 <- grep(position_list2, quant_files2, value = TRUE)
  sample_files<-c(sample_files1,sample_files2)
}




if (Species == "mouse") {
  ensembl_mart= useMart('ensembl', dataset="mmusculus_gene_ensembl") # connecting to ENSEMBL
  ensemble_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version","external_gene_name"), mart=ensembl_mart)
}

if (Species == "human") {
  ensembl_mart= useMart('ensembl', dataset="hsapiens_gene_ensembl") # connecting to ENSEMBL
  ensemble_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version","external_gene_name"), mart=ensembl_mart)
}



tx2gene <- ensemble_df %>% 
  as_tibble() %>% #here we convert data frame into tibble because we were getting an error without this line
  mutate(noVersion = sub('\\..*','', ensembl_transcript_id_version)) %>% # here we are creating a new column called noVersion which will replace ensembl_transcript_id_version because we are trying to remove the .1 at the end of the transcript name (because this is just the version of the gene naming) soo \\. means weverything that comes before the period we ignore, then .* means that we will be changing "mutating" the period and what comes after it, and will be replaced by what comes in the single paranth after the comma in this case we want to remove the .1 so we replace it with nothing '' 
  dplyr::select (noVersion, external_gene_name) # here we are selecting the 2 columns we want, the new ensembl_transcript_id_version without the version (noVersion) and  the gene names


txi <- tximport(sample_files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)  
txicounts<-as.matrix(txi$counts)
mode(txicounts) <- 'integer'

write.xlsx(as.data.frame(as.matrix(txi$counts)), paste("results/CountsTxi_",filenames2,"_comparison.xlsx",sep= "")
           , rowNames= T, overwrite = T)

if (Batch_effect == "Yes") {
  count_batch <- sample_data$batch
  count_adjusted <- ComBat_seq(txicounts, batch = count_batch, group = NULL)
  tpm_adjusted <- ComBat_seq(txi$abundance, batch = count_batch, group = NULL)
  txicounts<- count_adjusted
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "No") {
      dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA )  # condition * disease means condition + disease + condition:disease
    }
  }
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "Yes") {
      loadError=F
      a=try({dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*gender )})  # condition * disease means condition + disease + condition:disease
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+gender )  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "No") {
      loadError=F
      a=try({dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*variableB)})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+variableB)
      }
    }
  }
  
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "Yes") {
      loadError=F
      a=try({ dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*variableB+ variableA*gender+ variableB*gender )})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+variableB+gender)  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "No") {
      loadError=F
      a=try({ dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*variableB+ variableA*variableC + variableB*variableC )})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+variableB+variableC )  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "Yes") {
      loadError=F
      a=try({ dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*variableB+ variableA*variableC + variableB*variableC+ variableA*gender+ variableB*gender + variableC*gender )})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+variableB+variableC+gender)  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  
}

#count_batch <- sample_data$plate
#count_adjusted <- ComBat_seq(txicounts, batch = count_batch, group = NULL)
#tpm_adjusted <- ComBat_seq(txi$abundance, batch = count_batch, group = NULL)
if (Batch_effect == "No") {
  
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "No") {
      dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA )  # condition * disease means condition + disease + condition:disease
    }
  }
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "Yes") {
      loadError=F
      a=try({dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA*gender )})  # condition * disease means condition + disease + condition:disease
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA+gender )  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "No") {
      loadError=F
      a=try({dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA*variableB)})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA+variableB)
      }
    }
  }
  
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "Yes") {
      loadError=F
      a=try({ dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA*variableB+ variableA*gender+ variableB*gender )})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA+variableB+gender)  # condition * disease means condition + disease + condition:disease
      }
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "No") {
      loadError=F
      a=try({ dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA*variableB+ variableA*variableC + variableB*variableC )})
      loadError <- (is(a, 'try-error')|is(a,'error'))  
      if(loadError==TRUE){
        dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design =~ variableA+variableB+variableC )  # condition * disease means condition + disease + condition:disease
      }
    }
  }
}
if (Number_of_variables_to_compare == 3) {
  if (Genders_in_Data == "Yes") {
    loadError=F
    a=try({ dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA*variableB+ variableA*variableC + variableB*variableC+ variableA*gender+ variableB*gender + variableC*gender )})
    loadError <- (is(a, 'try-error')|is(a,'error'))  
    if(loadError==TRUE){
      dds <- DESeqDataSetFromMatrix(txicounts, colData = sample_data, design =~ variableA+variableB+variableC+gender)  # condition * disease means condition + disease + condition:disease
    }
  }
}


dds <- estimateSizeFactors(dds)
idx <- rowMeans(counts(dds, normalized= TRUE)) >= 5
dds <- dds[idx,]
if (Genders_in_Data == "Yes") {
  dds_LRT <- DESeq(dds, test = "LRT", reduced= ~ gender) #include reduced when doing LRT,reduced= ~ 1 means do not reduce/remove anything from analysis 
}
if (Genders_in_Data == "No") {
  dds_LRT <- DESeq(dds, test = "LRT", reduced= ~ 1) #include reduced when doing LRT,reduced= ~ 1 means do not reduce/remove anything from analysis 
}

DS_norm_counts <- counts(dds_LRT, normalized = TRUE)
mat1 <- DS_norm_counts
sample <- sample_data

if (Is_Order_in_excelfile == "No"){
  
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "No") {
      df0<- sample$variableA
      df <- sample$position
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      
      colnames(mat1) <- paste(df0,df,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$gender
      df2<- sample$position
      
      
      colnames(mat1) <- paste(df,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = G1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      colnames(mat1) <- paste(df0,df,df2,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0,df)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,2],decreasing = G1)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,1]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "No") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$position
      
      
      colnames(mat1) <- paste(df,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z2)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      colnames(mat1) <- paste(df0,df,df2,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0,df)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,2],decreasing = Z2)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,1]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$gender
      df3<- sample$position
      
      
      colnames(mat1) <- paste(df2,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = G1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z2)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df0,df,df2,df3,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0,df,df2)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,3],decreasing = G1)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,2],decreasing = Z2)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,1]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
      
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "No") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$variableC
      df3<- sample$position
      
      colnames(mat1) <- paste(df2,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z3)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z2)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      colnames(mat1) <- paste(df0,df,df2,df3,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0,df,df2)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,3],decreasing = Z3)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,2],decreasing = Z2)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,1]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
      
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$variableC
      df3<- sample$gender
      df4<- sample$position
      
      colnames(mat1) <- paste(df3,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = G1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      df4 <- df4[c(ix)]
      
      
      colnames(mat1) <- paste(df2,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z3)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      df4 <- df4[c(ix)]
      
      colnames(mat1) <- paste(df,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z2)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      df4 <- df4[c(ix)]
      
      colnames(mat1) <- paste(df0,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = Z1)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      df4 <- df4[c(ix)]
      
      colnames(mat1) <- paste(df0,df,df2,df3,df4,sep="-")
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df0,df,df2,df3)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,4],decreasing = G1)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,3],decreasing = Z3)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,2],decreasing = Z2)
      gaps <- gaps[c(ix)]
      ix<- order(gapss[,1],decreasing = Z1)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,1]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
}



if (Is_Order_in_excelfile == "Yes"){
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "No") {
      df1 <- sample$Order
      df0<- sample$variableA
      df <- sample$position
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      
      colnames(mat1) <- paste(df0,df,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  if (Number_of_variables_to_compare == 1) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$gender
      df2<- sample$position
      df1 <- sample$Order
      
      
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      
      colnames(mat1) <- paste(df0,df,df2,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  
  
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "No") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$position
      df1 <- sample$Order
      
      
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      
      
      colnames(mat1) <- paste(df0,df,df2,sep="-")
      
      cutoff<- heatmap_number_of_groups
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
  if (Number_of_variables_to_compare == 2) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$gender
      df3<- sample$position
      df1 <- sample$Order
      
      
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      
      colnames(mat1) <- paste(df0,df,df2,df3,sep="-")
      
      cutoff<- heatmap_number_of_groups
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
      
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "No") { 
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$variableC
      df3<- sample$position
      df1 <- sample$Order
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      
      
      colnames(mat1) <- paste(df0,df,df2,df3,sep="-")
      
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
      
    }
  }
  
  if (Number_of_variables_to_compare == 3) {
    if (Genders_in_Data == "Yes") {
      df0<- sample$variableA
      df <- sample$variableB
      df2<- sample$variableC
      df3<- sample$gender
      df4<- sample$position
      df1 <- sample$Order
      
      
      colnames(mat1) <- paste(df1,sep="-")
      
      colnamess<- colnames(mat1)
      ix<- order(colnamess,decreasing = FALSE)
      mat1 <- mat1[,c(ix)]
      
      df0 <- df0[c(ix)]
      df <- df[c(ix)]
      df2 <- df2[c(ix)]
      df3 <- df3[c(ix)]
      df4 <- df4[c(ix)]
      
      
      colnames(mat1) <- paste(df0,df,df2,df3,df4,sep="-")
      cutoff<- heatmap_number_of_groups
      
      gaps<- table(df1)
      gapss<- as.data.frame(gaps)
      ix<- order(gapss[,1],decreasing = FALSE)
      gaps <- gaps[c(ix)]
      
      gaps <- as.data.frame(gaps)
      gaps<-gaps[,2]
      
      i <- 0 
      
      for (var in gaps) {
        i <- i+1
        if (i==1)
          gaps[i]<-gaps[i]
        else
          gaps[i]<-gaps[i]+gaps[i-1]
      }
    }
  }
}




#### ___________________________________________________Statistics _____________________________________________________ ####

if (Number_of_variables_to_compare == 1) {
  if (Genders_in_Data == "No") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,sep= "_")
    filenames2<-Metadata1_Column_Name
    
  }
}
if (Number_of_variables_to_compare == 1) {
  if (Genders_in_Data == "Yes") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,"reduced4gender",sep= "_")
    filenames2<-paste(Metadata1_Column_Name,"reduced4gender",sep= "_")
    
  }
}


if (Number_of_variables_to_compare == 2) {
  if (Genders_in_Data == "No") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,Metadata2_Column_Name,sep= "_")
    filenames2<- paste(Metadata1_Column_Name,Metadata2_Column_Name,sep= "_")
    
  }
}
if (Number_of_variables_to_compare == 2) {
  if (Genders_in_Data == "Yes") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,Metadata2_Column_Name,"reduced4gender",sep= "_")
    filenames2<- paste(Metadata1_Column_Name,Metadata2_Column_Name,"reduced4gender",sep= "_")
    
  }
}

if (Number_of_variables_to_compare == 3) {
  if (Genders_in_Data == "No") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,Metadata2_Column_Name,Metadata3_Column_Name,sep= "_")
    filenames2<- paste(Metadata1_Column_Name,Metadata2_Column_Name,Metadata3_Column_Name,sep= "_")
    
  }
}
if (Number_of_variables_to_compare == 3) {
  if (Genders_in_Data == "Yes") {
    filenames<- paste(paste(Statistics,alpha,sep= ""),paste('_FC',logFC_cutoff,sep= ""),Metadata1_Column_Name,Metadata2_Column_Name,Metadata3_Column_Name,"reduced4gender",sep= "_")
    filenames2<- paste(Metadata1_Column_Name,Metadata2_Column_Name,Metadata3_Column_Name,"reduced4gender",sep= "_")
    
  }
}

resLRT <- results(dds_LRT, pAdjustMethod = 'BH', independentFiltering = T)

if (Statistics == "pvalue") {
  sig_res <- resLRT %>%
    data.frame() %>%
    rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
    as_tibble %>%
    filter(pvalue< alpha) %>%
    filter(!between(log2FoldChange, (-1*logFC_cutoff), logFC_cutoff)) %>%
    arrange(pvalue)
  
}


if (Statistics == "padj") {
  sig_res <- resLRT %>%
    data.frame() %>%
    rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
    as_tibble %>%
    filter(padj< alpha) %>%
    filter(!between(log2FoldChange, (-1*logFC_cutoff), logFC_cutoff)) %>%
    arrange(padj)
}
sorted_DEGenes <- sig_res$gene
resOrdered <- resLRT[order(resLRT$pvalue),]

if (Do_you_want_results_files == "Yes") {
  
  write.xlsx(as.data.frame(mat1), paste("results/CountsDS_",filenames2,"_comparison.xlsx",sep= "")
             , rowNames= T, overwrite = T)
  
  write.xlsx(as.data.frame(resOrdered), 
             paste("results/stats_allGenes_",filenames2,"_comparison.xlsx",sep= "")
             , rowNames= T, overwrite = T)
  
  write.xlsx(as.data.frame(sig_res), 
             paste("results/DifExpGenes_",filenames,"_comparison.xlsx",sep= ""),
             overwrite = T)
}

name_list <- c('gene', colnames(mat1))
name_list

sig_export <- mat1 %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene))

write.xlsx(as.data.frame(sig_export), 
           paste("results/CountsSig_",filenames,"_comparison.xlsx",sep= "")
           , overwrite = T)
#### ___________________________________________________Creating Heatmaps _____________________________________________________ ####
if (Do_you_want_heatmaps == "Yes") {
  miRNA_File <- paste("results/CountsSig_",filenames,"_comparison.xlsx",sep= "")
  cts <- read.xlsx(miRNA_File, rowNames = T)
  sample <- sample_data
  
  
  unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps','^Olf'), collapse = '|')
  
  cts_final <- cts %>%
    rownames_to_column('gene') %>%
    filter(!str_detect(gene, unwanted_genes)) %>%
    filter(!str_detect(.$gene, "Rik"))%>%
    column_to_rownames('gene')
  cts_subset <- cts_final %>% head(n = number_of_top_genes)
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  data_norm <- t(apply(round(cts_final), 1, cal_z_score))
  phm_full <- pheatmap(data_norm,
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(48),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = T,
                       cutree_row = cutoff,
                       cluster_cols = F,
                       gaps_col= gaps,
                       border_color = NA,
                       legend = TRUE,
                       show_rownames = F,
                       treeheight_col = 0,
                       treeheight_row = 0,
                       fontsize = 7,
                       scale = 'row')
  pdf(file = paste0("plots/heatmap_",filenames,"_comparison","_allDEGs.pdf"), pointsize = 10, width = 7, height =10)
  print(phm_full)
  dev.off()
  
  gene_clusters <- data.frame(sort(cutree(phm_full$tree_row, k=cutoff))) %>%
    rownames_to_column(var="gene") 
  
  write.xlsx(gene_clusters, 
             file= paste("results/GeneClusters_heatmap_",filenames,"_comparison.xlsx",sep= ""),overwrite = T)
  
  
  topGenes<- t(apply(round(cts_subset), 1, cal_z_score))
  heatmap_topgenes<- pheatmap(topGenes,
                              color = colorRampPalette(c("navy", "white", "firebrick3"))(48),
                              breaks = seq(-2, 2, by = 0.1),
                              kmeans_k = NA,
                              cluster_rows = T,
                              cutree_row = cutoff,
                              cluster_cols = F,
                              gaps_col= gaps,
                              border_color = NA,
                              legend = TRUE,
                              show_rownames = T,
                              treeheight_col = 0,
                              treeheight_row = 50,
                              fontsize = 7,
                              scale = 'row'
  )
  
  if (number_of_top_genes < 50){
    h = 5
  }
  if (number_of_top_genes < 120){
    h = 11
  }
  if (number_of_top_genes > 120){
    h = 20
  }
  if (number_of_top_genes > 200){
    h = 25
  }
  if (number_of_top_genes > 250){
    h = 30
  }
  
  pdf(file = paste0("plots/heatmap_",filenames,"_comparison","_top",number_of_top_genes,".pdf"), pointsize = 10, width = 7, height =h)
  print(heatmap_topgenes)
  dev.off()
}
#### ___________________________________________________Creating specific Gene bar plots _____________________________________________________ ####
if (Select_Genes_in_Data == "Yes") {
  
  metadata <- read_xlsx(Name_of_Metadata_File)
  genes_data <- data.frame(genes=eval(parse(text = Genes)))
  genes_data <- na.omit(genes_data) 
  genes_data<-as.character(genes_data[,1])
  
  genes_data2 <- data.frame(genes=eval(parse(text = Genes)))
  genes_data2 <- na.omit(genes_data2) 
  sig_export2 <- mat1 %>%
    data.frame() %>%
    rownames_to_column(var = 'genes') %>%
    as_tibble() %>%
    semi_join(genes_data2)
  
  write.xlsx(as.data.frame(sig_export2), 
             paste("results/CountsSig_Genes_of_Interest",filenames,"_comparison.xlsx",sep= "")
             , overwrite = T)
  
  miRNA_File_GenesofInteres <- paste("results/CountsSig_Genes_of_Interest",filenames,"_comparison.xlsx",sep= "")
  cts_SpecificGenes <- read.xlsx(miRNA_File_GenesofInteres, rowNames = T)
  sample <- sample_data
  
  
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  data_norm2 <- t(apply(round(cts_SpecificGenes), 1, cal_z_score))
  
  heatmap_genesofinterest<- pheatmap(data_norm2,
                                     color = colorRampPalette(c("navy", "white", "firebrick3"))(48),
                                     breaks = seq(-2, 2, by = 0.1),
                                     kmeans_k = NA,
                                     cluster_rows = T,
                                     cutree_row = 1,
                                     cluster_cols = F,
                                     gaps_col= gaps,
                                     border_color = NA,
                                     legend = TRUE,
                                     show_rownames = T,
                                     treeheight_col = 0,
                                     treeheight_row = 50,
                                     fontsize = 7,
                                     scale = 'row'
  )
  
  pdf(file = paste0("plots/heatmap_Genes_of_Interest",filenames,"_comparison.pdf"), pointsize = 10, width = 5, height =10)
  print(heatmap_genesofinterest)
  dev.off()
  
  for (i in genes_data){
    plot_list = list()
    
    GeneOfInterest <- as.data.frame(mat1) %>%
      rownames_to_column('gene') %>%
      filter(str_detect(gene, i)) %>%
      column_to_rownames('gene')
    
    
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "No") {
        colnames(GeneOfInterest)<- paste(df0)
        variables<- paste(df0)    }
    }
    if (Number_of_variables_to_compare == 1) {
      if (Genders_in_Data == "Yes") {
        colnames(GeneOfInterest)<- paste(df0,df)
        variables<- paste(df0,df)     }
    }
    
    
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "No") {
        colnames(GeneOfInterest)<- paste(df0,df)
        variables<- paste(df0,df)     }
    }
    if (Number_of_variables_to_compare == 2) {
      if (Genders_in_Data == "Yes") {
        colnames(GeneOfInterest)<- paste(df0,df,df2)
        variables<- paste(df0,df,df2)     }
    }
    
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "No") {
        colnames(GeneOfInterest)<- paste(df0,df,df2)
        variables<- paste(df0,df,df2)      }
    }
    if (Number_of_variables_to_compare == 3) {
      if (Genders_in_Data == "Yes") {
        colnames(GeneOfInterest)<- paste(df0,df,df2,df3)
        variables<- paste(df0,df,df2,df3)      }
    }
    
    
    GeneOfInterest<- t(GeneOfInterest)
    rownames(GeneOfInterest)<- NULL
    variables<- as.data.frame(variables)
    GeneOfInterest<- cbind(GeneOfInterest,variables)
    
    colors<- c(color1,color2,color3,color4,color5,color6,color7,color8,color9,color10)
    
    if (ErrorBars == 'SE') {
      barplot<- ggbarplot(
        GeneOfInterest, x = "variables", y =i , 
        add = c("mean_se", "jitter"),
        fill = "variables",
        palette  = colors
      )
    }
    
    if (ErrorBars == 'SD') {
      barplot<- ggbarplot(
        GeneOfInterest, x = "variables", y ==i , 
        add = c("mean_sd", "jitter"),
        fill = "variables",
        palette  = colors
      )
    }
    
    
    pdf(file = paste0("plots/BarPlot",Metadata1_Column_Name,i,"_comparison.pdf"), pointsize = 10, width = 8, height =6)
    print(barplot)
    dev.off()
    
    
    
  }
  
  
}




#### ___________________________________________________Creating PCA Plots _____________________________________________________ ##


rld <- vst(dds_LRT, blind=FALSE)
if (Number_of_variables_to_compare == 1) {
  if (Genders_in_Data == "No") {
    z <- plotPCA(rld, intgroup=c('variableA'))
  }
}
if (Number_of_variables_to_compare == 1) {
  if (Genders_in_Data == "Yes") {
    z <- plotPCA(rld, intgroup=c('variableA','gender'))
  }
}


if (Number_of_variables_to_compare == 2) {
  if (Genders_in_Data == "No") {
    z <- plotPCA(rld, intgroup=c('variableA','variableB'))
  }
}
if (Number_of_variables_to_compare == 2) {
  if (Genders_in_Data == "Yes") {
    z <- plotPCA(rld, intgroup=c('variableA','variableB','gender'))
  }
}

if (Number_of_variables_to_compare == 3) {
  if (Genders_in_Data == "No") {
    z <- plotPCA(rld, intgroup=c('variableA','variableB','variableC'))
  }
}
if (Number_of_variables_to_compare == 3) {
  if (Genders_in_Data == "Yes") {
    z <- plotPCA(rld, intgroup=c('variableA','variableB','variableC','gender'))  }
}

dev.off()
pdf(file = paste("plots/PCAplot_",filenames2,"_comparison","_labeled_Positions.pdf"), pointsize = 10, width =8, height =8)
z + geom_label(aes(label = sample_data$position))
dev.off()

dev.off()
pdf(file = paste("plots/PCAplot_",filenames2,"_comparison.pdf"), pointsize = 10, width =8, height =8)
z 
dev.off()

