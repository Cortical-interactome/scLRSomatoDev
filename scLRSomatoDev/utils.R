#######################
#Load libraries
######################
library(circlize)
library(reshape2)
library(ComplexHeatmap)
library(scales)
library(graphics)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(plotly)
library(data.table)
library(shiny)
library(DT)
library(esquisse)
library(shinyWidgets)
library(scrattch.io)
library(bslib)
library(heatmaply)
library(scrattch.vis)
library(rhdf5)
library(HDF5Array)
library(SummarizedExperiment)
library(paletteer)
library(scico)
library(grDevices)
library(ggthemes)
library(cowplot)
library(stringr)
library(shinycssloaders)




#######################################################
# Function to read H5 seu obj in Summarized experiment
#######################################################

read_from_h5_seu <- function(h5_file) {
  require(HDF5Array)
  require(SummarizedExperiment)
  
  
  meta <- h5read(h5_file, "/metadata")
  rownames(meta) <- meta$sample_name
  x <- HDF5Array(h5_file,"/counts",type="integer")
  colnames(x) <- as.vector(h5read(h5_file,"/sample_names"))
  rownames(x) <- as.vector(h5read(h5_file,"/genes"))
  e <- SummarizedExperiment(
    SimpleList(counts=x),
    colData = meta[colnames(x),,drop=FALSE],
  )
}


#####################################################
# Function read h5 landscape
####################################################

read_from_h5_DataFrameList <- function(h5_file, group_read="P0_wt/comm_res_sig") {
  
  require(rhdf5)
  require(SummarizedExperiment)
  
 meta <-h5read(h5_file,group_read)
  EdgesSE <- SummarizedExperiment(
    colData = meta
  )
}


###################
#Read Table
##################



order_celltype <- read.csv("Data/order_celltype.csv")
rownames(order_celltype) <- order_celltype$celltype_label




names_sig.list <- as.list(c("E18.5-P0","P1-P2","P4-P5","P8","P16","P30","Adult","Cortex"))
names(names_sig.list) <- c("E18.5-P0","P1-P2","P4-P5","P8","P16","P30","Adult","Cortex")

Cortex_edges_sig.list <- lapply(names_sig.list, function( name){
  
  edges.df <- as.data.table(colData(loadHDF5SummarizedExperiment( "Data/", paste0(name,"_comm_res_sig_"))))
  
  return(edges.df)
  
})






#load LRDB
LRintercellNetworkDB <- read.csv("Data/LRintercellNetworkDB.csv") 


#Load Gene Gen e transcriptional landscape
Gene_IT.TL <- readRDS("Data/Gene_IT.landscape.rds")
Gene_ET.TL <- readRDS("Data/Gene_ET.landscape.rds")
Gene_Sst.TL <- readRDS("Data/Gene_Sst.landscape.rds")
Gene_Pvalb.TL <- readRDS("Data/Gene_Pvalb.landscape.rds")
Gene_Vip.TL <- readRDS("Data/Gene_Vip.landscape.rds")
Gene_Lamp5.TL <- readRDS("Data/Gene_Lamp5.landscape.rds")
Gene_OtherGlutNs.TL <- readRDS("Data/Gene_Other GlutNs.landscape.rds")
Gene_OtherGABANs.TL <- readRDS("Data/Gene_Other GABANs.landscape.rds")



Gene_TL <- unique(c(Gene_IT.TL, Gene_ET.TL, Gene_Sst.TL, Gene_Pvalb.TL,Gene_Vip.TL,Gene_Lamp5.TL,Gene_OtherGlutNs.TL,Gene_OtherGABANs.TL))


#Load thr sig gene name list
thr_name.family.list <- readRDS("Data/thr_name.family.list.rds")

########################
# Load object
########################
GECortex_PostMsub <-loadHDF5SummarizedExperiment("Data", "GECortex_PostMsub_")


#Load pM object
Res_pt.list <- readRDS("Data/Res_pt.list.rds")

#Load pL object
Res_pL.list <- readRDS("Data/Res_pL.list.rds")








#########################
# Theme use for ggplot
#########################


mytheme_classic <- theme_classic()+ 
  theme(axis.text = element_text(size=14, color="black"), 
        axis.ticks = element_line(size=1), 
        axis.title = element_text(size=18, color="black"),
        axis.line = element_line(size=1),
        axis.ticks.length= unit(7, "points"),   
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))

mytheme_bw <- theme_bw()+ 
  theme(axis.text = element_text(size=14, color="black"), 
        axis.title = element_text(size=18, color="black"),
        axis.ticks.length= unit(3.5, "points"),   
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))

mytheme_gray <-theme_gray()+
  theme(axis.text = element_text(size=14, color="black"), 
        axis.title = element_text(size=18, color="black"),
        axis.ticks.length= unit(7, "points"),   
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14),
        panel.grid = element_line(colour = NA),
        panel.background = element_rect(fill = alpha("black",0.2), 
                                        colour = "black"),
        strip.background = element_rect(fill = "white", 
                                        colour = "black"))

mytheme_minimal <- theme_minimal()+ 
  theme(axis.text = element_text(size=14, color="black"), 
        axis.ticks.length= unit(7, "points"),   
        legend.text = element_text(size=12), 
        legend.title = element_text(size=14))

##################################
# Plot LRintercellNetworkDB family
##################################


plot_LigLRDB_family=function(LRDB= LRintercellNetworkDB){
  
  require(ggplot2)
  require(paletteer)
  require(dplyr)
  require(plotly)

  
  LRDB=LRDB[,c(1,7)] %>% unique()
  
  ######### Create data.frame
  family.df=as.data.frame(matrix(ncol=2,nrow=length(sort(unique(unlist(strsplit(LRDB$ligand_family, split=";")))))))
  colnames(family.df)= c("family", "Freq")
  family.df$family <- sort(unique(unlist(strsplit(LRDB$ligand_family, split=";"))))
  
  for (i in 1:nrow(family.df)) {
    family.df$Freq[i] = nrow(LRDB[unlist(sapply(family.df$family[i], grep, LRDB$ligand_family, fixed=T, USE.NAMES = F)),])
  }
  
  family.df= family.df[order(family.df$Freq),]
  family.df$family= factor(family.df$family, levels=family.df$family)
  
  order_color= colorRampPalette(as.character(paletteer_d("ggsci::default_igv")))(nrow(family.df))
  


  
  
  ############# Plot
  p <- ggplot(data=family.df, aes(x=family, y=Freq))+
    geom_col(aes(fill=family,
                   text = paste('Ligand family:', family,
                                '<br> Number of ligands:', Freq)), position="identity", color=NA)+
    scale_fill_manual(values = order_color, guide=guide_legend(title="Ligand family", order=1 ,override.aes = list(size=5)))+
    scale_x_discrete(position = "top", guide = guide_axis(angle=90))+
    xlab(NULL) + ylab("Number of ligands")+
    coord_flip()+
    mytheme_bw
  
  ########### Interactive plot
  plot_height= 950 + 50*length(unique(family.df$family))

  
  gp <-  ggplotly(p, tooltip=c( "text"), height=plot_height) 
  

  
  gp <-gp %>% layout(showlegend=F,
                     annotations = list(
                       visible = F),
                     legend= list(itemsizing='constant', y=1), margin=list(t=90), xaxis = list(side ="top" ))
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  return(gp)

  
  
}

plot_RecLRDB_family=function(LRDB= LRintercellNetworkDB){
  
  require(ggplot2)
  require(paletteer)
  require(dplyr)
  require(plotly)
  
  
  LRDB=LRDB[,c(2,8)] %>% unique()
  
  ######### Create data.frame
  family.df=as.data.frame(matrix(ncol=2,nrow=length(sort(unique(unlist(strsplit(LRDB$receptor_family, split=";")))))))
  colnames(family.df)= c("family", "Freq")
  family.df$family <- sort(unique(unlist(strsplit(LRDB$receptor_family, split=";"))))
  
  for (i in 1:nrow(family.df)) {
    family.df$Freq[i] = nrow(LRDB[unlist(sapply(family.df$family[i], grep, LRDB$receptor_family, fixed=T, USE.NAMES = F)),])
  }
  
  family.df= family.df[order(family.df$Freq),]
  family.df$family= factor(family.df$family, levels=family.df$family)
  
  order_color= colorRampPalette(as.character(paletteer_d("ggsci::default_igv")))(nrow(family.df))
  
  
  
  
  
  ############# Plot
  p <- ggplot(data=family.df, aes(x=family, y=Freq))+
    geom_col(aes(fill=family,
                 text = paste('Receptor family:', family,
                              '<br> Number of receptors:', Freq)), position="identity", color=NA)+
    scale_fill_manual(values = order_color, guide=guide_legend(title="receptor family", order=1 ,override.aes = list(size=5)))+
    scale_x_discrete(position = "top", guide = guide_axis(angle=90))+
    xlab(NULL) + ylab("Number of receptors")+
    coord_flip()+
    mytheme_bw
  
  ########### Interactive plot
  plot_height= 950 + 50*length(unique(family.df$family))
  
  
  gp <-  ggplotly(p, tooltip=c( "text"), height=plot_height) 
  
  
  
  gp <-gp %>% layout(showlegend=F,
                     annotations = list(
                       visible = F),
                     legend= list(itemsizing='constant', y=1), margin=list(t=90), xaxis = list(side ="top" ))
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  
  return(gp)
  
  
  
}

#####################
#Plot Dim_red 3D
#####################


#Grouping level=cluster
plot_Dim.red=function(se_obj,grouping.level="celltype",dim.red.type="UMAP_3D", ident.class=c("Glutamatergic"), ident.family=c("IT") ,ident.subclass=c("Lamp5"), 
                      ident.supertype=c("Lamp5 Egln3"), ident.celltype=c("Lamp5|L1 A7C/CNC") , ident.celltype_original=c("Lamp5|Egln3|L1 A7C/CNC"),
                      ident.age=c("P0"), ident.study=c("This study"), ident.region=c("SSp;SSs"), ident.RNAseq.method=c("scRNA-seq"), ident.platform=c("10X"),
                      GeneList=c("Gad1","Neurod2"),height="950px"){

  
ident.class <- intersect(se_obj$class_label, ident.class) 
ident.family <- intersect(se_obj$family_label, ident.family) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 


######Subset class
se_obj <- subset(se_obj, select= colData(se_obj)$class_label %in% ident.class)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.family <- intersect(se_obj$family_label, ident.family) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

######Subset family
se_obj <- subset(se_obj, select= colData(se_obj)$family_label %in% ident.family)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.family <- intersect(se_obj$family_label, ident.family) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 


######Subset subclass
se_obj <- subset(se_obj, select= colData(se_obj)$subclass_label %in% ident.subclass)


ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset celltype
se_obj <- subset(se_obj, select= colData(se_obj)$celltype_label %in% ident.celltype)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset supertype
se_obj <- subset(se_obj, select= colData(se_obj)$supertype_label %in% ident.supertype)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset celltype_original
se_obj <- subset(se_obj, select= colData(se_obj)$celltype_original_label %in% ident.celltype_original)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset age
se_obj <- subset(se_obj, select= colData(se_obj)$age.at.collection.grp %in% ident.age)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset region
se_obj <- subset(se_obj, select= colData(se_obj)$region_label %in% ident.region)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset platform
se_obj <- subset(se_obj, select= colData(se_obj)$platform_label %in% ident.platform)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset RNAseq.method
se_obj <- subset(se_obj, select= colData(se_obj)$RNAseq.method_label %in% ident.RNAseq.method)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 

#Subset platform
se_obj <- subset(se_obj, select= colData(se_obj)$study_label %in% ident.study)

ident.class <- intersect(se_obj$class_label, ident.class) 
ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
ident.region <- intersect(se_obj$region_label, ident.region) 
ident.platform <- intersect(se_obj$platform_label, ident.platform) 
ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
ident.study <- intersect(se_obj$study_label, ident.study) 


#Preapare dataframe for factor
celltype.df <- as.data.frame(colData(se_obj)) %>% 
  dplyr::select(celltype_id, celltype_label, celltype_color) %>% 
  unique() %>% 
  dplyr::arrange(celltype_id)

celltype_original.df <- as.data.frame(colData(GECortex_PostMsub)) %>% 
                        dplyr::select(subclass_id,subclass_label,subclass_color,celltype_original_id,celltype_original_label,celltype_original_color) %>% 
                        unique() %>%  arrange(as.numeric(celltype_original_id)) %>% arrange(as.numeric(subclass_id))

supertype.df <- as.data.frame(colData(GECortex_PostMsub)) %>% 
                dplyr::select(subclass_id,subclass_label,subclass_color,supertype_id,supertype_label,supertype_color) %>% 
                unique()  %>%  arrange(as.numeric(supertype_id)) %>% arrange(as.numeric(subclass_id))
supertype.df <- supertype.df[!(supertype.df$supertype_label%in%c("Undetermined")),]
supertype.df <- supertype.df %>% rbind(c(NA,NA,NA,500,"Undetermined","#BEBEBE"))

subclass.df <- as.data.frame(colData(se_obj)) %>% 
  dplyr::select(subclass_id, subclass_label, subclass_color) %>% 
  unique() %>% 
  dplyr::arrange(subclass_id)

family.df <- as.data.frame(colData(se_obj)) %>% 
  dplyr::select(family_id, family_label, family_color) %>% 
  unique() %>% 
  dplyr::arrange(family_id)

class.df <- as.data.frame(colData(se_obj)) %>% 
  dplyr::select(class_id, class_label, class_color) %>% 
  unique() %>% 
  dplyr::arrange(class_id)


#Prepare data for 3D
Dim.data_filter_3D <- as.data.frame(colData(se_obj)[,c("UMAP_1_3D", "UMAP_2_3D", "UMAP_3_3D", "tSNE_1_3D", "tSNE_2_3D", "tSNE_3_3D", 
                                                      "PC_1", "PC_2","PC_3", "IC_1", "IC_2", "IC_3",
                                                      "celltype_original_id","celltype_original_label","celltype_original_color",
                                                      "celltype_id","celltype_label", "celltype_color",
                                                      "supertype_id","supertype_label", "supertype_color",
                                                      "subclass_id","subclass_label","subclass_color",
                                                      "family_id","family_label", "family_color",
                                                      "class_id","class_label", "class_color",
                                                      "region_id","region_label", "region_color",
                                                      "study_id","study_label", "study_color",
                                                      "platform_id","platform_label", "platform_color",
                                                      "RNAseq.method_id","RNAseq.method_label", "RNAseq.method_color",
                                                      "age.at.collection.grp","nFeature_RNA", "nCount_RNA")])

colnames(Dim.data_filter_3D )[1:6] <- c("UMAP_1", "UMAP_2", "UMAP_3", "tSNE_1", "tSNE_2", "tSNE_3")

Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$platform_id),]
Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$RNAseq.method_id),]
Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$region_id),]
Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$study_id),]


Dim.data_filter_3D$platform_label <- factor(Dim.data_filter_3D$platform_label, levels=unique(Dim.data_filter_3D$platform_label))
Dim.data_filter_3D$RNAseq.method_label <- factor(Dim.data_filter_3D$RNAseq.method_label, levels=unique(Dim.data_filter_3D$RNAseq.method_label))
Dim.data_filter_3D$region_label <- factor(Dim.data_filter_3D$region_label, levels=unique(Dim.data_filter_3D$region_label))
Dim.data_filter_3D$study_label <- factor(Dim.data_filter_3D$study_label, levels=unique(Dim.data_filter_3D$study_label))

Dim.data_filter_3D$celltype_original_label <- factor(Dim.data_filter_3D$celltype_original_label, levels=celltype_original.df$celltype_original_label)
Dim.data_filter_3D$supertype_label <- factor(Dim.data_filter_3D$supertype_label, levels=supertype.df$supertype_label)
Dim.data_filter_3D$celltype_label <- factor(Dim.data_filter_3D$celltype_label, levels=celltype.df$celltype_label)
Dim.data_filter_3D$subclass_label <- factor(Dim.data_filter_3D$subclass_label, levels=subclass.df$subclass_label)
Dim.data_filter_3D$family_label <- factor(Dim.data_filter_3D$family_label, levels=family.df$family_label)
Dim.data_filter_3D$class_label <- factor(Dim.data_filter_3D$class_label, levels=class.df$class_label)
Dim.data_filter_3D$age.at.collection.grp <- factor(Dim.data_filter_3D$age.at.collection.grp, levels=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5",
                                                                                                      "E17.5", "E18.5", "P0","P1","P2","P4","P5","P8","P16",
                                                                                                      "P30","Adult"))

#Prepare data for 2D


Dim.data_filter_2D <- as.data.frame(colData(se_obj)[,c("UMAP_1_2D", "UMAP_2_2D",  "tSNE_1_2D", "tSNE_2_2D",
                                                       "PC_1", "PC_2", "IC_1", "IC_2",
                                                       "celltype_original_id","celltype_original_label","celltype_original_color",
                                                      "celltype_id","celltype_label", "celltype_color",
                                                      "supertype_id","supertype_label", "supertype_color",
                                                      "subclass_id","subclass_label","subclass_color",
                                                      "family_id","family_label", "family_color",
                                                      "class_id","class_label", "class_color",
                                                      "region_id","region_label", "region_color",
                                                      "study_id","study_label", "study_color",
                                                      "platform_id","platform_label", "platform_color",
                                                      "RNAseq.method_id","RNAseq.method_label", "RNAseq.method_color",
                                                      "age.at.collection.grp","nFeature_RNA", "nCount_RNA")])

colnames(Dim.data_filter_2D )[1:4] <- c("UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2")

Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$platform_id),]
Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$RNAseq.method_id),]
Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$region_id),]
Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$study_id),]


Dim.data_filter_2D$platform_label <- factor(Dim.data_filter_2D$platform_label, levels=unique(Dim.data_filter_2D$platform_label))
Dim.data_filter_2D$RNAseq.method_label <- factor(Dim.data_filter_2D$RNAseq.method_label, levels=unique(Dim.data_filter_2D$RNAseq.method_label))
Dim.data_filter_2D$region_label <- factor(Dim.data_filter_2D$region_label, levels=unique(Dim.data_filter_2D$region_label))
Dim.data_filter_2D$study_label <- factor(Dim.data_filter_2D$study_label, levels=unique(Dim.data_filter_2D$study_label))

Dim.data_filter_2D$celltype_original_label <- factor(Dim.data_filter_2D$celltype_original_label, levels=celltype_original.df$celltype_original_label)
Dim.data_filter_2D$supertype_label <- factor(Dim.data_filter_2D$supertype_label, levels=supertype.df$supertype_label)
Dim.data_filter_2D$celltype_label <- factor(Dim.data_filter_2D$celltype_label, levels=celltype.df$celltype_label)
Dim.data_filter_2D$subclass_label <- factor(Dim.data_filter_2D$subclass_label, levels=subclass.df$subclass_label)
Dim.data_filter_2D$family_label <- factor(Dim.data_filter_2D$family_label, levels=family.df$family_label)
Dim.data_filter_2D$class_label <- factor(Dim.data_filter_2D$class_label, levels=class.df$class_label)
Dim.data_filter_2D$age.at.collection.grp <- factor(Dim.data_filter_2D$age.at.collection.grp, levels=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5",
                                                                                                      "E17.5", "E18.5", "P0","P1","P2","P4","P5","P8","P16",
                                                                                                      "P30","Adult"))

#Plot

if(grouping.level=="original cell-type"){

if(dim.red.type=="UMAP_3D"){
  
gp <-  plot_ly(data = Dim.data_filter_3D, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~celltype_original_label, 
        colors = celltype_original.df$celltype_original_color,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 2, opacity=0.8), # controls size of points
        text=~paste('Cell:', rownames(Dim.data_filter_3D),
                    '<br>Class:', class_label,
                    '<br>Family:', family_label,
                    '<br>Subclass:', subclass_label,
                    '<br>Supertype:', supertype_label,
                    '<br>Cell-type:', celltype_label,
                    '<br>Original celltype:', celltype_original_label,
                    '<br>Age:', age.at.collection.grp,
                    '<br>Study:', study_label,
                    '<br>Region:', region_label,
                    '<br>RNA-seq method:', RNAseq.method_label,
                    '<br>Platform:', platform_label,
                    '<br>Number of genes:', nFeature_RNA,
                    '<br>Number of UMIs:', nCount_RNA), 
        hoverinfo="text", height=height) %>% 
  layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))

}


else if(dim.red.type=="tSNE_3D"){
  
  gp <- plot_ly(data = Dim.data_filter_3D, 
                x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 2, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_3D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA),  
                hoverinfo="text",height=height) %>% 
    layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
  
}

else if(dim.red.type=="PCA_3D"){
  
  gp <- plot_ly(data = Dim.data_filter_3D, 
                x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 2, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_3D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA),  
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
  
}

else if(dim.red.type=="ICA_3D"){
  
  gp <- plot_ly(data = Dim.data_filter_3D, 
                x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter3d", 
                mode = "markers", 
                marker = list(size = 2, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_3D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA), 
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
  
}

else if(dim.red.type=="UMAP_2D"){
  
  gp <- plot_ly(data = Dim.data_filter_2D, 
                x = ~UMAP_1, y = ~UMAP_2,
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter", 
                mode = "markers", 
                marker = list(size = 4, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_2D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA), 
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
  
}

else if(dim.red.type=="tSNE_2D"){
  
  gp <- plot_ly(data = Dim.data_filter_2D, 
                x = ~tSNE_1, y = ~tSNE_2,
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter", 
                mode = "markers", 
                marker = list(size = 4, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_2D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA),
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
  
}

else if(dim.red.type=="PCA_2D"){
  
  gp <- plot_ly(data = Dim.data_filter_2D, 
                x = ~PC_1, y = ~PC_2,
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter", 
                mode = "markers", 
                marker = list(size = 4, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_2D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA), 
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
  
}

else if(dim.red.type=="ICA_2D"){
  
  gp <- plot_ly(data = Dim.data_filter_2D, 
                x = ~IC_1, y = ~IC_2,
                color = ~celltype_original_label, 
                colors = celltype_original.df$celltype_original_color,
                type = "scatter", 
                mode = "markers", 
                marker = list(size = 4, opacity=0.8), # controls size of points
                text=~paste('Cell:', rownames(Dim.data_filter_2D),
                            '<br>Class:', class_label,
                            '<br>Family:', family_label,
                            '<br>Subclass:', subclass_label,
                            '<br>Supertype:', supertype_label,
                            '<br>Cell-type:', celltype_label,
                            '<br>Original celltype:', celltype_original_label,
                            '<br>Age:', age.at.collection.grp,
                            '<br>Study:', study_label,
                            '<br>Region:', region_label,
                            '<br>RNA-seq method:', RNAseq.method_label,
                            '<br>Platform:', platform_label,
                            '<br>Number of genes:', nFeature_RNA,
                            '<br>Number of UMIs:', nCount_RNA),
                hoverinfo="text", height=height) %>% 
    layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
  
}
  
}

else if(grouping.level=="cell-type"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~celltype_label, 
                   colors = celltype.df$celltype_color,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA),
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~celltype_label, 
                  colors = celltype.df$celltype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}

else if(grouping.level=="supertype"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~supertype_label, 
                   colors = supertype.df$supertype_color,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~supertype_label, 
                  colors = supertype.df$supertype_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}



else if(grouping.level=="subclass"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~subclass_label, 
                   colors = subclass.df$subclass_color,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~subclass_label, 
                  colors = subclass.df$subclass_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}



else if(grouping.level=="family"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~family_label, 
                   colors = family.df$family_color,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~family_label, 
                  colors = family.df$family_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}





else if(grouping.level=="class"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~class_label, 
                   colors = class.df$class_color,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA),  
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),   
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),   
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~class_label, 
                  colors = class.df$class_color,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}


else if(grouping.level=="age"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~age.at.collection.grp, 
                   colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA),   
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),  
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~age.at.collection.grp, 
                  colors = colorRampPalette(paletteer_d("grDevices::blues9"))(length(unique(Dim.data_filter_3D$age.at.collection.grp))),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA),
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}

else if(grouping.level=="study"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~study_label, 
                   colors = unique( Dim.data_filter_3D$study_color),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~study_label, 
                  colors = unique( Dim.data_filter_3D$study_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}

else if(grouping.level=="region"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~region_label, 
                   colors = unique( Dim.data_filter_3D$region_color),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~region_label, 
                  colors = unique( Dim.data_filter_3D$region_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}

else if(grouping.level=="RNA-seq method"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~RNAseq.method_label, 
                   colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~RNAseq.method_label, 
                  colors = unique( Dim.data_filter_3D$RNAseq.method_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}

else if(grouping.level=="platform"){
  
  if(dim.red.type=="UMAP_3D"){
    
    gp <-  plot_ly(data = Dim.data_filter_3D, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~platform_label, 
                   colors = unique( Dim.data_filter_3D$platform_color),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, opacity=0.8), # controls size of points
                   text=~paste('Cell:', rownames(Dim.data_filter_3D),
                               '<br>Class:', class_label,
                               '<br>Family:', family_label,
                               '<br>Subclass:', subclass_label,
                               '<br>Supertype:', supertype_label,
                               '<br>Cell-type:', celltype_label,
                               '<br>Original celltype:', celltype_original_label,
                               '<br>Age:', age.at.collection.grp,
                               '<br>Study:', study_label,
                               '<br>Region:', region_label,
                               '<br>RNA-seq method:', RNAseq.method_label,
                               '<br>Platform:', platform_label,
                               '<br>Number of genes:', nFeature_RNA,
                               '<br>Number of UMIs:', nCount_RNA), 
                   hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  
  else if(dim.red.type=="tSNE_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text",height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="PCA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube"))
    
  }
  
  else if(dim.red.type=="ICA_3D"){
    
    gp <- plot_ly(data = Dim.data_filter_3D, 
                  x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter3d", 
                  mode = "markers", 
                  marker = list(size = 2, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_3D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"),scene=list(aspectmode="cube") )
    
  }
  
  else if(dim.red.type=="UMAP_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~UMAP_1, y = ~UMAP_2,
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="tSNE_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~tSNE_1, y = ~tSNE_2,
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="PCA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~PC_1, y = ~PC_2,
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
  else if(dim.red.type=="ICA_2D"){
    
    gp <- plot_ly(data = Dim.data_filter_2D, 
                  x = ~IC_1, y = ~IC_2,
                  color = ~platform_label, 
                  colors = unique( Dim.data_filter_3D$platform_color),
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(size = 4, opacity=0.8), # controls size of points
                  text=~paste('Cell:', rownames(Dim.data_filter_2D),
                              '<br>Class:', class_label,
                              '<br>Family:', family_label,
                              '<br>Subclass:', subclass_label,
                              '<br>Supertype:', supertype_label,
                              '<br>Cell-type:', celltype_label,
                              '<br>Original celltype:', celltype_original_label,
                              '<br>Age:', age.at.collection.grp,
                              '<br>Study:', study_label,
                              '<br>Region:', region_label,
                              '<br>RNA-seq method:', RNAseq.method_label,
                              '<br>Platform:', platform_label,
                              '<br>Number of genes:', nFeature_RNA,
                              '<br>Number of UMIs:', nCount_RNA), 
                  hoverinfo="text", height=height) %>% 
      layout(legend = list(itemsizing="constant"), xaxis=list(zeroline=F), yaxis=list(zeroline=F))
    
  }
  
}




gp <- gp %>% config(
  toImageButtonOptions = list(
    format = "svg"
  )
)

return(gp)



}





###########################################
# Plot Dimension reduction Gene Expression
###########################################



plot_Dim.red_ExprGene=function(se_obj, tome_file,grouping.level="celltype",dim.red.type="UMAP_3D", ident.class=c("Glutamatergic"), ident.subclass=c("Lamp5"), 
                               ident.family=c("Lamp5"),ident.supertype=c("Lamp5 Egln3"), ident.celltype=c("Lamp5|L1 A7C/CNC") , ident.celltype_original=c("Lamp5|Egln3|L1 A7C/CNC"),
                      ident.age=c("P0"), ident.study=c("This study"), ident.region=c("SSp;SSs"), ident.RNAseq.method=c("scRNA-seq"), ident.platform=c("10X"),
                       GeneList=c("Gad1","Neurod2"), height="850px", color_gdt1="Reds", 
                      color_gdt2="Greens"){
  
  #require(Seurat)
  require(plotly)
  
  ########## Color for gene1
  if(color_gdt1=="Reds"){
    col1 =  brewer_pal(palette = "Reds")(8)
  }
  else if(color_gdt1=="Greens"){
    col1 =  brewer_pal(palette = "Greens")(8)
  }
  
  else if(color_gdt1=="Blues"){
    col1 = brewer_pal(palette = "Blues")(8)
  }
  else if(color_gdt1=="Purples"){
    col1 =brewer_pal(palette = "Purples")(8)
  }
  else if(color_gdt1=="Oranges"){
    col1 =brewer_pal(palette = "Oranges")(8)
  }
  else if(color_gdt1=="OrRd"){
    col1 =brewer_pal(palette = "OrRd")(8)
  }
  else if(color_gdt1=="YlOrRd"){
    col1 =brewer_pal(palette = "YlOrRd")(8)
  }
  else if(color_gdt1=="YlOrBr"){
    col1 =brewer_pal(palette = "YlOrBr")(8)
  }
  else if(color_gdt1=="BuYlRd"){
    col1 =brewer_pal(palette = "RdYlBu", direction = -1)(8)
  }

  ########## Color for gene2
  if(color_gdt2=="Reds"){
    col2 =  brewer_pal(palette = "Reds")(8)
  }
  else if(color_gdt2=="Greens"){
    col2 =  brewer_pal(palette = "Greens")(8)
  }
  
  else if(color_gdt2=="Blues"){
    col2 = brewer_pal(palette = "Blues")(8)
  }
  else if(color_gdt2=="Purples"){
    col2 =brewer_pal(palette = "Purples")(8)
  }
  else if(color_gdt2=="Oranges"){
    col2 =brewer_pal(palette = "Oranges")(8)
  }
  else if(color_gdt2=="OrRd"){
    col2 =brewer_pal(palette = "OrRd")(8)
  }
  else if(color_gdt2=="YlOrRd"){
    col2 =brewer_pal(palette = "YlOrRd")(8)
  }
  else if(color_gdt2=="YlOrBr"){
    col2 =brewer_pal(palette = "YlOrBr")(8)
  }
  else if(color_gdt2=="BuYlRd"){
    col2 =brewer_pal(palette = "RdYlBu", direction = -1)(8)
  }
  
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ######Subset class
  se_obj <- subset(se_obj, select= colData(se_obj)$class_label %in% ident.class)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ######Subset family
  se_obj <- subset(se_obj, select= colData(se_obj)$family_label %in% ident.family)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  ######Subset subclass
  se_obj <- subset(se_obj, select= colData(se_obj)$subclass_label %in% ident.subclass)
  
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset celltype
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_label %in% ident.celltype)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset supertype
  se_obj <- subset(se_obj, select= colData(se_obj)$supertype_label %in% ident.supertype)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset celltype_original
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_original_label %in% ident.celltype_original)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset age
  se_obj <- subset(se_obj, select= colData(se_obj)$age.at.collection.grp %in% ident.age)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset region
  se_obj <- subset(se_obj, select= colData(se_obj)$region_label %in% ident.region)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$platform_label %in% ident.platform)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset RNAseq.method
  se_obj <- subset(se_obj, select= colData(se_obj)$RNAseq.method_label %in% ident.RNAseq.method)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$study_label %in% ident.study)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.supertype <- intersect(se_obj$supertype_label, ident.supertype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  #Preapare dataframe for factor
  
  #Preapare dataframe for factor
  celltype.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(celltype_id, celltype_label, celltype_color) %>% 
    unique() %>% 
    dplyr::arrange(celltype_id)
  
  celltype_original.df <- as.data.frame(colData(GECortex_PostMsub)) %>% 
    dplyr::select(subclass_id,subclass_label,subclass_color,celltype_original_id,celltype_original_label,celltype_original_color) %>% 
    unique() %>%  arrange(as.numeric(celltype_original_id)) %>% arrange(as.numeric(subclass_id))
  
  supertype.df <- as.data.frame(colData(GECortex_PostMsub)) %>% 
    dplyr::select(subclass_id,subclass_label,subclass_color,supertype_id,supertype_label,supertype_color) %>% 
    unique()  %>%  arrange(as.numeric(supertype_id)) %>% arrange(as.numeric(subclass_id))
  supertype.df <- supertype.df[!(supertype.df$supertype_label%in%c("Undetermined")),]
  supertype.df <- supertype.df %>% rbind(c(NA,NA,NA,500,"Undetermined","#BEBEBE"))
  
  subclass.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(subclass_id, subclass_label, subclass_color) %>% 
    unique() %>% 
    dplyr::arrange(subclass_id)
  
  family.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(family_id, family_label, family_color) %>% 
    unique() %>% 
    dplyr::arrange(family_id)
  
  class.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(class_id, class_label, class_color) %>% 
    unique() %>% 
    dplyr::arrange(class_id)
  
  
  
  #Load Expr data
  GeneList<- GeneList[which(GeneList%in%rownames(se_obj))]
  
  norm.se_obj <- read_tome_gene_data(tome_file, genes = GeneList, format = "matrix")
  norm.se_obj <- norm.se_obj[rownames(norm.se_obj)%in%colnames(se_obj),,drop=F]
  norm.se_obj <- mat_to_data_df(norm.se_obj)
  
  
  #Prepare data for 3D
  Dim.data_filter_3D <- as.data.frame(colData(se_obj)[,c("UMAP_1_3D", "UMAP_2_3D", "UMAP_3_3D", "tSNE_1_3D", "tSNE_2_3D", "tSNE_3_3D", 
                                                         "PC_1", "PC_2","PC_3", "IC_1", "IC_2", "IC_3",
                                                         "celltype_original_id","celltype_original_label","celltype_original_color",
                                                         "celltype_id","celltype_label", "celltype_color",
                                                         "supertype_id","supertype_label", "supertype_color",
                                                         "subclass_id","subclass_label","subclass_color",
                                                         "class_id","class_label", "class_color",
                                                         "family_id","family_label", "family_color",
                                                         "region_id","region_label", "region_color",
                                                         "study_id","study_label", "study_color",
                                                         "platform_id","platform_label", "platform_color",
                                                         "RNAseq.method_id","RNAseq.method_label", "RNAseq.method_color",
                                                         "age.at.collection.grp")]) %>% cbind(norm.se_obj)
  
  
  colnames(Dim.data_filter_3D )[1:6] <- c("UMAP_1", "UMAP_2", "UMAP_3", "tSNE_1", "tSNE_2", "tSNE_3")
  
  Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$platform_id),]
  Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$RNAseq.method_id),]
  Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$region_id),]
  Dim.data_filter_3D <- Dim.data_filter_3D[order(Dim.data_filter_3D$study_id),]

  
  Dim.data_filter_3D$platform_label <- factor(Dim.data_filter_3D$platform_label, levels=unique(Dim.data_filter_3D$platform_label))
  Dim.data_filter_3D$RNAseq.method_label <- factor(Dim.data_filter_3D$RNAseq.method_label, levels=unique(Dim.data_filter_3D$RNAseq.method_label))
  Dim.data_filter_3D$region_label <- factor(Dim.data_filter_3D$region_label, levels=unique(Dim.data_filter_3D$region_label))
  Dim.data_filter_3D$study_label <- factor(Dim.data_filter_3D$study_label, levels=unique(Dim.data_filter_3D$study_label))

  Dim.data_filter_3D$celltype_original_label <- factor(Dim.data_filter_3D$celltype_original_label, levels=celltype_original.df$celltype_original_label)
  Dim.data_filter_3D$supertype_label <- factor(Dim.data_filter_3D$supertype_label, levels=supertype.df$supertype_label)
  Dim.data_filter_3D$celltype_label <- factor(Dim.data_filter_3D$celltype_label, levels=celltype.df$celltype_label)
  Dim.data_filter_3D$subclass_label <- factor(Dim.data_filter_3D$subclass_label, levels=subclass.df$subclass_label)
  Dim.data_filter_3D$family_label <- factor(Dim.data_filter_3D$family_label, levels=family.df$family_label)
  Dim.data_filter_3D$class_label <- factor(Dim.data_filter_3D$class_label, levels=class.df$class_label)
  Dim.data_filter_3D$age.at.collection.grp <- factor(Dim.data_filter_3D$age.at.collection.grp, levels=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5",
                                                                                                        "E17.5", "E18.5", "P0","P1","P2","P4","P5","P8","P16",
                                                                                                        "P30","Adult"))

  
  
  
  #Prepare data for 2D
  Dim.data_filter_2D <- as.data.frame(colData(se_obj)[,c("UMAP_1_2D", "UMAP_2_2D",  "tSNE_1_2D", "tSNE_2_2D", 
                                                         "PC_1", "PC_2", "IC_1", "IC_2", 
                                                         "celltype_original_id","celltype_original_label","celltype_original_color",
                                                         "celltype_id","celltype_label", "celltype_color",
                                                         "supertype_id","supertype_label", "supertype_color",
                                                         "subclass_id","subclass_label","subclass_color",
                                                         "family_id","family_label", "family_color",
                                                         "class_id","class_label", "class_color",
                                                         "region_id","region_label", "region_color",
                                                         "study_id","study_label", "study_color",
                                                         "platform_id","platform_label", "platform_color",
                                                         "RNAseq.method_id","RNAseq.method_label", "RNAseq.method_color",
                                                         "age.at.collection.grp")]) %>% cbind(norm.se_obj)
  
  colnames(Dim.data_filter_2D)[1:4] <- c("UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2")
  
  Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$platform_id),]
  Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$RNAseq.method_id),]
  Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$region_id),]
  Dim.data_filter_2D <- Dim.data_filter_2D[order(Dim.data_filter_2D$study_id),]

  
  Dim.data_filter_2D$platform_label <- factor(Dim.data_filter_2D$platform_label, levels=unique(Dim.data_filter_2D$platform_label))
  Dim.data_filter_2D$RNAseq.method_label <- factor(Dim.data_filter_2D$RNAseq.method_label, levels=unique(Dim.data_filter_2D$RNAseq.method_label))
  Dim.data_filter_2D$region_label <- factor(Dim.data_filter_2D$region_label, levels=unique(Dim.data_filter_2D$region_label))
  Dim.data_filter_2D$study_label <- factor(Dim.data_filter_2D$study_label, levels=unique(Dim.data_filter_2D$study_label))

  Dim.data_filter_2D$celltype_original_label <- factor(Dim.data_filter_2D$celltype_original_label, levels=celltype_original.df$celltype_original_label)
  Dim.data_filter_2D$supertype_label <- factor(Dim.data_filter_2D$supertype_label, levels=supertype.df$supertype_label)
  Dim.data_filter_2D$celltype_label <- factor(Dim.data_filter_2D$celltype_label, levels=celltype.df$celltype_label)
  Dim.data_filter_2D$subclass_label <- factor(Dim.data_filter_2D$subclass_label, levels=subclass.df$subclass_label)
  Dim.data_filter_2D$family_label <- factor(Dim.data_filter_2D$family_label, levels=family.df$family_label)
  Dim.data_filter_2D$class_label <- factor(Dim.data_filter_2D$class_label, levels=class.df$class_label)
  Dim.data_filter_2D$age.at.collection.grp <- factor(Dim.data_filter_2D$age.at.collection.grp, levels=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5",
                                                                                                        "E17.5", "E18.5", "P0","P1","P2","P4","P5","P8","P16",
                                                                                                        "P30","Adult"))

  
  #Plot
  genes=colnames(Dim.data_filter_3D)[45:ncol(Dim.data_filter_3D)]
  
  
  
  ###########UMAP_3D
    if(dim.red.type=="UMAP_3D"){
      
      if(ncol(Dim.data_filter_3D)==46){

      gp1 <- plot_ly(data = Dim.data_filter_3D, 
              x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
              color = ~Dim.data_filter_3D[,45], 
              colors = col1,
              type = "scatter3d", 
              mode = "markers", 
              marker = list(size = 2, opacity=0.8), # controls size of points
              text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                '<br>Class:', class_label,
                                '<br>Family:', family_label,
                                '<br>Subclass:', subclass_label,
                                '<br>Supertype:', supertype_label,
                                '<br>Cell-type:', celltype_label,
                                '<br>Original cell-type:', celltype_original_label,
                                '<br>Age:', age.at.collection.grp,
                                '<br>Study:', study_label,
                                '<br>Region:', region_label,
                                '<br>RNA-seq method:', RNAseq.method_label,
                                '<br>Platform:', platform_label,
                                '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
              hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
        colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                font=list(color="white")), outlinecolor="white",
                tickcolor="white",tickfont=list(color="white"),
                tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
        layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                  xref='paper', yref='paper', font=list(size=20, color="white")),
               paper_bgcolor=toRGB("black"))
      
      gp2 <- plot_ly(data = Dim.data_filter_3D, 
                     x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                     color = ~Dim.data_filter_3D[,46], 
                     colors = col2,
                     type = "scatter3d", 
                     mode = "markers", 
                     marker = list(size = 2, opacity=0.8), # controls size of points
                     text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                 '<br>Class:', class_label,
                                 '<br>Family:', family_label,
                                 '<br>Subclass:', subclass_label,
                                 '<br>Supertype:', supertype_label,
                                 '<br>Cell-type:', celltype_label,
                                 '<br>Original cell-type:', celltype_original_label,
                                 '<br>Age:', age.at.collection.grp,
                                 '<br>Study:', study_label,
                                 '<br>Region:', region_label,
                                 '<br>RNA-seq method:', RNAseq.method_label,
                                 '<br>Platform:', platform_label,
                                 '<br>CPM:', as.integer(2^Dim.data_filter_3D[,46]-1)), 
                     hoverinfo="text", showlegend=F,scene="scene2", height=height)%>% 
        colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                              font=list(color="white")), outlinecolor="white",
                 tickcolor="white",tickfont=list(color="white"),
                 tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)),
                 ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)))-1))) %>% 
        layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                  xref='paper', yref='paper', font=list(size=20, color="white")))

      
     gp <-  subplot(gp1, gp2)%>% 
       layout(scene=list(domain=list(x=c(0,0.48),y=c(0,1)),
                         aspectmode="cube", 
                         xaxis=list(color="white"),
                         yaxis=list(color="white"),
                         zaxis=list( color="white")),
              scene2=list(domain=list(x=c(0.52,1),y=c(0,1)),
                          aspectmode="cube",
                          xaxis=list(title="UMAP_1", color="white"),
                          yaxis=list(title="UMAP_2", color="white"),
                          zaxis=list(title="UMAP_3", color="white")))
     
      }
      
      else if(ncol(Dim.data_filter_3D)==45){
        gp <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                       color = ~Dim.data_filter_3D[,45], 
                       colors = col1,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"),scene=list(aspectmode="cube",
                 xaxis=list(title="UMAP_1", color="white"),
                 yaxis=list(title="UMAP_2", color="white"),
                 zaxis=list(title="UMAP_3", color="white")))
        
      }
      
      
    }
    
    ######### tSNE_3D
    else if(dim.red.type=="tSNE_3D"){
      
      
      if(ncol(Dim.data_filter_3D)==46){
        
        gp1 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                       color = ~Dim.data_filter_3D[,45], 
                       colors = col1,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"))
        
        gp2 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                       color = ~Dim.data_filter_3D[,46], 
                       colors = col2,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,46]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene2", height=height)%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2)%>% 
          layout(scene=list(domain=list(x=c(0,0.48),y=c(0,1)),
                            aspectmode="cube", 
                            xaxis=list(color="white"),
                            yaxis=list(color="white"),
                            zaxis=list( color="white")),
                 scene2=list(domain=list(x=c(0.52,1),y=c(0,1)),
                             aspectmode="cube",
                             xaxis=list(title="tSNE_1", color="white"),
                             yaxis=list(title="tSNE_2", color="white"),
                             zaxis=list(title="tSNE_3", color="white")))
        
      }
      
      else if(ncol(Dim.data_filter_3D)==45){
        gp <- plot_ly(data = Dim.data_filter_3D, 
                      x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
                      color = ~Dim.data_filter_3D[,45], 
                      colors = col1,
                      type = "scatter3d", 
                      mode = "markers", 
                      marker = list(size = 2, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                      hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"),scene=list(aspectmode="cube",
                                                         xaxis=list(title="tSNE_1", color="white"),
                                                         yaxis=list(title="tSNE_2", color="white"),
                                                         zaxis=list(title="tSNE_3", color="white")))
        
      }
      
    }
    
  
  ########## PCA_3D
    else if(dim.red.type=="PCA_3D"){
      
      if(ncol(Dim.data_filter_3D)==46){
        
        gp1 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                       color = ~Dim.data_filter_3D[,45], 
                       colors = col1,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"))
        
        gp2 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                       color = ~Dim.data_filter_3D[,46], 
                       colors = col2,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,46]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene2", height=height)%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2)%>% 
          layout(scene=list(domain=list(x=c(0,0.48),y=c(0,1)),
                            aspectmode="cube", 
                            xaxis=list(color="white"),
                            yaxis=list(color="white"),
                            zaxis=list( color="white")),
                 scene2=list(domain=list(x=c(0.52,1),y=c(0,1)),
                             aspectmode="cube",
                             xaxis=list(title="PC_1", color="white"),
                             yaxis=list(title="PC_2", color="white"),
                             zaxis=list(title="PC_3", color="white")))
        
      }
      
      else if(ncol(Dim.data_filter_3D)==45){
        gp <- plot_ly(data = Dim.data_filter_3D, 
                      x = ~PC_1, y = ~PC_2, z = ~PC_3, 
                      color = ~Dim.data_filter_3D[,45], 
                      colors = col1,
                      type = "scatter3d", 
                      mode = "markers", 
                      marker = list(size = 2, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                      hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"),scene=list(aspectmode="cube",
                                                         xaxis=list(title="PC_1", color="white"),
                                                         yaxis=list(title="PC_2", color="white"),
                                                         zaxis=list(title="PC_3", color="white")))
        
      }
      
    }
    
  
  ######## ICA_3D
    else if(dim.red.type=="ICA_3D"){
      
      if(ncol(Dim.data_filter_3D)==46){
        
        gp1 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                       color = ~Dim.data_filter_3D[,45], 
                       colors = col1,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"))
        
        gp2 <- plot_ly(data = Dim.data_filter_3D, 
                       x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                       color = ~Dim.data_filter_3D[,46], 
                       colors = col2,
                       type = "scatter3d", 
                       mode = "markers", 
                       marker = list(size = 2, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_3D[,46]-1)), 
                       hoverinfo="text", showlegend=F,scene="scene2", height=height)%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,46]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2)%>% 
          layout(scene=list(domain=list(x=c(0,0.48),y=c(0,1)),
                            aspectmode="cube", 
                            xaxis=list(color="white"),
                            yaxis=list(color="white"),
                            zaxis=list( color="white")),
                 scene2=list(domain=list(x=c(0.52,1),y=c(0,1)),
                             aspectmode="cube",
                             xaxis=list(title="IC_1", color="white"),
                             yaxis=list(title="IC_2", color="white"),
                             zaxis=list(title="IC_3", color="white")))
        
      }
      
      else if(ncol(Dim.data_filter_3D)==45){
        gp <- plot_ly(data = Dim.data_filter_3D, 
                      x = ~IC_1, y = ~IC_2, z = ~IC_3, 
                      color = ~Dim.data_filter_3D[,45], 
                      colors = col1,
                      type = "scatter3d", 
                      mode = "markers", 
                      marker = list(size = 2, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_3D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_3D[,45]-1)), 
                      hoverinfo="text", showlegend=F,scene="scene1", height=height)%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_3D[,45]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"),scene=list(aspectmode="cube",
                                                         xaxis=list(title="IC_1", color="white"),
                                                         yaxis=list(title="IC_2", color="white"),
                                                         zaxis=list(title="IC_3", color="white")))
        
      }
      
    }
    
  ######### UMAP_2D
    else if(dim.red.type=="UMAP_2D"){
      
      if(ncol(Dim.data_filter_2D)==41){
      gp1 <- plot_ly(data = Dim.data_filter_2D, 
                     x = ~UMAP_1, y = ~UMAP_2,
                     color = ~Dim.data_filter_2D[,40], 
                     colors = col1,
                     type = "scatter", 
                     mode = "markers", 
                     marker = list(size = 4, opacity=0.8), # controls size of points
                     text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                 '<br>Class:', class_label,
                                 '<br>Family:', family_label,
                                 '<br>Subclass:', subclass_label,
                                 '<br>Supertype:', supertype_label,
                                 '<br>Cell-type:', celltype_label,
                                 '<br>Original cell-type:', celltype_original_label,
                                 '<br>Age:', age.at.collection.grp,
                                 '<br>Study:', study_label,
                                 '<br>Region:', region_label,
                                 '<br>RNA-seq method:', RNAseq.method_label,
                                 '<br>Platform:', platform_label,
                                 '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                     hoverinfo="text", showlegend=F, height=height)%>% 
        layout(xaxis=list(title="UMAP_1", color="white", zeroline=F),
               yaxis=list(title="UMAP_2",color="white", zeroline=F))%>% 
        colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                              font=list(color="white")), outlinecolor="white",
                 tickcolor="white",tickfont=list(color="white"),
                 tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                 ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
        layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                  xref='paper', yref='paper', font=list(size=20, color="white")),
               paper_bgcolor=toRGB("black"), plot_bgcolor="black")
      
      gp2 <- plot_ly(data = Dim.data_filter_2D, 
                     x = ~UMAP_1, y = ~UMAP_2, 
                     color = ~Dim.data_filter_2D[,41], 
                     colors = col2,
                     type = "scatter", 
                     mode = "markers", 
                     marker = list(size = 4, opacity=0.8), # controls size of points
                     text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                 '<br>Class:', class_label,
                                 '<br>Family:', family_label,
                                 '<br>Subclass:', subclass_label,
                                 '<br>Supertype:', supertype_label,
                                 '<br>Cell-type:', celltype_label,
                                 '<br>Original cell-type:', celltype_original_label,
                                 '<br>Age:', age.at.collection.grp,
                                 '<br>Study:', study_label,
                                 '<br>Region:', region_label,
                                 '<br>RNA-seq method:', RNAseq.method_label,
                                 '<br>Platform:', platform_label,
                                 '<br>CPM:', as.integer(2^Dim.data_filter_2D[,41]-1)), 
                     hoverinfo="text", showlegend=F, height=height) %>% 
    layout(xaxis=list(title="UMAP_1",color="white", zeroline=F),
           yaxis=list(title="UMAP_2",color="white", zeroline=F))%>% 
        colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                              font=list(color="white")), outlinecolor="white",
                 tickcolor="white",tickfont=list(color="white"),
                 tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)),
                 ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)))-1))) %>%  
        layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                  xref='paper', yref='paper', font=list(size=20, color="white")))
      
      
      gp <-  subplot(gp1, gp2, titleX = T, titleY = T)
                     
      
      }
      
      if(ncol(Dim.data_filter_2D)==40){
        gp <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~UMAP_1, y = ~UMAP_2,
                       color = ~Dim.data_filter_2D[,40], 
                       colors = col1,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                       hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="UMAP_1", color="white", zeroline=F),
                 yaxis=list(title="UMAP_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
      }
      

      
    }
    
  ############### tSNE_2D
    else if(dim.red.type=="tSNE_2D"){
      
      if(ncol(Dim.data_filter_2D)==41){
        gp1 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~tSNE_1, y = ~tSNE_2,
                       color = ~Dim.data_filter_2D[,40], 
                       colors = col1,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                       hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="tSNE_1", color="white", zeroline=F),
                 yaxis=list(title="tSNE_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
        gp2 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~tSNE_1, y = ~tSNE_2, 
                       color = ~Dim.data_filter_2D[,41], 
                       colors = col2,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,41]-1)), 
                       hoverinfo="text", showlegend=F, height=height) %>% 
          layout(xaxis=list(title="tSNE_1",color="white", zeroline=F),
                 yaxis=list(title="tSNE_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)))-1))) %>%  
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2, titleX = T, titleY = T)
        
        
      }
      
      if(ncol(Dim.data_filter_2D)==40){
        gp <- plot_ly(data = Dim.data_filter_2D, 
                      x = ~tSNE_1, y = ~tSNE_2,
                      color = ~Dim.data_filter_2D[,40], 
                      colors = col1,
                      type = "scatter", 
                      mode = "markers", 
                      marker = list(size = 4, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                      hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="tSNE_1", color="white", zeroline=F),
                 yaxis=list(title="tSNE_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
      }
    }
    
  ########### PCA_2D
    else if(dim.red.type=="PCA_2D"){
      
      if(ncol(Dim.data_filter_2D)==41){
        gp1 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~PC_1, y = ~PC_2,
                       color = ~Dim.data_filter_2D[,40], 
                       colors = col1,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                       hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="PC_1", color="white", zeroline=F),
                 yaxis=list(title="PC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
        gp2 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~PC_1, y = ~PC_2, 
                       color = ~Dim.data_filter_2D[,41], 
                       colors = col2,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,41]-1)), 
                       hoverinfo="text", showlegend=F, height=height) %>% 
          layout(xaxis=list(title="PC_1",color="white", zeroline=F),
                 yaxis=list(title="PC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)))-1))) %>%  
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2, titleX = T, titleY = T)
        
        
      }
      
      if(ncol(Dim.data_filter_2D)==40){
        gp <- plot_ly(data = Dim.data_filter_2D, 
                      x = ~PC_1, y = ~PC_2,
                      color = ~Dim.data_filter_2D[,40], 
                      colors = col1,
                      type = "scatter", 
                      mode = "markers", 
                      marker = list(size = 4, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                      hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="PC_1", color="white", zeroline=F),
                 yaxis=list(title="PC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
      }
    }
    
  ############ ICA_2D
    else if(dim.red.type=="ICA_2D"){
      
      if(ncol(Dim.data_filter_2D)==41){
        gp1 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~IC_1, y = ~IC_2,
                       color = ~Dim.data_filter_2D[,40], 
                       colors = col1,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                       hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="IC_1", color="white", zeroline=F),
                 yaxis=list(title="IC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
        gp2 <- plot_ly(data = Dim.data_filter_2D, 
                       x = ~IC_1, y = ~IC_2, 
                       color = ~Dim.data_filter_2D[,41], 
                       colors = col2,
                       type = "scatter", 
                       mode = "markers", 
                       marker = list(size = 4, opacity=0.8), # controls size of points
                       text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                   '<br>Class:', class_label,
                                   '<br>Family:', family_label,
                                   '<br>Subclass:', subclass_label,
                                   '<br>Supertype:', supertype_label,
                                   '<br>Cell-type:', celltype_label,
                                   '<br>Original cell-type:', celltype_original_label,
                                   '<br>Age:', age.at.collection.grp,
                                   '<br>Study:', study_label,
                                   '<br>Region:', region_label,
                                   '<br>RNA-seq method:', RNAseq.method_label,
                                   '<br>Platform:', platform_label,
                                   '<br>CPM:', as.integer(2^Dim.data_filter_2D[,41]-1)), 
                       hoverinfo="text", showlegend=F, height=height) %>% 
          layout(xaxis=list(title="IC_1",color="white", zeroline=F),
                 yaxis=list(title="IC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[2], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,41]), length.out=5)))-1))) %>%  
          layout(annotations = list(x = 0.5 , y = 1, text = genes[2], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")))
        
        
        gp <-  subplot(gp1, gp2, titleX = T, titleY = T)
        
        
      }
      
      if(ncol(Dim.data_filter_2D)==40){
        gp <- plot_ly(data = Dim.data_filter_2D, 
                      x = ~IC_1, y = ~IC_2,
                      color = ~Dim.data_filter_2D[,40], 
                      colors = col1,
                      type = "scatter", 
                      mode = "markers", 
                      marker = list(size = 4, opacity=0.8), # controls size of points
                      text=~paste('Cell:', rownames(Dim.data_filter_2D),
                                  '<br>Class:', class_label,
                                  '<br>Family:', family_label,
                                  '<br>Subclass:', subclass_label,
                                  '<br>Supertype:', supertype_label,
                                  '<br>Cell-type:', celltype_label,
                                  '<br>Original cell-type:', celltype_original_label,
                                  '<br>Age:', age.at.collection.grp,
                                  '<br>Study:', study_label,
                                  '<br>Region:', region_label,
                                  '<br>RNA-seq method:', RNAseq.method_label,
                                  '<br>Platform:', platform_label,
                                  '<br>CPM:', as.integer(2^Dim.data_filter_2D[,40]-1)), 
                      hoverinfo="text", showlegend=F, height=height)%>% 
          layout(xaxis=list(title="IC_1", color="white", zeroline=F),
                 yaxis=list(title="IC_2",color="white", zeroline=F))%>% 
          colorbar(title = list(text=paste(genes[1], "(CPM)", sep=" "),
                                font=list(color="white")), outlinecolor="white",
                   tickcolor="white",tickfont=list(color="white"),
                   tickmode="array", tickvals=as.integer(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)),
                   ticktext=as.character(as.integer(as.integer(2^(seq(0,max(Dim.data_filter_2D[,40]), length.out=5)))-1))) %>% 
          layout(annotations = list(x = 0.5 , y = 1, text = genes[1], showarrow = F, 
                                    xref='paper', yref='paper', font=list(size=20, color="white")),
                 paper_bgcolor=toRGB("black"), plot_bgcolor="black")
        
      }
    }
  
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  return(gp)

      
  
    
  }


######################################
#Plot heatmap markers Expression Data
######################################

###### group_heatmap_plotly
group_heatmap_plotly=function (data, anno, genes, grouping,  group_order = NULL, stat = "median", 
                               log_scale = TRUE,normalize_rows = FALSE, colorset = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"), 
                               label_height = 0.05, height=height, fontsize=10) {
  
  
  require(ggplot2)
  require(scrattch.vis)
  require(dplyr)
  require(plotly)
  
  genes <- rev(genes)
  group_cols <- group_columns(grouping)
  gene_data <- filter_gene_data(data, genes, anno, group_cols, 
                                group_order, "sample_name")
  if (!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order, 
    ]
  }
  gene_stats1 <- group_stats(gene_data, value_cols = genes, 
                             anno = anno, grouping = group_cols$label, stat = stat)
  max_vals_unscaled <- max_gene_vals(gene_stats1, genes)
  if (log_scale) {
    gene_stats1 <- scale_gene_data(gene_stats1, genes, scale_type = "log10")
  }
  gene_stats <- data_df_to_colors(gene_stats1, value_cols = genes, 
                                  per_col = normalize_rows, colorset = colorset)
  plot_anno <- anno %>% dplyr::select(dplyr::one_of(group_cols$id, 
                                                    group_cols$label, group_cols$color)) %>% unique()
  group_counts <- anno %>% dplyr::group_by_(group_cols$id) %>% 
    dplyr::summarise(group_n = n())
  
  ######## Prepare plot.data
  colnames(gene_stats)[2:ncol(gene_stats)] <- rep("Gene_color", ncol(gene_stats)-1)
  
df_data_list <- vector(mode = 'list', length = length(genes))
for (i in 1:length(genes)) { # plot_data
  df_data_list[[i]] <- dplyr::left_join(
    plot_anno,
    gene_stats[, c(1, i + 1)],
    by = group_cols$label
  )
}
plot_data <- do.call(rbind, df_data_list)
    plot_data <- do.call(rbind,df_data_list)

  
  plot_data <- dplyr::left_join(plot_data, group_counts, by = group_cols$id)
  plot_data <- add_group_xpos(plot_data, group_cols = group_cols, 
                              group_order = group_order)
  
  ######## gene expression value
  colnames(gene_stats1)[2:ncol(gene_stats1)] <- rep("Gene_value", ncol(gene_stats1)-1)
  
  df_value_list <- vector(mode = 'list', length=length(genes))
  
  

  for (i in 1:length(genes)) {#plot_value
    df_value_list[[i]] <-   gene_stats1[,c(1,i+1)]
  }
  
  plot_value <- do.call(rbind,df_value_list)
  

  

  
  
  
  ####### Add expression and gene name to plot_data
  
  plot_data <- plot_data %>% cbind(Gene_value=plot_value[,2], Gene=rep(genes, each=length(unique(plot_data[,2]))))
  
  
  ######### Preapre header_label
  header_labels <- build_header_labels(data = plot_data, grouping = grouping, 
                                       group_order = group_order, ymin = length(genes) + 1, 
                                       label_height = label_height, label_type = "simple")
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  header_labels <- header_labels %>% cbind(group_n=plot_data$group_n[1:length(unique(plot_data[,2]))])
  header_labels$ymin <- length(genes)+0.5
  header_labels$ymax <-length(genes)+0.5+label_height
  
  ###### Plot
  label_grouping <- paste(grouping,":", sep="")
  
  plot_data <- plot_data[order(plot_data$Gene_value),]
  plot_data <- plot_data[order(plot_data[,2]),]
  
  plot_data$Gene <- factor(plot_data$Gene, levels = unique(plot_data$Gene))
  
  if(log_scale==T){
    
    if(length(unique(plot_data$Gene_value))>1){
    p <- ggplot()+
      geom_tile(data=plot_data, aes(x=plot_data[,2], y=Gene,fill=Gene_value,
                                    text=paste(label_grouping ,plot_data[,2],
                                               "<br>Number of cells:", group_n,
                                               "<br>Gene:", Gene,
                                               "<br>25% tmean log2(CPM+1):", 10^Gene_value-1)), color="white")+
      scale_fill_gradientn(colors=colorset,breaks=seq(0, max(plot_data$Gene_value), length.out = 5),
                           labels=format(round(10^seq(0, max(plot_data$Gene_value), length.out = 5)-1,2),nsmall=2))+ labs(fill="25% tmean log2(CPM+1)")+
      scale_y_discrete("",  expand = c(0, 0)) + scale_x_discrete("", expand = c(0, 0), position="top") + 
      theme_classic(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                               face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                      axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      geom_rect(data = header_labels, aes(xmin = xmin, 
                                          xmax = xmax, ymin = ymin, ymax = ymax, 
                                          text=paste(label_grouping ,label,
                                                     "<br>Number of cells:",group_n)), fill = header_labels$color)+
      geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                 size = 0.2, color="white")+
      geom_vline(ggplot2::aes(xintercept = 1.5:length(unique(plot_data[,2]))), 
                 size = 0.2, color="white")
    }
    
    else if(length(unique(plot_data$Gene_value))==1){
      p <- ggplot()+
        geom_tile(data=plot_data, aes(x=plot_data[,2], y=Gene,fill=Gene_value,
                                      text=paste(label_grouping ,plot_data[,2],
                                                 "<br>Number of cells:", group_n,
                                                 "<br>Gene:", Gene,
                                                 "<br>25% tmean log2(CPM+1):", 10^Gene_value-1)), color="white")+
        scale_fill_gradientn(colors=colorset)+ labs(fill="25% tmean log2(CPM+1)")+
        scale_y_discrete("",  expand = c(0, 0)) + scale_x_discrete("", expand = c(0, 0), position="top") + 
        theme_classic(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                                 face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                        axis.ticks.x = element_blank(), axis.line.x = element_blank())+
        geom_rect(data = header_labels, aes(xmin = xmin, 
                                            xmax = xmax, ymin = ymin, ymax = ymax, 
                                            text=paste(label_grouping ,label,
                                                       "<br>Number of cells:",group_n)), fill = header_labels$color)+
        geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                   size = 0.2, color="white")+
        geom_vline(ggplot2::aes(xintercept = 1.5:length(unique(plot_data[,2]))), 
                   size = 0.2, color="white")
    }
    
  }
  
  if(log_scale==F){
    
    if(length(unique(plot_data$Gene_value))>1){
    p <- ggplot()+
      geom_tile(data=plot_data, aes(x=plot_data[,2], y=Gene,fill=Gene_value,
                                    text=paste(label_grouping ,plot_data[,2],
                                               "<br>Number of cells:", group_n,
                                               "<br>Gene:", Gene,
                                               "<br>25% tmean log2(CPM+1):", Gene_value)), color="white")+
      scale_fill_gradientn(colors=colorset,breaks=seq(0, max(plot_data$Gene_value), length.out = 5),
                           labels=format(round(seq(0, max(plot_data$Gene_value), length.out = 5),2),nsmall=2))+ labs(fill="25% tmean log2(CPM+1)")+
      scale_y_discrete("",  expand = c(0, 0)) + scale_x_discrete("", expand = c(0, 0), position="top") + 
      theme_classic(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                               face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                      axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      geom_rect(data = header_labels, aes(xmin = xmin, 
                                          xmax = xmax, ymin = ymin, ymax = ymax, 
                                          text=paste(label_grouping ,label,
                                                     "<br>Number of cells:",group_n)), fill = header_labels$color)+
      geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                 size = 0.2, color="white")+
      geom_vline(ggplot2::aes(xintercept = 1.5:length(unique(plot_data[,2]))), 
                 size = 0.2, color="white")
    }
    
    if(length(unique(plot_data$Gene_value))==1){
      p <- ggplot()+
        geom_tile(data=plot_data, aes(x=plot_data[,2], y=Gene,fill=Gene_value,
                                      text=paste(label_grouping ,plot_data[,2],
                                                 "<br>Number of cells:", group_n,
                                                 "<br>Gene:", Gene,
                                                 "<br>25% tmean log2(CPM+1):", Gene_value)), color="white")+
        scale_fill_gradientn(colors=colorset)+ labs(fill="25% tmean log2(CPM+1)")+
        scale_y_discrete("",  expand = c(0, 0)) + scale_x_discrete("", expand = c(0, 0), position="top") + 
        theme_classic(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                                 face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                        axis.ticks.x = element_blank(), axis.line.x = element_blank())+
        geom_rect(data = header_labels, aes(xmin = xmin, 
                                            xmax = xmax, ymin = ymin, ymax = ymax, 
                                            text=paste(label_grouping ,label,
                                                       "<br>Number of cells:",group_n)), fill = header_labels$color)+
        geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                   size = 0.2, color="white")+
        geom_vline(ggplot2::aes(xintercept = 1.5:length(unique(plot_data[,2]))), 
                   size = 0.2, color="white")
    }
    
  }
  
  ###### Interactive plot
  gp <- ggplotly(p, tooltip = c("text"), height = height) %>% plotly::layout(xaxis=list(side="top"))
  
  return(gp)
  
  
}

Global.markers <- sample(c("Slc17a7","Trp73", "Prox1", "Slc30a3", "Cacna2d3","Rtn4rl1","Cux2", "Cxcl14", "Ndst4", "Grik1", "Cbln4", "Plch1", "Cfap58",
                           "Lef1", "Dcn", "Fign", "Otof", "Rorb", "Whrn", "Fezf2", "Sulf2", "Rspo1", "Scnn1a","Osr1", "Sulf1", "Tshz2", "Bcl6", "Pou3f1","Npr3", 
                           "Fam84b","Ptgfr", "Car3", "Fn1", "Nts", "Fibcd1","Il16", "Ntf3", "Ramp3", "Nxph3","Rai14", "Foxp2","Sla2", "Ly6g6e","Syt6",
                           "Clic5","Nxph4", "Cplx3", "Ctgf", "Adarb2","Prox1","Lamp5","Pax6","Ndnf","Sncg","Vip","Lhx6","Rxfp3","Ntf3","Pdlim5","Rxfp1","Dock5","Lsp1","Slc35d3",
                           "Jam2","Egln3","Fam19a1","Krt73","Serpinf1","Slc17a8","Ntng1","Pthlh","Pcdh11x","Cp","Mybpc1","Gpc3","Slc5a7","Cbln4",
                           "Chat","Rspo1","Lmo1","Tmem176a","Qrfpr","Calcb",
                           "Sst","Pvalb","Sox6","Rbp4","Chodl","Chrna2","Crh","Ptprk","Th", "Nts", "Myh8","Etv1","Calb2","Nmbr","Hpse", "Sfrp2","Necab1","Ctsc",
                           "Id3","Npffr1","Adamtsl1","Cxcr4","Sln","Cryba2","Pde3a","Npy2r","Grem1","Lpl","Vipr2","Sntb1"),10)




# plotly heatmap markers
plotly_heatmap_markers = function(se_obj, tome_file ,grouping="celltype",dim.red.type="UMAP_3D", ident.class=c("Glutamatergic"), 
                                  ident.subclass=c("Lamp5"), ident.family=c("Lamp5"),
                                  ident.celltype=c("Lamp5|L1 A7C/CNC") , ident.celltype_original=c("Lamp5|Egln3|L1 A7C/CNC"),
                                  ident.age=c("P0"), ident.study=c("This study"), ident.region=c("SSp;SSs"), ident.RNAseq.method=c("scRNA-seq"), ident.platform=c("10X"),
                                  GeneList=c("Gad1","Neurod2"), genes=genes, stat="tmean", colorset= "BuGyOrRd",
                                  log_scale=T, height="850px", fontsize=10){
  
  require(scrattch.io)
  #require(scrattch.hicat)
  require(scrattch.vis)
  require(scales)
  #require(conflicted)
  
  
  ########## Color for gene1
  if (colorset == "PuOr") {
  colorset = brewer_pal(palette = "PuOr", direction = -1)(11)
} else if (colorset == "GBBr") {
  colorset = brewer_pal(palette = "BrBG", direction = -1)(11)
}
  
  else if(colorset=="BuYlRd"){
    colorset = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(colorset=="BuPhRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="GnPhRd"){
    colorset =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="BuGyOrRd"){
    colorset =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(colorset=="BuGyRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="Jet"){
    colorset =jet.colors(10)
  }
  else if(colorset=="Rainbow"){
    colorset =rainbow(5,rev=T)
  }
  else if(colorset=="Heat"){
    colorset =heat.colors(10, rev=T)
  }
  else if(colorset=="viridis"){
    colorset =viridis_pal(option = "viridis")(10)
  }
  else if(colorset=="magma"){
    colorset =viridis_pal(option = "magma")(10)
  }
  else if(colorset=="inferno"){
    colorset =viridis_pal(option = "inferno")(10)
  }
  else if(colorset=="plasma"){
    colorset =viridis_pal(option = "plasma")(10)
  }
  else if(colorset=="cividis"){
    colorset =viridis_pal(option = "cividis")(10)
  }
  else if (colorset == "Purples_cb") {
  colorset = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ######Subset class
  se_obj <- subset(se_obj, select= colData(se_obj)$class_label %in% ident.class)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  
  ######Subset family
  se_obj <- subset(se_obj, select= colData(se_obj)$family_label %in% ident.family)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  ######Subset subclass
  se_obj <- subset(se_obj, select= colData(se_obj)$subclass_label %in% ident.subclass)
  
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family) 
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset celltype
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_label %in% ident.celltype)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  

  
  #Subset celltype_original
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_original_label %in% ident.celltype_original)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset age
  se_obj <- subset(se_obj, select= colData(se_obj)$age.at.collection.grp %in% ident.age)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset region
  se_obj <- subset(se_obj, select= colData(se_obj)$region_label %in% ident.region)
  
  ident.class <- intersect(se_obj$class_label, ident.class)
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$platform_label %in% ident.platform)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset RNAseq.method
  se_obj <- subset(se_obj, select= colData(se_obj)$RNAseq.method_label %in% ident.RNAseq.method)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$study_label %in% ident.study)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  #Prepare grouping in factor
  celltype.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(celltype_id, celltype_label, celltype_color) %>% 
    unique() %>% 
    dplyr::arrange(celltype_id)
  
  celltype_original.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(celltype_original_id, celltype_original_label, celltype_original_color) %>% 
    unique() %>% 
    dplyr::arrange(celltype_original_id)
  
  subclass.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(subclass_id, subclass_label, subclass_color) %>% 
    unique() %>% 
    dplyr::arrange(subclass_id)
  
  family.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(family_id, family_label, family_color) %>% 
    unique() %>% 
    dplyr::arrange(family_id)
  
  class.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(class_id, class_label, class_color) %>% 
    unique() %>% 
    dplyr::arrange(class_id)
  
  
  se_obj$celltype_original_label <- factor(se_obj$celltype_original_label, levels = celltype_original.df$celltype_original_label)
  se_obj$celltype_label <- factor(se_obj$celltype_label, levels = celltype.df$celltype_label)
  se_obj$subclass_label <- factor(se_obj$subclass_label, levels = subclass.df$subclass_label)
  se_obj$family_label <- factor(se_obj$family_label, levels = family.df$family_label)
  se_obj$class_label <- factor(se_obj$class_label, levels = class.df$class_label)

  
  ########## prepare data.frame

  genes<- genes[which(genes%in%rownames(se_obj))]
  
  norm.se_obj <- read_tome_gene_data(tome_file, genes = genes, format = "matrix")
  norm.se_obj <- norm.se_obj[rownames(norm.se_obj)%in%colnames(se_obj),,drop=F]
  norm.se_obj <- mat_to_data_df(norm.se_obj)
  
  


  genes<- make.names(genes)
  
  colnames(norm.se_obj) <- make.names(colnames(norm.se_obj) )
  

  
  ####### group_heatmap_plot
  plot_height= height + 20*length(unique(genes))
  
  label_height=plot_height/(7.5*height)
  
 
  
  gp <-  group_heatmap_plotly(norm.se_obj,
                              as.data.frame(colData(se_obj)), 
                              colorset = colorset,
                              genes = genes, 
                              grouping = grouping, 
                              stat = stat,
                              label_height = label_height,
                              fontsize = fontsize,
                              log_scale = log_scale,
                              height=plot_height)
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  return(gp)
  
  
}  


#######################
#group_dot_plotly
#######################

group_dot_plotly=function (data, anno, genes, grouping,  group_order = NULL, stat = "tmean", 
                           log_scale = TRUE,normalize_rows = FALSE, colorset = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"), 
                           label_height = 0.05, height=height, fontsize=10, max_size=6, size_stat="prop_gt0") {
  
  
  require(ggplot2)
  require(scrattch.vis)
  require(dplyr)
  require(plotly)
  
  genes <- rev(genes)
  group_cols <- group_columns(grouping)
  gene_data <- filter_gene_data(data, genes, anno, group_cols, 
                                group_order, "sample_name")
  if (!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order, 
    ]
  }
  gene_stats1 <- group_stats(gene_data, value_cols = genes, 
                             anno = anno, grouping = group_cols$label, stat = stat)
  max_vals_unscaled <- max_gene_vals(gene_stats1, genes)
  if (log_scale) {
    gene_stats1 <- scale_gene_data(gene_stats1, genes, scale_type = "log10")
  }
  gene_stats <- data_df_to_colors(gene_stats1, value_cols = genes, 
                                  per_col = normalize_rows, colorset = colorset)
                                  colnames(gene_stats)[2:ncol(gene_stats)] <- rep("Gene_color", ncol(gene_stats) - 1)
  plot_anno <- anno %>% dplyr::select(dplyr::one_of(group_cols$id, 
                                                    group_cols$label, group_cols$color)) %>% unique()
  group_counts <- anno %>% dplyr::group_by_(group_cols$id) %>% 
    dplyr::summarise(group_n = n())
  
  gene_size_stats <- group_stats(gene_data, value_cols = genes, 
                                 anno = anno, grouping = group_cols$label, stat = size_stat)
  
  ######## Prepare plot.data
df_data_list <- vector(mode = 'list', length = length(genes))
for (i in 1:length(genes)) { # plot_data
  df_data_list[[i]] <- dplyr::left_join(
    plot_anno,
    gene_stats[, c(1, i + 1)],
    by = group_cols$label
  )
}
  
  plot_data <- do.call(rbind,df_data_list)
  
  
  plot_data <- dplyr::left_join(plot_data, group_counts, by = group_cols$id)
  plot_data <- add_group_xpos(plot_data, group_cols = group_cols, 
                              group_order = group_order)
  
  ######## gene expression value
  colnames(gene_stats1)[2:ncol(gene_stats1)] <- rep("Gene_value", ncol(gene_stats1)-1)
  
  df_value_list <- vector(mode = 'list', length=length(genes))
  
  for (i in 1:length(genes)) {#plot_value
    df_value_list[[i]] <-   gene_stats1[,c(1,i+1)]
  }
  
  plot_value <- do.call(rbind,df_value_list)
  
  ######## gene prortion of samples >0
  colnames(gene_size_stats)[2:ncol(gene_size_stats)] <- rep("Gene_size", ncol(gene_size_stats)-1)
  
  df_size_list <- vector(mode = 'list', length=length(genes))
  
  for (i in 1:length(genes)) {#plot_size
    df_size_list[[i]] <-   gene_size_stats[,c(1,i+1)]
  }
  
  plot_size <- do.call(rbind,df_size_list)
  
  ####### Add expression and gene name to plot_data
  
  plot_data <- plot_data %>% cbind(Gene_value=plot_value[,2], Gene=rep(genes, each=length(unique(plot_data[,2]))), Gene_size=plot_size[,2])
  
  
  ######### Prepare header_label
  header_labels <- build_header_labels(data = plot_data, grouping = grouping, 
                                       group_order = group_order, ymin = length(genes) + 1, 
                                       label_height = label_height, label_type = "simple")
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  header_labels <- header_labels %>% cbind(group_n=plot_data$group_n[1:length(unique(plot_data[,2]))])
  header_labels$ymin <- length(genes)+0.5
  header_labels$ymax <-length(genes)+0.5+label_height
  
  ###### Plot
  label_grouping <- paste(grouping,":", sep="")
  
  plot_data <- plot_data[order(plot_data$Gene_value),]
  plot_data <- plot_data[order(plot_data[,2]),]
  
  plot_data$Gene <- factor(plot_data$Gene, levels = unique(plot_data$Gene))
  
  if(log_scale==T){
    
    if(length(unique(plot_data$Gene_value))>1){
    p <- ggplot()+
      geom_point(data=plot_data, aes(x=plot_data[,2], y=Gene,color=Gene_value, size=Gene_size,
                                     text=paste(label_grouping ,plot_data[,2],
                                                "<br>Number of cells:", group_n,
                                                "<br>Gene:", Gene,
                                                "<br>25% tmean log2(CPM+1):", 10^Gene_value-1,
                                                "<br>Proportion of samples > 0:", Gene_size)))+
      scale_color_gradientn(colors=colorset,breaks=seq(0, max(plot_data$Gene_value), length.out = 5),
                            labels=format(round(10^seq(0, max(plot_data$Gene_value), length.out = 5)-1,2),nsmall=2))+ labs(color="25% tmean log2(CPM+1)")+
      scale_size_area(max_size = max_size)+
      scale_y_discrete("") + scale_x_discrete("", expand = c(0, 0), position="top") + 
      theme_bw(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                          face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                 axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      geom_rect(data = header_labels, aes(xmin = xmin, 
                                          xmax = xmax, ymin = ymin, ymax = ymax, 
                                          text=paste(label_grouping ,label,
                                                     "<br>Number of cells:",group_n)), fill = header_labels$color, color="white")+
      geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                 size = 0.2, color="black")

    }
 
  
  else if(length(unique(plot_data$Gene_value))==1){
    p <- ggplot()+
      geom_point(data=plot_data, aes(x=plot_data[,2], y=Gene,color=Gene_value, size=Gene_size,
                                     text=paste(label_grouping ,plot_data[,2],
                                                "<br>Number of cells:", group_n,
                                                "<br>Gene:", Gene,
                                                "<br>25% tmean log2(CPM+1):", 10^Gene_value-1,
                                                "<br>Proportion of samples > 0:", Gene_size)))+
      scale_color_gradientn(colors=colorset)+ labs(color="25% tmean log2(CPM+1)")+
      scale_size_area(max_size = max_size)+
      scale_y_discrete("") + scale_x_discrete("", expand = c(0, 0), position="top") + 
      theme_bw(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                          face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                 axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      geom_rect(data = header_labels, aes(xmin = xmin, 
                                          xmax = xmax, ymin = ymin, ymax = ymax, 
                                          text=paste(label_grouping ,label,
                                                     "<br>Number of cells:",group_n)), fill = header_labels$color, color="white")+
      geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                 size = 0.2, color="black")
    
  }
}
    
  
  if(log_scale==F){
    
    if(length(unique(plot_data$Gene_value))>1){
    p <- ggplot()+
      geom_point(data=plot_data, aes(x=plot_data[,2], y=Gene,color=Gene_value, size=Gene_size,
                                     text=paste(label_grouping ,plot_data[,2],
                                                "<br>Number of cells:", group_n,
                                                "<br>Gene:", Gene,
                                                "<br>25% tmean log2(CPM+1):", Gene_value,
                                                "<br>Proportion of samples > 0:", Gene_size)))+
      scale_color_gradientn(colors=colorset,breaks=seq(0, max(plot_data$Gene_value), length.out = 5),
                            labels=format(round(seq(0, max(plot_data$Gene_value), length.out = 5),2),nsmall=2))+ labs(color="25% tmean log2(CPM+1)")+
      scale_size_area(max_size = max_size)+
      scale_y_discrete("") + scale_x_discrete("", expand = c(0, 0), position="top") + 
      theme_bw(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                          face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                 axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      geom_rect(data = header_labels, aes(xmin = xmin, 
                                          xmax = xmax, ymin = ymin, ymax = ymax, 
                                          text=paste(label_grouping ,label,
                                                     "<br>Number of cells:",group_n)), fill = header_labels$color)+
      geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                 size = 0.2, color="black")
    }
    
    else if(length(unique(plot_data$Gene_value))==1){
      p <- ggplot()+
        geom_point(data=plot_data, aes(x=plot_data[,2], y=Gene,color=Gene_value, size=Gene_size,
                                       text=paste(label_grouping ,plot_data[,2],
                                                  "<br>Number of cells:", group_n,
                                                  "<br>Gene:", Gene,
                                                  "<br>25% tmean log2(CPM+1):", Gene_value,
                                                  "<br>Proportion of samples > 0:", Gene_size)))+
        scale_color_gradientn(colors=colorset)+ labs(color="25% tmean log2(CPM+1)")+
        scale_size_area(max_size = max_size)+
        scale_y_discrete("") + scale_x_discrete("", expand = c(0, 0), position="top") + 
        theme_bw(fontsize) + theme(axis.text = element_text(size = rel(1), 
                                                            face = "italic"), axis.text.x = element_text(angle=90, hjust=0), 
                                   axis.ticks.x = element_blank(), axis.line.x = element_blank())+
        geom_rect(data = header_labels, aes(xmin = xmin, 
                                            xmax = xmax, ymin = ymin, ymax = ymax, 
                                            text=paste(label_grouping ,label,
                                                       "<br>Number of cells:",group_n)), fill = header_labels$color)+
        geom_hline(ggplot2::aes(yintercept = 1.5:length(unique(genes))), 
                   size = 0.2, color="black")
      
    }
  }
  
  ###### Interactive plot
  gp <- ggplotly(p, tooltip = c("text"), height = height) %>% plotly::layout(xaxis=list(side="top"))
  
  
  return(gp)
  
  
}



# plotly dot markers
plotly_dot_markers = function(se_obj, tome_file ,grouping="celltype",dim.red.type="UMAP_3D", ident.class=c("Glutamatergic"), ident.subclass=c("Lamp5"), 
                              ident.family=c("Lamp5"), 
                              ident.celltype=c("Lamp5|L1 A7C/CNC") , ident.celltype_original=c("Lamp5|Egln3|L1 A7C/CNC"),
                              ident.age=c("P0"), ident.study=c("This study"), ident.region=c("SSp;SSs"), ident.RNAseq.method=c("scRNA-seq"),
                              ident.platform=c("10X"),
                              GeneList=c("Gad1","Neurod2"), genes=genes, stat="tmean", colorset= "BuGyOrRd",
                              log_scale=T, height="850px", fontsize=10, max_size=6, size_stat="prop_gt0"){
  
  require(scrattch.io)
  require(scrattch.vis)
  require(scales)
  
  
  ########## Color for gene1
  if(colorset=="OrPu"){
    colorset =  brewer_pal(palette = "PuOr")(8)
  }
  else if(colorset=="BrBG "){
    colorset =  brewer_pal(palette = "BrBG")(8)
  }
  
  else if(colorset=="BuYlRd"){
    colorset = brewer_pal(palette = "RdYlBu", direction = -1)(8)
  }
  else if(colorset=="BuPhRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="GnPhRd"){
    colorset =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="BuGyOrRd"){
    colorset =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(colorset=="BuGyRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="Jet"){
    colorset =jet.colors(10)
  }
  else if(colorset=="Rainbow"){
    colorset =rainbow(5,rev=T)
  }
  else if(colorset=="Heat"){
    colorset =heat.colors(10, rev=T)
  }
  else if(colorset=="viridis"){
    colorset =viridis_pal(option = "viridis")(10)
  }
  else if(colorset=="magma"){
    colorset =viridis_pal(option = "magma")(10)
  }
  else if(colorset=="inferno"){
    colorset =viridis_pal(option = "inferno")(10)
  }
  else if(colorset=="plasma"){
    colorset =viridis_pal(option = "plasma")(10)
  }
  else if(colorset=="cividis"){
    colorset =viridis_pal(option = "cividis")(10)
  }
  else if (colorset == "Purples_cb") {
  colorset = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
} else if (colorset == "PuOr") {
  colorset = brewer_pal(palette = "PuOr", direction = -1)(10)
} else if (colorset == "GBBr") {
  colorset = brewer_pal(palette = "BrBG", direction = -1)(10)
}
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ######Subset class
  se_obj <- subset(se_obj, select= colData(se_obj)$class_label %in% ident.class)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study)
  
  ######Subset family
  se_obj <- subset(se_obj, select= colData(se_obj)$family_label %in% ident.family)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  ######Subset subclass
  se_obj <- subset(se_obj, select= colData(se_obj)$subclass_label %in% ident.subclass)
  
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset celltype
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_label %in% ident.celltype)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  

  
  #Subset celltype_original
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_original_label %in% ident.celltype_original)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset age
  se_obj <- subset(se_obj, select= colData(se_obj)$age.at.collection.grp %in% ident.age)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset region
  se_obj <- subset(se_obj, select= colData(se_obj)$region_label %in% ident.region)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$platform_label %in% ident.platform)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset RNAseq.method
  se_obj <- subset(se_obj, select= colData(se_obj)$RNAseq.method_label %in% ident.RNAseq.method)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$study_label %in% ident.study)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.celltype_original <- intersect(se_obj$celltype_original_label, ident.celltype_original) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Prepare grouping in factor
  celltype.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(celltype_id, celltype_label, celltype_color) %>% 
    unique() %>% 
    dplyr::arrange(celltype_id)
  celltype_original.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(celltype_original_id, celltype_original_label, celltype_original_color) %>% 
    unique() %>% 
    dplyr::arrange(celltype_original_id)

  subclass.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(subclass_id, subclass_label, subclass_color) %>% 
    unique() %>% 
    dplyr::arrange(subclass_id)
  
  family.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(family_id, family_label, family_color) %>% 
    unique() %>% 
    dplyr::arrange(family_id)
  
  class.df <- as.data.frame(colData(se_obj)) %>% 
    dplyr::select(class_id, class_label, class_color) %>% 
    unique() %>% 
    dplyr::arrange(class_id)
  
  
  se_obj$celltype_original_label <- factor(se_obj$celltype_original_label, levels = celltype_original.df$celltype_original_label)
  se_obj$celltype_label <- factor(se_obj$celltype_label, levels = celltype.df$celltype_label)
  se_obj$subclass_label <- factor(se_obj$subclass_label, levels = subclass.df$subclass_label)
  se_obj$family_label <- factor(se_obj$family_label, levels = family.df$family_label)
  se_obj$class_label <- factor(se_obj$class_label, levels = class.df$class_label)
  
  
  ########## prepare data.frame
  genes<- genes[which(genes%in%rownames(se_obj))]
  
  norm.se_obj <- read_tome_gene_data(tome_file, genes = genes, format = "matrix")
  norm.se_obj <- norm.se_obj[rownames(norm.se_obj)%in%colnames(se_obj),,drop=F]
  norm.se_obj <- mat_to_data_df(norm.se_obj)
  
  
  genes<- make.names(genes)
  
  colnames(norm.se_obj) <- make.names(colnames(norm.se_obj) )
  

  ####### group_dot_plot
  plot_height= height + 20*length(unique(genes))
 
  
  label_height=plot_height/(7.5*height)
  
 
  
  
  gp <-  group_dot_plotly(norm.se_obj,
                              as.data.frame(colData(se_obj)), 
                              colorset = colorset,
                              genes = genes, 
                              grouping = grouping, 
                              stat = stat,
                              label_height = label_height,
                              fontsize = fontsize,
                              log_scale = log_scale,
                              max_size = max_size,
                              size_stat = size_stat,
                              height=plot_height)
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  
  return(gp)
  
  
}  

#############################
#  Plot Heatmap by ct*Age
############################

plotlyHm_Exprgene_byAge <- function(se_obj, tome_file ,grouping="celltype", gene="Cdh13", ident.class=c("Glutamatergic"), 
                                    ident.family=c("Lamp5"), ident.subclass=c("Lamp5"), 
                                    ident.celltype=c("Lamp5|L1 A7C/CNC") , colorset="PuOr",
                                    ident.age=c("P0"), ident.study=c("This study"), ident.region=c("SSp;SSs"), ident.RNAseq.method=c("scRNA-seq"), 
                                    ident.platform=c("10X"),
                                    height="850px"){
  
  require(scrattch.io)
  require(scrattch.vis)
  require(scales)

  
  
  ########## Color gradient
  if(colorset=="PuOr"){
    colorset =  brewer_pal(palette = "PuOr")(11)
  }
  else if(colorset=="GBBr "){
    colorset =  brewer_pal(palette = "BrBG")(11)
  }
  
  else if(colorset=="BuYlRd"){
    colorset = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(colorset=="BuPhRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="GnPhRd"){
    colorset =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="BuGyOrRd"){
    colorset =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(colorset=="BuGyRd"){
    colorset =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(colorset=="Jet"){
    colorset =jet.colors(10)
  }
  else if(colorset=="Rainbow"){
    colorset =rainbow(5,rev=T)
  }
  else if(colorset=="Heat"){
    colorset =heat.colors(10, rev=T)
  }
  else if(colorset=="viridis"){
    colorset =viridis_pal(option = "viridis")(10)
  }
  else if(colorset=="magma"){
    colorset =viridis_pal(option = "magma")(10)
  }
  else if(colorset=="inferno"){
    colorset =viridis_pal(option = "inferno")(10)
  }
  else if(colorset=="plasma"){
    colorset =viridis_pal(option = "plasma")(10)
  }
  else if(colorset=="cividis"){
    colorset =viridis_pal(option = "cividis")(10)
  }
       else if(colorset=="Purples_cb"){
  colorset = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ######Subset class
  se_obj <- subset(se_obj, select= colData(se_obj)$class_label %in% ident.class)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  ######Subset family
  se_obj <- subset(se_obj, select= colData(se_obj)$family_label %in% ident.family)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  ######Subset subclass
  se_obj <- subset(se_obj, select= colData(se_obj)$subclass_label %in% ident.subclass)
  
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset celltype
  se_obj <- subset(se_obj, select= colData(se_obj)$celltype_label %in% ident.celltype)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  

  

  
  #Subset age
  se_obj <- subset(se_obj, select= colData(se_obj)$age.at.collection.grp %in% ident.age)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset region
  se_obj <- subset(se_obj, select= colData(se_obj)$region_label %in% ident.region)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$platform_label %in% ident.platform)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset RNAseq.method
  se_obj <- subset(se_obj, select= colData(se_obj)$RNAseq.method_label %in% ident.RNAseq.method)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  #Subset platform
  se_obj <- subset(se_obj, select= colData(se_obj)$study_label %in% ident.study)
  
  ident.class <- intersect(se_obj$class_label, ident.class) 
  ident.family <- intersect(se_obj$family_label, ident.family)
  ident.subclass <- intersect(se_obj$subclass_label, ident.subclass) 
  ident.celltype <- intersect(se_obj$celltype_label, ident.celltype) 
  ident.age <- intersect(se_obj$age.at.collection.grp, ident.age) 
  ident.region <- intersect(se_obj$region_label, ident.region) 
  ident.platform <- intersect(se_obj$platform_label, ident.platform) 
  ident.RNAseq.method <- intersect(se_obj$RNAseq.method_label, ident.RNAseq.method) 
  ident.study <- intersect(se_obj$study_label, ident.study) 
  
  
  ########## prepare data.frame
  gene<- gene[which(gene%in%rownames(se_obj))]
  
  norm.se_obj <- read_tome_gene_data(tome_file, gene = gene, format = "matrix")
  norm.se_obj <- norm.se_obj[rownames(norm.se_obj)%in%colnames(se_obj),,drop=F]
  norm.se_obj <- mat_to_data_df(norm.se_obj)
  
  gene<- make.names(gene)
  
  colnames(norm.se_obj) <- make.names(colnames(norm.se_obj) )
  
  #Prepare grouping in factor


  if(grouping=="celltype"){
  
  celltype.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(celltype_id,
                                                                  celltype_label, 
                                                                  celltype_color) %>% 
    unique() %>% arrange(as.numeric(celltype_id)) 
  
  
  
  Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","celltype_label","celltype_color","age.at.collection.grp")])
  Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))


  ####### group_dot_plot

  colnames(Hm.df)[5] <- "expression"
  
  
  age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
  Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
  
  Hm.df <-  aggregate(expression ~ age.at.collection.grp+celltype_label , data =Hm.df, FUN= "mean",trim=0.25 )
  
  
  Hm.df$celltype_label <- factor(Hm.df$celltype_label, levels=rev(celltype.df$celltype_label[celltype.df$celltype_label%in%Hm.df$celltype_label]))

  plot_height= height + 20*length(unique(celltype.df$celltype_label))
  
  }
  
  
  if(grouping=="celltype"){
    
    celltype.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(celltype_id,
                                                                    celltype_label, 
                                                                    celltype_color) %>% 
      unique() %>% arrange(as.numeric(celltype_id)) 
    
    
    
    Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","celltype_label","celltype_color","age.at.collection.grp")])
    Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))
    
    
    ####### group_dot_plot
    
    colnames(Hm.df)[5] <- "expression"
    
    
    age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
    Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
    
    Hm.df <-  aggregate(expression ~ age.at.collection.grp+celltype_label , data =Hm.df, FUN= "mean",trim=0.25 )
    
    
    Hm.df$celltype_label <- factor(Hm.df$celltype_label, levels=rev(celltype.df$celltype_label[celltype.df$celltype_label%in%Hm.df$celltype_label]))
    
    plot_height= height + 20*length(unique(celltype.df$celltype_label))
    
  }
  
  
  else if(grouping=="class"){
    
    class.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(class_id,
                                                                    class_label, 
                                                                    class_color) %>% 
      unique() %>% arrange(as.numeric(class_id)) 
    
    
    
    Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","class_label","class_color","age.at.collection.grp")])
    Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))
    
    
    ####### group_dot_plot
    
    colnames(Hm.df)[5] <- "expression"
    
    
    age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
    Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
    
    Hm.df <-  aggregate(expression ~ age.at.collection.grp+class_label , data =Hm.df, FUN= "mean",trim=0.25 )
    
    
    Hm.df$class_label <- factor(Hm.df$class_label, levels=rev(class.df$class_label[class.df$class_label%in%Hm.df$class_label]))
    
    plot_height= height + 20*length(unique(class.df$class_label))
    
  }
  
  
  else if(grouping=="subclass"){
    
    subclass.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(subclass_id,
                                                                 subclass_label, 
                                                                 subclass_color) %>% 
      unique() %>% arrange(as.numeric(subclass_id)) 
    
    
    
    Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","subclass_label","subclass_color","age.at.collection.grp")])
    Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))
    
    
    ####### group_dot_plot
    
    colnames(Hm.df)[5] <- "expression"
    
    
    age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
    Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
    
    Hm.df <-  aggregate(expression ~ age.at.collection.grp+subclass_label , data =Hm.df, FUN= "mean",trim=0.25 )
    
    
    Hm.df$subclass_label <- factor(Hm.df$subclass_label, levels=rev(subclass.df$subclass_label[subclass.df$subclass_label%in%Hm.df$subclass_label]))
    
    plot_height= height + 20*length(unique(subclass.df$subclass_label))
    
  }
  
  
  else if(grouping=="family"){
    
    family.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(family_id,
                                                                 family_label, 
                                                                 family_color) %>% 
      unique() %>% arrange(as.numeric(family_id)) 
    
    
    
    Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","family_label","family_color","age.at.collection.grp")])
    Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))
    
    
    ####### group_dot_plot
    
    colnames(Hm.df)[5] <- "expression"
    
    
    age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
    Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
    
    Hm.df <-  aggregate(expression ~ age.at.collection.grp+family_label , data =Hm.df, FUN= "mean",trim=0.25 )
    
    
    Hm.df$family_label <- factor(Hm.df$family_label, levels=rev(family.df$family_label[family.df$family_label%in%Hm.df$family_label]))
    
    plot_height= height + 20*length(unique(family.df$family_label))
    
  }
  
  
  
  else if(grouping=="celltype_original"){
    
    celltype_original.df <- as.data.frame(colData(se_obj)) %>% dplyr::select(celltype_original_id,
                                                                 celltype_original_label, 
                                                                 celltype_original_color) %>% 
      unique() %>% arrange(as.numeric(celltype_original_id)) 
    
    
    
    Hm.df <- as.data.frame(colData(se_obj)[,c("sample_name","celltype_original_label","celltype_original_color","age.at.collection.grp")])
    Hm.df <- dplyr::left_join(Hm.df, norm.se_obj, by=c("sample_name"))
    
    
    ####### group_dot_plot
    
    colnames(Hm.df)[5] <- "expression"
    
    
    age.var=c("E11.5","E12.5","E13.5","E14.5","E15.5","E16.5","E17.5","E18.5","P0","P1","P2","P4","P5","P8","P16","P30","Adult")
    Hm.df$age.at.collection.grp <- factor(Hm.df$age.at.collection.grp, levels=age.var[age.var%in%Hm.df$age.at.collection.grp])
    
    Hm.df <-  aggregate(expression ~ age.at.collection.grp+celltype_original_label , data =Hm.df, FUN= "mean",trim=0.25 )
    
    
    Hm.df$celltype_original_label <- factor(Hm.df$celltype_original_label, levels=rev(celltype_original.df$celltype_original_label[celltype_original.df$celltype_original_label%in%Hm.df$celltype_original_label]))
    
    plot_height= height + 20*length(unique(celltype_original.df$celltype_original_label))
    
  }
  
  
  
  ########## Plot


  label_grouping <- paste0("<br>",grouping,":")

  p <- ggplot(Hm.df, aes(x=age.at.collection.grp,y=Hm.df[,2]))+
    geom_tile(aes(fill=expression,
                  text=paste("<br>Age:", age.at.collection.grp,
                             label_grouping, Hm.df[,2],
                             "<br>25% tmean log2(CPM+1):", expression)))+
    scale_fill_gradientn(colors = colorset, na.value = "black", limits=c(0,max(Hm.df$expression)))+
    geom_hline(ggplot2::aes(yintercept = 0.5:(nrow(Hm.df)-0.5)), 
               size = 0.2, color="white")+
    geom_vline(ggplot2::aes(xintercept = 0.5:(nrow(Hm.df)-0.5)), 
               size = 0.2, color="white")+
    xlab("Age")+ylab("Cell-type")+labs(fill="25% tmean log2(CPM+1)")+
    mytheme_classic+theme(axis.line = element_blank(),
                          axis.ticks =  element_blank(),
                          axis.text.x = element_text(angle=45, hjust=1),
                          axis.ticks.length = unit(1,"points"),
                          plot.title = element_text(size = 14),
                          plot.margin = margin(c(0,0,0,10,"cm")) )+
  ggtitle(gene)+
  ggeasy::easy_center_title()
  
 

  gp <-ggplotly(p, tooltip = c("text"), height = plot_height) 
  

  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"))
  
  return(gp)
  
}

#############################
# Plot Curve pseudomaturation
#############################

Gene.pM <- sort(unique(unlist(lapply(Res_pt.list, function(Res){rownames(Res$smooth.scale)}))))

plotlyGene_smoothSpline_pt= function(Res.list,celltype=c("Lamp5|L1 A7C/CNC"), gene = NULL, free.scale = TRUE,  ncol = NA, 
                                     line.size = 1,  showlegend=F, height="850px",width="850px"){
  
  require(ggplot2)
  require(dplyr)
  require(plotly)
  
  
  #Set facet_wrap option
  if (is.na(ncol)) 
    nrow = round(sqrt(length(gene)))
  else nrow = NA
  a <- if (free.scale) 
    "free_y"
  else "fixed"  
  
  Res.list <- Res.list[celltype]
  
  name.list <- as.list(names(Res.list))
  
  wave.pattern.list <- mapply(function(Res, name){
    
    wave.pattern.df <- reshape2::melt(Res$smooth.scale[rownames(Res$smooth.scale)%in%gene,,drop=F])
    colnames(wave.pattern.df)[1:3] <- c("gene","n.point","expression")
    wave.pattern.df$n.point <- as.character(wave.pattern.df$n.point)
    
    pt.df <- reshape2::melt(Res$pseudotime)
    pt.df <- pt.df %>% cbind(n.point=rownames(pt.df))
    colnames(pt.df)[1] <- "pseudotime"
    
    
    
    wave.pattern.df <- left_join(wave.pattern.df, pt.df, by=c("n.point"))
    
    wave.pattern.df <-wave.pattern.df %>% cbind(celltype=name)
    return(wave.pattern.df)
    
  }, Res=Res.list, name=name.list, SIMPLIFY = F)
  
  wave.pattern <- do.call(rbind, wave.pattern.list)
  
  wave.pattern$celltype <- factor(wave.pattern$celltype, levels = unique(wave.pattern$celltype))
  wave.pattern$gene <- factor(wave.pattern$gene, levels = gene)
  
  color.celltype <- unlist(lapply(Res.list,function(x){x$color.celltype}))
  
  p <- ggplot(wave.pattern, aes(y=expression, x=pseudotime))+
    geom_line(aes(y = expression, x = pseudotime,  color=celltype),show.legend = showlegend, size=line.size)+
    scale_color_manual(values = color.celltype)+
    scale_x_continuous(breaks= c(0, 1))+xlab("Pseudo-maturation")+ylab("Expression")+
    mytheme_gray
  
  
  
  

    p <- p + facet_wrap(~gene, ncol = round(sqrt(length(gene))), scales = a) +theme(strip.text = element_text(size = 14))

  
  plot_height=  height + 10*length(unique(gene))
  plot_width=  width + 10*length(unique(gene))
  
    
  gp <-  ggplotly(p,  height=plot_height, width=plot_width) 
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  print(gp)
  
  
}



plotGene_Hm_pt= function(Res.list, celltype=c("Lamp5|L1 A7C/CNC"), gene = NULL, free.scale = TRUE,  ncol = NA, 
                         col_gdt=c("blue3", "skyblue","white", "pink", "red3")){
  
  require(ggplot2)
  require(dplyr)
  require(plotly)
  
  ########## Color gradient
  if(col_gdt=="PuOr"){
    col_gdt =  brewer_pal(palette = "PuOr",direction = -1)(11)
  }
  else if(col_gdt=="GBBr"){
    col_gdt =  brewer_pal(palette = "BrBG",direction = -1)(11)
  }
  
  else if(col_gdt=="BuYlRd"){
    col_gdt = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(col_gdt=="BuRd"){
    col_gdt = brewer_pal(palette = "RdBu", direction = -1)(11)
  }
  else if(col_gdt=="BuPhRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="GnPhRd"){
    col_gdt =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="BuGyOrRd"){
    col_gdt =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(col_gdt=="BuGyRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="Jet"){
    col_gdt =jet.colors(10)
  }
  else if(col_gdt=="BuSkyWhPkRd"){
    col_gdt =c("blue3", "skyblue","white", "pink", "red3")
  }
  else if(col_gdt=="Vik"){
    col_gdt =scico(11, palette="vik")
  }
     else if(col_gdt=="Purples_cb"){
  col_gdt = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}

  
  
  #Set facet_wrap option
  if (is.na(ncol)) 
    nrow = round(sqrt(length(gene)))
  else nrow = NA
  a <- if (free.scale) 
    "free_y"
  else "fixed"  
  
  name.list <- as.list(names(Res.list))
  names(name.list) <- names(Res.list)
  
  
  wave.pattern.list <- mapply(function(Res, name){
    
    wave.pattern.df <- reshape2::melt(Res$smooth.scale[rownames(Res$smooth.scale)%in%gene,,drop=F])
    colnames(wave.pattern.df)[1:3] <- c("gene","n.point","expression")
    wave.pattern.df$n.point <- as.character(wave.pattern.df$n.point)
    
    pt.df <- reshape2::melt(Res$pseudotime)
    pt.df <- pt.df %>% cbind(n.point=rownames(pt.df))
    colnames(pt.df)[1] <- "pseudotime"
    
    
    
    wave.pattern.df <- left_join(wave.pattern.df, pt.df, by=c("n.point"))
    
    wave.pattern.df <-wave.pattern.df %>% cbind(celltype=name)
    return(wave.pattern.df)
    
  }, Res=Res.list[names(Res.list)%in%celltype], name=name.list[names(name.list)%in%celltype], SIMPLIFY = F)
  
  wave.pattern <- do.call(rbind, wave.pattern.list)
  
  wave.pattern$celltype <- factor(wave.pattern$celltype, levels = unique(wave.pattern$celltype))
  wave.pattern$gene <- factor(wave.pattern$gene, levels = gene)
  wave.pattern$pseudotime <- as.character(wave.pattern$pseudotime)

  p <- ggplot(wave.pattern, aes(y=celltype, x=pseudotime))+
    geom_tile(aes(fill = expression),show.legend = NA, color=NA,width=5)+
    scale_fill_gradientn(colors = colorRampPalette(col_gdt)(100))+
    scale_x_discrete(breaks= c(0, 1), expand = c(0,0))+xlab("Pseudo-maturation")+ylab("Cell-type")+labs(fill="Expression")+
    scale_y_discrete(limits=rev)+
    mytheme_classic+theme(axis.line=element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(0.5, "points"))
  
  p <- p + facet_wrap(~gene, ncol=round(sqrt(length(gene))), scales = "fixed") +theme(strip.text = element_text(size = 14)) 
  
  
  
  print(p)
  
  
}




######################
#Plot pseudo-layer
######################



Gene_IT.pL <- sort(unique(unlist(lapply(Res_pL.list$IT, function(Res){rownames(Res$smooth.scale)}))))
Gene_ET.pL <- sort(unique(unlist(lapply(Res_pL.list$ET, function(Res){rownames(Res$smooth.scale)}))))
Gene_Sst.pL <- sort(unique(unlist(lapply(Res_pL.list$Sst, function(Res){rownames(Res$smooth.scale)}))))
Gene_Pvalb.pL <- sort(unique(unlist(lapply(Res_pL.list$Pvalb, function(Res){rownames(Res$smooth.scale)}))))
Gene_Vip.pL <- sort(unique(unlist(lapply(Res_pL.list$Vip, function(Res){rownames(Res$smooth.scale)}))))
Gene_Lamp5.pL <- sort(unique(unlist(lapply(Res_pL.list$Lamp5, function(Res){rownames(Res$smooth.scale)}))))

Gene.pL<- unique(c(Gene_IT.pL, Gene_ET.pL, Gene_Sst.pL, Gene_Pvalb.pL, Gene_Vip.pL, Gene_Lamp5.pL))



plotlyGene_smoothSpline_pL= function(Res.list, age=c("P8"), gene = NULL, free.scale = TRUE,  ncol = NA, 
                                     col_gdt=paletteer_d("grDevices::blues9"),
                                   line.size = 1, showlegend=F, height="850px", width="850px"){
  
  ########## Color gradient
  if(col_gdt=="Blues"){
    col_gdt =  blues9
  }
  else if(col_gdt=="Reds"){
    col_gdt =  brewer_pal(palette = "Reds")(9)
  }
  
  
  
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  
  #Set facet_wrap option
  if (is.na(ncol)) 
    nrow = round(sqrt(length(gene)))
  else nrow = NA
  a <- if (free.scale) 
    "free_x"
  else "fixed"  
  
  
  Res.list <- Res.list[names(Res.list)%in%age]
  
  name.list <- as.list(names(Res.list))
  
  wave.pattern.list <- mapply(function(Res, name){
    
    wave.pattern.df <- reshape2::melt(Res$smooth.scale[rownames(Res$smooth.scale)%in%gene,,drop=F])
    colnames(wave.pattern.df)[1:3] <- c("gene","n.point","expression")
    wave.pattern.df$n.point <- as.character(wave.pattern.df$n.point)
    
    pL.df <- reshape2::melt(Res$pseudolayer)
    pL.df <- pL.df %>% cbind(n.point=rownames(pL.df))
    colnames(pL.df)[1] <- "pseudolayer"

    
    wave.pattern.df <- left_join(wave.pattern.df, pL.df, by=c("n.point"))

    
    wave.pattern.df <-wave.pattern.df %>% cbind(age=name)
    return(wave.pattern.df)
    
  }, Res=Res.list, name=name.list, SIMPLIFY = F)
  
  wave.pattern <- do.call(rbind, wave.pattern.list)
  
  wave.pattern$age <- factor(wave.pattern$age, levels = unique(wave.pattern$age))
  wave.pattern$gene <- factor(wave.pattern$gene, levels = gene)


    color.age=colorRampPalette(col_gdt)(length(Res.list))


  
  p <- ggplot(wave.pattern, aes(x=pseudolayer, y=expression))+
    geom_line(aes(y = expression, x = pseudolayer,  color=age),orientation = "y", show.legend = showlegend, size=line.size)+
    scale_color_manual(values = color.age)+
    scale_x_reverse(breaks= c(0, 1))+ylab("Pseudo-layer")+xlab("Expression")+
    coord_flip()+
    mytheme_gray
  

    p <- p + facet_wrap(~gene, ncol = round(sqrt(length(gene))), scales = a) +theme(strip.text = element_text(size = 14))

  
    plot_height=  height + 50*length(unique(gene))
    plot_width=  width + 10*length(unique(gene))
    
    
    gp <-  ggplotly(p,  height=plot_height, width=plot_width) 
    
    gp <- gp %>% config(
      toImageButtonOptions = list(
        format = "svg"
      )
    )
    
    
  print(gp)
  
}



plotGene_Hm_pL= function(Res.list, age=c("P8"), gene = NULL, free.scale = TRUE,  ncol = NA, 
                                     col_gdt=paletteer_d("grDevices::blues9"),
                                     line.size = 1, showlegend=F){
  
  require(ggplot2)
  require(ggnewscale)
  require(dplyr)
  
  ########## Color gradient
  if(col_gdt=="PuOr"){
    col_gdt =  brewer_pal(palette = "PuOr", direction=-1)(11)
  }
  else if(col_gdt=="GBBr"){
    col_gdt =  brewer_pal(palette = "BrBG", direction=-1)(11)
  }
  
  else if(col_gdt=="BuYlRd"){
    col_gdt = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(col_gdt=="BuRd"){
    col_gdt = brewer_pal(palette = "RdBu", direction = -1)(11)
  }
  else if(col_gdt=="BuPhRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="GnPhRd"){
    col_gdt =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="BuGyOrRd"){
    col_gdt =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(col_gdt=="BuGyRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="Jet"){
    col_gdt =jet.colors(10)
  }
  else if(col_gdt=="BuSkyWhPkRd"){
    col_gdt =c("blue3", "skyblue","white", "pink", "red3")
  }
  else if(col_gdt=="Vik"){
    col_gdt =scico(11, palette="vik")
  }
     else if(col_gdt=="Purples_cb"){
  col_gdt = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  #Set facet_wrap option
  if (is.na(ncol)) 
    nrow = round(sqrt(length(gene)))
  else nrow = NA
  a <- if (free.scale) 
    "free_x"
  else "fixed"  
  
  
  Res.list <- Res.list[names(Res.list)%in%age]
  
  name.list <- as.list(names(Res.list))
  
  wave.pattern.list <- mapply(function(Res, name){
    
    wave.pattern.df <- reshape2::melt(Res$smooth.scale[rownames(Res$smooth.scale)%in%gene,,drop=F])
    colnames(wave.pattern.df)[1:3] <- c("gene","n.point","expression")
    wave.pattern.df$n.point <- as.character(wave.pattern.df$n.point)
    
    pL.df <- reshape2::melt(Res$pseudolayer)
    pL.df <- pL.df %>% cbind(n.point=rownames(pL.df))
    colnames(pL.df)[1] <- "pseudolayer"
    
    
    wave.pattern.df <- left_join(wave.pattern.df, pL.df, by=c("n.point"))
    
    
    wave.pattern.df <-wave.pattern.df %>% cbind(age=name)
    return(wave.pattern.df)
    
  }, Res=Res.list, name=name.list, SIMPLIFY = F)
  
  wave.pattern <- do.call(rbind, wave.pattern.list)
  
  wave.pattern$age <- factor(wave.pattern$age, levels = unique(wave.pattern$age))
  wave.pattern$gene <- factor(wave.pattern$gene, levels = gene)
  wave.pattern$pseudolayer <- as.character(wave.pattern$pseudolayer)
  
  

  
  
  p <- ggplot(wave.pattern, aes(y=pseudolayer, x=age))+
    geom_tile(aes(fill = expression),show.legend = NA, color=NA,height=5)+
    scale_fill_gradientn(colors = colorRampPalette(col_gdt)(100))+
    scale_y_discrete(breaks= c(min(as.numeric(wave.pattern$pseudolayer)),
                               max(as.numeric(wave.pattern$pseudolayer))),
                     limits=rev,
                     label=c(round(min(as.numeric(wave.pattern$pseudolayer)),2),
                             round(max(as.numeric(wave.pattern$pseudolayer)),2)))+xlab("Age")+ylab("Pseudo-layer")+labs(fill="Expression")+
    mytheme_classic+theme(axis.line=element_blank(), axis.ticks.x = element_blank(), axis.ticks.length = unit(0.5, "points"),
                          axis.text.x = element_text(angle=90))
  
  p <- p + facet_wrap(~gene, nrow=nrow, scales = "fixed") +theme(strip.text = element_text(size = 14)) 
  
  
  
  print(p)
  
  
}




#################################
#Plot Transcriptional landscape
#################################

Plot.Tlandscape = function(celltype, col_gdt="BuRd", gene=NULL, label_pM=c("E18.5","P30"), label_pL=c("L1","L6"), 
                           aspect.ratio=1, sig.gene=thr_name.family.list[[c("IT")]]){
  
  ########## Color gradient
  if(col_gdt=="PuOr"){
    col_gdt =  brewer_pal(palette = "PuOr", direction=-1)(11)
  }
  else if(col_gdt=="GBBr"){
    col_gdt =  brewer_pal(palette = "BrBG", direction=-1)(11)
  }
  
  else if(col_gdt=="BuYlRd"){
    col_gdt = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(col_gdt=="BuRd"){
    col_gdt = brewer_pal(palette = "RdBu", direction = -1)(11)
  }
  else if(col_gdt=="BuPhRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="GnPhRd"){
    col_gdt =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="BuGyOrRd"){
    col_gdt =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(col_gdt=="BuGyRd"){
    col_gdt =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col_gdt=="Jet"){
    col_gdt =jet.colors(10)
  }
  else if(col_gdt=="BuSkyWhPkRd"){
    col_gdt =c("blue3", "skyblue","white", "pink", "red3")
  }
  else if(col_gdt=="Vik"){
    col_gdt =scico(11, palette="vik")
  }
     else if(col_gdt=="Purples_cb"){
  col_gdt = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  
  ######## Gene selection
  if(celltype=="IT"){
    gene=gene[gene%in%Gene_IT.TL]
  }
  else if(celltype=="ET"){
    gene=gene[gene%in%Gene_ET.TL]
  }
  else if(celltype=="Sst"){
    gene=gene[gene%in%Gene_Sst.TL]
  }
  else if(celltype=="Pvalb"){
    gene=gene[gene%in%Gene_Pvalb.TL]
  }
  else if(celltype=="Vip"){
    gene=gene[gene%in%Gene_Vip.TL]
  }
  else if(celltype=="Lamp5"){
    gene=gene[gene%in%Gene_Lamp5.TL]
  }
  else if(celltype=="Other GlutNs"){
    gene=gene[gene%in%Gene_OtherGlutNs.TL]
  }
  else if(celltype=="Other GABANs"){
    gene=gene[gene%in%Gene_OtherGABANs.TL]
  }
  
  
  
  
  
  names.list <- as.list(gene)
  gene.list <- c("pMpLaxis_1","pMpLaxis_2",gene)
  names(gene.list) <- c("pMpLaxis_1", "pMpLaxis_2", gene)
  
  test <- lapply(gene.list, function(x){
    
    
    gene.landscape <- read_from_h5_DataFrameList("Data/landscape.h5", paste0(celltype,"/gene.landscape","/",x))
    gene.landscape <- as.data.frame(colData(gene.landscape))
    colnames(gene.landscape) <- x
    return(gene.landscape)
  })
  
  pList <- mapply(function(gene, name){
    
    pMpL <- do.call(cbind, test[c(1:2)])
    landscape.df <- cbind(pMpL, gene)
    colnames(landscape.df)[3] <- "expression"
    
    if(any(sig.gene==name)){
    p <- ggplot(landscape.df, aes(x=pMpLaxis_1, y=pMpLaxis_2))+
      geom_tile(aes(fill = ifelse(expression>quantile(expression,0.99),quantile(expression,0.99),expression) ),show.legend = NA,size=0)+
      scale_fill_gradientn(colors = col_gdt, 
                           limits=c(min(landscape.df$expression),max(landscape.df$expression)),
                           breaks=c(min(landscape.df$expression),max(landscape.df$expression)),
                           labels=c("low","high"))+
      scale_x_continuous(limits=c(0,1),breaks= c(0, 1), labels = label_pM)+xlab("Pseudo-maturation")+
      scale_y_reverse(limits=c(1,0),breaks= c(0, 1), labels = label_pL)+
      ylab("Pseudo-layer")+labs(fill="Expression")+
      geom_rangeframe()+
      theme_tufte(base_size = 14, base_family = "") +
      theme( aspect.ratio = aspect.ratio,legend.key.height = unit(0.5, "cm"),title = element_text( face = "bold.italic"),
             axis.title = element_text( face = "plain"),legend.title = element_text( face = "plain"),)+
      ggtitle(paste(colnames(gene),"_",celltype)) +
      ggeasy::easy_center_title()
    }
    
    else  {
      p <- ggplot(landscape.df, aes(x=pMpLaxis_1, y=pMpLaxis_2))+
        geom_tile(aes(fill = ifelse(expression>quantile(expression,0.99),quantile(expression,0.99),expression) ),show.legend = NA, size=0)+
        scale_fill_gradientn(colors = colorRampPalette("grey")(11), 
                             limits=c(min(landscape.df$expression),max(landscape.df$expression)),
                             breaks=c(min(landscape.df$expression),max(landscape.df$expression)),
                             labels=c("low","high"))+
        scale_x_continuous(limits=c(0,1),breaks= c(0, 1), labels = label_pM)+xlab("Pseudo-maturation")+
        scale_y_reverse(limits=c(1,0),breaks= c(0, 1), labels = label_pL)+
        annotate("text",label="NS Expr.",x=0.5,y=0.5,size=NA)+
        ylab("Pseudo-layer")+labs(fill="Expression")+
        geom_rangeframe()+
        theme_tufte(base_size = 14, base_family = "") +
        theme( aspect.ratio = aspect.ratio,legend.key.height = unit(0.5, "cm"), title = element_text(face = "bold.italic"),
               axis.title = element_text( face = "plain"),legend.title = element_text( face = "plain"))+
      ggtitle(colnames(gene)) +
        ggeasy::easy_center_title()
      
      
    }
    
    return(p) 
  },gene=test[c(3:length(test))],name=names.list, SIMPLIFY=F)
  
  print(x = cowplot::plot_grid(plotlist = pList, ncol = round(sqrt(length(gene)))))
  
  
}



###########################################
#Plot number of interactions
###########################################


filter_function <- function(data, field, operator, type, input_item, selected_item, NA_included){
  
  if(type == "string"){
    selected_tmp <- data[,field] %in% input_item
  }
  if(type == "numeric"){
    selected_tmp <- data[,field] > input_item
  }
  if(type=="split"){
    
    selected_tmp <- rownames(data) %in% rownames(data[unlist(sapply(input_item, str_detect, data[,field], USE.NAMES = F)),])

  }
  
  selected_tmp[is.na(selected_tmp)] <- NA_included
  
  if(tolower(operator) == "and"){
    selected_item <- selected_item & selected_tmp
  }
  if(tolower(operator) == "or"){
    selected_item <- selected_item | selected_tmp
  }
  
  return(selected_item)
}



scSeqComm_select_bis <- function(data, operator = "and",
                                 ligand = NULL, receptor = NULL, LR_pair = NULL,
                                 cluster_L = NULL, cluster_R = NULL, interaction = NULL, 
                                 L_score_S_lr = NULL, R_score_S_lr = NULL, S_inter = NULL, SOURCE.class=NULL, TARGET.class=NULL,
                                 mean.count_L = NULL, mean.count_R = NULL, mean.product = NULL, regularized.product = NULL, detection.rate_L=NULL,
                                 detection.rate_R=NULL,numdup.LR_pair.interaction=NULL, numdup.LR_pair.cluster_L=NULL, numdup.LR_pair.cluster_R=NULL,
                                 Criterion.Thr.mean.product.bis=NULL,
                                 Zhou_score= NULL, Skelly_score= NULL, S_inter_diff= NULL,Criterion.Thr.mean.product=NULL,Criterion.Thr.mean.product.w=NULL,
                                 Criterion.Thr.S_inter.w=NULL,
                                 S_intra = NULL, pathway = NULL,Ligand.category=NULL, Receptor.category=NULL,Ligand.family=NULL, Receptor.family=NULL,
                                 Criterion.Thr=NULL,Criterion.Thr.ct_pair=NULL,Criterion.Thr.ct_pair.w_mean=NULL,Criterion.Thr.mean.product.w_mean=NULL,
                                 Criterion.Thr.mean.product.nct=NULL,Criterion.Thr.S_inter.w_mean=NULL,
                                 NA_included = TRUE, output_field = NULL){
  
  # if no filter is selected, return all rows
  if(is.null(ligand)&is.null(receptor)&is.null(LR_pair)&is.null(cluster_L)&is.null(cluster_R)&is.null(interaction)&
     is.null(L_score_S_lr)&is.null(R_score_S_lr)&is.null(S_inter)&is.null(Criterion.Thr)&is.null(Criterion.Thr.ct_pair)&is.null(Criterion.Thr.ct_pair.w_mean)&
     is.null(Ligand.category)&is.null(Receptor.category)&is.null(Ligand.family)&is.null(Receptor.family)&is.null(Criterion.Thr.mean.product)&
     is.null(SOURCE.class)&is.null(TARGET.class)&is.null(Criterion.Thr.mean.product.w)&is.null(Criterion.Thr.mean.product.bis)&
     is.null(Criterion.Thr.mean.product.w_mean)&is.null(Criterion.Thr.mean.product.nct)&is.null(Criterion.Thr.S_inter.w_mean)&
     is.null(numdup.LR_pair.interaction)&is.null(numdup.LR_pair.cluster_L)&is.null(numdup.LR_pair.cluster_R)&
     is.null(mean.count_L)&is.null(mean.count_R)&is.null(mean.product)&is.null(regularized.product)&is.null(S_inter_diff)&is.null(Zhou_score)&
     is.null(Skelly_score)&is.null(Criterion.Thr.S_inter.w)&is.null(detection.rate_L)&is.null(detection.rate_R)&
     is.null(S_intra)&is.null(pathway)){
    selected_row <- rep(TRUE, nrow(data))
  }else{# if at least one filter is selected
    if(tolower(operator) == "and"){
      selected_row <- rep(TRUE, nrow(data))
    }
    if(tolower(operator) == "or"){
      selected_row <- rep(FALSE, nrow(data))
    }
    
    # string field
    if(!is.null(ligand)){
      selected_row <- filter_function(data = data, field = "ligand", operator = operator, type = "string", input_item = ligand, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(receptor)){
      selected_row <- filter_function(data = data, field = "receptor", operator = operator, type = "string", input_item = receptor, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(LR_pair)){
      selected_row <- filter_function(data = data, field = "LR_pair", operator = operator, type = "string", input_item = LR_pair, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(cluster_L)){
      selected_row <- filter_function(data = data, field = "cluster_L", operator = operator, type = "string", input_item = cluster_L, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(cluster_R)){
      selected_row <- filter_function(data = data, field = "cluster_R", operator = operator, type = "string", input_item = cluster_R, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(interaction)){
      selected_row <- filter_function(data = data, field = "interaction", operator = operator, type = "string", input_item = interaction, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(pathway)){
      selected_row <- filter_function(data = data, field = "pathway", operator = operator, type = "string", input_item = pathway, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr", operator = operator, type = "string", input_item = Criterion.Thr, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.ct_pair)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.ct_pair", operator = operator, type = "string", input_item = Criterion.Thr.ct_pair, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.ct_pair.w_mean)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.ct_pair.w_mean", operator = operator, type = "string", input_item = Criterion.Thr.ct_pair.w_mean, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(SOURCE.class)){
      selected_row <- filter_function(data = data, field = "SOURCE.class", operator = operator, type = "string", input_item = SOURCE.class, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(TARGET.class)){
      selected_row <- filter_function(data = data, field = "TARGET.class", operator = operator, type = "string", input_item = TARGET.class, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Zhou_score)){
      selected_row <- filter_function(data = data, field = "Zhou_score", operator = operator, type = "string", input_item = Zhou_score, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Skelly_score)){
      selected_row <- filter_function(data = data, field = "Skelly_score", operator = operator, type = "string", input_item = Skelly_score, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.mean.product)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.mean.product", operator = operator, type = "string", input_item = Criterion.Thr.mean.product, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.mean.product.w)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.mean.product.w", operator = operator, type = "string", input_item = Criterion.Thr.mean.product.w, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.S_inter.w)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.S_inter.w", operator = operator, type = "string", input_item = Criterion.Thr.S_inter.w, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.mean.product.bis)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.mean.product.bis", operator = operator, type = "string", input_item = Criterion.Thr.mean.product.bis, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.mean.product.w_mean)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.mean.product.w_mean", operator = operator, type = "string", input_item = Criterion.Thr.mean.product.w_mean, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.mean.product.nct)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.mean.product.nct", operator = operator, type = "string", input_item = Criterion.Thr.mean.product.nct, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Criterion.Thr.S_inter.w_mean)){
      selected_row <- filter_function(data = data, field = "Criterion.Thr.S_inter.w_mean", operator = operator, type = "string", input_item = Criterion.Thr.S_inter.w_mean, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(numdup.LR_pair.interaction)){
      selected_row <- filter_function(data = data, field = "numdup.LR_pair.interaction", operator = operator, type = "string", input_item = numdup.LR_pair.interaction, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(numdup.LR_pair.cluster_L)){
      selected_row <- filter_function(data = data, field = "numdup.LR_pair.cluster_L", operator = operator, type = "string", input_item = numdup.LR_pair.cluster_L, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(numdup.LR_pair.cluster_R)){
      selected_row <- filter_function(data = data, field = "numdup.LR_pair.cluster_R", operator = operator, type = "string", input_item = numdup.LR_pair.cluster_R, selected_item = selected_row, NA_included = NA_included)
    }
    
    #split field
    if(!is.null(Ligand.category)){
      selected_row <- filter_function(data = data, field = "Ligand.category", operator = operator, type = "split", input_item = Ligand.category, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Receptor.category)){
      selected_row <- filter_function(data = data, field = "Receptor.category", operator = operator, type = "split", input_item = Receptor.category, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Ligand.family)){
      selected_row <- filter_function(data = data, field = "Ligand.family", operator = operator, type = "split", input_item = Ligand.family, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(Receptor.family)){
      selected_row <- filter_function(data = data, field = "Receptor.family", operator = operator, type = "split", input_item = Receptor.family, selected_item = selected_row, NA_included = NA_included)
    }
    
    
    
    # numeric field
    if(!is.null(S_inter)){
      selected_row <- filter_function(data = data, field = "S_inter", operator = operator, type = "numeric", input_item = S_inter, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(S_intra)){
      selected_row <- filter_function(data = data, field = "S_intra", operator = operator, type = "numeric", input_item = S_intra, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(L_score_S_lr)){
      selected_row <- filter_function(data = data, field = "L_score_S_lr", operator = operator, type = "numeric", input_item = L_score_S_lr, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(R_score_S_lr)){
      selected_row <- filter_function(data = data, field = "R_score_S_lr", operator = operator, type = "numeric", input_item = R_score_S_lr, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(mean.count_L)){
      selected_row <- filter_function(data = data, field = "mean.count_L", operator = operator, type = "numeric", input_item = mean.count_L, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(mean.count_R)){
      selected_row <- filter_function(data = data, field = "mean.count_R", operator = operator, type = "numeric", input_item = mean.count_R, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(mean.product)){
      selected_row <- filter_function(data = data, field = "mean.product", operator = operator, type = "numeric", input_item = mean.product, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(regularized.product)){
      selected_row <- filter_function(data = data, field = "regularized.product", operator = operator, type = "numeric", input_item = regularized.product, selected_item = selected_row, NA_included = NA_included)
    }
    
    if(!is.null(S_inter_diff)){
      selected_row <- filter_function(data = data, field = "S_inter_diff", operator = operator, type = "numeric", input_item = S_inter_diff, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(detection.rate_L)){
      selected_row <- filter_function(data = data, field = "detection.rate_L", operator = operator, type = "numeric", input_item = detection.rate_L, selected_item = selected_row, NA_included = NA_included)
    }
    if(!is.null(detection.rate_R)){
      selected_row <- filter_function(data = data, field = "detection.rate_R", operator = operator, type = "numeric", input_item = detection.rate_R, selected_item = selected_row, NA_included = NA_included)
    }
    
  }
  
  # select output columns
  output_col <- colnames(data)
  if(!is.null(output_field)){
    output_col <- output_field
  }
  
  return(data[selected_row,output_col])
}



scSeqComm_heatmaply_cardinality <- function(full_data,data, y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                            ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",   cl.order=order_celltype,
                                            key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,30,0), sizefont=10, 
                                            height="850px", main="E18.5-P0", scale="none",
                                            Sending.cl=unique(P0wt_edges_sig$Sending.cluster),
                                            Target.cl=unique(P0wt_edges_sig$Target.cluster),SOURCE.CT=unique(P0wt_edges_sig$SOURCE.Cell.type), 
                                            TARGET.CT=unique(P0wt_edges_sig$TARGET.Cell.type),
                                            col.heatmap=c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"))
{
  #full data filter only what you don't want to be displayed
  require(heatmaply)
  require(reshape2)
  require(dplyr)
  

  
  ########### color gradient
  if(col.heatmap=="Sevilla"){
    col.heatmap =  c("#FFFFFF", "#D6D6D7", "#C7B0B0","#C98F8F", "#CA6E6E", "#C64D4D", "#9F2D2D", "#761A1A", "#4C1111", "#270A0A", "#000000" )
  }
  else if(col.heatmap=="Oslo"){
    col.heatmap =  scico(11,palette = "oslo",direction=-1)
  }
  else if(col.heatmap=="PuOr"){
    col.heatmap =  brewer_pal(palette = "PuOr", direction=-1)(11)
  }
  else if(col.heatmap=="GBBr"){
    col.heatmap =  brewer_pal(palette = "BrBG", direction=-1)(11)
  }
  
  else if(col.heatmap=="BuYlRd"){
    col.heatmap = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(col.heatmap=="BuPhRd"){
    col.heatmap =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="GnPhRd"){
    col.heatmap =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="BuGyOrRd"){
    col.heatmap =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(col.heatmap=="BuGyRd"){
    col.heatmap =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="Jet"){
    col.heatmap =jet.colors(10)
  }
  else if(col.heatmap=="Rainbow"){
    col.heatmap =rainbow(5, rev=T)
  }
  else if(col.heatmap=="Heat"){
    col.heatmap =heat.colors(10, rev = T)
  }
  else if(col.heatmap=="viridis"){
    col.heatmap =viridis_pal(option = "viridis")(11)
  }
  else if(col.heatmap=="magma"){
    col.heatmap =viridis_pal(option = "magma")(11)
  }
  else if(col.heatmap=="inferno"){
    col.heatmap =viridis_pal(option = "inferno")(11)
  }
  else if(col.heatmap=="plasma"){
    col.heatmap =viridis_pal(option = "plasma")(11)
  }
  else if(col.heatmap=="cividis"){
    col.heatmap =viridis_pal(option = "cividis")(11)
  }
    else if(col.heatmap=="Purples_cb"){
  col.heatmap = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  #select variables from data
  data <- distinct(data, ligand, receptor, interaction, .keep_all = T)
  data1 <- data[,c(x,y)]
  
  #count the number of observations of x-y pair
  df <- data1 %>% dplyr::count(data1[[x]],data1[[y]])
  colnames(df) <- c(x,y,"cardinality")
  
  
  #Whole dataset
  full_data <- full_data[,c("cluster_L","cluster_R","SOURCE.class","TARGET.class")]
  
  full_data <- left_join(full_data, df, by=c("cluster_L","cluster_R"))
  full_data <-full_data %>% unique()
  
  #Convert into matrix
  
  count_matrix <- acast(full_data, cluster_L ~ cluster_R, value.var = "cardinality")
  count_matrix[is.na(count_matrix)==T] <- 0
  count_matrix= count_matrix[cl.order[cl.order$class_label%in%unique(full_data$SOURCE.class),1], cl.order[cl.order$class_label%in%unique(full_data$TARGET.class),1]]
  
  #Subset cluster and cell type
  SOURCE.CT=cl.order$celltype_label[cl.order$class_label%in%SOURCE.CT]
  TARGET.CT=cl.order$celltype_label[cl.order$class_label%in%TARGET.CT]
  Sending.cl=intersect(Sending.cl, SOURCE.CT)
  Target.cl=intersect(Target.cl, TARGET.CT)
  count_matrix=count_matrix[rownames(count_matrix)%in%Sending.cl, colnames(count_matrix)%in%Target.cl, drop=F]
  
  #plotheatmap
  col.heatmap <- colorRampPalette(col.heatmap)(300)
  
  
  ht <- heatmaply(count_matrix, dendrogram = "none", grid_gap = 2, column_text_angle = 45, label_names = c("SOURCE", "TARGET", hover_anno),
                  colors = col.heatmap, main=main, scale=scale, 
                  col_side_colors = cl.order[cl.order$celltype_label%in%colnames(count_matrix),2:3],
                  row_side_colors = cl.order[cl.order$celltype_label%in%rownames(count_matrix),2:3], key.title = key.title,
                  hide_colorbar = hide_colorbar,  margins = margins, fontsize_col = sizefont, fontsize_row = sizefont,
                  xlab= "Cell-type expressing Receptor (TARGET)", ylab= "Cell-type expressing Ligand (SOURCE)",
                  heatmap_layers = list( theme(axis.line = element_blank())),
                  col_side_palette = colorRampPalette(c("grey95","black")), row_side_palette = colorRampPalette(c("grey95","black"))) %>% 
    plotly::layout( showlegend=F,
                    annotations = list(
                      visible = FALSE
                    ),
                    autosize=T, 
                    height=height)
  
  #Arrange the annotations
  ht$x$layout$xaxis2$showticklabels <- "" #Discard legend of annotations
  ht$x$layout$yaxis$showticklabels <- "" #Discard legend of annotations
  ht$x$layout$xaxis$ticks <- "" #discard ticks of x
  ht$x$layout$yaxis2$ticks <- "" #discard ticks of y 
  ht$x$layout$xaxis2$ticktext <- c("","") #Discard names of row annotations
  ht$x$layout$yaxis$ticktext <- c("","") #Discard names of col annotations
  
  #Change the cell type color and the text display on hover for side annotations
  for (i in 1:length(ht$x$data)) {
    
    ht$x$data[[i]]$fillcolor[ht$x$data[[i]]$text%like%c("GABAergic")] <- "rgba(216,13,0,1)"
    ht$x$data[[i]]$fillcolor[ht$x$data[[i]]$text%like%c("Glutamatergic")] <- "rgba(17,68,227,1)"
    
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] <- gsub("value", "Normalized soma depth",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] <- gsub("value", "Class",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("row")] <- gsub("row", "SOURCE",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("row")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("column")] <- gsub("column", "TARGET",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("column")] )
  }
  
  for (i in 1:length(ht$x$data)) {
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] <- gsub("<br />variable: Normalized.soma.depth","",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")])
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] <- gsub("<br />variable: class_label","",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")])
    
  }
  
  #Change the hover order for heatmap
  

  #Adjust the colorscale in case of only one value in heatmap
  colorscale.df=as.data.frame(matrix(c(seq(0,max(count_matrix), length.out = 5), rep(col.heatmap[length(col.heatmap)],5)), nrow=5, ncol=2))
  colnames(colorscale.df)=NULL
  
  
  if(length(unique(ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+2]]$colorscale[,2]))==1){
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+2]]$colorscale=colorscale.df
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorscale[,1] <- seq(0,max(count_matrix),length.out = 300) #range of values of colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorscale[,2] <- col.heatmap #range of color of colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorbar$tickvals=seq(0,max(count_matrix), length.out = 5) #ticks of the colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorbar$ticktext="0" #Set the first tick value of the colorscale to 0
  }
  
  ht <- ht %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  return(ht)
  
}






###### contiguous Hm from E18.5-P0 to P30
prepare_contigHm <- function(comm){
  
  comm <- dplyr::distinct(comm, ligand, receptor,interaction,age.at.collection,.keep_all = T)
  
  comm1 <- comm[comm$age.at.collection%in%c("E18.5-P0","P1-P2"),]
  intersect_numdup1 <- plyr::ddply(comm1, c("interaction","LR_pair"),nrow)
  
  
  comm1 <-left_join(comm1,intersect_numdup1,by=c("interaction","LR_pair")) 
  
  
  
  comm2 <- comm[comm$age.at.collection%in%c("P1-P2", "P4-P5"),]
  intersect_numdup2 <- plyr::ddply(comm2, c("interaction","LR_pair"),nrow)
  
  comm2 <-left_join(comm2,intersect_numdup2,by=c("interaction","LR_pair")) 
  
  
  comm3 <- comm[comm$age.at.collection%in%c("P4-P5", "P8"),]
  intersect_numdup3 <- plyr::ddply(comm3, c("interaction","LR_pair"),nrow)
  
  comm3 <-left_join(comm3,intersect_numdup3,by=c("interaction","LR_pair")) 
  
  
  
  comm4 <- comm[comm$age.at.collection%in%c("P8", "P16"),]
  intersect_numdup4 <- plyr::ddply(comm4, c("interaction","LR_pair"),nrow)
  
  comm4 <-left_join(comm4,intersect_numdup4,by=c("interaction","LR_pair")) 
  
  
  
  comm5 <- comm[comm$age.at.collection%in%c( "P16", "P30"),]
  intersect_numdup5 <- plyr::ddply(comm5, c("interaction","LR_pair"),nrow)
  
  
  comm5 <-left_join(comm5,intersect_numdup5,by=c("interaction","LR_pair")) 
  
  
  comm6 <- comm[comm$age.at.collection%in%c( "P30", "Adult"),]
  intersect_numdup6 <- plyr::ddply(comm6, c("interaction","LR_pair"),nrow)
  
  
  comm6 <-left_join(comm6,intersect_numdup6,by=c("interaction","LR_pair")) 
  
  
  
  comm_contig.unique <- comm1 %>% rbind(comm2) %>% rbind(comm3) %>% 
    rbind(comm4) %>% rbind(comm5) %>% 
    rbind(comm6)
  comm_contig.unique <- comm_contig.unique[comm_contig.unique$V1%in%c("2") & !(comm_contig.unique$age.at.collection%in%c("Adult")),]
  comm_contig.unique <- dplyr::distinct(comm_contig.unique, ligand, receptor,interaction,.keep_all = T)
  
  return(comm_contig.unique)
  
  
}




scSeqComm_contigheatmaply_cardinality <- function(full_data,data, y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                            ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",   cl.order=order_celltype,
                                            key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,30,0), sizefont=10, 
                                            height="850px", main="E18.5-P0", scale="none",
                                            Sending.cl=unique(P0wt_edges_sig$Sending.cluster),
                                            Target.cl=unique(P0wt_edges_sig$Target.cluster),SOURCE.CT=unique(P0wt_edges_sig$SOURCE.Cell.type), 
                                            TARGET.CT=unique(P0wt_edges_sig$TARGET.Cell.type),
                                            col.heatmap=c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"))
{
  #full data filter only what you don't want to be displayed
  require(heatmaply)
  require(reshape2)
  require(dplyr)
  
  
  
  ########### color gradient
  if(col.heatmap=="Sevilla"){
    col.heatmap =  c("#FFFFFF", "#D6D6D7", "#C7B0B0","#C98F8F", "#CA6E6E", "#C64D4D", "#9F2D2D", "#761A1A", "#4C1111", "#270A0A", "#000000" )
  }
  else if(col.heatmap=="Oslo"){
    col.heatmap =  scico(11,palette = "oslo",direction=-1)
  }
  else if(col.heatmap=="PuOr"){
    col.heatmap =  brewer_pal(palette = "PuOr", direction=-1)(11)
  }
  else if(col.heatmap=="GBBr"){
    col.heatmap =  brewer_pal(palette = "BrBG", direction=-1)(11)
  }
  
  else if(col.heatmap=="BuYlRd"){
    col.heatmap = brewer_pal(palette = "RdYlBu", direction = -1)(11)
  }
  else if(col.heatmap=="BuPhRd"){
    col.heatmap =c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="GnPhRd"){
    col.heatmap =c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="BuGyOrRd"){
    col.heatmap =c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500")
  }
  else if(col.heatmap=="BuGyRd"){
    col.heatmap =c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF")
  }
  else if(col.heatmap=="Jet"){
    col.heatmap =jet.colors(10)
  }
  else if(col.heatmap=="Rainbow"){
    col.heatmap =rainbow(5, rev=T)
  }
  else if(col.heatmap=="Heat"){
    col.heatmap =heat.colors(10, rev = T)
  }
  else if(col.heatmap=="viridis"){
    col.heatmap =viridis_pal(option = "viridis")(11)
  }
  else if(col.heatmap=="magma"){
    col.heatmap =viridis_pal(option = "magma")(11)
  }
  else if(col.heatmap=="inferno"){
    col.heatmap =viridis_pal(option = "inferno")(11)
  }
  else if(col.heatmap=="plasma"){
    col.heatmap =viridis_pal(option = "plasma")(11)
  }
  else if(col.heatmap=="cividis"){
    col.heatmap =viridis_pal(option = "cividis")(11)
  }
  else if(col.heatmap=="Purples_cb"){
  col.heatmap = c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d", "#2d004b")
}
  
  #select variables from data
  data <- prepare_contigHm(data)
  data1 <- data[,c(x,y)]
  
  #count the number of observations of x-y pair
  df <- data1 %>% dplyr::count(data1[[x]],data1[[y]])
  colnames(df) <- c(x,y,"cardinality")
  
  
  #Whole dataset
  full_data <- full_data[,c("cluster_L","cluster_R","SOURCE.class","TARGET.class")]
  
  full_data <- left_join(full_data, df, by=c("cluster_L","cluster_R"))
  full_data <-full_data %>% unique()
  
  #Convert into matrix
  
  count_matrix <- acast(full_data, cluster_L ~ cluster_R, value.var = "cardinality")
  count_matrix[is.na(count_matrix)==T] <- 0
  count_matrix= count_matrix[cl.order[cl.order$class_label%in%unique(full_data$SOURCE.class),1], cl.order[cl.order$class_label%in%unique(full_data$TARGET.class),1]]
  
  #Subset cluster and cell type
  SOURCE.CT=cl.order$celltype_label[cl.order$class_label%in%SOURCE.CT]
  TARGET.CT=cl.order$celltype_label[cl.order$class_label%in%TARGET.CT]
  Sending.cl=intersect(Sending.cl, SOURCE.CT)
  Target.cl=intersect(Target.cl, TARGET.CT)
  count_matrix=count_matrix[rownames(count_matrix)%in%Sending.cl, colnames(count_matrix)%in%Target.cl, drop=F]
  
  #plotheatmap
  col.heatmap <- colorRampPalette(col.heatmap)(300)
  
  
  ht <- heatmaply(count_matrix, dendrogram = "none", grid_gap = 2, column_text_angle = 45, label_names = c("SOURCE", "TARGET", hover_anno),
                  colors = col.heatmap, main=main, scale=scale, 
                  col_side_colors = cl.order[cl.order$celltype_label%in%colnames(count_matrix),2:3],
                  row_side_colors = cl.order[cl.order$celltype_label%in%rownames(count_matrix),2:3], key.title = key.title,
                  hide_colorbar = hide_colorbar,  margins = margins, fontsize_col = sizefont, fontsize_row = sizefont,
                  xlab= "Cell-type expressing Receptor (TARGET)", ylab= "Cell-type expressing Ligand (SOURCE)",
                  heatmap_layers = list( theme(axis.line = element_blank())),
                  col_side_palette = colorRampPalette(c("grey95","black")), row_side_palette = colorRampPalette(c("grey95","black"))) %>% 
    plotly::layout( showlegend=F,
                    annotations = list(
                      visible = FALSE
                    ),
                    autosize=T, 
                    height=height)
  
  #Arrange the annotations
  ht$x$layout$xaxis2$showticklabels <- "" #Discard legend of annotations
  ht$x$layout$yaxis$showticklabels <- "" #Discard legend of annotations
  ht$x$layout$xaxis$ticks <- "" #discard ticks of x
  ht$x$layout$yaxis2$ticks <- "" #discard ticks of y 
  ht$x$layout$xaxis2$ticktext <- c("","") #Discard names of row annotations
  ht$x$layout$yaxis$ticktext <- c("","") #Discard names of col annotations
  
  #Change the cell type color and the text display on hover for side annotations
  for (i in 1:length(ht$x$data)) {
    
    ht$x$data[[i]]$fillcolor[ht$x$data[[i]]$text%like%c("GABAergic")] <- "rgba(216,13,0,1)"
    ht$x$data[[i]]$fillcolor[ht$x$data[[i]]$text%like%c("Glutamatergic")] <- "rgba(17,68,227,1)"
    
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] <- gsub("value", "Normalized soma depth",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] <- gsub("value", "Class",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("row")] <- gsub("row", "SOURCE",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("row")] )
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("column")] <- gsub("column", "TARGET",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("column")] )
  }
  
  for (i in 1:length(ht$x$data)) {
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")] <- gsub("<br />variable: Normalized.soma.depth","",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("Normalized.soma.depth")])
    ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")] <- gsub("<br />variable: class_label","",ht$x$data[[i]]$text[ht$x$data[[i]]$text%like%c("class_label")])
    
  }
  
  #Change the hover order for heatmap
  

  
  #Adjust the colorscale in case of only one value in heatmap
  colorscale.df=as.data.frame(matrix(c(seq(0,max(count_matrix), length.out = 5), rep(col.heatmap[length(col.heatmap)],5)), nrow=5, ncol=2))
  colnames(colorscale.df)=NULL
  
  
  if(length(unique(ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+2]]$colorscale[,2]))==1){
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+2]]$colorscale=colorscale.df
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorscale[,1] <- seq(0,max(count_matrix),length.out = 300) #range of values of colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorscale[,2] <- col.heatmap #range of color of colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorbar$tickvals=seq(0,max(count_matrix), length.out = 5) #ticks of the colorscale
    ht$x$data[[2*nrow(cl.order[which(cl.order$celltype_label%in%colnames(count_matrix) | cl.order$celltype_label%in%rownames(count_matrix)),])+3]]$marker$colorbar$ticktext="0" #Set the first tick value of the colorscale to 0
  }
  
  ht <- ht %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  return(ht)
  
}



################################################################################
#
#                             Plot DotPlot LR
#
################################################################################


summaryze_S_intra=function (data) 
{
  return(as.data.frame(data %>% dplyr::group_by(interaction, 
                                                LR_pair) %>% dplyr::arrange(dplyr::desc(S_intra)) %>% 
                         dplyr::slice(1) %>% dplyr::ungroup()))
}

scSeqComm_heatmaply=function (data, x, y, fill = "S_intra", size = "S_inter",  height,width,
                              xlab="", ylab = "",  gdt_pal="lajolla", gap_width=40, gap_height=130,
                              limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                          0.5, 0.8), 
                              facet_grid_x = NULL, facet_grid_y = NULL, 
                              annotation_GO = NULL, cutoff = 0.05, topGO = 5, GO_ncol = 1) 
{
  
  if(gdt_pal=="lajolla"){
    gdt_pal =  scico(11,palette = "lajolla")
  }
  else if(gdt_pal=="bilbao"){
    gdt_pal = scico(11,palette = "bilbao")
  }
  else if(gdt_pal=="acton"){
    gdt_pal = scico(11,palette = "acton", direction = -1)
  }
  else if(gdt_pal=="batlowK"){
    gdt_pal = scico(11,palette = "batlowK", direction = -1)
  }
  else if(gdt_pal=="davos"){
    gdt_pal = scico(11, palette="davos", direction=-1)
  }
  else if(gdt_pal=="hawaii"){
    gdt_pal =scico(11, palette="hawaii", direction=-1)
  }
  else if(gdt_pal=="imola"){
    gdt_pal =scico(11,palette="imola", direction=-1)
  }
  else if(gdt_pal=="lapaz"){
    gdt_pal =scico(11,palette="lapaz", direction=-1)
  }
  else if(gdt_pal=="nuuk"){
    gdt_pal =scico(11,palette="nuuk", direction=-1)
  }
  else if(gdt_pal=="tokyo"){
    gdt_pal =scico(11,palette="tokyo", direction=-1)
  }
  else if(gdt_pal=="turku"){
    gdt_pal =scico(11,palette="turku", direction=-1)
  }
  
  if(facet_grid_x=="NULL"){
    facet_grid_x= NULL
  }
  
  
  if(facet_grid_y=="NULL"){
    facet_grid_y= NULL
  }
  
    
  
  if((x=="cluster_L" & y=="cluster_R") | (x=="cluster_R" & y=="cluster_L")){
    
    data <- summaryze_S_intra(data)
    

    data <- data %>% dplyr::arrange(S_inter)
    
    plot_width= width + gap_width/2*length(unique(data$cluster_L))
    plot_height= height+ gap_height/2*length(unique(data$cluster_R))
    
    
  }
  
  else{
    
    
    data$cluster_L <- factor(data$cluster_L, levels=order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_L)])
    data$cluster_R <- factor(data$cluster_R, levels=order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_L)])
    
    data <- data %>% dplyr::arrange(S_inter)
    
    plot_height= height + gap_height*length(unique(data$LR_pair))
    plot_width= width + gap_width*length(unique(data$interaction))
    
  }
  
 
  if(x=="cluster_R" & y=="cluster_L"){
    
    data$cluster_L <- factor(data$cluster_L, levels=rev(order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_L)]))
    data$cluster_R <- factor(data$cluster_R, levels=order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_R)])
  }
  if(x=="cluster_L" & y=="cluster_R"){

    data$cluster_L <- factor(data$cluster_L, levels=order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_L)])
    data$cluster_R <- factor(data$cluster_R, levels=rev(order_celltype$celltype_label[order_celltype$celltype_label%in%unique(data$cluster_R)]))
  }
  
  if (is.null(limit_fill)) {
    limit_fill <- c(min(data[[fill]]), max(data[[fill]]))
  }
  if (is.null(limit_size) & !is.null(size)) {
    limit_size <- c(min(data[[size]]), max(data[[size]]))
  }
  if (!is.null(annotation_GO)) {
    data$LR_pair <- factor(data$LR_pair, levels = unique(data$LR_pair[order(data$receptor)]))
  }
  if (!is.null(size)) {
    g <- ggplot(data, aes_string(x = x, y = y, col = fill, 
                                 size = size)) + 
      geom_point(aes(text = paste('SOURCE --> TARGET:', interaction,
                                  '<br> Ligand - Receptor:', LR_pair,
                                  '<br> S_inter:', format(round(S_inter,3), nsmall = 3),
                                  '<br> S_intra:', format(round(S_intra,3), nsmall = 3),
                                  '<br> S_inter_diff:', format(round(S_inter_diff,3), nsmall = 3),
                                  '<br> pathway:', pathway,
                                  '<br> Ligand category:', Ligand.category,
                                  '<br> Receptor category:', Receptor.category,
                                  '<br> Ligand family:', Ligand.family,
                                  '<br> Receptor family:', Receptor.family))) + 
      labs(x = xlab, y = ylab) + 
      scale_size_binned(range = c(0.01, 5), breaks = breaks_size, 
                        trans = "exp", limits = limit_size) + 
      add2ggplot::theme_white() + 
      scale_color_gradientn(colors =  gdt_pal, limits = limit_fill, 
                               na.value = "grey")
  }
  else {
    g <- ggplot(data, aes_string(x = x, y = y, fill = fill)) + 
      labs(x = xlab, y = ylab, fill = fill) + 
      geom_tile(colour = "white", size = 0.2) + scale_fill_gradientn(colors = gdt_pal, 
                                                                        limits = limit_fill, na.value = "grey")+
      scale_x_discrete(position="top")
  }
  g <- g + add2ggplot::theme_white() + theme(legend.position = "right", 
                                             legend.direction = "vertical", legend.title = element_text(colour = "black"), 
                                             legend.text = element_text(colour = "black", 
                                                                        size = 7), legend.key.height = grid::unit(0.6, "cm"),                                                
                                             legend.key.width = grid::unit(0.5, "cm"), strip.text.x = element_text(size = 7, 
                                                                                                                   color = "black", face = "bold.italic"), strip.text.y = element_text(size = 7, 
                                                                                                                                                                                       color = "black", face = "bold.italic", angle = 270), 
                                             axis.text.x = element_text(size = 8, colour = "black", 
                                                                        angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(vjust = 0.2, 
                                                                                                                                      colour = "black", size = 8), axis.ticks = element_line(size = 0.4), 
                                             plot.background = element_blank(), panel.border = element_blank(), 
                                             panel.grid = element_line(colour = "grey92"), plot.margin = margin(0, 0, 0, 4, "cm"))
                                                                                                                                                                                                        
  if (!is.null(facet_grid_x) & !is.null(facet_grid_y)) {
    g <- g + facet_grid(data[[facet_grid_y]] ~ data[[facet_grid_x]], 
                        switch = "x", scales = "free", space = "free")
  }
  if (is.null(facet_grid_x) & !is.null(facet_grid_y)) {
    g <- g + facet_grid(data[[facet_grid_y]] ~ ., scales = "free_y",  
                        space = "free_y")
  }
  if (!is.null(facet_grid_x) & is.null(facet_grid_y)) {
    g <- g + facet_grid(~data[[facet_grid_x]], switch = "x", 
                        scales = "free_x", space = "free_x")
  }
  if (!is.null(annotation_GO)) {
    if (class(annotation_GO) == "list") {
      gtables <- sapply(names(annotation_GO), function(y) create_gtable(annotation_GO[[y]],
                                                                        y, cutoff, topGO))
      g1 <- gridExtra::arrangeGrob(grobs = gtables, ncol = GO_ncol)
    }
    else {
      g1 <- create_gtable(annotation_GO, " ", cutoff, topGO)
    }
    g <- gridExtra::grid.arrange(g, g1, ncol = 2)
  }
  
  
  
  gp <- ggplotly(g,tooltip=c("text"), width=plot_width, height=plot_height)
  
  gp <- gp %>% config(
    toImageButtonOptions = list(
      format = "svg"
    )
  )
  
  
  return(gp)
}
