#### te.pipeline 2025-10-19
#### main code in comparing across tissue/age
#packages----
if(T){
  library(stringr)
  library(openxlsx)
  library(readr)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(progress)
  library(tidyr)
  library(reshape2)
  library(stats)
  library(tidyverse)
  library(ggrepel)
  library(ggalluvial)
  library(data.table)
  library(ggsci)
  library(cowplot)
  library(scales)
  library(Cairo)
  library(ComplexHeatmap)
  library(ggforce)
}
#path----
tissue ='tissue'
if(T){
  path <- 'limma_result/'
  path2 <- 'limma_pdf/'
  rds1 <- paste0(path,tissue,"_lcpm_ori.rds")
  rds2 <- paste0(path,tissue,"_list_te_deg.rds")
  rds3 <- paste0(path,tissue,"_te_degtem.rds")
  rds4 <- paste0(path,tissue,"_te_lcpm_updown_list.rds")
  rds5 <- paste0(path,tissue,"_gene_lcpmtem.rds")
  rds6 <- paste0(path,tissue,"_gene_degtem.rds")
  rds7 <- paste0(path,tissue,"_protein_lcpmtem.rds")
  rds8 <- paste0(path,tissue,"_protein_degtem.rds")
  uprds <- paste0(path,tissue,"_up_temer.rds")
  downrds <- paste0(path,tissue,"_down_temer.rds")
  plotrds1 <- paste0(path,tissue,"_degplot1.rds")
  plotrds2 <- paste0(path,tissue,"_degplot2.rds")
  plotrds2n <- paste0(path,tissue,"_degplot2n.rds")
  genelist_addrds <- paste0(path,tissue,"_genelist_add.rds")
  bpdf1 <- paste0(path2,tissue,'_plot1.pdf')
  bpdf2 <- paste0(path2,tissue,'_plot2.pdf')
  bpdfflow <- paste0(path2,tissue,'_plot_flow.pdf')
  bpdf3 <- paste0(path2,tissue,'_plot3.pdf')
  bpdfgo <- paste0(path2,tissue,'_plotgo.pdf')
  bpdfpie <- paste0(path2,tissue,'_plotpie.pdf')
  bpdf4 <- paste0(path2,tissue,'_plot4.pdf')
  bpdf5 <- paste0(path2,tissue,'_plot5.pdf')
  bpdf6 <- paste0(path2,tissue,'_plot6.pdf')
  bpdf7 <- paste0(path2,tissue,'_plot7.pdf')
  bpdf8 <- paste0(path2,tissue,'_plot8.pdf')
  bpdf9 <- paste0(path2,tissue,'_plot1229_01.pdf')
  bpdf10 <- paste0(path2,tissue,'_plot1229_02.pdf')
  bpdf11 <- paste0(path2,tissue,'_plot1229_03.pdf')
}

#Step1 Normolization----
gene_all <- read.table('te_all.txt',header = T) #loading te maritx
if(T){
  group <- factor(c(rep("L3",20),rep("L1",20),rep("L2",20),rep("L4",20)))
  cor_list <- c('L2_1','L3_1','L4_1','L3_2','L4_2','L4_3')
  color_red <- c("#DD9898", "#CD6868", "#B43C3C", "#561B1B")
  color_org <- c("#F4E6AE", "#F3D767", "#DDB309", "#937C01")
  color_gre <- c("#D1F7C8", "#9AF386", "#2DD307", "#167500")
  color_blu <- c("#B0DBFE", "#4EACFC", "#0479DC", "#023764")
  color_pur <- c("#D4BAFF", "#A26AFF", "#5700E5", "#260065")
  corlist <- list(color_red,color_org,color_gre,color_blu,color_pur)
}
cor_list_sub <- cor_list[c(13:18)]
cor5 <- corlist[[3]]
if(T){
  x <- gene_all
  m <- as.character(group)
  y <- DGEList(counts=x,group=m)
  ##normalization
  cpm <- cpm(y)
  lcpm <- cpm(y, log=TRUE)
  table(rowSums(y$counts==0)==5)
  keep.exprs <- filterByExpr(y, group=group,min.count = 5, 
                             min.total.count = 5, 
                             large.n = 5, min.prop = 0.05)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log=TRUE)
}
saveRDS(lcpm,rds1)
#Step2 DEGs----
##deg
colnames(lcpm)
logC <- as.data.frame(lcpm)
if(T){
  cor_list <- c('L2_1','L3_1','L4_1','L3_2','L4_2','L4_3')
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(logC)
  fit <- lmFit(logC,design)
  contr.matrix <- makeContrasts(L2_1=L2-L1,L3_1=L3-L1,L4_1=L4-L1,
                                L3_2=L3-L2,L4_2=L4-L2,L4_3=L4-L3,
                                levels = colnames(design))
  fit <- contrasts.fit(fit,contr.matrix)
  fit <- eBayes(fit)
  list_deg <- list()#cor_list
}
for(i in 1:length(cor_list)){
  cor <- cor_list[i]
  A <- topTable(fit,coef = i, 
                adjust.method = "BH",
                number = Inf)
  colnames(A) <- paste0(cor,"_",colnames(A))
  A$id <- rownames(A)
  if(i==1){
    DEG <- A
  }else{
    DEG <- merge(DEG,A,by.x="id",by.y="id")
  }
  list_deg[[cor]] <- A
  print(i)
  print(cor)
}
saveRDS(list_deg, rds2)
#Step3 MDS+pointplot----
list_deg <- readRDS(rds2)
lcpm <- readRDS(rds1)
##mds plot
for (i in 1:length(list_deg)) {
  A <- list_deg[[i]]
  logFC <- A[,1]
  top_genes <- A$id[order(abs(logFC),decreasing = T)[1:200]]
  if(i==1){gene <- top_genes}else{gene <- c(gene,top_genes)}
  print(i)
}
gene <- unique(gene)
lcpm2 <- lcpm[gene,]
brain_ll <- c('L1','L2','L3','L4')
if(T){
  lcpm <- as.data.frame(lcpm)
  iris = as.data.frame(t(lcpm2))
  dis_iris = dist(iris,p=2)
  mds_x = cmdscale(dis_iris)
  mds_x = data.frame(mds_x)
  xy = cbind(mds_x, m)
  colnames(xy) <- c('Dim1','Dim2','samples')
  xy$samples <- factor(xy$samples,levels=brain_ll)
  
}
CairoPDF(bpdf1,height = 4,width = 3.2,family = 'Times New Roman')
if(T){
  xy$samples2 <- 'Adolescent'
  xy$samples2[which(xy$samples=='L2')] <- 'Adult'
  xy$samples2[which(xy$samples=='L3')] <- 'Middle-aged'
  xy$samples2[which(xy$samples=='L4')] <- 'Old'
  xy$samples2 <- factor(xy$samples2,levels=c('Adolescent','Adult','Middle-aged','Old'))
  ggplot(xy,aes(x = Dim1,y= Dim2, color =samples2, 
                fill = samples2,shape = samples2))+
    geom_point(aes(colour = factor(samples2)), size = 3) +
    #geom_point(colour = "grey90", size = 0.5)+
    scale_color_manual(values = c('#cd3232','#e59f00','#cd7ba9','#019e74'))+
    scale_fill_manual(values =c('#cd3232','#e59f00','#cd7ba9','#019e74'))+
    labs(title = tissue)+
    scale_shape_manual(values=c(21,22,24,25)) + theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12,family = "Times New Roman",
                                    hjust=0.5,face = 'bold'),
          axis.text.x = element_text(size = 9,family = "Times New Roman"),
          axis.text.y = element_text(size = 9,family = "Times New Roman"),
          axis.title.x = element_text(size = 10,family = "Times New Roman"),
          axis.title.y = element_text(size = 10,family = "Times New Roman"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10,family = "Times New Roman"),
          legend.position = 'bottom' ) + 
    guides(color = guide_legend(nrow = 2)) +
    theme(legend.box = "horizontal",legend.key.spacing.y = unit(0, 'cm'),  # 垂直间距
          legend.key.spacing.x = unit(1, 'cm'),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")  # 图例边缘间距
    )
  
}
dev.off()
##pointplot
if(T){
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(logC)
  fit <- lmFit(logC,design)
  contr.matrix <- makeContrasts(L2_1=L2-L1,L3_1=L3-L1,L4_1=L4-L1,levels = colnames(design))
  fit <- contrasts.fit(fit,contr.matrix)
  fit <- eBayes(fit)
  T2_1 <- topTable(fit,coef = 'L2_1', adjust.method = "BH",number = Inf)
  T2_1$ns <- 'NS'
  T2_1$ns[which(T2_1$adj.P.Val<0.05&T2_1$logFC>1)] <- 'Up-regulated'
  T2_1$ns[which(T2_1$adj.P.Val<0.05&T2_1$logFC< -1)] <- 'Down-regulated'
  table(T2_1$ns)
  T3_1 <- topTable(fit,coef = 'L3_1', adjust.method = "BH",number = Inf)
  T3_1$ns <- 'NS'
  T3_1$ns[which(T3_1$adj.P.Val<0.05&T3_1$logFC>1)] <- 'Up-regulated'
  T3_1$ns[which(T3_1$adj.P.Val<0.05&T3_1$logFC< -1)] <- 'Down-regulated'
  table(T3_1$ns)
  T4_1 <- topTable(fit,coef = 'L4_1', adjust.method = "BH",number = Inf)
  T4_1$ns <- 'NS'
  T4_1$ns[which(T4_1$adj.P.Val<0.05&T4_1$logFC>1)] <- 'Up-regulated'
  T4_1$ns[which(T4_1$adj.P.Val<0.05&T4_1$logFC< -1)] <- 'Down-regulated'
  table(T4_1$ns)
  T2_1$cluster <- 'Adult'
  T3_1$cluster <- 'Middle-Aged'
  T4_1$cluster <- 'Old'
  T2_1$gene <- rownames(T2_1)
  T3_1$gene <- rownames(T3_1)
  T4_1$gene <- rownames(T4_1)
  df <- rbind(T2_1,T3_1)
  df <- rbind(df,T4_1)
  df$cluster <- factor(df$cluster,levels = c('Adult','Middle-Aged','Old'))
}
saveRDS(df,rds3)
df <- readRDS(rds3)
table(df[which(df$cluster=='Adult'),]$ns)
table(df[which(df$cluster=='Middle-Aged'),]$ns)
table(df[which(df$cluster=='Old'),]$ns)
table(df$cluster)
if(T){
  ##
  df <- df[which(df$ns!='NS'),]
  colnames(df)
  df_s <- as.data.frame(table(df$cluster,df$ns))
  p <- ggplot()+
    geom_jitter(data = df,
                aes(x = cluster, y = logFC, color = ns),
                size = 0.85, width =0.4)
  dfbar<-data.frame(x=c('Adult','Middle-Aged','Old'), y=c(4, 6.5, 7))
  dfbar1<-data.frame(x=c('Adult','Middle-Aged','Old'),y=c(-4, -6.5, -6.5))
  p2 <-ggplot()+
    geom_col(data = dfbar,mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar1,mapping = aes(x = x,y = y),
             fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = df,aes(x = cluster, y = logFC, color = ns),
                size = 0.85,width =0.4)
  dfcol<-data.frame(x=c(1:3),y=0,label=c('Adult','Middle-Aged','Old'))
  dfcol2<-data.frame(x=c(1:3),y=c(-5, -7.5, -7.5),label=df_s$Freq[c(4:6)])
  dfcol3<-data.frame(x=c(1:3),y=c(5, 7.5, 8),label=df_s$Freq[c(1:3)])
  mycol <- c('#e59f00','#cd7ba9','#019e74')
  p3 <- p2 + geom_tile(data = dfcol,aes(x=x,y=y), height=1,
                       color = "#dcdcdc",fill = mycol,alpha = 0.6,show.legend = F)
  p5 <- p3 + scale_color_manual(name=NULL,values = c("blue","red"))
  p6 <- p5 + labs(x="",y="average logFC")+
    geom_text(data=dfcol,aes(x=x,y=y,label=c('Adult','Middle-Aged','Old')),
              size =4,color ="white",fontface = "bold",family = "Times New Roman")+
    geom_text(data=dfcol2,aes(x=x,y=y,label=c('6547','10155','29149')),
              size =4,color ="black",fontface = "bold",family = "Times New Roman")+
    geom_text(data=dfcol3,aes(x=x,y=y,label=c('6108','13444','7821')),
              size =4,color ="black",fontface = "bold",family = "Times New Roman")
  
  p7 <- p6+theme_minimal()+
    theme(axis.title = element_text(size = 13,color = "black",
                                    face = "bold",family = "Times New Roman"),
          axis.line.y = element_line(color = "black", size = 1),
          axis.line.x = element_blank(),axis.text.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.text = element_text(size = 12,,family = "Times New Roman"))
}  
CairoPDF(bpdf2,height = 5,width = 4,family = 'Times New Roman')
p7
dev.off()
#Step4 Data for heatmap(TG)----
##data
updown_lcpm <- readRDS(rds4)
if(T){
  up1_lcpm <- updown_lcpm[[1]]
  up1_lcpm <- as.data.frame(up1_lcpm)
  up1_lcpm$gene <- rownames(up1_lcpm)
  up1_lcpm$cluster <- 'Up1'
  up2_lcpm <- updown_lcpm[[2]]
  up2_lcpm <- as.data.frame(up2_lcpm)
  up2_lcpm$gene <- rownames(up2_lcpm)
  up2_lcpm$cluster <- 'Up2'
  up3_lcpm <- updown_lcpm[[3]]
  up3_lcpm <- as.data.frame(up3_lcpm)
  up3_lcpm$gene <- rownames(up3_lcpm)
  up3_lcpm$cluster <- 'Up3'
  down1_lcpm <- updown_lcpm[[5]]
  down1_lcpm <- as.data.frame(down1_lcpm)
  down1_lcpm$gene <- rownames(down1_lcpm)
  down1_lcpm$cluster <- 'Down1'
  down2_lcpm <- updown_lcpm[[6]]
  down2_lcpm <- as.data.frame(down2_lcpm)
  down2_lcpm$gene <- rownames(down2_lcpm)
  down2_lcpm$cluster <- 'Down2'
  down3_lcpm <- updown_lcpm[[7]]
  down3_lcpm <- as.data.frame(down3_lcpm)
  down3_lcpm$gene <- rownames(down3_lcpm)
  down3_lcpm$cluster <- 'Down3'
}
re_list <- list(up1_lcpm,up2_lcpm,up3_lcpm, down1_lcpm,down2_lcpm,down3_lcpm)
data1 <- rbind(rbind(up1_lcpm,up2_lcpm),up3_lcpm)
data2 <- rbind(rbind(down1_lcpm,down2_lcpm),down3_lcpm)
data <- rbind(data1,data2)
#gene
lcpm <- readRDS(rds5)
genelist <- rownames(lcpm)
#position #select
position_unique <- read.table('position_unique_m.txt',header = T) ##TE-Gene pairs
anno1_tem <- list()
for (i in 1:6) {
  k <- re_list[[i]]
  anno1_tem1 <- position_unique[which(position_unique$TE_transcript_id %in% k$gene),]
  anno1_tem1 <- anno1_tem1[which(anno1_tem1$Gene_name%in% genelist),]
  anno1_tem1 <- anno1_tem1[order(anno1_tem1$Gene_name),]
  anno1_tem[[i]] <- anno1_tem1
  if(i==1){anno1 <- anno1_tem[[i]]}else{
    anno1 <- rbind(anno1,anno1_tem[[i]])}
}
####data
group3 <- factor(c(rep("Middle_Aged",20),rep("Adolescent",20), rep("Adult",20),rep("Old",20)))
sample <- c('Adolescent','Adult','Middle_Aged','Old')
te_s_list <- anno1$TE_transcript_id
gene_s_list <- anno1$Gene_name
#te ##te_select_data ##te_select_anno
data_te <- data[te_s_list,]
#gene #gene_select_anno
lcpm <- as.data.frame(lcpm)
lcpm$gene <- rownames(lcpm)
lcpm2 <- lcpm[gene_s_list,]
deg_plot2_list <- list(anno1,data_te,lcpm2)
saveRDS(deg_plot2_list,plotrds2)

#Step5 plot----
##select-heatmap####
degplot2 <- readRDS(plotrds2)
anno1 <- degplot2[[1]]
df <- readRDS(rds6)
#Down-regulated NS Up-regulated
table(df$cluster)
k2 <- df[which(df$cluster=='Adult' & df$ns=='NS'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='NS'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Up-regulated'),]$gene
upgene <- intersect(k2,intersect(k3,k4))
k2 <- df[which(df$cluster=='Adult' & df$ns=='NS'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='Up-regulated'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Up-regulated'),]$gene
upgene2 <- intersect(k2,intersect(k3,k4))
k2 <- df[which(df$cluster=='Adult' & df$ns=='Up-regulated'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='Up-regulated'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Up-regulated'),]$gene
upgene3 <- intersect(k2,intersect(k3,k4))
upgene <- c(upgene,upgene2,upgene3)
upgene <- unique(upgene)
length(upgene)#1352
k2 <- df[which(df$cluster=='Adult' & df$ns=='NS'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='NS'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Down-regulated'),]$gene
downgene <- intersect(k2,intersect(k3,k4))
k2 <- df[which(df$cluster=='Adult' & df$ns=='NS'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='Down-regulated'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Down-regulated'),]$gene
downgene2 <- intersect(k2,intersect(k3,k4))
k2 <- df[which(df$cluster=='Adult' & df$ns=='Down-regulated'),]$gene
k3 <- df[which(df$cluster=='Middle-Aged' & df$ns=='Down-regulated'),]$gene
k4 <- df[which(df$cluster=='Old' & df$ns=='Down-regulated'),]$gene
downgene3 <- intersect(k2,intersect(k3,k4))
downgene <- c(downgene,downgene2,downgene3)
downgene <- unique(downgene)
length(downgene)#2664
k1 <- anno1[c(1:22),]
k1 <- unique(k1)#19
k1u <- k1[which(k1$Gene_name %in% upgene),]#14
k1d <- k1[which(k1$Gene_name %in% downgene),]#0
anno1.new <- rbind(k1u,k1d)
genelist_add <- list()
genelist_add[['uu']] <- as.data.frame(table(k1u$Gene_name))
genelist_add[['ud']] <- as.data.frame(table(k1d$Gene_name))
k1 <- anno1[c(23:113),]
k1 <- unique(k1)#91
k1u <- k1[which(k1$Gene_name %in% upgene),]#21
k1d <- k1[which(k1$Gene_name %in% downgene),]#18
anno1.new2 <- rbind(k1u,k1d)
genelist_add[['du']] <- as.data.frame(table(k1u$Gene_name))
genelist_add[['dd']] <- as.data.frame(table(k1d$Gene_name))
anno1.new <- rbind(anno1.new,anno1.new2)
anno1.new <- unique(anno1.new)#53
dim(table(anno1.new$Gene_name))#14
saveRDS(genelist_add,genelist_addrds)

te_s_list <- anno1.new$TE_transcript_id
gene_s_list <- anno1.new$Gene_name
lcpm <- readRDS(rds1)
lcpm <- as.data.frame(lcpm)
lcpm_select <- lcpm[te_s_list,]
lcpm <- readRDS(rds5)
lcpm <- as.data.frame(lcpm)
lcpm$gene <- rownames(lcpm)
lcpm2 <- lcpm[gene_s_list,]
dim(lcpm2)
deg_plot2_list <- list(anno1.new,lcpm_select,lcpm2)
saveRDS(deg_plot2_list,plotrds2n)

##heatmap####
degplot2 <- readRDS(plotrds2n)
#color9 <- c('#002C5B','#195696','#4D9AC7','#B6D8E8','#F7F7F7','#F9C5AB','#DC6E56','#9C1129','#630019')
color6 <- c('#e04532','#f2764d','#fddf91','#F7F7F7','#c2dfed','#73a0cc','#4c7cbb')
if(T){
  anno1 <- degplot2[[1]]
  data_te <- degplot2[[2]]
  lcpm2 <- degplot2[[3]]
  ##left te heatmap
  data_plot <- data_te[,c(21:60,1:20,61:80)]
  #table(data_te$cluster)
  anno4 <- c(rep('A', 14),rep('B', 39))
  #row anno
  annotation_row <- as.data.frame(c(rep('Up', 14), rep('Down', 39)))
  colnames(annotation_row) <- 'regulate'
  annotation_row$regulate2 <-  c(rep('UpUp', 14), rep('UpDown', 0),
                                 rep('DownUp',21),rep('DownDown', 18))
  annotation_row$regulate <- factor(annotation_row$regulate,
                                    levels = c('Up','Down'))
  rownames(annotation_row) <- rownames(data_te)
  cor6 = list(Cluster = c('Up' = '#e04532', 'Down' ='#4c7cbb'))
  row_anno <- HeatmapAnnotation(Cluster = annotation_row$regulate,
                                col = cor6,show_annotation_name = F,
                                which = c("row"),
                                annotation_legend_param = list(Cluster = list(
                                  title='Cluster',
                                  title_gp=gpar(fontsize=9,fontface='bold',
                                                fontfamily='Times New Roman'),
                                  labels_gp=gpar(fontsize=8,
                                                 fontfamily='Times New Roman'))),
                                simple_anno_size = unit(2, "mm")) 
  #col anno
  annotation_col <- data.frame(c(rep("Adolescent",20),rep("Adult",20),
                                 rep("Middle_Aged",20),rep("Old",20)))
  colnames(annotation_col) <- 'Group'
  rownames(annotation_col) = colnames(data_plot)
  cor4 = list(Group = c('Adolescent' = '#cd3232', 'Adult' ='#e59f00',
                        'Middle_Aged' = '#cd7ba9', 'Old' ='#019e74'))
  col_anno <- HeatmapAnnotation(Group = annotation_col$Group,
                                col = cor4,show_annotation_name = F,
                                which = c("column"),
                                gap = unit(2,'mm'),show_legend = TRUE,
                                annotation_legend_param = list(Group = list(
                                  title='Group',
                                  title_gp=gpar(fontsize=9,fontface='bold',
                                                fontfamily='Times New Roman'),
                                  labels_gp=gpar(fontsize=8,
                                                 fontfamily='Times New Roman'))),
                                simple_anno_size = unit(2, "mm")) 
  #te
  data_tem <- data_plot
  data_tem <- as.data.frame(t(scale(t(data_tem))))
  data_tem <- as.matrix(data_tem)
  col_fun = circlize::colorRamp2(c(3,2,1,0,-1,-2,-3),
                                 color6)
  h3 <- Heatmap(data_tem,col = col_fun,
                width = 8,
                cluster_rows = F, show_row_names = F,
                row_title = NULL,row_order = rownames(data_tem),  
                row_split =anno4,row_gap = unit(c(1), "mm"),
                cluster_columns = F,show_column_names = F,
                column_title = "TE",
                column_title_gp = gpar(fontsize = 10,fontface='bold',
                                       fontfamily='Times New Roman'),
                top_annotation = col_anno,
                left_annotation = row_anno,
                show_heatmap_legend = F)
  #gene
  data_tem <- lcpm2[,-81]
  data_tem <- data_tem[,c(21:60,1:20,61:80)]
  data_tem <- as.data.frame(t(scale(t(data_tem))))
  data_tem <- as.matrix(data_tem)
  h4 <- Heatmap(data_tem,col = col_fun,width = 8,
                cluster_rows = F,show_row_names = F,
                row_title = NULL,row_order = rownames(data_tem),
                row_split =anno4,row_gap = unit(c(1), "mm"),
                cluster_columns = F,show_column_names = F,
                column_title = "mRNA",
                column_title_gp = gpar(fontsize = 10,fontface='bold',
                                       fontfamily='Times New Roman'),
                top_annotation = col_anno,
                heatmap_legend_param = list(
                  title = 'Expression',
                  title_gp=gpar(fontsize=9,fontface='bold',
                                fontfamily='Times New Roman'),
                  labels_gp=gpar(fontsize=8,fontfamily='Times New Roman')
                  #direction='vertical',
                  #legend_position='left',
                  # title_position='leftcenter-rot'
                ))
  #position
  data_tem <- anno1[,c(5,15)]
  data_tem <- as.matrix(data_tem[,2])
  colnames(data_tem) <- 'Position'
  annotation_col <- data.frame("Position")
  colnames(annotation_col) <- 'Position'
  rownames(annotation_col) = 'Position'
  cor1 = list(Position = c('Position' = "#AAE1C4"))
  col_anno2 <- HeatmapAnnotation(Position = annotation_col$Position,
                                 col = cor1,show_annotation_name = F,
                                 which = c("column"),gap = unit(2,'mm'),
                                 show_legend = F,
                                 simple_anno_size = unit(2, "mm")) 
  col_fun2 = circlize::colorRamp2(c(40000,20000,1000,0,-10000,-20000,-40000),
                                  color6)
  h5 <- Heatmap(data_tem,col = col_fun2,
                width = 0.3,
                cluster_rows = F,show_row_names = F,
                row_title = NULL,row_order = rownames(data_tem),
                row_split =anno4,row_gap = unit(c(1), "mm"),
                cluster_columns = F,show_column_names = F,
                column_title = "Position",
                column_title_gp = gpar(fontsize = 10,fontface='bold',
                                       fontfamily='Times New Roman'),
                top_annotation = col_anno2,
                heatmap_legend_param = list(
                  title = 'Position',
                  title_gp=gpar(fontsize=9,fontface='bold',
                                fontfamily='Times New Roman'),
                  labels_gp=gpar(fontsize=8,fontfamily='Times New Roman')))
  #class
  data_tem <- anno1[,c(5,7)]
  data_tem <- as.matrix(data_tem[,2])
  colnames(data_tem) <- 'Class'
  annotation_col <- data.frame("Class")
  colnames(annotation_col) <- 'Class'
  rownames(annotation_col) = 'Class'
  cor1 = list(Class = c('Class' = "#c3add3"))
  col_anno2 <- HeatmapAnnotation(Class = annotation_col$Class,
                                 col = cor1,show_annotation_name = F,
                                 which = c("column"),gap = unit(2,'mm'),
                                 show_legend = F,
                                 simple_anno_size = unit(2, "mm")) 
  col_fun2 =  c('#019e74',#D
                '#e59f00',#LI
                '#cd7ba9',#LT
                '#3888cb',#U
                '#cd3232'#SI
                
  )
  h6 <- Heatmap(data_tem,col =  col_fun2 ,
                width = 0.3,
                cluster_rows = F,show_row_names = F,
                row_title = NULL,row_order = rownames(data_tem),
                row_split =anno4,row_gap = unit(c(1), "mm"),
                cluster_columns = F,show_column_names = F,
                column_title = "Class",
                column_title_gp = gpar(fontsize = 10,fontface='bold',
                                       fontfamily='Times New Roman'),
                top_annotation = col_anno2,
                heatmap_legend_param = list(
                  title = 'Class',
                  at=c('SINE','LINE','LTR','DNA','Unknown'),
                  
                  title_gp=gpar(fontsize=9,fontface='bold',
                                fontfamily='Times New Roman'),
                  labels_gp=gpar(fontsize=8,fontfamily='Times New Roman')))
  #blank
  data_tem <-  annotation_row
  data_tem <- as.matrix(data_tem[,2])
  colnames(data_tem) <- 'Cla'
  annotation_col <- data.frame("Cla")
  colnames(annotation_col) <- 'Cla'
  rownames(annotation_col) = 'Cla'
  col_fun2 =  c('#e7efcd',#D
                '#c1e6f0',#LI
                #'#fff2e2',#LT
                '#fecccb'
  )
  h7 <- Heatmap(data_tem,col =  col_fun2 ,
                width = 4,
                cluster_rows = F,show_row_names = F,
                row_title = NULL,row_order = rownames(data_tem),
                row_split =anno4,row_gap = unit(c(1), "mm"),
                cluster_columns = F,show_column_names = F,
                column_title_gp = gpar(fontsize = 10,fontface='bold',
                                       fontfamily='Times New Roman'),
                show_heatmap_legend = F)
  
  ht_list <- h6+h3+h5+h4+h7
}
CairoPDF(bpdf7,height = 3,width = 10,family = 'Times New Roman')
draw(ht_list,ht_gap = unit(c(2,2,2,2), "mm"),
     heatmap_legend_side='bottom',merge_legend=T,
     annotation_legend_side='bottom',row_title = tissue,
     row_title_gp = gpar(fontsize = 12,fontface='bold',fontfamily='Times New Roman'))
dev.off()

##heatmap_sub#####
tissue ='Lung'
if(T){
  path <- 'limma_result/'
  path2 <- 'limma_pdf/'
  hscsv1 <- paste0(path2,tissue,'_anno1.csv')
  hspdf1 <- paste0(path2,tissue,'_plot0121_01.pdf')
  hspdf2 <- paste0(path2,tissue,'_plot0121_02.pdf')
  hspdf3 <- paste0(path2,tissue,'_plot0121_03.pdf')
  hspdf4 <- paste0(path2,tissue,'_plot0121_04.pdf')
  hspdf5 <- paste0(path2,tissue,'_plot0121_05.pdf')
  hspdf6 <- paste0(path2,tissue,'_plot0121_06.pdf')
  hspdf7 <- paste0(path2,tissue,'_plot0121_07.pdf')
  hspdf8 <- paste0(path2,tissue,'_plot0121_08.pdf')
  hspdf9 <- paste0(path2,tissue,'_plot0121_09.pdf')
  hspdf10 <- paste0(path2,tissue,'_plot0121_010.pdf')
  hspdf11 <- paste0(path2,tissue,'_plot0121_011.pdf')
  hspdf12 <- paste0(path2,tissue,'_plot0121_012.pdf')
}
degplot2 <- readRDS(plotrds2)
anno1 <- degplot2[[1]]
write.csv(anno1,hscsv1)
df <- readRDS(rds8)
genelist <- unique(df$gene)
genelist <- intersect(genelist,unique(anno1$Gene_name))
df <- df[which(df$gene %in% genelist),]
df <- df[which(df$cluster=='Old'),]
df <- df[which(df$ns != 'NS'),]
df2 <- readRDS(rds6)
df2 <- df2[which(df2$gene %in% genelist),]
df2 <- df2[which(df2$cluster=='Old'),]
df2 <- df2[which(df2$ns != 'NS'),]
genelist <- unique(c(df$gene,df2$gene))
#lcpm+plot
theme_dlqpect <- function(...){
  theme(
    panel.grid = element_blank(),
    axis.text =  element_text(size = 9, family = 'Times New Roman'),
    axis.title = element_text(size = 10, family = 'Times New Roman'),
    plot.title = element_text(hjust = 0.5, vjust = 1,size = 10,
                              margin = margin(t = 10),face = "bold",
                              family = 'Times New Roman'),
    legend.text = element_text(size = 9,,family = "Times New Roman"),
    legend.title = element_blank() )
}
library(reshape2)
library(ggplot2)
library(patchwork)
if(T){
  group <- factor(c(rep('Middle_Aged',20),rep('Adolescent',20),
                    rep('Adult',20),rep('Old',20)))
}
m <- as.character(group)
genelist
trpplot_pipeline <- function(j) {
  telist <- anno1[which(anno1$Gene_name==genelist[j]),]
  telist <- telist$TE_transcript_id
  #mrna_plot
  lcpm <- readRDS(rds5)
  lcpm.count <- lcpm[genelist[j],]
  lcpm.count <- as.data.frame(t(lcpm.count))
  lcpm.count <- 2^(lcpm.count)
  lcpm.count$group <- m
  lcpm.count <- lcpm.count[which(lcpm.count$group=='Adolescent'| lcpm.count$group=='Old'),]
  data_plot =  reshape2::melt(lcpm.count)
  colnames(data_plot) = c('group','te','lcpm')
  data_plot$group <- factor(data_plot$group)
  data_plot$te <- factor(data_plot$te)
  k <- max(data_plot$lcpm)
  plota <- ggplot(data_plot, aes(x = group, y = lcpm,fill = group))+ 
    geom_bar(stat = "summary", fun = mean, width = 0.7,color='black')+
    geom_jitter(width = 0.25,size=1)+
    scale_fill_manual("expression",values = c('#f2764d','#73a0cc'))+
    stat_compare_means(comparisons = list(c("Old","Adolescent")),paired = F,
                       method = "wilcox.test",label = "p.signif",label.y = k )+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color='black', width=0.2) +
    stat_summary(fun.y=mean, geom="point", color='black')+
    labs(title = paste0(genelist[j],'_mRNA'),y='Expression',x='')+
    scale_y_continuous(limits = c(0,k*1.2),
                       expand = c(0,0))+
    theme_bw()+theme_dlqpect()
  plot_mrna <- plota
  #protein_plot
  lcpm <- readRDS(rds7)
  lcpm.count <- lcpm[genelist[j],]
  lcpm.count <- as.data.frame(t(lcpm.count))
  lcpm.count <- 2^(lcpm.count)
  lcpm.count$group <- m
  lcpm.count <- lcpm.count[which(lcpm.count$group=='Adolescent'|
                                   lcpm.count$group=='Old'),]
  data_plot =  reshape2::melt(lcpm.count)
  colnames(data_plot) = c('group','te','lcpm')
  data_plot$group <- factor(data_plot$group)
  data_plot$te <- factor(data_plot$te)
  k <- max(data_plot$lcpm)
  plota <- ggplot(data_plot, aes(x = group, y = lcpm,fill = group))+ 
    geom_bar(stat = "summary", fun = mean, width = 0.7,color='black')+
    geom_jitter(width = 0.25,size=1)+
    scale_fill_manual("expression",values = c('#fddf91','#c2dfed'))+
    stat_compare_means(comparisons = list(c("Old","Adolescent")),paired = F,
                       method = "wilcox.test",label = "p.signif",label.y = k )+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color='black', width=0.2) +
    stat_summary(fun.y=mean, geom="point", color='black')+
    labs(title = paste0(genelist[j],'_Protein'),y='Expression',x='')+
    scale_y_continuous(limits = c(0,k*1.2),
                       expand = c(0,0))+
    theme_bw()+theme_dlqpect()
  plot_pro <- plota
  plot <- plot_mrna+plot_pro
  #te_plot
  lcpm <- readRDS(rds1)
  lcpm.count <- lcpm[telist,]
  if(class(lcpm.count)[1]=="numeric"){
    lcpm.count <- as.data.frame(lcpm[telist,],ncol=1)
    colnames(lcpm.count) <- telist
  }else{
    lcpm.count <- as.data.frame(t(lcpm.count))
  }
  lcpm.count <- 2^(lcpm.count)
  lcpm.count$group <- m
  lcpm.count <- lcpm.count[which(lcpm.count$group=='Adolescent'|
                                   lcpm.count$group=='Old'),]
  numberofte <- length(lcpm.count)-1
  for (i in 1:numberofte) {
    tename <- as.character(colnames(lcpm.count))[i]
    data_plot =  reshape2::melt(lcpm.count[,c(i,length(lcpm.count))])
    colnames(data_plot) = c('group','te','lcpm')
    data_plot$group <- factor(data_plot$group)
    data_plot$te <- factor(data_plot$te)
    k <- max(data_plot$lcpm)
    plota <- ggplot(data_plot, aes(x = group, y = lcpm,fill = group))+ 
      geom_bar(stat = "summary", fun = mean, width = 0.7,color='black')+
      geom_jitter(width = 0.25,size=1)+
      scale_fill_manual("expression",values = c('#e04532','#4c7cbb'))+
      stat_compare_means(comparisons = list(c("Old","Adolescent")),paired = F,
                         method = "t.test",label = "p.signif",label.y = k )+
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                   geom="errorbar", color='black', width=0.2) +
      stat_summary(fun.y=mean, geom="point", color='black')+
      labs(title = tename,y='Expression',x='')+
      scale_y_continuous(limits = c(0,k*1.2),
                         expand = c(0,0))+
      theme_bw()+theme_dlqpect()
    plot <- plot+plota
  }
  return(plot)
  
}
length(genelist)
plot <- trpplot_pipeline(1) #1-12
CairoPDF(hspdf1 ,height = 5,width = 6,family = 'Times New Roman')
plot+plot_layout(nrow = 2,byrow =F)
dev.off()

###plot_scatter####
library(ggplot2)
library(patchwork)
lcpm_updown_list <- readRDS(rds4)
#position
position_unique <- read.table('te_gtf.txt', header = T)

if(T){
  up_temer <- position_unique[which(
    position_unique$transcript_id %in% rownames(lcpm_updown_list[[4]])),]#up
  table(up_temer$family_id)
  dim(up_temer)
  down_temer <- position_unique[which(
    position_unique$transcript_id %in% rownames(lcpm_updown_list[[8]])),]#down
  table(up_temer$family_id)
  dim(down_temer)
}
up_temer <- rbind(up_temer,down_temer)
##up
#plot1
telist <- up_temer
table(telist$class_id)
telist$class_id[which(telist$class_id =='RC'|
                        telist$class_id =='Retroposon'|
                        telist$class_id =='Satellite'|
                        telist$class_id =='Unknown')] <- 'Others'
table(telist$class_id)

dflist <- as.data.frame(table(telist$class_id))
sum(dflist$Freq)
colnames(dflist) <- c('Class','Count')
dflist$percentage <- dflist$Count/59.58
dflist$percentage2 <- paste0(round(dflist$percentage,digits = 2),'%')
dflist$Class <- factor(dflist$Class,levels=c('SINE','LINE','LTR','DNA','Others') )
#cor5 <-  c('#cd3232','#e59f00','#cd7ba9','#019e74','#0070C0')
cor5 <- c("#F0C1C1","#F7E2B2","#F0D7E5","#B2E2D5",'#B2D4EC')

CairoPDF(bpdfpie,height = 3,width = 3,family = 'Times New Roman')
if(T){
  ggplot()+
    geom_arc_bar(data=dflist,stat="pie",
                 aes(x0=0,y0=0,r0=0,r=1,
                     amount=percentage,fill=Class),
                 show.legend = FALSE,linewidth=0.5)+theme_bw()+
    theme(panel.grid = element_blank(),panel.border = element_blank(),
          axis.title = element_blank(),axis.text = element_blank(),
          axis.ticks = element_blank())+coord_equal()+
    annotate(geom = "text", x= -0.2, y=  0.65, label="DNA",   size=4,family = "Times New Roman")+
    annotate(geom = "text", x= -0.2, y=  0.5, label="5.05%", size=4,family = "Times New Roman")+
    annotate("segment", x = 0.1, y = 0.6, xend = -0.07, yend = 0.6,color = "black")+
    annotate(geom = "text", x=  0.5, y=  0.5, label="LINE",  size=4,family = "Times New Roman")+
    annotate(geom = "text", x=  0.5, y=  0.35, label="21.4%",size=4,family = "Times New Roman")+
    annotate(geom = "text", x=  0.3, y= -0.3, label="LTR",   size=4,family = "Times New Roman")+
    annotate(geom = "text", x=  0.3, y= -0.45, label="29.12%",size=4,family = "Times New Roman")+
    annotate(geom = "text", x= 0, y= -0.65, label="Others",size=4,family = "Times New Roman")+
    annotate(geom = "text", x= 0, y= -0.8, label="1.46%", size=4,family = "Times New Roman")+
    annotate("segment", x = -0.15, y = -0.7, xend = -0.3, yend = -0.7,color = "black")+
    annotate(geom = "text", x= -0.5, y=  0.05, label="SINE",  size=4,family = "Times New Roman")+
    annotate(geom = "text", x= -0.5, y= -0.1, label="42.97%",size=4,family = "Times New Roman")+
    scale_fill_manual(values = cor5[c(1:5)])
}
dev.off()

#c('SINE','LINE','LTR','DNA','Others')
dflist <- as.data.frame(table(telist$class_id,telist$family_id))
dflist <- dflist[which(dflist$Freq != 0),]
dflist <- dflist[order(dflist$Freq,decreasing = T),]
sum(dflist$Freq)
colnames(dflist) <- c('Class','Family','Count')
dflist$percentage <- dflist$Count/59.58
dflist$percentage2 <- paste0(round(dflist$percentage,digits = 2),'%')
if(T){
  dflist2 <- dflist[which(dflist$Class=='SINE'),]
  dflist2 <- dflist2[c(1:5),]
  dflist2$Family <- factor(dflist2$Family,levels = dflist2$Family)
  p1 <- ggplot(dflist2, aes(x = Family, y = percentage, fill = Class)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percentage2), vjust = -0.5, size = 3,family = "Times New Roman",) +
    labs(title = '',x = '', y = "Percentage") +
    theme_bw()+
    scale_fill_manual(values =c('#cd3232'))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          #axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12,family = "Times New Roman",
                                    hjust=0.5,face = 'bold'),
          axis.text.x = element_text(color = 'black',angle = 45, hjust = 1,,
                                     size = 9,family = "Times New Roman"),
          axis.text.y = element_text(size = 9,family = "Times New Roman"),
          axis.title.x = element_text(size = 9,family = "Times New Roman"),
          axis.title.y = element_text(size = 9,family = "Times New Roman"),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey90"),
          axis.line.y = element_line(colour = 'black',size = 0.25),
          axis.line.x = element_line(colour = 'black',size = 0.25),
          legend.text = element_text(size = 9,family = "Times New Roman")
    ) +
    scale_y_continuous(limits = c(0,20),expand = c(0,0))
}
if(T){
  dflist2 <- dflist[which(dflist$Class=='LINE'),]
  dflist2 <- dflist2[c(1:5),]
  dflist2$Family <- factor(dflist2$Family,levels = dflist2$Family)
  p2 <- ggplot(dflist2, aes(x = Family, y = percentage, fill = Class)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percentage2), vjust = -0.5, size = 3,family = "Times New Roman",) +
    labs(title = '',x = '', y = "Percentage") +
    theme_bw()+
    scale_fill_manual(values =c('#e59f00'))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          #axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12,family = "Times New Roman",
                                    hjust=0.5,face = 'bold'),
          axis.text.x = element_text(color = 'black',angle = 45, hjust = 1,,
                                     size = 9,family = "Times New Roman"),
          axis.text.y = element_text(size = 9,family = "Times New Roman"),
          axis.title.x = element_text(size = 9,family = "Times New Roman"),
          axis.title.y = element_text(size = 9,family = "Times New Roman"),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey90"),
          axis.line.y = element_line(colour = 'black',size = 0.25),
          axis.line.x = element_line(colour = 'black',size = 0.25),
          legend.text = element_text(size = 9,family = "Times New Roman")
    ) +
    scale_y_continuous(limits = c(0,25),expand = c(0,0))
}
if(T){
  dflist2 <- dflist[which(dflist$Class=='LTR'),]
  dflist2 <- dflist2[c(1:5),]
  dflist2$Family <- factor(dflist2$Family,levels = dflist2$Family)
  p3 <- ggplot(dflist2, aes(x = Family, y = percentage, fill = Class)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percentage2), vjust = -0.5, size = 3,family = "Times New Roman",) +
    labs(x = '', y = "Percentage") +
    theme_bw()+
    scale_fill_manual(values =c('#cd7ba9'))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          #axis.ticks.x = element_blank(),
          plot.title = element_text(size = 12,family = "Times New Roman",
                                    hjust=0.5,face = 'bold'),
          axis.text.x = element_text(color = 'black',angle = 45, hjust = 1,,
                                     size = 9,family = "Times New Roman"),
          axis.text.y = element_text(size = 9,family = "Times New Roman"),
          axis.title.x = element_text(size = 9,family = "Times New Roman"),
          axis.title.y = element_text(size = 9,family = "Times New Roman"),
          legend.title = element_blank(),
          #axis.text.x = element_blank(),
          #panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_line(color = "grey90"),
          axis.line.y = element_line(colour = 'black',size = 0.25),
          axis.line.x = element_line(colour = 'black',size = 0.25),
          legend.text = element_text(size = 9,family = "Times New Roman")
    ) +
    scale_y_continuous(limits = c(0,20),expand = c(0,0))
}
CairoPDF(bpdf9,height = 6,width = 5,family = 'Times New Roman')
p1/p2/p3 + plot_layout(heights = c(2,2,2))
dev.off()


##te_family/class####
library(ggpubr)
library(ggplot2)
library(reshape2)
library(plyr)
library(Cairo)
tissue ='Lung'
if(T){
  path <- 'limma_result/'
  path2 <- 'limma_pdf/'
  rds11 <- paste0(path,tissue,'_family.rds')
  rds12 <- paste0(path,tissue,'_class.rds')
  rds13 <- paste0(path,tissue,'_score.rds')
  bpdf9 <- paste0(path2,tissue,'_plot1229_01.pdf')
  bpdf10 <- paste0(path2,tissue,'_plot1229_02.pdf')
  bpdf11 <- paste0(path2,tissue,'_plot1229_03.pdf')
  
  bpdfy1 <- paste0(path2,tissue,'_plot0103_01.pdf')
  bpdfy2 <- paste0(path2,tissue,'_plot0103_02.pdf')
  bpdfy3 <- paste0(path2,tissue,'_plot0103_03.pdf')
  bpdfy4 <- paste0(path2,tissue,'_plot0103_04.pdf')
  bpdfy5 <- paste0(path2,tissue,'_plot0103_05.pdf')
  bpdfy6 <- paste0(path2,tissue,'_plot0103_06.pdf')
  bpdfy7 <- paste0(path2,tissue,'_plot0103_07.pdf')
  bpdfy8 <- paste0(path2,tissue,'_plot0103_08.pdf')
  bpdfy9 <- paste0(path2,tissue,'_plot0103_09.pdf')
}
getwd()#"F:/te_0305/mouse"#D:/R/Rdata/TE/te_0305/mous e
gene_all <- read.table('te_all.txt',header = T)
table(rowSums(gene_all)==0)
###Lung
if(T){
  group <- factor(c(rep("B3",20),rep("B1",20),rep("B2",20),rep("B4",20),
                    rep("H3",20),rep("H1",20),rep("H2",20),rep("H4",20),
                    rep("L3",20),rep("L1",20),rep("L2",20),rep("L4",20),
                    rep("M3",20),rep("M1",20),rep("M2",20),rep("M4",20),
                    rep("S3",20),rep("S1",20),rep("S2",20),rep("S4",20)))
  cor_list <- c('B2_1','B3_1','B4_1','B3_2','B4_2','B4_3',#6
                'H2_1','H3_1','H4_1','H3_2','H4_2','H4_3',
                'L2_1','L3_1','L4_1','L3_2','L4_2','L4_3',
                'M2_1','M3_1','M4_1','M3_2','M4_2','M4_3',
                'S2_1','S3_1','S4_1','S3_2','S4_2','S4_3',
                'BH1','BL1','BM1','BS1','HL1','HM1','HS1','LM1','LS1','MS1',
                'BH2','BL2','BM2','BS2','HL2','HM2','HS2','LM2','LS2','MS2',
                'BH3','BL3','BM3','BS3','HL3','HM3','HS3','LM3','LS3','MS3',
                'BH4','BL4','BM4','BS4','HL4','HM4','HS4','LM4','LS4','MS4')
  color_red <- c("#DD9898", "#CD6868", "#B43C3C", "#561B1B")
  color_org <- c("#F4E6AE", "#F3D767", "#DDB309", "#937C01")
  color_gre <- c("#D1F7C8", "#9AF386", "#2DD307", "#167500")
  color_blu <- c("#B0DBFE", "#4EACFC", "#0479DC", "#023764")
  color_pur <- c("#D4BAFF", "#A26AFF", "#5700E5", "#260065")
  corlist <- list(color_red,color_org,color_gre,color_blu,color_pur)
}
cor_list_sub <- cor_list[c(13:18)]
cor5 <- corlist[[3]]
x <- gene_all[,c(161:240)]
m <- as.character(group[161:240])

te.gtf <- read.table('te_gtf.txt', header = T)#TE annotation file
head(x)
colnames(te.gtf)
te.transcript_id <- as.data.frame(rownames(x))
colnames(te.transcript_id) <- 'transcript_id'
te.gtf2 <- merge(te.transcript_id,te.gtf,
                 by.x = 'transcript_id',
                 by.y = 'transcript_id',all.x = T)
colnames(te.gtf2)
head(te.gtf2)
te.gtf2 <- te.gtf2[,c(1,12,13)]
x$transcript_id <- rownames(x)
x2 <- merge(x,te.gtf2,
            by.x = 'transcript_id',
            by.y = 'transcript_id',all.x = T)
colnames(x2)

theme_dlqdot<- function(...){
  theme(text = element_text(size = 12),
        panel.grid=element_blank(),
        #panel.border = element_rect(fill=NA,colour="black",size=0.8),
        panel.border = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.title =  element_text(family = 'Times New Roman',size = 12,face='bold'),
        axis.line = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(family = 'Times New Roman',size = 12),
        axis.title.y = element_text(family = 'Times New Roman',size = 12),
        axis.text = element_text(family = 'Times New Roman'),
        legend.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(family = 'Times New Roman',size = 10))
}
count <- x2[,-c(1)]
result <- count %>% group_by(class_id,family_id) %>% 
  summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result)
saveRDS(result.k,rds11)
result.k <- readRDS(rds11)
for (i in 3:length(result.k)) {
  result.k[,i] <- result.k[,i]/sum(result.k[,i])
}
head(result.k)
if(T){
  result.k <- result.k[which(result.k$class_id %in% 
                               c('SINE','LINE','LTR','DNA')),]
  color6 <- c('#019e74','#5a9a89','#be98ad','#3888cb','#f0e643',
              "#F0C1C1","#F7E2B2","#F0D7E5","#B2E2D5",'#B2D4EC')
  #p1-p2
  result.k2 <- result.k[which(result.k$class_id =='SINE'),]
  if(T){
    result.k2 <- result.k2[,-1]
    rownames(result.k2) <- result.k2[,1]
    result.k2 <- result.k2[,-1]
    #result.k <- log(result.k+1)
    #result.k2 <- result.k2/1000000
    result.k2 <- as.data.frame(t(result.k2))
    result.k2$sample <- m
    result.k2 <- result.k2 %>% 
      group_by(sample) %>% 
      summarise(across(everything(), mean,na.rm = TRUE))
    result.k2 <- as.data.frame(result.k2)
    result.k2$sample <- c('Adolescent',"Adult", "Middle_Aged", "Old")
    #rownames(result.k2) <- result.k2[,1]
    #result.k2 <- result.k2[,-1]
    plotdata = reshape2::melt(result.k2)
    colnames(plotdata) <- c('Sample','Family','Count')
    plotdata <- plotdata[order(plotdata$Count),]
    table(plotdata$Family)
    plotdata$Family <- factor(plotdata$Family,
                              levels = c('Alu','B4','B2','MIR',
                                         'ID','5S-Deu-L2'))
    plotdata$Sample <- factor(plotdata$Sample,
                              levels = c('Adolescent',"Adult", "Middle_Aged", "Old"))
  }
  p2 <- ggline(plotdata, x = "Sample", y = "Count",
               color= 'Family',palette = "jco")+
    scale_color_manual(values = color6)+
    scale_y_break(breaks = c(0.12, 0.16), scales = 'fixed', space = 0.5)+
    scale_y_continuous(limits = c(0,0.2),
                       breaks = seq(0,0.2,by=0.02),
                       expand = c(0,0))+
    labs(x='',y='Pecentage',title = 'SINE')+
    theme_dlqdot()
  #p3-p4
  result.k2 <- result.k[which(result.k$class_id =='LINE'),]
  if(T){
    result.k2 <- result.k2[,-1]
    rownames(result.k2) <- result.k2[,1]
    result.k2 <- result.k2[,-1]
    #result.k <- log(result.k+1)
    #result.k2 <- result.k2/1000000
    result.k2 <- as.data.frame(t(result.k2))
    result.k2$sample <- m
    result.k2 <- result.k2 %>% 
      group_by(sample) %>% 
      summarise(across(everything(), mean,na.rm = TRUE))
    result.k2 <- as.data.frame(result.k2)
    result.k2$sample <- c('Adolescent',"Adult", "Middle_Aged", "Old")
    #rownames(result.k2) <- result.k2[,1]
    #result.k2 <- result.k2[,-1]
    plotdata = reshape2::melt(result.k2)
    colnames(plotdata) <- c('Sample','Family','Count')
    plotdata <- plotdata[order(plotdata$Count),]
    table(plotdata$Family)
    plotdata <- plotdata[which(plotdata$Family %in% 
                                 c('L1','L2','CR1','RTE-Bo.B',
                                   'Penelope','I-Jockey') ),]
    plotdata$Family <- factor(plotdata$Family,
                              levels = c('L1','L2','CR1','RTE-Bo.B',
                                         'Penelope','I-Jockey'))
    plotdata$Sample <- factor(plotdata$Sample,
                              levels = c('Adolescent',"Adult", "Middle_Aged", "Old"))
  }
  p4 <- ggline(plotdata, x = "Sample", y = "Count",
               color= 'Family',palette = "jco")+
    scale_color_manual(values = color6)+
    scale_y_break(breaks = c(0.04, 0.1), scales = 'fixed', space = 0.5)+
    scale_y_continuous(limits = c(0,0.12),
                       breaks = seq(0,0.12,by=0.01),
                       expand = c(0,0))+
    labs(x='',y='Pecentage',title = 'LINE')+
    theme_dlqdot()
  #p5-p6
  result.k2 <- result.k[which(result.k$class_id =='LTR'),]
  if(T){
    result.k2 <- result.k2[,-1]
    rownames(result.k2) <- result.k2[,1]
    result.k2 <- result.k2[,-1]
    #result.k <- log(result.k+1)
    #result.k2 <- result.k2/1000000
    result.k2 <- as.data.frame(t(result.k2))
    result.k2$sample <- m
    result.k2 <- result.k2 %>% 
      group_by(sample) %>% 
      summarise(across(everything(), mean,na.rm = TRUE))
    result.k2 <- as.data.frame(result.k2)
    result.k2$sample <- c('Adolescent',"Adult", "Middle_Aged", "Old")
    #rownames(result.k2) <- result.k2[,1]
    #result.k2 <- result.k2[,-1]
    plotdata = reshape2::melt(result.k2)
    colnames(plotdata) <- c('Sample','Family','Count')
    plotdata <- plotdata[order(plotdata$Count),]
    table(plotdata$Family)
    plotdata$Family <- factor(plotdata$Family,
                              levels = c('ERVK','ERVL-MaLR',
                                         'ERVL','ERV1',
                                         'LTR','Gypsy'))
    plotdata$Sample <- factor(plotdata$Sample,
                              levels = c('Adolescent',"Adult", "Middle_Aged", "Old"))
  }
  p6 <- ggline(plotdata, x = "Sample", y = "Count",
               color= 'Family',palette = "jco")+
    scale_color_manual(values =color6)+
    scale_y_break(breaks = c(0.04, 0.08), 
                  scales = 'fixed', space = 0.5)+
    scale_y_continuous(limits = c(0,0.11),
                       breaks = seq(0,0.11,by=0.01),
                       expand = c(0,0))+
    labs(x='',y='Pecentage',title = 'LTR')+
    theme_dlqdot()
}

count <- x2[,-c(1,82)]
result2 <- count %>% group_by(class_id) %>% 
  summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
saveRDS(result.k,rds12)
result.k <- readRDS(rds12)
for (i in 2:length(result.k)) {
  result.k[,i] <- result.k[,i]/sum(result.k[,i])
}
head(result.k)
if(T){
  rownames(result.k) <- result.k[,1]
  result.k <- result.k[,-1]
  #result.k <- log(result.k+1)
  #result.k <- result.k/1000000
  result.k <- as.data.frame(t(result.k))
  result.k$sample <- m
  result.k2 <- result.k %>% group_by(sample) %>% 
    summarise(across(everything(), mean,na.rm = TRUE))
  result.k2 <- as.data.frame(result.k2)
  result.k2 <- result.k2[,c('sample','SINE','LINE','LTR','DNA')]
  result.k2$sample <- c('Adolescent',"Adult", "Middle_Aged", "Old")
  #rownames(result.k2) <- result.k2[,1]
  #result.k2 <- result.k2[,-1]
  plotdata = reshape2::melt(result.k2)
  colnames(plotdata) <- c('Sample','Class','Count')
}
p9 <- ggline(plotdata, x = "Sample", y = "Count",
             color= 'Class',palette = "jco")+
  scale_color_manual(values = c('#cd3232','#e59f00','#cd7ba9','#019e74'))+
  scale_y_break(breaks = c(0.28, 0.42), 
                scales = 'fixed', space = 0.5)+
  scale_y_continuous(limits = c(0,0.5),
                     breaks = seq(0,0.5,by=0.05),
                     expand = c(0,0))+
  labs(x='',y='Pecentage',title = tissue)+
  theme_dlqdot()

library(patchwork)
CairoPDF(bpdf10, width =6, height =4,family = 'Times New Roman')
p9 
p2
p4
p6
dev.off()


##score####
#library(BiocManager)
#BiocManager::install('GSVA')
library(GSVA)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
tissue ='tissue'
if(T){
  path <- 'limma_result/'
  path2 <- 'limma_pdf/'
  rds11 <- paste0(path,tissue,'_family.rds')
  rds12 <- paste0(path,tissue,'_class.rds')
  rds13 <- paste0(path,tissue,'_score.rds')
  bpdf9 <- paste0(path2,tissue,'_plot1229_01.pdf')
  bpdf10 <- paste0(path2,tissue,'_plot1229_02.pdf')
  bpdf11 <- paste0(path2,tissue,'_plot1229_03.pdf')
  
  bpdfy1 <- paste0(path2,tissue,'_plot0103_01.pdf')
  bpdfy2 <- paste0(path2,tissue,'_plot0103_02.pdf')
  bpdfy3 <- paste0(path2,tissue,'_plot0103_03.pdf')
  bpdfy4 <- paste0(path2,tissue,'_plot0103_04.pdf')
  bpdfy5 <- paste0(path2,tissue,'_plot0103_05.pdf')
  bpdfy6 <- paste0(path2,tissue,'_plot0103_06.pdf')
  bpdfy7 <- paste0(path2,tissue,'_plot0103_07.pdf')
  bpdfy8 <- paste0(path2,tissue,'_plot0103_08.pdf')
  bpdfy9 <- paste0(path2,tissue,'_plot0103_09.pdf')
}
lcpm <- readRDS(rds5)
m <- as.character(group)
raw_counts = lcpm
raw_counts = as.data.frame(raw_counts)
head(raw_counts)
genelist = read.csv('macro score.csv',header = T)#supplement.table
genelist <- genelist[,c(1,2)]
genelist$Metagene2 <- str_to_title(genelist$Metagene)
list <- list()
for(i in 1:8){
  # i=1
  list[[i]] <- genelist$Metagene2[genelist$Cell.type==(unique(genelist$Cell.type)[i])]
}
names(list)<- unique(genelist$Cell.type)[c(1:8)]
expr=as.matrix(raw_counts)
gsvaPar <- ssgseaParam(exprData = expr, 
                       geneSets = list,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)
gsva_data[1:3,1:3]#head(gsva_data)
gsva_mat = as.data.frame(gsva_data)
gsva_mat_t = gsva_mat %>% t() %>% as.data.frame() %>% 
  mutate(sample_id = rownames(.))
gsva_mat_t$group <- m
colnames(gsva_mat_t) <- str_to_title(colnames(gsva_mat_t))
pathwaylist <- colnames(gsva_mat_t)
ppl1 <- c('Adolescent',"Adult", "Middle_Aged", "Old")
gsva_mat_t$Group[which(gsva_mat_t$Group=='B1')] <-'Adolescent' 
gsva_mat_t$Group[which(gsva_mat_t$Group=='B2')] <-'Adult' 
gsva_mat_t$Group[which(gsva_mat_t$Group=='B3')] <-'Middle_Aged' 
gsva_mat_t$Group[which(gsva_mat_t$Group=='B4')] <-'Old' 
saveRDS(gsva_mat_t,rds13)

#cell cycle #DNA Replication #Chromatin
geneset1 <- msigdbr(species="Mus musculus", collection="H") %>%
  dplyr::filter(gs_name=="HALLMARK_G2M_CHECKPOINT")
geneset2 <- msigdbr(species="Mus musculus", collection="C5", 
                  subcollection="GO:BP") %>%
  dplyr::filter(grepl("GOBP_MITOTIC_CELL_CYCLE_CHECKPOINT", gs_name))
geneset3 <- msigdbr(species="Mus musculus",  collection="C2", 
                    subcollection="CP:REACTOME") %>%
  dplyr::filter(grepl("REACTOME_CHROMATIN_ORGANIZATION", gs_name))
geneset4 <- msigdbr(species="Mus musculus", collection="C5", 
                    subcollection="GO:BP") %>%
  dplyr::filter(gs_name=="GOBP_HISTONE_MRNA_METABOLIC_PROCESS")
geneset5 <- msigdbr(species="Mus musculus", collection="C5", 
                    subcollection="GO:BP") %>%
  dplyr::filter(grepl("GOBP_NUCLEAR_ENVELOPE_ORGANIZATION", gs_name))
geneset6 <- msigdbr(species="Mus musculus", collection="C5", 
                    subcollection="GO:CC") %>%
  dplyr::filter(gs_name=="GOCC_NUCLEAR_LAMINA")
geneset7 <- msigdbr(species="Mus musculus", collection="C2", 
                    subcollection="CP:REACTOME") %>%
  dplyr::filter(grepl("REACTOME_DNA_REPLICATION", gs_name))
geneset8 <- msigdbr(species="Mus musculus", collection="C5", 
                    subcollection="GO:BP") %>%
  dplyr::filter(grepl("GOBP_MRNA_PROCESSING", gs_name))
geneset9 <- msigdbr(species="Mus musculus", collection="C2", 
                    subcollection="CP:REACTOME") %>%
  dplyr::filter(gs_name=="REACTOME_NONSENSE_MEDIATED_DECAY_NMD")


#IMMUNE
geneset10 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:KEGG_MEDICUS") %>% #MEDICUS/LEGACY
  dplyr::filter(gs_name=="KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY")
geneset11 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:REACTOME") %>%
  dplyr::filter(grepl("REACTOME_CELLULAR_SENESCENCE", gs_name))
geneset12 <- msigdbr(species="Mus musculus", collection="H") %>%
  dplyr::filter(grepl("HALLMARK_P53_PATHWAY", gs_name))
geneset13 <- msigdbr(species="Mus musculus", collection="H") %>%
  dplyr::filter(gs_name=="HALLMARK_APOPTOSIS")
geneset14 <- msigdbr(species="Mus musculus", collection="H") %>%
  dplyr::filter(gs_name=="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
geneset15 <- msigdbr(species="Mus musculus", collection="C5", 
                     subcollection="GO:BP") %>%
  dplyr::filter(grepl("GOBP_AUTOPHAGY", gs_name))

#
geneset16 <- msigdbr(species="Mus musculus", collection="H") %>%
  dplyr::filter(grepl("HALLMARK_MTORC1_SIGNALING", gs_name))
geneset17 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:REACTOME") %>%
  dplyr::filter(gs_name=="REACTOME_PI3K_AKT_ACTIVATION")
geneset18 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:KEGG_LEGACY") %>%#MEDICUS/LEGACY
  dplyr::filter(gs_name=="KEGG_INSULIN_SIGNALING_PATHWAY")
geneset19 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:REACTOME") %>%
  dplyr::filter(gs_name=="REACTOME_ACTIVATION_OF_AMPK_DOWNSTREAM_OF_NMDARS")
geneset20 <- msigdbr(species="Mus musculus", collection="C2") %>%
  dplyr::filter(grepl("WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1", gs_name))
geneset21 <- msigdbr(species="Mus musculus", collection="C2", 
                     subcollection="CP:REACTOME") %>%
  dplyr::filter(gs_name=="REACTOME_FOXO_MEDIATED_TRANSCRIPTION")
geneset22 <- msigdbr(species="Mus musculus", collection="C5", 
                     subcollection="GO:BP") %>%
  dplyr::filter(gs_name=="GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION")

if(T){
list <- list(geneset1,geneset2,geneset3,geneset4,geneset5,
             geneset6,geneset7,geneset8,geneset9,geneset10,
             geneset11,geneset12,geneset13,geneset14,geneset15,
             geneset16,geneset17,geneset18,geneset19,geneset20,
             geneset21,geneset22)
pathwayname <- c('HALLMARK_G2M_CHECKPOINT',
                 'GOBP_MITOTIC_CELL_CYCLE_CHECKPOINT',
                 'REACTOME_CHROMATIN_ORGANIZATION',
                 'GOBP_HISTONE_MRNA_METABOLIC_PROCESS',
                 'GOBP_NUCLEAR_ENVELOPE_ORGANIZATION',
                 'GOCC_NUCLEAR_LAMINA', 
                 'REACTOME_DNA_REPLICATION',
                 'GOBP_MRNA_PROCESSING',
                 'REACTOME_NONSENSE_MEDIATED_DECAY_NMD',
                 'KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY',
                 'REACTOME_CELLULAR_SENESCENCE',
                 'HALLMARK_P53_PATHWAY',
                 'HALLMARK_APOPTOSIS',
                 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                 'GOBP_AUTOPHAGY',
                 'HALLMARK_MTORC1_SIGNALING',
                 'REACTOME_PI3K_AKT_ACTIVATION',
                 'KEGG_INSULIN_SIGNALING_PATHWAY',
                 'REACTOME_ACTIVATION_OF_AMPK_DOWNSTREAM_OF_NMDARS',
                 'WANG_ADIPOGENIC_GENES_REPRESSED_BY_SIRT1',
                 'REACTOME_FOXO_MEDIATED_TRANSCRIPTION',
                 'GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION'
                 )
}
pathwaylist <- list()
for (i in c(1:22)) {
  print(i)
  geneset <- list[[i]]$gene_symbol
  a <- pathwayname[i]
  pathwaylist[[a]] <- geneset
  print(a)
  print(length(geneset))
}
saveRDS(pathwaylist,'pathwaylist.rds')

pathwaylist <- readRDS('pathwaylist.rds')
gsvaPar <- ssgseaParam(exprData = expr, 
                       geneSets = pathwaylist,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)
gsva_data[1:3,1:3]#head(gsva_data)
gsva_mat = as.data.frame(gsva_data)
gsva_mat_t = gsva_mat %>% t() %>% as.data.frame() %>% 
  mutate(sample_id = rownames(.))
gsva_mat_t$group <- m
#colnames(gsva_mat_t) <- str_to_title(colnames(gsva_mat_t))
pathwaylist <- colnames(gsva_mat_t)
ppl1 <- c('Adolescent',"Adult", "Middle_Aged", "Old")
gsva_mat_t$group[which(gsva_mat_t$group=='B1')] <-'Adolescent' 
gsva_mat_t$group[which(gsva_mat_t$group=='B2')] <-'Adult' 
gsva_mat_t$group[which(gsva_mat_t$group=='B3')] <-'Middle_Aged' 
gsva_mat_t$group[which(gsva_mat_t$group=='B4')] <-'Old' 
saveRDS(gsva_mat_t,rds14)
##cor####
library(corrplot)
result.k <- readRDS(rds11)
for (i in 3:length(result.k)) {
  result.k[,i] <- result.k[,i]/sum(result.k[,i])
}
head(result.k)
gsva_mat_t <- readRDS(rds13)
gsva_mat_t2 <- readRDS(rds14)
result.k <- result.k[which(result.k$class_id %in% 
                             c('SINE','LINE','LTR')),]
result.k <- result.k[,-1]
table(result.k$family_id)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k <-result.k[,c('Alu','B4','B2','MIR',
                       'ID','5S-Deu-L2',
                       'L1','L2','CR1','RTE-Bo.B',
                       'Penelope','I-Jockey',
                       'ERVK','ERVL-MaLR',
                       'ERVL','ERV1',
                       'LTR','Gypsy')]
table(rownames(result.k)==rownames(gsva_mat_t))
ym <- cbind(result.k,gsva_mat_t[,c(2:4)])
length(ym)
#ym <- ym[,-c(22,23)]
ym <- cbind(ym,gsva_mat_t2)
length(ym)
ym <- ym[,-c(44,45)]
for (i in 1:length(ym)) {
  ym[,i] <- as.numeric(ym[,i])
}
tdc <- cor (ym, method="spearman")
addcol <- colorRampPalette(c("#0000ff", "#4040ff", "white","#ff4040", "#ff0000"))
testRes = cor.mtest(ym, method="spearman",conf.level = 0.95)
#corrplot(tdc, method = "ellipse", 
#         type = "upper",
#         tl.col = "black", tl.cex = 1.2, tl.srt = 45)
CairoPDF(bpdfy8, width = 12, height = 10,family = 'Times New Roman')
corrplot(tdc[c(19:43),c(1:18)],  col = addcol(100), #method = "ellipse",
         tl.col = "black", tl.cex = 1, #tl.srt = 45,
         p.mat = testRes$p[c(19:43),c(1:18)], diag = T, #type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex =0.8,
         insig = 'label_sig', pch.col = 'grey20')
dev.off()

