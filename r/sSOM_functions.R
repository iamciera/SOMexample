### Visualize by Cluster

# Read in data used for GO enrichment analysis

geneLength <- read.csv("../data//normalized_genes_length.csv")
cate <- read.table("../data/melted.GOTable.txt", header=TRUE)


#clusterVis Function
#displays transformed data in a box plot 
clusterVis <- function(clustNum){
  
  sub_cluster <- subset(plot.data, ssom.unit.classif==clustNum)
  sub_data <- sub_cluster[,c(1, 9:14)] # just the sample types
  m.data <- melt(sub_data) 
  p <- ggplot(m.data, aes(x=variable, y=value, color = genotype))
  p + geom_point(alpha=0.5,position="jitter", size=1) + 
    geom_boxplot(alpha=0.75, outlier.size=0) + 
    theme_bw()
}

clusterGO <- function(clustNum){
  ##GO Enrichment on the catergories
  dev.off()
  plot.new()
  
  #we need to first get the data in the right format.
  #First get the list of ITAG
  
  #sub_cluster
  sub_cluster <- subset(plot.data, ssom.unit.classif==clustNum)
  
  itag.sc <- as.data.frame(sub_cluster$gene) 
  colnames(itag.sc)[1] <- "itag"
  itag.sc$sc <- 1    
  
  #Since each orthologue between tf2 and wt are represented twice in this set, 
  #we have to keep only the unique ITAGs.
  
  itag.sc <- unique(itag.sc) #Check. 
  #Should cut the list in half. # dim(itag.sc) before and after
  
  #Merge all by itag
  matrixGO <- merge(itag.sc, geneLength, by = "itag", all = TRUE)
  matrixGO[is.na(matrixGO)] <- 0
  pat <- matrixGO
  
  #Now that we have the data in the right format we can proceed with GO enrichment.
  
  genes = as.integer(pat[,"sc"])
  names(genes) = pat$itag
  table(genes)
  length(genes)
  
  pwf = nullp(genes,bias.data=pat$length)
  
  GO.wall = goseq(pwf,gene2cat = cate)
  head(GO.wall)
  
  #This is going to correct for multiple testing.  
  #You can specify the p-value cut-off of GO categories you are interested.
  
  enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < 0.05]
  enriched.GO
  
  my.GO <- as.character(enriched.GO)
  my.GO.table <- Term(my.GO)
  my.GO.table
  t <- as.matrix(my.GO.table)
  
  print(t) #this is for the knitr document
}

clusterVis_line <- function(clustNum) {
  sub_cluster <- subset(plot.data, ssom.unit.classif==clustNum)
  sub_data <- sub_cluster[,c(1, 2, 9:14)] # just the sample types
  sub_data <- melt(sub_data)
  sub_data <- within(sub_data, lineGroup <- paste(genotype, gene,sep='.'))
  ggplot(sub_data, aes(variable, value, group = lineGroup, color =  genotype )) + 
    geom_line(alpha = .1, (aes(color = factor(genotype)))) + 
    geom_point(alpha = .0) + 
    theme_bw()
}

#Prereq annotation files for function

annotation1<- read.delim("../data/ITAG2.3_all_Arabidopsis_ITAG_annotations.tsv", header=FALSE)  #Changed to the SGN human readable annotation
colnames(annotation1) <- c("ITAG", "SGN_annotation")
annotation2 <- read.delim("../data/ITAG2.3_all_Arabidopsis_annotated.tsv")
annotation <- merge(annotation1,annotation2, by = "ITAG")

#Only Gene Name and ITAG
annotation <- annotation[,c(1,5)]

genesInClust <- function(clustNum, plot.data, annotation) {
  sub_cluster <- subset(plot.data, ssom.unit.classif==clustNum)
  sub_data <- as.data.frame(sub_cluster[,2])
  colnames(sub_data) <- "ITAG"
  resultsTable <- merge(sub_data,annotation,by = "ITAG", all.x=TRUE)
  print(nrow(unique(resultsTable)))
  return(unique(resultsTable))
}

clusterVis_PCA <- function(clustNum) {
  
  #make dataset for visualization
  data.val3 <- plot.data
  names(data.val3)
  data.val3$cluster[data.val3[,21] == clustNum] <- "subcluster"
  data.val3$cluster[data.val3[,21] != clustNum] <- "other"
  
  #plot
  
  p <- ggplot(data.val3, aes(PC1, PC2, color = cluster)) 
  p + geom_point(size=I(2), alpha = 0.6) +
    scale_colour_manual(values=c("#cccccc", "#000000")) + 
    theme_bw() + 
    theme(legend.text = element_text(
      size = 16, 
      face = "bold")) + 
    facet_grid(. ~ genotype)
}

clusterVis_PCAsub <- function(clustNum) {
  
  
  #make dataset for visualization
  plot.data <- subset(plot.data, ssom.unit.classif==clustNum)
  data.val3 <- plot.data
  
  #plot
  
  p <- ggplot(data.val3, aes(PC1, PC2, color = genotype)) 
  p + geom_point(size=I(2), alpha = 0.6) +
    scale_colour_manual(values=c("#ef8a62", "#67a9cf")) + 
    theme_bw() + 
    theme(legend.text = element_text(
      size = 16, 
      face = "bold"))
}



