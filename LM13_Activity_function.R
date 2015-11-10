##### FUNCTIONS DENEFLAB GITHUB ########



##### Normalization #######

# Better rounding function than R's base round
matround <- function(x){trunc(x+0.5)}


# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}


##### ADONIS ###########

# Function to run adonis test on a phyloseq object and a variable from metadata
# Make sure OTU data is standardized/normalized before 
phyloseq_to_adonis <- function(physeq, distmat = NULL, dist = "bray", formula) {
  
  if(!is.null(distmat)){
    phydist <- distmat
  } else {
    phydist <- phyloseq::distance(physeq, dist)
  }
  
  metadata <- as(sample_data(physeq), "data.frame")
  
  # Adonis test
  f <- reformulate(formula, response = "phydist")
  adonis.test <- adonis(f, data = metadata)
  print(adonis.test)
  
  # Run homogeneity of dispersion test if there is only 1 variable
  if (length(formula) == 1) {
    
    group <- metadata[,formula]
    beta <- betadisper(phydist, group)
    disper.test = permutest(beta)
    print(disper.test)
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test, 
      disper = disper.test
    )
    
  } else {
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test
    )
  }
  return (l)
}

########## Bar Plots #################

# This function takes a phyloseq object, agglomerates OTUs to the desired taxonomic rank, 
# prunes out OTUs below a certain relative proportion in a sample (ie 1% ) 
# and melts the phyloseq object into long format which is suitable for ggplot stacked barplots.
taxglom_and_melt <- function(physeq, taxrank, prune){
  
  # Agglomerate all otu's by given taxonomic level
  pglom <- tax_glom(physeq, taxrank = taxrank)
  
  # Create a new phyloseq object which removes taxa from each sample below the prune parameter
  pglom_prune <- transform_sample_counts(pglom,function(x) {x/sum(x)})
  otu_table(pglom_prune)[otu_table(pglom_prune) < prune] <- 0
  pglom_prune <- prune_taxa(taxa_sums(pglom_prune) > 0, pglom_prune)
  
  # Melt into long format and sort by taxonomy
  physeq_long <- psmelt(pglom_prune)
  physeq_long <- physeq_long[order(physeq_long[ ,taxrank]), ]
  
  # Return long data frame
  return(physeq_long)
}


###### Merge functions ############

# Merge samples by averaging OTU counts instead of summing
merge_samples_mean <- function(physeq, group, round){
  # Calculate the number of samples in each group
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  
  # Merge samples by summing
  merged <- merge_samples(physeq, group)
  
  # Divide summed OTU counts by number of samples in each group to get mean
  # Calculation is done while taxa are columns
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  
  # Pick the rounding functions
  if (round == "floor"){
    out <- floor(t(x/group_sums))
  } else if (round == "round"){
    out <- matround(t(x/group_sums))
  }
  
  # Return new phyloseq object with taxa as rows
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

# Merge samples, just including OTUs that were present in all merged samples
# Call this function before running merge_samples()
merge_OTU_intersect <- function(physeq, group){
  
  # Make sure we're not starting with more taxa than we need 
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  
  s <- data.frame(sample_data(physeq))
  l <- levels(s[,group])
  o <- otu_table(physeq)
  
  # Loop through each category
  for (cat in 1:length(l)) {
    
    # Get the index of all samples in that category
    w <- which(s[,group]==l[cat])
    
    # subset to just those columns of OTU table
    cat.sub<-o[,w]
    print(dim(cat.sub))
    
    # Find the indices of 0's in the OTU table
    zeros <- apply(cat.sub, 1, function(r) any(r == 0))
    
    # If an OTU had a 0 in at least one sample, change all samples to 0
    cat.sub[zeros,] <- 0
  }
  
  o[,w] <- cat.sub
  otu_table(physeq) <- o
  
  return(physeq)
  
}

#################################################################################### 4 + 5
#OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU OTU 
## I wrote these 2 functions off of the following tutorial from the Phyloseq GitHub page:
#http://joey711.github.io/phyloseq-extensions/DESeq2.html

deSEQ <- function(data, valuetest){
  data_pruned=prune_taxa(taxa_sums(data)>(147*nrow(sample_data(data))),data)
  de_data = phyloseq_to_deseq2(data_pruned, valuetest)
  de_data2 = DESeq(de_data, test="Wald", fitType="parametric")
  res_data = results(de_data2, cooksCutoff = FALSE)
  alpha = 0.01
  sig_data = res_data[which(res_data$padj < alpha), ]
  sigtab_sherm = cbind(as(sig_data, "data.frame"), as(tax_table(data_pruned)[rownames(sig_data), ], "matrix"))
} 

deSEQ_noprune <- function(data, valuetest){
  data=prune_taxa(taxa_sums(data)>0,data)
  de_data = phyloseq_to_deseq2(data, valuetest)
  de_data2 = DESeq(de_data, test="Wald", fitType="parametric")
  res_data = results(de_data2, cooksCutoff = FALSE, contrast=c("DNA","cD","D"))
  plotMA(res_data)
  alpha = 1
  sig_data = res_data[which(res_data$padj < alpha), ]
  sigtab_sherm = cbind(as(sig_data, "data.frame"), as(tax_table(data)[rownames(sig_data), ], "matrix"))
} 

plot_deSEQ <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Species, function(x) max(x))
  y = sort(y, TRUE)
  deSEQdata$Species = factor(as.character(deSEQdata$Species), levels=names(y))
  ggplot(deSEQdata, aes(x=Phylum, y=log2FoldChange, color=Phylum, size=plotvalue, shape=sig)) + 
    geom_point(alpha=0.9) + theme_bw() +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    #    scale_color_manual(name="p-value",  breaks = c("0", "1"), labels = c("p>0.05", "p<0.05"), values = c("0" = "grey", "1"="black")) +
    scale_shape_manual(name = "p-value", breaks = c("0", "1"), 
                       labels = c("p>0.05", "p<0.05"),
                       values = c("0" = 21, "1"= 19)) +
    scale_size_continuous("abundance") +
    scale_y_continuous(breaks=seq(-15,15,1),limits=c(-7,7)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}


#################################################################################### 4 + 5
#################################################################################### 4 + 5

plot_deSEQ_combo <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Species, function(x) max(x))
  y = sort(y, TRUE)
  deSEQdata$Species = factor(as.character(deSEQdata$Species), levels=names(y))
  ggplot(deSEQdata, aes(x=Phylum_plot, y=log2FoldChange, color=Phylum, size=plotvalue, shape=sig)) + 
    geom_point(alpha=0.9) + theme_bw() +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    #    scale_color_manual(name="p-value",  breaks = c("0", "1"), labels = c("p>0.05", "p<0.05"), values = c("0" = "grey", "1"="black")) +
    scale_shape_manual(name = "p-value", breaks = c("0", "1"), 
                       labels = c("p>0.05", "p<0.05"),
                       values = c("0" = 21, "1"= 19)) +
    scale_size_continuous("abundance") +
    scale_y_continuous(breaks=seq(-15,15,1),limits=c(-7,7)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}
#################################################################################### 6
#Class 
plot_class_deSEQ <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Class, function(x) max(x))
  #  y = sort(y, TRUE)
  deSEQdata$Class = factor(as.character(deSEQdata$Class), levels=names(y))
  ggplot(deSEQdata, aes(x=Class, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=6) +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    scale_y_continuous(breaks=seq(-15,15,1)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}

#Order 
plot_order_deSEQ <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Order, function(x) max(x))
  #  y = sort(y, TRUE)
  deSEQdata$Order = factor(as.character(deSEQdata$Order), levels=names(y))
  ggplot(deSEQdata, aes(x=Order, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=6) +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    scale_y_continuous(breaks=seq(-15,15,1)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}


#################################################################################### 7
#################################################################################### 7
#PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM PHYLUM 
plot_phylum_deSEQ <- function(deSEQdata, title){
  x = tapply(deSEQdata$log2FoldChange, deSEQdata$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  deSEQdata$Phylum = factor(as.character(deSEQdata$Phylum), levels=names(x))
  ggplot(deSEQdata, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=6) +  ggtitle(title) +  
    scale_y_continuous(breaks=seq(-15,15,1)) + 
    scale_color_manual(values = phylum.colors,name="Phylum") +
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}
#################################################################################### 7
#################################################################################### 7

# Multiple plot function is from the R-Cookbook website 
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#################################################################################### 1
#################################################################################### 1

######################### Functions for clustering #########################

## Wrapper function to cluster OTUs, cut the dendrogram into groups and plot a heatmap
# dist.method can be spearman, bray-curtis,  or any distance metric accepted by dist()
# k is the number of groups to split the dendrogram into
heat_wrapper <- function(physeq, taxrank, dist.method, k) {
  
  # Scale OTUs
  scaled.otu <- scale_otu(
    physeq = physeq, 
    taxrank = taxrank
  )
  
  # hierarchical clustering of OTUs using complete linkage
  clustered.otu <- cluster_otu(scaled.otu, dist.method)
  
  # cut the dendrogram into groups
  groups <- cutree(clustered.otu, k = k)[clustered.otu$order]
  
  # Melt to long format
  melted.otu <- melt(scaled.otu, varnames = c("Group_cDvsD", "OTU"))
  
  # Reorder levels of OTUs to match hierarchical clustering
  otu.order <- levels(melted.otu$OTU)[clustered.otu$order]
  melted.otu$OTU <- factor(melted.otu$OTU, levels = otu.order)
  
  return(list(clusters = groups, data = melted.otu))
  
}


# Scales OTU relative abundance and changes the names of the OTUs to taxrank + OTU number. 
# Inputs are a phyloseq object and the taxonomic rank to combine with OTU # for the label 
scale_otu <- function(physeq, taxrank) {
  
  # Scale relative abundance of each OTU
  # Assumes taxa are rows in original physeq object
  scaled.otu <- scale(t(otu_table(physeq)))
  
  # Change the taxa names to be taxrank + OTU number
  colnames(scaled.otu) <- paste(tax_table(physeq)[ ,taxrank], 
                                tax_table(physeq)[ ,"Species"])
  
  # Fix formatting of rownames and column names
  rownames(scaled.otu) <- sample_data(physeq)$Group_cDvsD
  colnames(scaled.otu) <- gsub("[.]", " ", colnames(scaled.otu))
  
  return(scaled.otu)
} 


# Computes dissimilarity between OTUs (spearman, bray, or input to dist() function) and performs 
# hierarchical clustering with complete linkage.
cluster_otu <- function(data, dist.method) {
  if (dist.method == "spearman") {
    data.sim <- cor(data, method = "spearman") 
    # Conversion from spearman correlation to dissimilarity. Problematic ...
    data.dis <- sqrt(2*(1-data.sim))
    data.clust <- hclust(d = as.dist(data.dis), method = "complete")
  } else if (dist.method == "bray") {
    data.dis <- vegdist(t(data), method = "bray")
    data.clust <- hclust(d = data.dis, method = "complete")
  } else {
    data.dis <- dist(t(data), method = dist.method)
    data.clust <- hclust(d = data.dis, method = "complete")
  }
  
  return(data.clust)
  
}


