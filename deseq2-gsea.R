######################################Deseq2 Bovine Embryo############################

getwd()
setwd("/mnt/data/home/sarahsczelecki/osm/output-files/final")
library(tidyverse)
library(ggplot2)
library(tidyr)
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(edgeR)
library(pheatmap)
library(ComplexHeatmap)

##create required variables
FDR <- 0.05 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis
group_colours <- c(
  "Control" = "#A09E9F", 
  "DNOSM"   = "#1E88E5",  
  "COCSM"   = "#FFC107",  
  "CO-CUL"  = "#004D40") 

## load count file
count_file <- read.csv("/mnt/data/home/sarahsczelecki/osm/osm_total_counts_final.csv", sep=',', header = TRUE, row.names = "Geneid")
metadata <- read.csv("/mnt/data/home/sarahsczelecki/osm/metadata.csv", sep=',', header = TRUE)

#reverse order of metadata
rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1]
metadata <- t(metadata)
metadata <- as.data.frame(metadata)

#check metadata file and count file match up
all(rownames(metadata) == colnames(count_file))

#use EdgeR's filter by  Expr to determine which genes should stay - or else you get a lot of DEGs downstream
# Create DGEList object
dge <- edgeR::DGEList(counts = count_file, group = metadata$group, min.count = 10)

# Use filterByExpr - automatically determines appropriate thresholds
keep <- edgeR::filterByExpr(dge, group = metadata$group)

dge <- dge[keep, , keep.lib.sizes = FALSE]

# Extract filtered counts
count_file <- dge$counts

print(paste("Genes retained:", sum(keep)))
print(paste("Genes filtered:", sum(!keep)))

count_file_cpm <- edgeR::cpm(dge)
print(count_file_cpm)

#log2 transformation of the cpm + c (c is the pseudo count so that log2 isn't performed on 0, so 0 = 2)
pseudo <- 1
count_file_cpm <- log2(count_file_cpm + pseudo)
rownames(count_file_cpm) <- rownames(count_file)

write.csv(count_file_cpm, "osm_logcpmc.csv", row.names = TRUE)

## Creation of the DESeqDataSet to include metadata information
count_file_dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_file,
                                                 colData = metadata,
                                                 design = ~ group)

#set the factor level, tell Deseq which level to compare against
count_file_dds$group <- relevel(count_file_dds$group, ref = "Control")
order_group <- c("Control", "DNOSM", "COCSM", "CO-CUL")
count_file_dds$group <- factor(count_file_dds$group, levels = order_group)

##data transformation for PCA and Corr Plots of data
#you can choose multiple options, but i will go with VST in the deseq package
count_file_dds_vst <- vst(count_file_dds, blind = TRUE)

#plot PCA with Deseq2, but you can't do two groups
#DESeq2::plotPCA(count_file_dds_vst, intgroup = c("condition"))

#so generate the PCA plot manually with ggplot
pcaData <- plotPCA(count_file_dds_vst, intgroup = c("group"), 
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3, alpha = 1.0) + 
  scale_color_manual(values = group_colours) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

print(pcaplot)

ggsave(filename = "pca_small.png", plot = pcaplot, width = 4, height = 3, dpi = 800)
ggsave(filename = "pca_big.png", plot = pcaplot, width = 8, height = 6, dpi = 800)


#plot correlation (non complex heatmap)
count_file_mat_vst <- assay(count_file_dds_vst) #extract the vst matrix from the object
corr_value <- cor(count_file_mat_vst) #compute pairwise correlation values

corrplot <- pheatmap(
  corr_value,
  annotation = dplyr::select(metadata, group),
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Correlation Values of Sample Transcriptomes")

print(corrplot)
ggsave(filename = "corr_small.png", plot = corrplot, width = 4, height = 3, dpi = 800)
ggsave(filename = "corr_big.png", plot = corrplot, width = 8, height = 6, dpi = 800)

#######################ComplexHeatmap Correlation Plot############################
####making a heatmap using the ComplexHeatmap function
#need to make the annotation bars for heatmap annotation
#top annotation
# Relevel the metadata to match the PCA plot
metadata$group <- factor(metadata$group, levels = order_group)

# Top annotation (columns)
ha_top <- HeatmapAnnotation(
  group = metadata$group,
  col = list(group = group_colours),
  annotation_height = unit(3, "mm"),
  border = TRUE,
  annotation_legend_param = list(
    group = list(
      nrow = 5,
      title = "Group",
      title_position = "topleft",
      legend_direction = "horizontal",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12, fontface = "plain")
    )
  )
)
# Row annotation
# Row annotation (only group)
ha_row <- HeatmapAnnotation(
  which = "row",
  group = metadata$group,
  col = list(group = group_colours),
  annotation_width = unit(0.8, 'cm'),
  border = TRUE,
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    group = list(
      nrow = 5,
      title = "Group",
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain')
    )
  )
)
# Create the heatmap
corrplot2 <- ComplexHeatmap::Heatmap(
  corr_value, 
  name = "Correlation of Expressed mRNAs", 
  top_annotation = ha_top,
  right_annotation = ha_row, 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  row_names_gp = gpar(fontsize = 10), 
  heatmap_legend_param = list(
    color_bar = "continuous", 
    legend_direction = "vertical", 
    title = "Correlation",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12, fontface = "plain"),  # Fixed typo: fonface -> fontface
    grid_width = unit(3, "mm"),
    grid_height = unit(3, "mm")
  ), 
  rect_gp = gpar(col = "grey10", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(fill = NA, col = "grey10", lwd = 0.5))
  }
)

# Save to file (optional - uncomment to use)
# png(file = "E:/paper-files/mlcm_small_corr.png", width = 1200, height = 800, res = 800)

# Draw the heatmap
draw(corrplot2,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

# Close the graphics device (if saving to file)
# dev.off()


########################################Deseq2 DEGS###############################
##then do the actual differential gene expression analysis
# Install
#BiocManager::install("org.Mm.eg.db")
# Install biomaRt from Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("biomaRt")

library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(DESeq2)
library(tibble)
library(dplyr)

# Annotation (unchanged)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
gene_list <- rownames(count_file_dds)

UCD2 <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id',
                 'gene_biotype', 'description'),
  filters = 'ensembl_gene_id',
  values = gene_list,
  mart = ensembl
)

colnames(UCD2) <- c("ensgene", "symbol", "entrez", "biotype", "description")

# DESeq2
count_file_dds <- DESeq2::estimateSizeFactors(count_file_dds)
dds <- DESeq2::DESeq(count_file_dds)

# Contrasts
contrasts_list <- list(
 DNOSM_vs_Control = c("group", "DNOSM", "Control"),
 COCSM_vs_Control = c("group", "COCSM", "Control"),
  COCUL_vs_Control = c("group", "CO-CUL", "Control"),
  COCUL_vs_COCSM = c("group", "CO-CUL", "COCSM")
)

# Loop through contrasts
for (comparison in names(contrasts_list)) {
  
  # Get ALL DESeq2 results
  res <- DESeq2::results(dds, contrast = contrasts_list[[comparison]],
                         alpha = alpha,
                         independentFiltering = TRUE)
  res_df <- rownames_to_column(as.data.frame(res), var = "ensgene")
  
  # Annotate all genes
  res_anno <- left_join(
    res_df,
    UCD2[, c("ensgene", "symbol", "entrez", "biotype", "description")],
    by = "ensgene"
  )
  
  # Save ONE CSV only
  write.csv(
    res_anno,
    file = paste0(comparison, "_all_DEGs_annotated.csv"),
    row.names = FALSE
  )
  
  cat("Saved:", paste0(comparison, "_all_DEGs_annotated.csv\n"))
}


#annotating logcpm_count file:
library(tibble)
count_file_cpm <- as.data.frame(count_file_cpm)
count_file_cpm <- rownames_to_column(count_file_cpm, var = "ensgene")

count_file_cpm_anno <- left_join(
  count_file_cpm,
  UCD2[, c("ensgene", "symbol", "entrez", "biotype", "description")],
  by = "ensgene"
)
View(count_file_cpm_anno)
write.csv(count_file_cpm_anno, "osm_logcpmc_anno.csv", row.names = TRUE)







#visualise genes - just to check because this is count, which I assume is not normalised
gene <- plotCounts(dds, gene="ENSMUSG00000074899", intgroup="group", returnData=TRUE)

ggplot(gene, aes(x = group, y = count, colour = group)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_y_log10() + 
  scale_colour_manual(values = group_colours) +
  theme_bw()


######visualise normalised data!!
vsd <- vst(dds)

gene_vst <- data.frame(
  sample = colnames(vsd),
  expression = assay(vsd)["ENSMUSG00000074899", ],
  group = dds$group
)

ggplot(gene_vst, aes(x = group, y = expression, colour = group, fill = group)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_colour_manual(values = group_colours) +
  scale_fill_manual(values = group_colours) +
  theme_bw()





###GSEA preparation - to use with local software
# Create .rnk files for all contrasts
for (comparison in names(contrasts_list)) {
  
  # Get DESeq2 results for this contrast
  res <- DESeq2::results(dds, contrast = contrasts_list[[comparison]],
                         alpha = alpha,
                         independentFiltering = TRUE)
  res_df <- as.data.frame(res)
  
  # Add gene ID as a column
  res_df <- tibble::rownames_to_column(res_df, var = "gene")
  
  # Drop genes without a ranking value
  res_df <- dplyr::filter(res_df, !is.na(stat))
  
  # Keep only the columns needed
  rnk <- dplyr::select(res_df, gene, stat)
  
  # Save as .rnk (tab-delimited, NO header)
  write.table(
    rnk,
    file = paste0(comparison, ".rnk"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  cat("Saved:", paste0(comparison, ".rnk\n"))
}

# ------------------------------
# GSEA Preparation
# ------------------------------
rnk_list <- list()
for (comparison in names(contrasts_list)) {
  
  # 1 Get DESeq2 results
  res <- DESeq2::results(
    dds,
    contrast = contrasts_list[[comparison]],
    alpha = alpha,
    independentFiltering = TRUE
  )
  
  # 2 Convert to data frame and add gene symbols
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    dplyr::filter(!is.na(stat))
  
  # 3 Map gene IDs
  gene_ids <- bitr(
    res_df$ENSEMBL,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  # 4 Join with DE results
  res_df <- res_df %>%
    dplyr::inner_join(gene_ids, by = "ENSEMBL")
  
  # 5 Remove duplicates by max |stat|
  res_df <- res_df %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # 6 Create ranked list
  gene_list <- res_df$stat
  names(gene_list) <- res_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Store for later use
  rnk_list[[comparison]] <- gene_list
  
  # Save .rnk file
  write.table(
    data.frame(ENTREZID = names(gene_list), stat = gene_list),
    file = paste0(comparison, ".rnk"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  
  cat("✓", comparison, ":", length(gene_list), "genes\n")
}

# VERIFY rnk_list was created properly
cat("\n=== Checking rnk_list ===\n")
cat("Number of comparisons in rnk_list:", length(rnk_list), "\n")
cat("Comparison names:", paste(names(rnk_list), collapse = ", "), "\n\n")

# --------------------------------------------------------------------
# KEGG + GO - Loop over all comparisons
# --------------------------------------------------------------------
output_dir <- "/mnt/data/home/sarahsczelecki/osm/output-files/final"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (comparison in names(rnk_list)) {
  
  cat("\n========================================\n")
  cat("Processing:", comparison, "\n")
  cat("========================================\n")
  
  gene_list <- rnk_list[[comparison]]
  cat("Gene list length:", length(gene_list), "\n")
  
  # ------------------------------
  # KEGG GSEA
  # ------------------------------
  cat("Running KEGG GSEA...\n")
  
  tryCatch({
    kegg_gsea <- gseKEGG(
      geneList     = gene_list,
      organism     = "mmu",
      minGSSize    = 10,
      pvalueCutoff = 0.25,
      verbose      = FALSE
    )
    
    if (!is.null(kegg_gsea) && nrow(kegg_gsea) > 0) {
      kegg_table <- as.data.frame(kegg_gsea)
      write.csv(kegg_table,
                file = file.path(output_dir, paste0(comparison, "_KEGG.csv")),
                row.names = FALSE)
      
      png(file.path(output_dir, paste0(comparison, "_KEGG_dotplot.png")),
          width = 1200, height = 800, res = 150)
      print(dotplot(kegg_gsea, showCategory = 20) +
              ggtitle(paste0(comparison, " KEGG")))
      dev.off()
      
      cat("✓ KEGG results saved:", nrow(kegg_table), "pathways\n")
    } else {
      cat("⚠ No significant KEGG results\n")
    }
  }, error = function(e) {
    cat("✗ KEGG GSEA failed:", conditionMessage(e), "\n")
  })
  
  # ------------------------------
  # GO GSEA (BP, CC, MF)
  # ------------------------------
  ontologies <- c("BP", "CC", "MF")
  
  # Initialize list to store all GO results for this comparison
  all_go_results <- list()
  
  for (ont in ontologies) {
    cat("Running GO", ont, "GSEA...\n")
    
    tryCatch({
      go_gsea <- gseGO(
        geneList     = gene_list,
        OrgDb        = org.Mm.eg.db,
        keyType      = "ENTREZID",
        ont          = ont,
        minGSSize    = 10,
        pvalueCutoff = 0.25,
        verbose      = FALSE,
        eps          = 0
      )
      
      if (!is.null(go_gsea) && nrow(go_gsea) > 0) {
        go_table <- as.data.frame(go_gsea)
        
        # Add ontology column to identify the source
        go_table$Ontology <- ont
        
        # Store in list for combining later
        all_go_results[[ont]] <- go_table
        
        # Still save individual files
        write.csv(
          go_table,
          file = file.path(output_dir,
                           paste0(comparison, "_GO_", ont, ".csv")),
          row.names = FALSE
        )
        
        png(file.path(output_dir,
                      paste0(comparison, "_GO_", ont, "_dotplot.png")),
            width = 1200, height = 800, res = 150)
        print(dotplot(go_gsea, showCategory = 20) +
                ggtitle(paste0(comparison, " GO ", ont)))
        dev.off()
        
        cat("✓ GO", ont, "results saved:", nrow(go_table), "terms\n")
      } else {
        cat("⚠ No significant GO", ont, "results\n")
      }
    }, error = function(e) {
      cat("✗ GO", ont, "GSEA failed:", conditionMessage(e), "\n")
    })
  }
  
  # ------------------------------
  # Combine all GO results into one file
  # ------------------------------
  if (length(all_go_results) > 0) {
    combined_go <- dplyr::bind_rows(all_go_results)
    
    # Reorder columns to put Ontology near the front
    col_order <- c("ID", "Description", "Ontology", 
                   setdiff(names(combined_go), c("ID", "Description", "Ontology")))
    combined_go <- combined_go[, col_order]
    
    # Save combined file
    write.csv(
      combined_go,
      file = file.path(output_dir, paste0(comparison, "_GO_ALL.csv")),
      row.names = FALSE
    )
    
    cat("✓ Combined GO file saved:", nrow(combined_go), "total terms\n")
  } else {
    cat("⚠ No GO results to combine\n")
  }
}

cat("\n=== All analyses complete ===\n")