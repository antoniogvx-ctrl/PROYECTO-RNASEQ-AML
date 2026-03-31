# Instalamos lo necesario si no está instalado
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")

# Cargamos las librerías
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# CCargamos datos
counts_data <- read.table("./gene_counts_matrix.txt", 
                          header=TRUE, skip=1, row.names=1)

# Extraemos columnas de las muestras 
count_matrix <- counts_data[, 6:ncol(counts_data)]
colnames(count_matrix) <- gsub(".bam", "", colnames(count_matrix))

#Analisis estadístico
# Crear los Metadatos (ColData)
samples <- colnames(count_matrix)
condition <- factor(c(rep("UNC1999", 3), rep("DMSO", 3)), levels = c("DMSO", "UNC1999"))
coldata <- data.frame(row.names=samples, condition=condition)

# Construir el objeto y ejecutar DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "UNC1999", "DMSO"))

#Anotamos genes
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res$plot_label <- ifelse(is.na(res$symbol), rownames(res), res$symbol)
write.csv(as.data.frame(res), file="Resultados_Completos_UNC1999_vs_DMSO.csv")

# VOLCANO PLOT 
pdf("VolcanoPlotNextFlow.pdf", width=12, height=10)

EnhancedVolcano(res,
                lab = res$plot_label,      
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Targeting Polycomb: UNC1999 vs DMSO',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,            
                FCcutoff = 1.0,            
                pointSize = 3.0,
                labSize = 4.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendLabels = c('NS', expression(Log[2] ~ FC), 'p-value', 
                                 expression(p-value ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                drawConnectors = TRUE,     
                widthConnectors = 0.5,
                colConnectors = 'grey30')

dev.off()

message(paste("Volcano Plot generado con éxito y guardado en:", getwd()))
