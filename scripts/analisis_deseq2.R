
# Instalamos lo necesario si no está instalado
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!require("knitr", quietly = TRUE)) install.packages("knitr")

# Cargamos las librerías
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(knitr)

# CARGAMOS DATOS Y LIMPIAMOS
counts_data <- read.table("C:/Users/anton/Desktop/Master/Segundo cuatri/ADO/Entrega/analisis_GSE276955_resultados/rnaseq_GSE276955/results/counts/gene_counts_matrix.txt", 
                          header=TRUE, skip=1, row.names=1)

# Extraemos columnas de las muestras (de la 6 en adelante)
count_matrix <- counts_data[, 6:ncol(counts_data)]
colnames(count_matrix) <- gsub(".bam", "", colnames(count_matrix))

#ANALISIS ESTADISTICO
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

#ANOTACION GENES
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
res$plot_label <- ifelse(is.na(res$symbol), rownames(res), res$symbol)
write.csv(as.data.frame(res), file="Resultados_Completos_UNC1999_vs_DMSO.csv")

# VOLCANO PLOT
pdf("VolcanoPlot_Final_UNC1999_Symbols.pdf", width=12, height=10)

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

# TOP 10 GENES UP y TOP 10 DOWN
res_df <- as.data.frame(res)
res_sig <- na.omit(res_df)
res_sig <- res_sig[res_sig$padj < 0.05, ] 

# Extraemos y unimos
top10_up <- res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ][1:10, ]
top10_down <- res_sig[order(res_sig$log2FoldChange, decreasing = FALSE), ][1:10, ]
top_20_genes <- rbind(top10_up, top10_down)

# Formateamos la tabla final
tabla_final <- top_20_genes[, c("plot_label", "baseMean", "log2FoldChange", "padj")]
colnames(tabla_final) <- c("Gen", "Expresión_Media", "Log2_FoldChange", "P_Valor_Ajustado")

# Imprimimos en consola, guardamos CSV 
print(kable(tabla_final, digits = 3, caption = "Top 10 genes UP y DOWN (UNC1999 vs DMSO)"))
write.csv(tabla_final, file = "Top20_Genes_UNC1999.csv", row.names = FALSE)

message(paste("Archivos guardados en:", getwd()))