#!/usr/bin/env Rscript

library(jsonlite)
library(tximport)
library(DESeq2)
library(data.table)
library(readr)

# Get inputs, outputs, and parameters from Snakemake
input_files <- snakemake@input[['quants']]
output_file <- snakemake@output[[1]]
conditionA <- snakemake@params[['conditionA']]
conditionB <- snakemake@params[['conditionB']]
tx2gene_file <- snakemake@input[['tx2gene']]

# Read the transcript-to-gene mapping file
tx2gene <- fread(tx2gene_file)

# Read in Salmon quant files
names(input_files) <- c(conditionA, conditionB)

# Ensure the names vector matches the length of files
sample_names <- c(paste0("conditionA_", seq_along(conditionA)), paste0("conditionB_", seq_along(conditionB)))
names(input_files) <- sample_names

txi <- tximport(input_files, type="salmon", tx2gene=tx2gene)

coldata <- data.frame(
  condition=factor(c(rep("conditionA", length(conditionA)), rep("conditionB", length(conditionB)))),
  row.names=names(input_files)
)

dds <- DESeqDataSetFromTximport(txi, colData=coldata, design= ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res_table <- data.table(as.data.frame(res))
res_table$gene_ID <- rownames(as.data.frame(res))

fwrite(res_table[ , .(gene_ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)], output_file, sep='\t')
