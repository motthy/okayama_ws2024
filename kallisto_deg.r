if (!require("BiocManager",  quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")

library(tximport)
library(tidyverse)

# count table
drr <- paste0("DRR", 357080:357084)
file_kallisto <- file.path("data/kallisto", drr, "abundance.tsv")
names(file_kallisto) <- drr
txi.kallisto <- tximport(file_kallisto,
                         type = "kallisto",
                         txOut = TRUE)
# transcripts-gene
tx2gene <- data.frame(
  TXNAME = rownames(txi.kallisto$counts),
  GENEID = sapply(str_split(rownames(txi.kallisto$counts), '\\|'), '[',  2)
)

write_tsv(tx2gene,
          file = "data/kallisto/tx2gene.tsv")

# geneレベルの発現量(raw counts)
gene.exp_raw <- summarizeToGene(txi.kallisto,
                                tx2gene)
#head(gene.exp_raw$counts)
rawcount <- data.frame(gene.exp_raw$counts) %>%
  rownames_to_column() %>%
  separate(rowname,  "\\.",  into = c("gene_id",  "postnum")) %>%
  dplyr::select(!postnum)

write_tsv(rawcount,
          file = "data/kallisto/kallisto_count.tsv")

# geneレベルの発現量(scaledTPM)
# for edgeR or DESeq2
gene.exp <- summarizeToGene(txi.kallisto,
                            tx2gene,
                            countsFromAbundance = "scaledTPM")
#head(gene.exp$counts)
tpm <- data.frame(gene.exp$counts) %>%
  rownames_to_column() %>%
  separate(rowname,  "\\.",  into = c("gene_id",  "postnum")) %>%
  dplyr::select(!postnum)

write_tsv(tpm,
          file = "data/kallisto/kallisto_scaledTPM.tsv")

# DEG抽出
# 参考　二群間比較（edgeR）
# https://bi.biopapyrus.jp/rnaseq/analysis/de-analysis/2g-edger.html
library(edgeR)

group <- factor(c("C", "C", "T", "T", "T"))
design <- model.matrix(~ group)

#edgeR を利用して解析
d <- DGEList(counts = gene.exp$counts, 
             group = group)

# TMM正規化
d <- calcNormFactors(d) 

#分散の推定を行うときデザイン行列を与える
d <- estimateDisp(d, design)

# 推測された分布に観測値を当てはめる
fit <- glmFit(d, design)

# 尤度比検定
lrt <- glmLRT(fit, coef = 2)

#検定結果のうち p-value がもっとも小さい10遺伝子
topTags(lrt)

top10 <- data.frame(topTags(lrt)) %>%
  rownames_to_column() %>%
  separate(rowname, "\\.", into = c("gene_id", "postnum")) %>%
  select(!postnum)

# ここからは付録
# 遺伝子アノテーションをつける
library(mygene)

query_list <- c(top10$gene_id)
annotation <- queryMany(query_list,scopes = "ensembl.gene")

result <- data.frame(
  annotation$query,
  annotation$symbol,
  annotation$name
)

lt_data <- as.data.frame(topTags(lrt, n=nrow(d)))
lt_dataFDR <- lt_data[lt_data$FDR < 0.05,] %>%
  rownames_to_column() %>%
  separate(rowname,"\\.",into = c("gene_id", "postnum")) %>%
  dplyr::select(!postnum)

# annotation
query_list2 <- lt_dataFDR$gene_id
annotation2 <- queryMany(query_list2,scopes = "ensembl.gene")

result2 <- data.frame(
  annotation2$query,
  annotation2$symbol,
  annotation2$name
) 
# カラム名揃える
colnames(result2) <- c("gene_id", "symbol", "description")

# join
lt_dataFDR_anno <- lt_dataFDR %>%
  left_join(result2, by="gene_id")

# annotation付きDEG結果
write_tsv(lt_dataFDR_anno,
          file = "data/kallisto/kallisto_edgeR_LRT_FDR0.05_withanno.tsv")

