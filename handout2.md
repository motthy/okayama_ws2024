# 岡山大　バイオインフォマティクスワークショップ
2024-02-06 ~ 2024-02-08

2日目（2024-02-07 10:00-16:00）

- RNA-seqデータの処理
1. rsem
2. kallisto

### RNA-seqデータの処理

結果を入れるディレクトリ作成

```{sh}
mkdir ref rsem kallisto
```

### 1. rsem

#### 1-1. index作成

公式ページ
https://github.com/deweylab/RSEM#built

GENCODE https://www.gencodegenes.org/human/ からリファレンスゲノムをダウンロードしてindexを作成する。

`ref`ディレクトリに移動
```{sh}
cd ref
```

Genome sequence primary assembly(GRCh38)fastaをダウンロード
```{sh}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
```

ファイルを解凍
```{sh}
gunzip GRCh38.primary_assembly.genome.fa.gz
```

アノテーションファイル(GTF)をダウンロード
```{sh}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
```

ファイルを解凍
```
gunzip gencode.v44.annotation.gtf.gz
```

index作成スクリプト
```{rsem_index.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G
#$ -l s_vmem=4G

rsem-prepare-reference --gtf gencode.v44.annotation.gtf --bowtie2 GRCh38.primary_assembly.genome.fa human_gencode

```

これを`rsem_index.sh`として`ref`ディレクトリに保存、`qsub`で実行
```{sh}
qsub rsem_index.sh
```


#### 1-2. マッピングとリードカウント

okayana_wsディレクトリに移動
```{sh}
cd ../okayama_ws
```

アレイジョブで実行するスクリプト
```{rsem_mapping.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G
#$ -l s_vmem=4G
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

rsem-calculate-expression --paired-end -p 4 --bowtie2 --strandedness reverse --estimate-rspd fastq/${DRR}_1.trim.fq.gz fastq/${DRR}_2.trim.fq.gz ref/human_gencode rsem/$DRR

```

`rsem_mapping.sh` として保存し、`qsub`で実行
```{sh}
qsub rsem_mapping.sh
```
（2時間半くらいかかる）

#### 1-3. カウントテーブル作成

実行スクリプト
```{rsem_countmtx.sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs

# gene level
rsem-generate-data-matrix rsem/DRR*genes.results > rsem/ds1_all.genes.results

# transcripts(isoform) level
rsem-generate-data-matrix rsem/DRR*isoforms.results > rsem/ds1_all.isoforms.results
```

`rsem_countmtx.sh`として保存、`qsub`で投入

数分で終わるのでshortノードでOK
```{sh}
qsub ーl short rsem_countmtx.sh
```

`rsem`ディレクトリに移動
```{sh}
cd rsem
```

遺伝子IDは`ENSG00000001084.13`のような形式になっていて、このままではiDEPで遺伝子IDを変換できないので、 ".13" の部分をワンライナーで削除する

```
 "rsem_idep/DRR357080.genes.results" "rsem_idep/DRR357081.genes.results" "rsem_idep/DRR357082.genes.results" "rsem_idep/DRR357083.genes.results" "rsem_idep/DRR357084.genes.results"
"ENSG00000000003.16" 469.00 494.00 353.00 411.00 489.26
"ENSG00000000005.6" 0.00 0.00 0.00 0.00 0.00
```

**ワンライナーによる処理**

gene  levelのカウントテーブル
```
more ds1_all.genes.results |sed -r 's/\.[0-9]+\"/\"/g' > ds1.genes.results.rename.tsv
```

transcript levelのカウントテーブル
```
more ds1.isoforms.results |sed -r 's/\.[0-9]+\"/\"/g' > ds1.isoforms.results.rename.tsv
```

**headerの編集**

テキストエディタでファイルを開いてディレクトリ名の部分を削除し`DRRxxxxxxx`だけ残す。

つぎに、実験条件のcsvをつくる
```
id DRR357080 DRR357081 DRR357082 DRR357083 DRR357084
treatment DMSO DMSO 5H4PB 5H4PB 5H4PB
```
`ds1_experiment.csv`として保存

`ds1_all.genes.results.rename.tsv`と`ds1_experiment.csv`をiDEPで使用する。

### 2. kallisto

以下の例では`~/tools/kallisto`にインストールしたものを利用する。

#### 2-1. index作成

GENCODE https://www.gencodegenes.org/human/ からトランスクリプトームをダウンロードしてindexを作成する。

`ref`ディレクトリに移動
```{sh}
cd ref
```

transcriptome配列をダウンロード
```{sh}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
```

indexを作成するスクリプト
```{kallisto_index.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G
#$ -l s_vmem=4G

~/tools/kallisto/kallisto index -i kallisto_index gencode.v44.transcripts.fa.gz -k 31
```

`kallisto_index.sh`として保存、実行
```
qsub kallisto_index.sh
```
10分以内で完了したのでshortノードでもOK

#### 2-2. リードカウント

結果をいれるディレクトリをつくる
```
mkdir kallisto
```

サンプルごとのディレクトリをまとめてつくる
```{kallisto_dir.sh}
#!/bin/bash

FQLIST="DRR357080 DRR357081 DRR357082 DRR357083 DRR357084"

for SRR in $FQLIST;do
　mkdir kallisto/$SRR
done
```

`kallisto_dir.sh`として保存、実行

確認
```{sh}
ls -l kallisto
```

準備ができたので、カウントを行う
```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

~/tools/kallisto/kallisto quant -i ref/kallisto_index -o kallisto/$DRR -t 4 -b 100 --rf-stranded fastq/${DRR}_1.trim.fq.gz fastq/${DRR}_2.trim.fq.gz

```

`kallisto_count.sh`として保存、`qsub`で実行（shortノードでOK）
```
qsub kallisto_count.sh
```
10分程度で終わる。各サンプル名のディレクトリに以下の3個のファイルができているはずである。

- abundance.h5
- abundance.tsv
- run_info.json

#### 2-3. tximportによるカウントテーブル作成

**参考URL**　tximport edgeR（二群間比較）https://bi.biopapyrus.jp/rnaseq/analysis/

はじめに、`kallisto/DRR*` をlocalにダウンロードする

Rで`tximport`を用いて以下のような処理を行う。講師の環境ではR 4.3.2を使用した。

（RコンソールでもRstudioでもお好みの環境で実行してください）
```{r}
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

BiocManager::install("tximport")

library(tximport)
library(tidyverse)

# count table
drr <- paste0("DRR", 357080:357084)
file_kallisto <- file.path("kallisto", drr, "abundance.tsv")
names(file_kallisto) <- drr
txi.kallisto <- tximport(file_kallisto,
                         type = "kallisto",
                         txOut = TRUE)
# transcripts-gene
tx2gene <- data.frame(
  TXNAME = rownames(txi.kallisto$counts)
  GENEID = sapply(str_split(rownames(txi.kallisto$counts), '\\|'), '[',  2)
)

write_tsv(tx2gene,
          file = "tx2gene.tsv")

# geneレベルの発現量(raw counts)
gene.exp_raw <- summarizeToGene(txi.kallisto,
                            tx2gene)
#head(gene.exp_raw$counts)
rawcount <- data.frame(gene.exp_raw$counts) %>%
  rownames_to_column() %>%
  separate(rowname,  "\\.",  into = c("gene_id",  "postnum")) %>%
  dplyr::select(!postnum)

write_tsv(rawcount,
          file = "kallisto_count.tsv")

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
          file = "kallisto_scaledTPM.tsv")

```
iDEPにはcountデータ`kallisto_count.tsv`を入力する。

edgeRやDEseq2にはscaledTPM`kallisto_scaledTPM.tsv`を使う。

GitHubには、つづけてedgeRで2群間比較を行うスクリプト`kallisto_deg.r`を置きました。

#### 参考：sleuthの利用
本講習では扱いませんが、transcript levelのカウントテーブルだけ共有します。

**参考URL**

https://pachterlab.github.io/sleuth/walkthroughs

https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

- geneレベルにまとめようとするとプロセス開きすぎで落ちる
- transcriptレベルでは解析可能、table作れた

### GitHubにおいたカウントテーブル

https://github.com/motthy/okayama_ws2024

rsemのカウントテーブル　-> data/rsem

kallistoのカウントテーブル
- tximportで作成したgene levelのカウント-> data/kallisto/tximport/
- sleuthで作成したtranscript levelのカウント -> data/kallisto/kallisto_table.csv.gz
