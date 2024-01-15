# 岡山大　バイオインフォマティクスワークショップ
2024-02-06 ~ 2024-02-08

1日目（2024-02-06　10:00-16:00）

- スパコン利用の基礎
- 解析ツールについて
- RNA-seqデータの取得
- RNA-seqデータのクオリティチェック

## 参考図書

レシピ本「改訂版RNA-seqデータ解析 WETラボのための超鉄板レシピ」坊農秀雅編

## 解析ツールについて

遺伝研スパコンを利用する場合、`/usr/bin/`に入っているツールを利用できる。PATHが通っているので、コマンドを入力するだけでよい。ただし、ツールのバージョン更新が頻繁に行われるとは限らない。最新バージョンを使いたい場合は自分でインストールする必要がある。

PATHが通っているか（コマンドの所在）は以下のようにする
```
which fastp
```
遺伝研スパコンの場合、`/usr/bin/fastp`と表示される（2024-01-15）

fastpのバージョン確認は以下のコマンドを入力
```
fastp --version
```
遺伝研スパコンの場合、`fastp 0.20.1`と表示（2024-01-15）

今回利用するコマンドのバージョン（遺伝研スパコン`/usr/bin/`以下）
```
"fasterq-dump" version 2.11.3

/usr/bin/bowtie2-align-s version 2.4.4

RSEM v1.3.1

kallisto version 0.46.2
```

Apptainer(singularity)コンテナを利用する場合、`/usr/local/biotools/ツールの頭文字/`以下から目的のツールのコンテナを探す。たとえば、rsemを探したいときは
```{sh}
ls /usr/local/biotools/r/
```

または正規表現を利用して
```{sh}
ls /usr/local/biotools/r/rsem*
```
とする。

### kallistoのインストール

kallistoは遺伝研スパコンの/usr/bin/kallistoでは、HDF5形式の`~.h5`が出力されないので、公式サイトから最新の実行可能ファイルをダウンロードしてインストールする。

公式サイト　https://github.com/pachterlab/kallisto

配布ページ　https://github.com/pachterlab/kallisto/releases

ツール用のディレクトリを作成(すでにツール用のディレクトリがある場合はそちらを使ってください)
```
mkdir tools
```

`tools`ディレクトリに移動
```{sh}
cd tools
```

リンクをコピーして`wget`コマンドでダウンロード
```{sh}
wget https://github.com/pachterlab/kallisto/releases/download/v0.50.1/kallisto_linux-v0.50.1.tar.gz
```

`tar`コマンドで解凍
```{sh}
tar xzvf kallisto_linux-v0.50.1.tar.gz
```

動くかどうか確認
```{sh}
kallisto/kallisto version
```

### localの環境構築
もし遺伝研スパコンを利用しない場合、condaを使って環境構築をすることをおすすめする。以下の例では、rnaseq_envという仮想環境をつくり、そのなかにsratoolkit、rsem、bowtie2、kallistoをインストールする方法を紹介する。

```
conda create -n rnaseq_env
conda activate rnaseq_env
conda install conda install bioconda::sra-tools
conda install bioconda::rsem
conda install bioconda::kallisto
```
kallistoは**Kallistoのインストール**の項で紹介したように、GitHubのReleaseのページからバイナリをダウンロードしてもよい。

## RNA-seqデータの取得

薬剤処理したHaCaT細胞のデータセット：薬剤処理3サンプルとコントロール（DMSO処理）2サンプル（PRJDB13297）

**論文タイトル** Fucosyltransferase 8 (FUT8) and core fucose expression in oxidative stress response

https://pubmed.ncbi.nlm.nih.gov/36780470/

SRA Run Selectorのリンク

https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRP009171&o=acc_s%3Aa

| DRA | BioSample | 処理 |
| ---- | ---- | ---- |
| DRR357080 | SAMD00451567 | DMSO処理(Control) |
| DRR357081 | SAMD00451568 | DMSO処理(Control) |
| DRR357082 | SAMD00451569 | 5H4PB処理 |
| DRR357083 | SAMD00451570 | 5H4PB処理 |
| DRR357084 | SAMD00451571 | 5H4PB処理 |

作業ディレクトリ `okayama_ws`をつくる
```{sh}
mkdir okayama_ws
```

`okayama_ws`に移動
```{sh}
cd oykayama_ws
```

fastqとlogを入れるディレクトリを作る

```{sh}
mkdir fastq logs
```

SRAToolkitの`fasterq-dump`でfastqを取得、**a~cのいずれか**で実行する。

基本のコマンド; `-O`オプションで出力先を指定、末尾に&をつけてバックグラウンド実行
```{sh}
fastrq-dump DRR357080 -O fastq　&
```

a. forループをつかったスクリプトにして、qsubで投入
```{get_fq.sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")

for FQ in ${FQLIST[@]};do
    fasterq-dump $FQ -O fastq
done
```

b. アレイジョブにする場合
```{get_fq.sh}
#$ -S /bin/bash
#$ -cwd
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")

fasterq-dump ${FQLIST[SGE_TASK_ID - 1]} -O fastq
```

c(参考). singularityコンテナを使い、アレイジョブにする場合
```{get_fq.sh}
#$ -S /bin/bash
#$ -cwd
#$ -t 1-5:1
#$ -o logs
#$ -e logs

SIMS='singularity exec /usr/local/biotools/s/sra-tools:2.11.0--pl5321ha49a11a_3'

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")

$SIMS fasterq-dump ${FQLIST[SGE_TASK_ID - 1]} -O fastq
```

a~c のいずれかを`get_fq.sh`として保存、qsubでジョブ投入
```
$ qsub get_fq.sh
```

（aまたはbの実行中）途中経過としてfastqディレクトリを見る。

```{sh}
ls -l fastq
```
`fasterq.tmp.xxxx`というディレクトリができている。そのなかに、fasterq-dumpの中間ファイルがある。

ジョブが完了すると、fastqディレクトリの中にDRR357080~DRR357084のfastqができる。
```{sh}
$ ls -l fastq
```

出力例
```{sh}
total 125509648
-rw-rw-r-- 1 youraccount yourgroup 13417485894 Dec 27 20:12 DRR357080_1.fastq
-rw-rw-r-- 1 youraccount yourgroup 13417485894 Dec 27 20:12 DRR357080_2.fastq
-rw-rw-r-- 1 youraccount yourgroup 13695918594 Dec 27 20:13 DRR357081_1.fastq
-rw-rw-r-- 1 youraccount yourgroup 13695918594 Dec 27 20:14 DRR357081_2.fastq
-rw-rw-r-- 1 youraccount yourgroup 11045294394 Dec 27 19:55 DRR357082_1.fastq
-rw-rw-r-- 1 youraccount yourgroup 11045294394 Dec 27 19:55 DRR357082_2.fastq
-rw-rw-r-- 1 youraccount yourgroup 11949574134 Dec 27 20:04 DRR357083_1.fastq
-rw-rw-r-- 1 youraccount yourgroup 11949574134 Dec 27 20:04 DRR357083_2.fastq
-rw-rw-r-- 1 youraccount yourgroup 14152601814 Dec 27 20:20 DRR357084_1.fastq
-rw-rw-r-- 1 youraccount yourgroup 14152601814 Dec 27 20:20 DRR357084_2.fastq
```
アレイジョブで実行して2時間くらいかかった。この時点のデータ容量は130GBなので、自分のパソコンで実行している人はストレージの空き容量に注意。

### fastqの圧縮

fastqのままではストレージを圧迫するのでgz圧縮しておく。最近のツールはgz圧縮したfastqを解凍せずに読み込んでくれるものが多い。今回利用するrsemやkallistoもgz圧縮したfastqを入力することができる。

以下の**aまたはb**で実行する。

a. forループをつかったスクリプトにして、qsubで投入
```{compress_fq.sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")

cd fastq

for FQ in ${FQLIST[@]};do
    gzip ${FQLIST[SGE_TASK_ID - 1]}_1.fastq
    gzip ${FQLIST[SGE_TASK_ID - 1]}_2.fastq
done
```

b. アレイジョブで実行
```{compress_fq.sh}
#$ -S /bin/bash
#$ -cwd
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")

cd fastq

gzip ${FQLIST[SGE_TASK_ID - 1]}_1.fastq
gzip ${FQLIST[SGE_TASK_ID - 1]}_2.fastq
```

aまたはbを`compress_fq.sh`として保存、qsubでジョブ投入

```
$ qsub compress_fq.sh
```

（aまたはbの実行が終わったら）圧縮できているかどうか確認する。
```{sh}
$ ls -l fastq
```

出力例
```
total 40242856
-rw-rw-r-- 1 youraccount yourgroup 4136763053 Dec 27 20:12 DRR357080_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4466632721 Dec 27 20:12 DRR357080_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4249129156 Dec 27 20:13 DRR357081_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4511478479 Dec 27 20:14 DRR357081_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3429993615 Dec 27 19:55 DRR357082_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3682755126 Dec 27 19:55 DRR357082_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3689479375 Dec 27 20:04 DRR357083_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3930767256 Dec 27 20:04 DRR357083_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4381069126 Dec 27 20:20 DRR357084_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4730558757 Dec 27 20:20 DRR357084_2.fastq.gz
```
アレイジョブで40分くらいかかった。

### 講習会用の共有データ

fasterq-dumpがうまくいかない　or なかなかジョブが入らない場合、`/home/ayanosatoh/Bioinfo2024/fastq`にgz圧縮済みのfastqファイルを置いてあるので、こちらのファイルを使ってください。
```{sh}
$ ls -l /home/ayanosatoh/Bioinfo2024/fastq
```

## RNA-seqデータのクオリティチェック

結果を出力するディレクトリ`fastp`を`okayama_ws`の下に作成する。
```{sh}
mkdir fastp
```

fastpを使用、**aまたはbのいずれか**で実行する。

a.アレイジョブで実行
```{qc_fq.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

fastp -i fastq/${DRR}_1.fastq.gz \
 -o fastq/${DRR}_1.trim.fq.gz \
 -I fastq/${DRR}_2.fastq.gz \
 -O fastq/${DRR}_2.trim.fq.gz \
 -h fastp/${DRR}_fastp_report.html \
 -j fastp/${DRR}_fastp_report.json \
 -w 4
```

b.singularityコンテナを利用し、アレイジョブで実行
```{qc_fq.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G
#$ -t 1-5:1
#$ -o logs
#$ -e logs

SIMS='singularity exec /usr/local/biotools/f/fastp:0.23.4--hadf994f_2'
FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

$SIMS fastp -i fastq/${DRR}_1.fastq.gz \
 -o fastq/${DRR}_1.trim.fq.gz \
 -I fastq/${DRR}_2.fastq.gz \
 -O fastq/${DRR}_2.trim.fq.gz \
 -h fastp/${DRR}_fastp_report.html \
 -j fastp/${DRR}_fastp_report.json \
 -w 4
```
~.jsonを作っておくと細かい情報が得られる

`qc_fq.sh`として保存、qsubで投入

```{sh}
$ qsub qc_fq.sh
```

ジョブが完了したら結果を確認する。
```{sh}
ls -l fastq
```

出力例
```
total 80079420
-rw-rw-r-- 1 youraccount yourgroup 4136763053 Dec 27 20:12 DRR357080_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4101069241 Dec 27 22:21 DRR357080_1.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 4466632721 Dec 27 20:12 DRR357080_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4412621028 Dec 27 22:21 DRR357080_2.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 4249129156 Dec 27 20:13 DRR357081_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4216015252 Dec 27 22:23 DRR357081_1.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 4511478479 Dec 27 20:14 DRR357081_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4466873955 Dec 27 22:23 DRR357081_2.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 3429993615 Dec 27 19:55 DRR357082_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3398888209 Dec 27 22:26 DRR357082_1.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 3682755126 Dec 27 19:55 DRR357082_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3637405775 Dec 27 22:26 DRR357082_2.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 3689479375 Dec 27 20:04 DRR357083_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3660869285 Dec 27 22:28 DRR357083_1.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 3930767256 Dec 27 20:04 DRR357083_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 3891971561 Dec 27 22:28 DRR357083_2.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 4381069126 Dec 27 20:20 DRR357084_1.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4338148545 Dec 27 22:31 DRR357084_1.trim.fq.gz
-rw-rw-r-- 1 youraccount yourgroup 4730558757 Dec 27 20:20 DRR357084_2.fastq.gz
-rw-rw-r-- 1 youraccount yourgroup 4668715035 Dec 27 22:31 DRR357084_2.trim.fq.gz
```
（実行時間3~4分）

fastpディレクトリを自分のパソコンにダウンロードして、~.htmlをダブルクリックで開くと、それぞれのQC結果をみることができる。主に注意するところは
- 読み込まれたリード数が極端に少ないもの <-fastqが壊れているかも？
- フィルタリング前後でリード数が極端に減っているもの <-fastqのクオリティがよくない

などである。ほかの記載事項は以下を参考にしてほしい。

fastpのGitHubページ　https://github.com/OpenGene/fastp


### mappingとカウント

#### dataset1

##### HISAT2-~~rsem~~ StringTie

rsemはSTARかbowtie2でないとうまく行かないかも?

結果を入れるディレクトリ作成

```{sh}
$ mkdir ref hisat2_rsem kallisto
```


###### HISAT2 index

/home/mikasaka/misc/ayano/okayama_ws/ref

HISAT公式からダウンロードする
http://daehwankimlab.github.io/hisat2/download/

```{sh}
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xzvf grch38_genome.tar.gz
```
6m 31s で終わる。通信の問題。。。



##### rsem(mappingはbowtie2)
HISAT2とrsemってどうなんでしょう、となったのでbowtie2を使うタイプにする
Dockerあった（HumanCellAtlas;たぶんこれはHISAT2を使うタイプ）
https://github.com/HumanCellAtlas/skylab/blob/master/docker/rsem/Dockerfile

###### 1 rsem index作成

GENCODEからダウンロードしてindex作成
https://www.gencodegenes.org/human/

Genome sequence primary assembly (GRCh38)fastaとGTFをダウンロード
```{sh}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

gunzip GRCh38.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

gunzip gencode.v44.annotation.gtf.gz

$ ls -l
total 8712912
-rw-rw-r-- 1 youraccount yourgroup 3151417447 Jul 18 05:32 GRCh38.primary_assembly.genome.fa
-rw-rw-r-- 1 youraccount yourgroup 1568837472 Jul 18 05:31 gencode.v44.annotation.gtf
-rw-rw-r-- 1 youraccount yourgroup 4210306865 Jun 26 2020 grch38_genome.tar.gz
```

index作成スクリプト https://github.com/deweylab/RSEM#built
rsem bowtie2は/usr/bin/にはいっている
Bowtie 2 version 2.4.4
rsem~~はわからない...スパコンの日付2021/9/9が最新、v1.3.3が2020/2/15~~
```
$ rsem-calculate-expression --version
Current version: RSEM v1.3.1
```

```
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G

rsem-prepare-reference --gtf gencode.v44.annotation.gtf \
 --bowtie2 \
 GRCh38.primary_assembly.genome.fa \
 ref/human_gencode
```
これをrsem_index.shとして保存、qsubで実行
```{sh}
$ qsub rsem_index.sh
```
Your job 24934887 ("rsem_index.sh") has been submitted #singularity bowtieのパスが通ってないよね？
いや、色々入ってる...fasterq-dump2.11.3、fastp(2020/4/12)

コマンドで動くよう修正したrsem_index.sh
```{sh}
$ qsub rsem_index.sh
Your job 24934948 ("rsem_index.sh") has been submitted
```
動くわ...

###### 2 mapping
実行スクリプト
```{sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G
#$ -t 1-5:1
#$ -o logs
#$ -e logs


FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

# /usr/bin/に入っているコマンド
rsem-calculate-expression --paired-end -p 4 \
 --bowtie2 \
 --strandedness reverse \
 --append-names \
 --estimate-rspd \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz \
 ref/human_gencode rsem/$DRR
```
実行
```{sh}
$ qsub rsem_mapping.sh
```
Your job-array 24934992.1-5:1 ("rsem_mapping.sh") has been submitted #optionのタイポ
Your job-array 24934999.1-5:1 ("rsem_mapping.sh") has been submitted
18:49-22:22 2時間半

--append-nameを入れるとiDEPでパスウェイ解析ができなくなる（カウントテーブルだけ見るには便利）


###### 3 カウントテーブル作成

スクリプト
```{sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs


#FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
#DRR=${FQLIST[SGE_TASK_ID - 1]}

#/usr/bin/にあるコマンドを使う
# gene level
rsem-generate-data-matrix rsem/DRR*genes.results > rsem/ds1_all.genes.results

#transcripts(isoform) level
rsem-generate-data-matrix rsem/DRR*isoforms.results > rsem/ds1_all.isoforms.results
```
rsem_countmtx.shとして保存、qsubで投入
```{sh}
$ qsub rsem_countmtx.sh
```
Your job 24935741 ("rsem_countmtx.sh") has been submitted

###### 4 iDEP

ENSG00000001084.13_GCLC ".13_GCLC" の部分を削除

`% more ds1_all.genes.results.tsv |sed -r 's/_[A-Z0-9a-z]+\"/\"/g'|sed -r 's/\.[0-9]+\"/\"/g' > ds1_all.genes.results.rename.tsv` # 遺伝子名に-が含まれているものが残ってしまう
```{sh}
% more ds1_all.genes.results.tsv |sed -r 's/_[A-Z0-9a-z\-].+\"/\"/g'|sed -r 's/\.[0-9]+\"/\"/g' > ds1_all.genes.results.rename.tsv
```
これならOK、ただしheaderが消えてしまうので元データからコピーしてやること

実験条件のcsvをつくる

```
 DRR357080 DRR357081 DRR357082 DRR357083 DRR357084
treatment DMSO DMSO 5H4PB 5H4PB 5H4PB
```
ds1_experiment.csvとして保存
ds1_all.genes.results.rename.tsv ds1_experiment.csvをiDEPに入力する

###### 5 iDEP用 マッピングとカウントテーブル作成

`--append-name`オプションを入れない

スクリプトと作業ディレクトリ
```
$ cp rsem_mapping.sh rsem_mapping_idep.sh
$ vi rsem_mapping_idep.sh
$ mkdir rsem_idep
```
実行
```
$ qsub rsem_mapping_idep.sh
```
Your job-array 24965706.1-5:1 ("rsem_mapping_idep.sh") has been submitted

```
mikasaka@at138:okayama_ws $ qstat
job-ID prior name user state submit/start at queue	jclass slots ja-task-ID
------------------------------------------------------------------------------------------------------------------------------------------------
 24965706 0.25164 rsem_mappi mikasaka r 01/09/2024 11:01:00 epyc.q@at143	4 1
 24965706 0.25164 rsem_mappi mikasaka r 01/09/2024 11:01:00 epyc.q@at161	4 2
 24965706 0.25164 rsem_mappi mikasaka r 01/09/2024 11:01:00 epyc.q@at150	4 3
 24965706 0.25164 rsem_mappi mikasaka r 01/09/2024 11:01:00 epyc.q@at162	4 4
 24965706 0.25164 rsem_mappi mikasaka r 01/09/2024 11:01:00 epyc.q@at155	4 5
 24965695 0.25016 QLOGIN mikasaka r 01/09/2024 10:54:44 login.q@at138	1
```
2024-01-09 11:01-14:33 3.5h

カウントテーブル作成 数分で終わるのでshortノードでOK
```
$ qsub -l short rsem_countmtx_idep.sh
Your job 24966457 ("rsem_countmtx_idep.sh") has been submitted
```
遺伝子IDだけになったが、下記のように.1などの枝番が残るので処理が必要
```
 "rsem_idep/DRR357080.genes.results" "rsem_idep/DRR357081.genes.results" "rsem_idep/DRR357082.genes.results" "rsem_idep/DRR357083.genes.results" "rsem_idep/DRR357084.genes.results"
"ENSG00000000003.16" 469.00 494.00 353.00 411.00 489.26
"ENSG00000000005.6" 0.00 0.00 0.00 0.00 0.00
```
処理
```
% more ds1_all.genes.results |sed -r 's/\.[0-9]+\"/\"/g' > ds1.genes.results.rename.tsv
% more ds1.isoforms.results |sed -r 's/\.[0-9]+\"/\"/g' > ds1.isoforms.results.rename.tsv
```
ヘッダを編集(DRRxxxxxxだけにする)


##### kallisto-tximport

###### 1 index作成

GENCODEからtranscriptomeをダウンロードしてindex作成
https://www.gencodegenes.org/human/

```
$ cd ref
$ wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
```
スクリプト
```
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G

kallisto index -i kallisto_index \
 gencode.v44.transcripts.fa.gz \
 -k 31
```
kallisto_index.shとして保存、実行
```
$ qsub kallisto_index.sh
Your job 24966034 ("kallisto_index.sh") has been submitted
```
2024-01-09 11:38-11:46 8min（shortノードでもいける）

###### 2 リードカウント
結果ディレクトリをつくる
```
mkdir kallisto
```
スクリプトでまとめてつくる
```
#!/bin/bash
FQLIST="DRR357080 DRR357081 DRR357082 DRR357083 DRR357084"
for SRR in $FQLIST;do
mkdir kallisto/$SRR
done
```
これをkallistoresult_dir.shとして保存、実行

スクリプト
```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

kallisto quant -i ref/kallisto_index \
 -o kallisto/$DRR \
 -t 4 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz
```
これをkallisto_count.shとして保存

実行;5分/1サンプルなのでshortノードで実行する
```
$ qsub -l short kallisto_count.sh
Your job-array 24966333.1-5:1 ("kallisto_count.sh") has been submitted
```
2024-01-09 13:19-13:24 5min!!はやい

###### 3 カウントテーブル作成
tximport使用 tximport_1.30.0

`kallisto/DRR*` をlocalにダウンロード

https://bi.biopapyrus.jp/rnaseq/analysis/
tximport edgeR（二群間比較）

###### 4 iDEP

dockerでやってみる
docker desktopを立ち上げてから
```
docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest
```
pull完了するまで数分かかる
ブラウザでlocalhost:3838を開く

データ読んだところで落ちる（やばい）
```
docker stop idep
docker rm idep
```
webでやろう

##### kallistoやりなおし

~.h5ができないのでsingularity containerでやりなおし # ~.h5ができるのは-b（ブートストラップオプションのあるとき）
/usr/local/biotools/k/kallisto:0.50.1--hc877fd6_0

WD: kallisto_sims/

```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

#singularity
#SIMS=/usr/local/biotools/k/kallisto:0.50.1--hc877fd6_0
SIMS="singularity exec /usr/local/biotools/k/kallisto:0.50.0--hc877fd6_0"

$SIMS kallisto quant -i ref/kallisto_index \
 -o kallisto_sims/$DRR \
 -t 4 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz

```
kallisto_count_sims.sh として保存

実行
```
$ qsub -l short kallisto_count_sims.sh
```
Your job-array 24966981.1-5:1 ("kallisto_count_sims.sh") has been submitted

コアダンプした
`#$ -l mem_req=8G s_vmem=8G`を追加

Your job-array 24966989.1-5:1 ("kallisto_count_sims.sh") has been submitted
またコアダンプ

8スレッドにする
Your job-array 24966990.1-5:1 ("kallisto_count_sims.sh") has been submitted

qloginしなおし
Your job-array 24966994.1-5:1 ("kallisto_count_sims.sh") has been submitted
やっぱりだめ

```
singularity exec /usr/local/biotools/k/kallisto:0.50.1--hc877fd6_0 kallisto quant -i ref/kallisto_index -o kallisto_sims/DRR357080 --rf-stranded --verbose fastq/DRR357080_1.trim.fq.gz fastq/DRR357080_2.trim.fq.gz
```
だめ

コンテナにするとメモリ足りないのかも
GitHubにはHDF5のサポートは将来無くすと書いてあった
Sluethつかわなくていいか。。。

##### kallisto with bootstrap（-bオプション)
v0.46.2（/usr/bin/にあるもの）
```
$ mkdir kallisto_b
$ cp kallisto_count.sh kallisto_count_b.sh
$ vi kallisto_count_b.sh
$ vi kallistoresult_dir.sh
$ ./kallistoresult_dir.sh
$ ls kallisto_b
DRR357080 DRR357081 DRR357082 DRR357083 DRR357084

$ more kallisto_count_b.sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

#/usr/bin/にあるコマンドを使う

kallisto quant -i ref/kallisto_index \
 -o kallisto_b/$DRR \
 -t 4 -b 100 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz
```
実行
```
$ qsub kallisto_count_b.sh
```
Your job-array 24967170.1-5:1 ("kallisto_count_b.sh") has been submitted

logに
```
Warning: kallisto was not compiled with HDF5 support so no bootstrapping
will be performed. Run quant with --plaintext option or recompile with
HDF5 support to obtain bootstrap estimates.
```
コンパイルの時にHDF5を出力できるようにしないといかんらしい
しかし動いているから結果を見よう
```
$ ls -l kallisto_b/DRR357080
total 34004
-rw-rw-r-- 1 youraccount yourgroup 34815064 Jan 9 22:49 abundance.tsv
-rw-rw-r-- 1 youraccount yourgroup 421 Jan 9 22:49 run_info.json
```
やっぱりできてない

###### バイナリを落としてきた
~/toolsにインストール
```
$ wget kallisto_linux-v0.50.1.tar.gz
$ tar xzvf kallisto_linux-v0.50.1.tar.gz
$ ~/tools/kallisto/kallisto version
kallisto version 0.50.1

# new WD
$ rm -rf kallisto_b/
$ mkdir kallisto_b

$ rm -rf kallisto_b/
$ mkdir kallisto_b
$ vi kallisto_count_b.sh
$ ./kallistoresult_dir.sh
$ ls -l kallisto_b
total 20
drwxrwxr-x 2 youraccount yourgroup 4096 Jan 9 22:59 DRR357080
drwxrwxr-x 2 youraccount yourgroup 4096 Jan 9 22:59 DRR357081
drwxrwxr-x 2 youraccount yourgroup 4096 Jan 9 22:59 DRR357082
drwxrwxr-x 2 youraccount yourgroup 4096 Jan 9 22:59 DRR357083
drwxrwxr-x 2 youraccount yourgroup 4096 Jan 9 22:59 DRR357084

$ more kallisto_count_b.sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

#~/toolsにあるコマンドを使う

~/tools/kallisto/kallisto quant -i ref/kallisto_index \
 -o kallisto_b/$DRR \
 -t 4 -b 100 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz
```
実行
```
$ qsub kallisto_count_b.sh
```
Your job-array 24967181.1-5:1 ("kallisto_count_b.sh") has been submitted

###### index作成 v0.50.0
インデックスをもう一度作らないといかん
```
$ qsub kallisto_index0_50_0.sh
```
Your job 24967185 ("kallisto_index0_50_0.sh") has been submitted
2024-01-09 23:06-23:13

###### count
bootstrapあり（-b 100）
スクリプト
```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-5:1
#$ -o logs
#$ -e logs

FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}

#~/toolsにあるコマンドを使う

~/tools/kallisto/kallisto quant -i ref/kallisto_index0_50_0 \
 -o kallisto_b/$DRR \
 -t 4 -b 100 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz
```
実行
```
$ qsub kallisto_count_b.sh
```
Your job-array 24968556.1-5:1 ("kallisto_count_b.sh") has been submitted
2024-01-10 13:50~14:01 10分くらいで終わる

###### sleuth
チュートリアル
https://pachterlab.github.io/sleuth/walkthroughs

https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

geneレベルにまとめようとするとプロセス開きすぎで落ちる
transcriptレベルでは解析可能、table作れた

###### tximport
カウントテーブル作っておく

#### githubにいれたもの

2024-01-10

https://github.com/motthy/okayama_ws2024

mappingとリードカウント（サンプルごと）
- rsem_idep（gene nameつけない）-> data/rsemにrename
- kallisto_bootstrap（v0.50.0）-> data/kallistoにrename

kallistoのカウントテーブル
- kallisto_edgeR.r（tximport) -> data/kallisto/tximport
- kallisto_sleuth.r (sleuth;transcript levelのみ) -> data/kallisto/kallisto_table.csv.gz



#### library type推定
MGIのキットなのでRFかFRか気になった
RF(dUTP)だと思うけどいちおう.

こちらのページにやり方
https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/

- If sequences of read 1 align to the RNA strand the library is “stranded”.
- If sequences of read 2 align to the RNA strand the library is “reversed stranded”.
- Sometimes sequences of read 1 align to the RNA strand; the other times sequences of read 2 align to the RNA strand. The library is “unstranded”.

```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -o logs
#$ -e logs

#Strandness in RNASeq by Hong Zheng / 2017-08-17
#https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/

DRR="DRR357080"

kallisto quant -i ref/kallisto_index \
 -o libtest/un \
 -t 4 \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz

kallisto quant -i ref/kallisto_index \
 -o libtest/rf \
 -t 4 \
 --rf-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz

kallisto quant -i ref/kallisto_index \
 -o libtest/fr \
 -t 4 \
 --fr-stranded \
 fastq/${DRR}_1.trim.fq.gz \
 fastq/${DRR}_2.trim.fq.gz

paste libtest/fr/abundance.tsv libtest/rf/abundance.tsv libtest/un/abundance.tsv | cut -f1 4 9 14 | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1 sum2 sum3}' > libtest/test.libtypetesting

less libtest/test.libtypetesting | awk '{print $2/$1 $3/$1 $3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' > libtest/test.libtype
```
kallisto_libtypetest.shとして保存、実行
2024-01-09 13:01-13:15

countは5分くらいで終わる（shortノードでもOK）

```
$ more test.libtype
reverse
```

##### 結論
–rf-strandedでOK(=dUTPを使用したキットだった)

#### dataset2

##### rsem

```{rsem_mapping_idep.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G
#$ -t 1-9:1
#$ -o logs
#$ -e logs

SAMPLEID=E326S_$SGE_TASK_ID
WORKDIR=/home/mikasaka/misc/ayano/okayama_ws

# /usr/bin/に入っているコマンド
rsem-calculate-expression --paired-end -p 4 \
 --bowtie2 \
 --strandedness reverse \
 --estimate-rspd \
 fastq/${SAMPLEID}_1.trim.fq.gz \
 fastq/${SAMPLEID}_2.trim.fq.gz \
 ${WORKDIR}/ref/human_gencode rsem_idep/$SAMPLEID
```

実行
```
$ qsub rsem_mapping_idep.sh
```
Your job-array 24983031.1-9:1 ("rsem_mapping_idep.sh") has been submitted
01/11/2024 00:17:08-2:26


```{rsem_countmtx_idep.sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs

#SAMPLEID=E326S_{1..6}

#/usr/bin/にあるコマンドを使う
# gene level
rsem-generate-data-matrix rsem_idep/E326S*genes.results > rsem_idep/ds2.genes.results

#transcripts(isoform) level
rsem-generate-data-matrix rsem_idep/E326S*isoforms.results > rsem_idep/ds2.isoforms.results
```
実行
```
$ qsub rsem_countmtx_idep.sh
```
Your job 24983416 ("rsem_countmtx_idep.sh") has been submitted
01/11/2024 10:31:44-10:32

##### kallisto with bootstrap(-bオプション)

```
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-9:1
#$ -o logs
#$ -e logs

SAMPLEID=E326S_$SGE_TASK_ID
WORKDIR=/home/mikasaka/misc/ayano/okayama_ws

#~/toolsにあるコマンドを使う

~/tools/kallisto/kallisto quant \
 -i ${WORKDIR}/ref/kallisto_index0_50_0 \
 -o kallisto_b/$SAMPLEID \
 -t 4 -b 100 \
 --rf-stranded \
 fastq/${SAMPLEID}_1.trim.fq.gz \
 fastq/${SAMPLEID}_2.trim.fq.gz
```
ディレクトリ作成
```
$ mkdir kallisto_b/E326S_{1..9}
```
実行
```
$ qsub kallisto_count_b.sh
```
Your job-array 24983032.1-9:1 ("kallisto_count_b.sh") has been submitted
01/11/2024 00:17:23-00:40

#### dataset3

##### rsem

```{rsem_mapping_idep.sh}
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -cwd
#$ -l mem_req=4G s_vmem=4G
#$ -t 1-10:1
#$ -o logs
#$ -e logs

SAMPLELIST=(
"A1_S11"
"A2_S12"
"A3_S13"
"A4_S14"
"A5_S15"
"A6_S16"
"A7_S17"
"A8_S18"
"A9_S19"
"A10_S20"
)

SAMPLE=${SAMPLELIST[SGE_TASK_ID - 1]}
WORKDIR=/home/mikasaka/misc/ayano/okayama_ws

# /usr/bin/に入っているコマンド
rsem-calculate-expression --paired-end -p 4 \
 --bowtie2 \
 --strandedness reverse \
 --estimate-rspd \
 fastq/${SAMPLE}_R1.trim.fq.gz \
 fastq/${SAMPLE}_R2.trim.fq.gz \
 ${WORKDIR}/ref/human_gencode rsem_idep/$SAMPLE
```

実行
```
$ qsub rsem_mapping_idep.sh
```
Your job-array 24983462.1-10:1 ("rsem_mapping_idep.sh") has been submitted
01/11/2024 10:50:38-15:19

```{rsem_countmtx_idep.sh}
#$ -S /bin/bash
#$ -cwd
#$ -o logs
#$ -e logs

#/usr/bin/にあるコマンドを使う
# gene level
rsem-generate-data-matrix rsem_idep/A*genes.results > rsem_idep/ds3.genes.results

#transcripts(isoform) level
rsem-generate-data-matrix rsem_idep/A*isoforms.results > rsem_idep/ds3.isoforms.results
```
実行
```
$ qsub rsem_countmtx_idep.sh
```
Your job 24984430 ("rsem_countmtx_idep.sh") has been submitted
01/11/2024 18:11:22-18:11

##### kallisto with bootstrap(-bオプション)

```kallisto_count_b.sh
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -t 1-10:1
#$ -o logs
#$ -e logs

SAMPLELIST=(
"A1_S11"
"A2_S12"
"A3_S13"
"A4_S14"
"A5_S15"
"A6_S16"
"A7_S17"
"A8_S18"
"A9_S19"
"A10_S20"
)

SAMPLE=${SAMPLELIST[SGE_TASK_ID - 1]}
WORKDIR=/home/mikasaka/misc/ayano/okayama_ws

#~/toolsにあるコマンドを使う

~/tools/kallisto/kallisto quant \
 -i ${WORKDIR}/ref/kallisto_index0_50_0 \
 -o kallisto_b/$SAMPLE \
 -t 4 -b 100 \
 --rf-stranded \
 fastq/${SAMPLE}_R1.trim.fq.gz \
 fastq/${SAMPLE}_R2.trim.fq.gz
```
ディレクトリ作成
```
$ cd kallisto_b/
$ mkdir A10_S20 A1_S11 A2_S12 A3_S13 A4_S14 A5_S15 A6_S16 A7_S17 A8_S18 A9_S19
```
実行
```
$ qsub kallisto_count_b.sh
```
Your job-array 24983459.1-10:1 ("kallisto_count_b.sh") has been submitted
01/11/2024 10:49:19-11:17

#### insert size

SRA/DRAにペアエンドリードを登録する時に必要

Picard CollectInsertSizeMetrics
https://gatk.broadinstitute.org/hc/en-us/articles/360037225252-CollectInsertSizeMetrics-Picard-

##### picard（うまくいかない）
```{insert_size.sh}
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -o logs
#$ -e logs

SIMS="apptainer exec /usr/local/biotools/p/picard:3.1.1--hdfd78af_0"

# dataset2

$SIMS picard CollectInsertSizeMetrics \
 I=dataset2/rsem_idep/E326S_1.transcript.bam \
 O=dataset2/E326S_1.insert_size_metrics.txt \
 H=dataset2/E326S_1.insert_size_histogram.pdf \
 M=0.5


# dataset3

$SIMS picard CollectInsertSizeMetrics \
 I=dataset3/rsem_idep/A1_S11.transcript.bam \
 O=dataset3/A1_S11.insert_size_metrics.txt \
 H=dataset3/A1_S11.insert_size_histogram.pdf \
 M=0.5
```
実行
```
$ qsub insert_size.sh
```
Your job 24984694 ("insert_size.sh") has been submitted
qw 01/11/2024 20:55:03

コマンド新しくなったって
```
CollectInsertSizeMetrics -I dataset2/rsem_idep/E326S_1.transcript.bam -O dataset2/E326
S_1.insert_size_metrics.txt -H dataset2/E326S_1.insert_size_histogram.pdf -M 0.5
```
pdfできてない、txtはある
sample MEDIAN_INSERT_SIZE MODE_INSERT_SIZE MEDIAN_ABSOLUTE_DEVIATION MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE
E326S_1 251 232 16.5 232 278 253
A1_S11 167 167 0 167 167 167

ゲノムリードちょっと短いのでほかのもみてみる
```
$ qsub insert_size_ds3.sh
```
Your job 24984747 ("insert_size_ds3.sh") has been submitted
qw 01/11/2024 21:54:54
おちた

Your job 24984752 ("insert_size_ds3.sh") has been submitted
qw 01/11/2024 22:02:23
20240111-22:23:14-22:24:41
まただめ？

Your job 24984895 ("insert_size_ds3.sh") has been submitted
だめー

/usr/local/biotools/p/picard:2.26.4--hdfd78af_0にする
```
$ qsub insert_size.sh
```
Your job 24984976 ("insert_size.sh") has been submitted
01/11/2024 23:31:05

もうひとつ memory増したもの
```
$ qsub insert_size_ds3.sh
```
Your job 24984987 ("insert_size_ds3.sh") has been submitted
01/11/2024 23:38:17

どちらもRがはいっていないとエラーになる（つかえねえ）
あとたぶんLANG=Cを書かないとエラーになる？

##### goleft
2024-01-12
bamの分析に使うバイオインフォマティクスのツールキット goleft
https://kazumaxneo.hatenablog.com/entry/2018/02/13/214230

コンテナ
/usr/local/biotools/g/goleft:0.2.4--h9ee0642_1

```
$ singularity exec /usr/local/biotools/g/goleft:0.2.4--h9ee0642_1 goleft -h
goleft Version: 0.2.4

covstats : coverage stats across bams by sampling
depth : parallelize calls to samtools in user-defined windows
depthwed : matricize output from depth to n-sites * n-samples
indexcov : quick coverage estimate using only the bam index
indexsplit : create regions of even coverage across bams/crams
samplename : report samplename(s) from a bam's SM tag
```
実行
```
$ singularity exec /usr/local/biotools/g/goleft:0.2.4--h9ee0642_1 goleft covstats rsem_idep/E326S_1.transcript.bam
coverage insert_mean insert_sd insert_5th insert_95th template_mean template_sd pct_unmapped pct_bad_readspct_duplicate pct_proper_pair read_length bam sample
panic: open rsem_idep/E326S_1.transcript.bai: no such file or directory

goroutine 1 [running]:
github.com/brentp/goleft/covstats.pcheck(...)
 /home/brentp/go/go/src/github.com/brentp/goleft/covstats/covstats.go:30
github.com/brentp/goleft/covstats.Main()
 /home/brentp/go/go/src/github.com/brentp/goleft/covstats/covstats.go:246 +0xed9
main.main()
 /home/brentp/go/go/src/github.com/brentp/goleft/cmd/goleft/goleft.go:68 +0x170
```
indexがないとだめらしい
```
# indexing
$ samtools index E326S_1.transcript.bam
samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libtinfow.so.6: no version information available (required by samtools)
samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)
samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)
[E::hts_idx_push] Unsorted positions on sequence #54395: 8806 followed by 8678
[E::sam_index] Read 'A00917:1083:HTTKMDSX3:3:1101:1741:1000' with ref_name='ENST00000687169.1' ref_length=9959 flags=163 pos=8678 cannot be indexed
samtools index: failed to create index for "E326S_1.transcript.bam"

```
indexできない...

##### kazumaxさんのスクリプト(成功)
2024-01-12

```{calc_insert_size.sh}
#!/bin/bash

# this program needs samtools seqkit minimap2 and goleft
# version1 2018/06/11
# modified 2024/01/12 using Appatainer container

CMDNAME=`basename $0`
if [ $# -ne 3 ]; then
 echo "Usage: $CMDNAME ref.fa pair1.fq pair2.fq" 1>&2
 exit 1
fi

# sampling reads
seqkit sample -n 100000 $2 > sampling1.fq 2>/dev/null
seqkit sample -n 100000 $3 > sampling2.fq 2>/dev/null

#mapping
minimap2 -t 8 -ax sr $1 sampling1.fq sampling2.fq 2>/dev/null | samtools sort -@ 8 -O BAM -o sampling.bam - 2>/dev/null
#index
samtools index sampling.bam

# output
echo
echo ********results********
singularity exec /usr/local/biotools/g/goleft:0.2.4--h9ee0642_1 goleft covstats sampling.bam |column -t
echo

# remove temporary files
rm sampling1.fq sampling2.fq sampling.bam*
exit
```
実行権限つける
```
$ chmod +x calc_insert_size.sh
```
実行
```
$ nohup ./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset2/fastq/E326S_1_1.trim.fq.gz dataset2/fastq/E326S_1_2.trim.fq.gz > ds2_E326S_1.insertsize.txt &

$ more ds2_E326S_1.insertsize.txt

samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libtinfow.so.6: no version information available (required by samtools)
samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)
samtools: /lustre7/home/mikasaka/tools/miniconda3/bin/../lib/libncursesw.so.6: no version information available (required by samtools)

********results********
chromosomes: GL000008.2 GL000009.2 GL000194.1 GL000208.1 GL000213.1 GL000214.1 GL000216.2 GL000221.1 GL000224.1 GL000225.1 GL000226.1 KI270311
.1 KI270317.1 KI270322.1 KI270435.1 KI270438.1 KI270442.1 KI270512.1 KI270519.1 KI270538.1 KI270579.1 KI270589.1 KI270706.1 KI270707.1 KI27070
8.1 KI270709.1 KI270710.1 KI270713.1 KI270714.1 KI270715.1 KI270716.1 KI270717.1 KI270718.1 KI270720.1 KI270722.1 KI270723.1 KI270724.1 KI2707
25.1 KI270726.1 KI270729.1 KI270730.1 KI270731.1 KI270732.1 KI270734.1 KI270735.1 KI270736.1 KI270737.1 KI270738.1 KI270739.1 KI270740.1 KI270
741.1 KI270743.1 KI270744.1 KI270745.1 KI270746.1 KI270747.1 KI270748.1 KI270749.1 KI270750.1 KI270751.1 KI270752.1 KI270753.1 KI270754.1 KI27
0755.1 KI270756.1 KI270757.1 not found in sampling.bam
coverage insert_mean insert_sd insert_5th insert_95th template_mean template_sd pct_unmapped pct_bad_reads pct_duplicate pct_proper_
pair read_length bam sample
0.01 -53.74 53.33 -131 37 243.43 53.14 0.25 0.0 0.0 71.4
 150 sampling.bam <no-read-groups>

$ nohup ./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset3/fastq/A1_S11_R1.trim.fq.gz dataset3/fastq/A1_S11_R2.trim.fq.gz > ds3_A1_S11.insertsize.txt &

coverage insert_mean insert_sd insert_5th insert_95th template_mean template_sd pct_unmapped pct_bad_reads pct_duplicate pct_proper_
pair read_length bam sample
0.01 -88.76 58.12 -146 29 209.10 59.35 0.49 0.0 0.0 78.6
 151 sampling.bam <no-read-groups>
```
ゲノムリードさんのインサート、ちょっと短い

batch処理にする
```{batch_calc_insertsize.sh}
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l mem_req=4G s_vmem=4G
#$ -o logs
#$ -e logs

for DS2 in $(seq 2 9);do
./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset2/fastq/E326S_${DS2}_1.trim.fq.gz dataset2/fastq/E326S_${DS2}_2.trim.fq.gz
> insertsize/ds2_E326S_${DS2}.insertsize.txt
done

DS3_LIST="A2_S12 A3_S13 A4_S14 A5_S15 A6_S16 A7_S17 A8_S18 A9_S19 A10_S20"

for DS3 in ${DS3_LIST};do

./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset3/fastq/${DS3}_R1.trim.fq.gz dataset3/fastq/${DS3}_R2.trim.fq.gz > insertsize/
ds3_${DS3}.insertsize.txt

done
```

実行
```
$ qsub batch_calc_insertsize.sh
```
Your job 24985541 ("batch_calc_insertsize.sh") has been submitted

うまくいかないのでさらに修正
- seqkit sample -p にする（-nは一度に全部読み込むので大きなfastqでは-pを使うように書いてあった）

```{calc_insertsize.sh}
#!/bin/bash

# this program needs samtools seqkit minimap2 and goleft
# version1 2018/06/11
# modified 2024/01/12 using Appatainer container
# modified 2024/01/12 for large fastq

CMDNAME=`basename $0`
if [ $# -ne 3 ]; then
 echo "Usage: $CMDNAME ref.fa pair1.fq pair2.fq" 1>&2
 exit 1
fi

# sampling reads
seqkit sample -p 0.1 $2 > sampling1.fq 2>/dev/null
seqkit sample -p 0.1 $3 > sampling2.fq 2>/dev/null

#mapping
minimap2 -t 8 -ax sr $1 sampling1.fq sampling2.fq 2>/dev/null | samtools sort -@ 8 -O BAM -o sampling.bam
 - 2>/dev/null
#index
samtools index sampling.bam

# output
echo
echo ********results********
singularity exec /usr/local/biotools/g/goleft:0.2.4--h9ee0642_1 goleft covstats sampling.bam |column -t
echo

# remove temporary files
rm sampling1.fq sampling2.fq sampling.bam*
exit
```
これをバッチ処理

```{batch_calc_insertsize.sh}
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe def_slot 8
#$ -l mem_req=4G s_vmem=4G
#$ -o logs
#$ -e logs

for DS2 in $(seq 2 9);do
./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset2/fastq/E326S_${DS2}_1.trim.fq.gz data
set2/fastq/E326S_${DS2}_2.trim.fq.gz > insertsize/ds2_E326S_${DS2}.insertsize.txt
done

DS3_LIST="A2_S12 A3_S13 A4_S14 A5_S15 A6_S16 A7_S17 A8_S18 A9_S19 A10_S20"

for DS3 in $DS3_LIST;do
./calc_insert_size.sh ref/GRCh38.primary_assembly.genome.fa dataset3/fastq/${DS3}_R1.trim.fq.gz dataset3/
fastq/${DS3}_R2.trim.fq.gz > insertsize/ds3_${DS3}.insertsize.txt
done
```
うまくいったっぽい



##### 参考
What is the difference between a Read and a Fragment in RNA-seq?
https://www.biostars.org/p/106291/

outer insert size = リード1 + リード2 + 間に挟まれた部分
inner insert size = 間に挟まれた部分