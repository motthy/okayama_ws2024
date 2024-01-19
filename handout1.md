# 岡山大　バイオインフォマティクスワークショップ
2024-02-06 ~ 2024-02-08

1日目（2024-02-06　10:00-16:00）

- スパコン利用の基礎
- 解析ツールについて
- RNA-seqデータの取得
- RNA-seqデータのクオリティチェック

## 参考図書

「改訂版RNA-seqデータ解析 WETラボのための超鉄板レシピ」坊農秀雅編、羊土社

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

- fasterq-dump version 2.11.3

- /usr/bin/bowtie2-align-s version 2.4.4

- RSEM v1.3.1

- kallisto version 0.46.2


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
kallistoは**kallistoのインストール**の項で紹介したように、GitHubのReleaseのページからバイナリをダウンロードしてもよい。

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

homeに移動
```{sh}
cd ~
```

作業ディレクトリ `okayama_ws`をつくる
```{sh}
mkdir okayama_ws
```

`okayama_ws`に移動
```{sh}
cd okayama_ws
```

### お試し

fastqの取得からQCまでの実行例（1サンプル分のみ）を以下に示す。

（スパコンを利用しているなど、ストレージに余裕がある場合は、次の「**fastq取得**」から始めてください）

お試し用ディレクトリの作成
```{sh}
mkdir test
```

1. fastq取得

fasterq-dumpで
```{sh}
fasterq-dump DRR357080 -O test &
```
`-O`オプションで出力先を指定、末尾に&をつけてバックグラウンド実行

2. fastq圧縮

testディレクトリへ移動
```{sh}
cd test
```

DRR357080_1.fastqを`gzip`で圧縮（&をつけてバックグラウンド実行）
```{sh}
gzip DRR357080_1.fastq &
```

同じようにDRR357080_2.fastqも圧縮（&をつけてバックグラウンド実行）
```
gzip DRR357080_2.fastq &
```

完了後、testディレクトリに以下のファイルがあることを確認
- DRR357080_1.fastq.gz
- DRR357080_2.fastq.gz


3. QC

fastpでクオリティチェックとトリミング（&をつけてバックグラウンド実行）
```{sh}
fastp -i DRR357080_1.fastq.gz -o DRR357080_1.trim.fq.gz -I DRR357080_2.fastq.gz -O DRR357080_2.trim.fq.gz -h DRR357080_fastp_report.html &
```
- -i、-I  入力fastq
- -o、-O　出力fastq(トリミング後)
- -h QC結果ファイル（html）

完了後、testディレクトリに以下のファイルがあることを確認
- DRR357080_1.fastq.gz
- DRR357080_2.fastq.gz
- DRR357080_1.trim.fq.gz
- DRR357080_2.trim.fq.gz
- DRR357080_fastp_report.html

### fastq取得
スパコンを利用している場合、バッチジョブまたはアレイジョブでデータセットをまとめて取得することができる。
(localではforループを使ったバッチジョブなら可能)

fastqとlogを入れるディレクトリを作る

```{sh}
mkdir fastq logs
```

SRAToolkitの`fasterq-dump`でfastqを取得、**a~cのいずれか**で実行する。

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

a~c のいずれかを`get_fq.sh`として保存、`qsub`でジョブ投入
```
qsub get_fq.sh
```

（aまたはbの実行中）途中経過としてカレントディレクトリを見る。

```{sh}
ls -l
```
`fasterq.tmp.xxxx`というディレクトリができている。そのなかに、`fasterq-dump`の中間ファイルがある。

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

a. forループをつかったスクリプトにして、`qsub`で投入
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

aまたはbを`compress_fq.sh`として保存、`qsub`でジョブ投入

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

`fasterq-dump`がうまくいかない　or なかなかジョブが入らない場合、`/home/ayanosatoh/Bioinfo2024/fastq`にgz圧縮済みのfastqファイルを置いてあるので、こちらのファイルを使ってください。
```{sh}
$ ls -l /home/ayanosatoh/Bioinfo2024/fastq
```

## RNA-seqデータのクオリティチェック

結果を出力するディレクトリ`fastp`を`okayama_ws`の下に作成する。
```{sh}
mkdir fastp
```

`fastpを使用`、**aまたはbのいずれか**で実行する。

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
`~.json`を作っておくと、より詳細な情報が得られる

`qc_fq.sh`として保存、`qsub`で投入

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

`fastp`ディレクトリを自分のパソコンにダウンロードして、`~.html`をダブルクリックで開くと、それぞれのQC結果をみることができる。主に注意するところは
- 読み込まれたリード数が極端に少ないもの <-fastqが壊れているかも？
- フィルタリング前後でリード数が極端に減っているもの <-fastqのクオリティがよくない

などである。ほかの記載事項は以下を参考にしてほしい。

fastpのGitHubページ　https://github.com/OpenGene/fastp

### 参考（講習会用の共有データの利用）

`/home/ayanosatoh/Bioinfo2024/fastq`にあるgz圧縮済みのfastqファイルを使用する場合のスクリプト

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
FQDIR=/home/ayanosatoh/Bioinfo2024/fastq

fastp -i ${FQDIR}/${DRR}_1.fastq.gz \
      -o fastq/${DRR}_1.trim.fq.gz \
      -I ${FQDIR}/${DRR}_2.fastq.gz \
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

SIMS="singularity exec -B /home/ayanosatoh/Bioinfo2024/fastq /usr/local/biotools/f/fastp:0.23.4--hadf994f_2"
FQLIST=("DRR357080" "DRR357081" "DRR357082" "DRR357083" "DRR357084")
DRR=${FQLIST[SGE_TASK_ID - 1]}
FQDIR=/home/ayanosatoh/Bioinfo2024/fastq

$SIMS fastp -i ${FQDIR}/${DRR}_1.fastq.gz \
            -o fastq/${DRR}_1.trim.fq.gz \
            -I ${FQDIR}/${DRR}_2.fastq.gz \
            -O fastq/${DRR}_2.trim.fq.gz \
            -h fastp/${DRR}_fastp_report.html \
            -j fastp/${DRR}_fastp_report.json \
            -w 4
```

## 参考（遺伝研スパコンのshared_data）

DDBJ SearchでDRA(SRA) accessionを検索、Downloadの項の「fastq」にリンクがあれば、スパコンの`/usr/local/shared_data/dra/fastq/`以下にfastq（bzip2圧縮）がある。

例: DRR357080 https://ddbj.nig.ac.jp/resource/sra-run/DRR357080

sra-submission(DRA013778)とsra-experiment(DRX342991)からディレクトリを探す

```{sh}
ls -l /usr/local/shared_data/dra/fastq/DRA013/DRA013778
total 60
-rw-r--r-- 1 tracesys tracesys 12989 Nov  2  2022 DRA013778.experiment.xml
-rw-r--r-- 1 tracesys tracesys  1257 Nov  2  2022 DRA013778.run.xml
-rw-r--r-- 1 tracesys tracesys  9734 Nov  2  2022 DRA013778.sample.xml
-rw-r--r-- 1 tracesys tracesys  1347 Nov  2  2022 DRA013778.study.xml
-rw-r--r-- 1 tracesys tracesys   243 Nov  2  2022 DRA013778.submission.xml
drwxr-xr-x 2 tracesys tracesys  4096 Nov  2  2022 DRX342991
drwxr-xr-x 2 tracesys tracesys  4096 Nov  2  2022 DRX342992
drwxr-xr-x 2 tracesys tracesys  4096 Nov  2  2022 DRX342993
drwxr-xr-x 2 tracesys tracesys  4096 Nov  2  2022 DRX342994
drwxr-xr-x 2 tracesys tracesys  4096 Nov  2  2022 DRX342995
```

fastqの場所
```{sh}
ls -l /usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342991
total 7169212
-rw-r--r-- 1 tracesys tracesys 3537298471 Apr  1  2022 DRR357080_1.fastq.bz2
-rw-r--r-- 1 tracesys tracesys 3803962633 Apr  1  2022 DRR357080_2.fastq.bz2
```
bz2をそのまま使えないツールの場合、自分の`$HOME`以下にコピーし、解凍して使う（さらにgz圧縮する場合もある）

今回使用したデータセットPRJDB13297のfastqのPATH
```{sh}
/usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342991
/usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342992
/usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342993
/usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342994
/usr/local/shared_data/dra/fastq/DRA013/DRA013778/DRX342995
```