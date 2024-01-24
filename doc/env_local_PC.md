# local PCの環境設定

condaの仮想環境を使って、RNA-seqの解析ツールをインストールします。

## 1. miniforgeのインストール

miniforge
https://github.com/conda-forge/miniforge

上記リンク先から、自分のOSに対応したインストーラをダウンロードし、説明に沿ってインストールしてください。

M1/M2 Macの場合は **arm64 (Apple Silicon)** を選んでください。

**MacでHomebrewを使っている場合はbrewでインストールしてください**
```{sh}
brew install miniforge
```

## 2. 仮想環境の作成

rnaseq_envという仮想環境をつくります。ターミナルapp（Windowsの場合はTeraPadなどの端末アプリ）から、

```{sh}
conda create -n rnaseq_env
```

仮想環境の有効化
```{sh}
conda activate rnaseq_env
```

SRA-toolkitのインストール
```{sh}
conda install -y bioconda::sra-tools
```

fastpのインストール
```{sh}
conda install -y bioconda::fastp
```

bowtie2のインストール
```{sh}
conda install -y bioconda::bowtie2
```

rsemのインストール
```{sh}
conda install -y bioconda::rsem
```

kallistoのインストール
```{sh}
conda install -y bioconda::kallisto
```

仮想環境を閉じる
```{sh}
conda deactivate
```

## 3. 解析を始める時

仮想環境をふたたび有効化
```{sh}
conda activate rnaseq_env
```
