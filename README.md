<h3 align=center>HmnTrimmer</h3>

[![Github Release](https://img.shields.io/github/v/release/guillaume-gricourt/HmnTrimmer)](version) [![Conda Release](https://img.shields.io/conda/vn/bioconda/hmntrimmer.svg)](https://anaconda.org/bioconda/hmntrimmer)  
[![DOI](https://zenodo.org/badge/270999070.svg)](https://zenodo.org/badge/latestdoi/270999070)  

# Introduction
A trimmer of reads produced by NGS dedicated for common applications like genomic, transcriptomic, targeted metagenomic and shotgun metagenomic.

# Install

## Conda (recommanded)

```shell
conda install -c bioconda hmntrimmer
```

## Docker

```shell
# From docker hub
docker pull hmntrimmer:<VERSION>
# From github
docker pull ghcr.io/guillaume-gricourt/hmntrimmer:<VERSION>
```

## Manual

### Prerequisites
Use software with debian systems :  
   * `yasm`  
   * `build-essential`  
   * `zlib1g-dev`  
GCC used for compilation must be > 4 and < 8.  

Test software :  
   * `python3`  

Create statistic report :  
With conda :  
   * `python3 django matplotlib seaborn packaging`  
With ubuntu/debian using pip :  
   * `python3-pip`
   * `django matplotlib seaborn packaging`  
   
### Compile

Install first `igzip`  
`hmndir=./HmnTrimmer`  
`cd ./lib/igzip-042/igzip && make slib0c`  
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD`

Then  
`make`

## Test

`make test`

## Getting started
Software is available by :  
`HmnTrimmer [OPTIONS] [TRIMMERS]`

Minimal example :  
```shell
./HmnTrimmer \
  --input-fastq-forward INPUT_FILE \
  --output-fastq-forward OUTPUT_FILE \
  --length-min 50
```

# Commands

## Input/Output

Files are indicated with these differents commands :  
```shell
  --input-fastq-forward INPUT_FILE
  --input-fastq-reverse INPUT_FILE
  --input-fastq-interleaved INPUT_FILE

  --output-fastq-forward OUTPUT_FILE
  --output-fastq-reverse OUTPUT_FILE
  --output-fastq-interleaved OUTPUT_FILE
```

Discarded sequences are optionnaly output with this command.
If sequencing is paired, file produced is interleaved.

```shell
  --output-fastq-discard OUTPUT_FILE
```

## Trimmers
Several categories : quality, length and information.  
Firstly trimmers based on information are applied, then based on quality finaly based on length.

### Quality Tail 
Based on a successive number of bases from end of read which are below a cut off.  
Two parameters : quality, optionaly the number of bases below the quality firstly indicated (default 1 base) and the length percent cut off request to keep read if it was truncated (default not removed).  
Format :  `<int>:<int>:<int>`  
```shell
  --quality-tail STRING
```

### Quality Sliding Window 
Based on a sliding window of bases from end of read which are below a minimal mean.  
Two parameters : mean quality and size of window.  
Format : `<int>:<int>`  
```shell
  --quality-sliding-window STRING
```

### Length Min
Minimal length to keep a read.
```shell
  --length-min INTEGER
```

### Information Dust
Based on Dust score.
```shell
  --information-dust INTEGER
```

## Performance/Other Options:

### Report
Optionaly save a report, with differents statistics. Format Json.
```shell
  --output-report OUTPUT_FILE
```
### Threads
Specify number of threads to use.
```shell
  --threads 1..8
```
### Reads batch
Reads are read in batch. Defined size of batch.
```shell
  --reads-batch 100..50000000
```

### Verbose
Log level to use.
```shell
  --verbose 1..6 (error..trace)
```

## Statistic report
To create html report :
```shell
# Clone the repository
git clone git@github.com:guillaume-gricourt/HmnTrimmer.git
# Run
./HmnTrimmer/script/RenderingReportFile.py \
  --template-file ./HmnTrimmer/script/template.html \
  --input-file JSON_FILE \
  --output-file HTML_FILE
```

## Docker

### Trimming
```shell
docker run \
    -it \
    --rm \
    -v $PWD:$PWD \
    hmntrimmer:<VERSION> \
    --input-fastq-forward $PWD/test/GoldInput/BIG.R1.fastq \
    --output-fastq-forward $PWD/test/DockerTest.R1.fastq.gz \
    --output-report $PWD/test/DockerTest.json \
    --length-min 50
```

### Statistic report
```shell
docker run \
    -it \
    --rm \
    -v $PWD:$PWD \
    --entrypoint /opt/HmnTrimmer/script/RenderingReportFile.py \
    hmntrimmer:<VERSION> \
    --input-file $PWD/test/DockerTest.json \
    --output-file $PWD/test/DockerTest.html \
    --template-file /opt/HmnTrimmer/script/template.html
```

# Built with these main libraries

* [SeqAn](https://seqan.readthedocs.io/en/master/index.html) - Essential library to work with HTS files, algorithms
* [rapidjson](https://rapidjson.org/) - Read/Write Json files efficiently
* [spdlog](https://github.com/gabime/spdlog) - Nice log manager
* [igzip](https://software.intel.com/content/www/us/en/develop/articles/igzip-a-high-performance-deflate-compressor-with-optimizations-for-genomic-data.html) - Very fast deflate algorithm

# Versioning

[SemVer](http://semver.org/) is used for versioning.

# Authors

* **Guillaume Gricourt**

# License

See the [LICENSE.md](LICENSE) file for details
