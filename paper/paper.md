---
title: 'HmnTrimmer: a fast-flexible trimmer used for several applications of next-generation sequencing'
tags:
  - Python
authors:
  - name: Guillaume Gricourt
    orcid: 0000-0003-0143-5535
    affiliation: 1
    corresponding: true
  - name: Christophe Rodriguez
    affiliation: "2, 3, 5"
  - name: Jean-Michel Pawlotsky
    affiliation: "2, 3, 5"
  - name: Abdel Aissat
    affiliation: "2, 4"
affiliations:
 - name: Paris-Saclay University, INRAE, AgroParisTech, Micalis Institute, Jouy-en-Josas, 78352, France
   index: 1
 - name: Genomics platform, Henri-Mondor Hospital, APHP, Creteil, 94010, France
   index: 2
 - name: Department of Prevention, Diagnosis and Treatment of Infections, Henri-Mondor Hospital, APHP, Creteil, 94010, France
   index: 3
 - name: Department of Genetics, Henri-Mondor Hospital, APHP, Creteil, 94010, France
   index: 4
 - name: Mondor Institute for Biomedical Research, Inserm U955, Creteil, 94010, France
   index: 5
date: 6 July 2022
bibliography: paper.bib
---


# Summary
Next-Generation Sequencing technologies have the capacity to generate millions of short DNA sequence reads by massive parallel sequencing.
While new technologies are continuously being developed, two dominate the market: one by Illumina, the other one by Ion Torrent, using dedicated devices.
Despite differences in the length of their reads and single or paired-end reading modes, both manufacturers produce sequences with a quality that decreases at the tail of the reads.
Moreover, sequence errors may arise as a result of sample/library preparation or amplification, with an error rate increasing with extreme nucleic acid contents such as GC-rich regions or homopolymers.
Even the spatial position along the flow-cell can produce single position errors or artefacts such as indels [@Laehnemann2016].
Additionally, sequencers themselves may add errors due to misinterpreted signals during base calling.
Contamination of adapters can also occur but, in most of protocols, they are removed by the manufacturer's software.
NGS applications, such as genomics, transcriptomics, shotgun or targeted metagenomics, require read trimming tools to discard low-quality sequences.

# Statement of need

Trimming tools are usually dedicated to only one application, while they partially cover others, and to reads generated exclusively by Illumina or Ion Torrent instruments.  
Popular softwares like fastp [@Chen2018], one of the latest and fastest developed tool, and Trimmomatic [@Bolger2014] use mainly algorithms based on the quality scores and the length of the reads.
Others like Reaper or Kraken tools [@Davis2013] trim sequences using DUST algorithm.
Rare softwares, like Prinseq [@Schmieder2011] apply all algorithms to trim the reads.
Surprinsigly, none of them optimize the step to deflate sequences and deal with single or paired-end FASTQ files with all these algorithms.


HmnTrimmer is a tool built to trim data coming from genomics, shotgun metagenomics and targeted metagenomics.
Each application requires one or several algorithms to discard or crop the read.
Two algorithms use quality scores and work in $O(n)$ complexity where $n$ is the length of the read:
The first one sweeps a score window if a minimum quality is not reached.
If a score window is found with a mean quality below the requested one, the tail of the read is trimmed in addition to the adjacent bases with a score below the requested mean [@Bolger2014].
The second one trims the tail of each read sequence if a group of bases below a specific score of quality is found.
This way of trimming is recommended to clean Operational Taxonomic Unit in targeted metagenomics [@Bokulich2013].  
Other algorithms trim sequences based on sequence information.
The DUST algorithm, implemented by BLAST to mask low-complexity regions, is suited to many applications such as shotgun metagenomics.
DUST works in $0(n+Kmer)$ complexity where $Kmer$ is the number of k-mers found in the read.  
The last algorithm discards sequences based on the number of uncalled bases.
It works in $0(n)$ complexity.
All implemented algorithms work in $0(1)$ in space.  


HmnTrimmer is a user-friendly, efficient trimming tool writtent in C++ templates using Seqan extended with igzip [@Guilford2014], spdlog and rapidJson libraries.
It can be used with Docker, images are available at GitHub and DockerHub.
HmnTrimmer accepts standard file formats, either compressed or not, as input or output: FASTQ or interleaved FASTQ, compressed or not.
A JSON report can be produced, given an evidence of multiple elements about the software, the command line, the running time or the number of sequences.
This file can be integrated into dashboard reports or a convenient script is provided to visualize these statistics.
HmnTrimmer uses igzip library to deflate data which is optimised for genomic data.  
HmnTrimmer is easy to use and may become a helpful new companion for analysis pipelines required to preprocess FASTQ files coming from different sequencing platform.  

# Acknowledgements

We thank all members of the genomic platform at hospital Henri Mondor, France.
We are grateful to the Genotoul bioinformatics platform Toulouse Occitanie for providing help and/or computing and/or storage resources.

# References
