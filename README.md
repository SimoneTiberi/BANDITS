# BANDITS: Bayesian ANalysis of DIfferenTial Splicing
`BANDITS` is a Bayesian hierarchical model for detecting differential splicing of genes and transcripts,
via differential transcript usage (DTU), 
between two or more conditions.
The method uses a Bayesian hierarchical framework, which allows for sample specific proportions
in a Dirichlet-Multinomial model, and samples the allocation of fragments to the transcripts.
Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a DTU test is performed 
via a multivariate Wald test on the posterior densities for the average relative abundance of transcripts.

## Bioconductor installation 
`BANDITS` is available on [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/BANDITS.html) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("BANDITS")
```

## Devel installation from Github
To install the latest development version of the package from github, use `devtools` (available [here](https://github.com/hadley/devtools)):
``` r
devtools::install_github("SimoneTiberi/BANDITS")
```

To install the package jointly with its vignette use \code{build_vignettes = TRUE}:
``` r
devtools::install_github("SimoneTiberi/BANDITS", build_vignettes = TRUE)
```

## Vignette
The vignette illustrating how to use the package can be accessed on the 
[Bioconductor website](https://www.bioconductor.org/packages/release/bioc/vignettes/BANDITS/inst/doc/BANDITS.pdf)
or from R via:
``` r
vignette("BANDITS")
```
or
``` r
browseVignettes("BANDITS")
```

## Alignment
The package inputs the equivalence classes and respective counts, representing what transcripts each read is compatible with.
These can be obtained by aligning reads either directly to a reference transcriptome with pseudo-alignmers, such as \software{salmon} or \software{kallisto}, or to a reference genome with splice-aware genome alignment algorithms, such as \software{STAR} or \software{TopHat2}, and checking the transcripts compatible with each genome alignment.

NOTE: currently \software{BANDITS} only inputs equivalence classes computed from \software{salmon}.
We are extending our current framework to allows reads to be aligned with \software{kallisto}.

Importantly, when using \software{salmon}, use the option `--dumpEq` to obtain the equivalence classes, and when using \software{STAR}, use the option `--quantMode TranscriptomeSAM` to obtain alignments translated into transcript coordinates.

Below we show a pipeline for aligning reads in both ways.

## Download the example input data
To obtain the example raw data, download or clone the ARMOR github repository:
``` bash
git clone https://github.com/csoneson/ARMOR.git

# set a base_dir variable to the downloaded repo
base_dir="~/ARMOR-master"

# input reads:
fastq_files=$base_dir/example_data/FASTQ
```

The example data consits of four paired-end RNA-seq reads of 63 base pairs.

## Align reads to the genome with STAR, and compute the equivalence classes with salmon
Create a variable for the fasta format reference genome (DNA) and gtf file:
``` bash
fasta=$base_dir/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.dna.chromosome.1.1.10M.fa
gtf=$base_dir/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.93.1.1.10M.gtf
```

Make the directory for the STAR output and the genome index
``` bash
# create a directory for STAR
mkdir $base_dir/STAR

# create a directory for the genome index
mkdir $base_dir/STAR/genome_index
GDIR=$base_dir/STAR/genome_index

# create a directory for the output of the alignment
mkdir $base_dir/STAR/alignment
outDir=$base_dir/STAR/alignment
```

Generate the Genome index:
``` bash
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $GDIR  \
	   --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 62
```
Note that sjdbOverhang ideally should be set to the lenght of the reads -1 (our reads are 63 bps).

Align reads with STAR:
``` bash
cd $outDir
STAR --runMode alignReads --runThreadN 4 --genomeDir $GDIR \
--readFilesIn <(zcat $fastq_files/SRR1039508_R1.fastq.gz) <(zcat $fastq_files/SRR1039508_R2.fastq.gz) \
--outFileNamePrefix sample1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
```
When running STAR \code{--quantMode TranscriptomeSAM} is essential to obtain the transcript alignments.

Use gffread to build a reference transcriptome (fasta format) compatible with the DNA fasta and gtf files used for STAR:
``` bash
gffread -w cDNA.fa -g $fasta $gtf
cdna=$outDir/cDNA.fa
```

Use salmon on the transcript alignments to compute the equivalence classes:
``` bash
salmon quant -t $cdna -l A -a sample1Aligned.toTranscriptome.out.bam -o sample1 -p 4 --dumpEq
```
The option \code{--dumpEq} is essential to obtain the equivalence classes from salmon.


## Align reads to the transcriptome with salmon
Create a variable for the fasta format reference transcriptome (cDNA):
``` bash
fasta_tr=$base_dir/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.cdna.all.1.1.10M.fa.gz
```

Make the directory for the salmon output and the genome index
``` bash
# create a directory for salmon
mkidr $base_dir/salmon

# create a directory for the genome index
mkdir $base_dir/salmon/Salmon_index
idx=$base_dir/salmon/Salmon_index

# create a directory for the output of the alignment
mkdir $base_dir/salmon/alignment
out_Salmon=$base_dir/salmon/alignment
```

Build salmon index
``` bash
salmon index -i $idx -t $fasta_tr -p 4 --type quasi -k 31
```

Align reads and quantify transcript abundance with salmon:
``` bash
salmon quant -i $idx -l A -1 $fastq_files/SRR1039508_R1.fastq.gz -2 $fastq_files/SRR1039508_R2.fastq.gz \
-p 4 -o $out_Salmon/sample1 --seqBias --gcBias --dumpEq
```
The option \code{--dumpEq} is essential to obtain the equivalence classes from salmon.
