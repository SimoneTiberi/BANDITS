# BANDITS: Bayesian ANalysis of DIfferenTial Splicing

<img src="inst/extdata/BANDITS.png" width="200" align="right"/> 

`BANDITS` is a Bayesian hierarchical model for detecting differential splicing of genes and transcripts,
via differential transcript usage (DTU), 
between two or more conditions.
The method uses a Bayesian hierarchical framework, which allows for sample specific proportions
in a Dirichlet-multinomial model, and samples the allocation of fragments to the transcripts.
Parameters are inferred via Markov chain Monte Carlo (MCMC) techniques and a DTU test is performed 
via a multivariate Wald test on the posterior densities for the average relative abundance of transcripts.

>  Simone Tiberi and Mark D Robinson (2020).
BANDITS: Bayesian differential splicing accounting for sample-to-sample variability and mapping uncertainty.
> 
> *Genome Biology* **21 (69)**.
doi: [10.1186/s13059-020-01967-8](https://doi.org/10.1186/s13059-020-01967-8)

## Bioconductor installation 
`BANDITS` is available on [Bioconductor](https://bioconductor.org/packages/BANDITS) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("BANDITS")
```

## Vignette
The vignette illustrating how to use the package can be accessed on 
[Bioconductor](https://bioconductor.org/packages/BANDITS)
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
These can be obtained by aligning reads either directly to a reference transcriptome with pseudo-alignmers, via `salmon` or `kallisto`, or to a reference genome with splice-aware genome alignment algorithms, via `STAR`, and checking the transcripts compatible with each genome alignment with `salmon`

NOTE: when using `salmon`, use the option `--dumpEq` to obtain the equivalence classes, when using `STAR`, use the option `--quantMode TranscriptomeSAM` to obtain alignments translated into transcript coordinates, and when using `kallisto`, run both the `quant` and `pseudo` modes to obtain the transcript estimated counts and equivalence classes, respectively.

Below we show three pipelines for aligning reads with `salmon`, `kallisto` and `STAR`.

### Download the example input data
To obtain the example raw data, download or clone the ARMOR github repository:
``` bash
git clone https://github.com/csoneson/ARMOR.git

# set a base_dir variable to the downloaded repo
base_dir="~/ARMOR"

# input reads:
fastq_files=$base_dir/example_data/FASTQ
```

The example data consits of four paired-end RNA-seq reads of 63 base pairs.

### Align reads to the transcriptome with salmon
Create a variable for the fasta format reference transcriptome (cDNA):
``` bash
fasta_tr=$base_dir/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.cdna.all.1.1.10M.fa.gz
```

Make the directory for the salmon output and the genome index
``` bash
# create a directory for salmon
mkdir $base_dir/salmon

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
The option `--dumpEq` is essential to obtain the equivalence classes from salmon.

In the output folder (`$out_Salmon/sample1`), the file `quant.sf` contains the estimated transcripts abundances, while the equivalence classes (and respective counts) are stored in `aux_info/eq_classes.txt`.

### Align reads to the transcriptome with kallisto
Create a variable for the fasta format reference transcriptome (cDNA):
``` bash
fasta_tr=$base_dir/example_data/reference/Ensembl.GRCh38.93/Homo_sapiens.GRCh38.cdna.all.1.1.10M.fa.gz
```

Make the directory for the kallisto output and the genome index
``` bash
# create a directory for kallisto
mkdir $base_dir/kallisto

# create a directory for the genome index
mkdir $base_dir/kallisto/kallisto_index
idx=$base_dir/kallisto/kallisto_index

# create a directory for the output of the alignment
mkdir $base_dir/kallisto/alignment
out_kallisto_alignment=$base_dir/kallisto/alignment

# create a directory for the output of the equivalence classes
mkdir $base_dir/kallisto/pseudo
out_kallisto_pseudo=$base_dir/kallisto/pseudo
```

Build kallisto index
``` bash
kallisto index -i $idx/transcripts.idx $fasta_tr
```

Align reads and quantify transcript abundance with kallisto:
``` bash
kallisto quant -i $idx/transcripts.idx -o $out_kallisto_alignment/sample1 --bias --threads 4 \
$fastq_files/SRR1039508_R1.fastq.gz $fastq_files/SRR1039508_R2.fastq.gz 
```

In the output folder (`$out_kallisto_alignment/sample1`), the file `abundance.tsv` contains the estimated transcripts abundances.

Compute the equivalence classes and respective counts with kallisto:
``` bash
kallisto pseudo -i $idx/transcripts.idx -o $out_kallisto_pseudo/sample1  \
$fastq_files/SRR1039508_R1.fastq.gz $fastq_files/SRR1039508_R2.fastq.gz 
```

In the output folder (`$out_kallisto_pseudo/sample1`), the file `pseudoalignments.ec` contains the transcripts forming each equivalence class, while `pseudoalignments.tsv` contains the equivalence classes counts.

### Align reads to the genome with STAR, and compute the equivalence classes with salmon
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
When running STAR `--quantMode TranscriptomeSAM` is essential to obtain the transcript alignments.

Use gffread to build a reference transcriptome (fasta format) compatible with the DNA fasta and gtf files used for STAR:
``` bash
gffread -w cDNA.fa -g $fasta $gtf
cdna=$outDir/cDNA.fa
```

Use salmon on the transcript alignments to compute the equivalence classes:
``` bash
salmon quant -t $cdna -l A -a sample1Aligned.toTranscriptome.out.bam -o sample1 -p 4 --dumpEq
```
The option `--dumpEq` is essential to obtain the equivalence classes from salmon.

In the output folder (`$outDir/sample1`), the file `quant.sf` contains the estimated transcripts abundances, while the equivalence classes (and respective counts) are stored in `aux_info/eq_classes.txt`.
