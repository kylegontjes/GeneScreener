# Gene Screener 

Snakemake for quickly screening gene presence in genes.

The current version is configured to run BLASTN using a predefined database on a collection of one or more isolates. However, this pipeline can be modified to run any version of BLAST.

# Create an index of the reference genome
```
module load Bioinformatics

module load bwa

bwa index [reference genome]
```

# Two purposes
## 1. Screen for the presence and position of genes in your database 
Example: Identify the presence and position of blaKPC genes in genomes

Query: Genome(s)

Database: Set of genes

Snakemake: GeneScreener.smk

Config: config/config.yaml

## 2. Identify genes in the genome that correspond to a query/queries
Example: Identify what genome genes correspond to those in an operon of interest

Query: Set of genes

Database: Genome

Snakemake: GeneScreener_2.smk

Config: config/config_2.yaml

## Installation
```
git clone https://github.com/kylegontjes/GeneScreener.git
```
# Easy run

## Set up sample directory
```
path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"
sample_id="sample_id"
sample_names=$(ls -1 $path | grep fasta | cut -d. -f1 | sort | uniq)
echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv
```

## Singularity and snakemake are necessary to run this pipeline
```
module load singularity
module load snakemake
```

## Dry run to check the ability to run
### 1. Screen for the presence and position of genes in your database 
```
snakemake -s GeneScreener.smk --dryrun -p
```
### 2. Identify genes in the genome that correspond to a query/queries
```
snakemake -s GeneScreener_2.smk --dryrun -p
```

## Full run

Sbat files corresponding to each implementation have been included to enable quick implementation. Please edit the sbat files accordingly.
