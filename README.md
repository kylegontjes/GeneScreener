# Gene Screener 

Snakemake to screen gene presence in genes quickly.

The current version (8/23) is set to run blastn using a set database on a collection of one or more isolates. However, this pipeline can altered to run any flavor of blast.

# Create index of reference genome
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
Example: Identify what genome genes correspond to that in an operon of interest

Query: Set of genes

Database: Genome

Snakemake: GeneScreener_2.smk

Config: config/config_2.yaml

## Installation
```
git clone https://github.com/kylegontjes/GeneScreener.git
```
# Easy run

## Set  up sample directory
```
path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"
sample_id="sample_id"
sample_names=$(ls -1 $path | grep fasta | cut -d. -f1 | sort | uniq)
echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv
```

## Singularity and snakemake are necessary
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
### 1. Screen for the presence and position of genes in your database 
```
snakemake -s GeneScreener.smk --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 30
```
### 2. Identify genes in the genome that correspond to a query/queries
```
snakemake -s GeneScreener_2.smk --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config_2.yaml --latency-wait 30
```
# To run many isolates at the same time (and possibly close the computer, try running this command)
```
srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash
```
