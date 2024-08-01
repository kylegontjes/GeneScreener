# GeneScreener
Snakemake to quickly screen gene presence in genes 

# Installation
git clone https://github.com/kylegontjes/GeneScreener.git

## Easy run
### Singularity and snakemake are necessary 
module load singularity
module load snakemake

### Do dry run to check ability to run
snakemake -s GeneScreener.smk --dryrun -p

### Full run
snakemake -s GeneScreener.smk--use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 30