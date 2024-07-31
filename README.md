# GeneScreener
Snakemake to quickly screen gene presence in genes 

# Installation
git clone https://github.com/kylegontjes/GeneScreener.git

## Easy run
### Singularity and snakemake are necessary 
module load singularity
module load snakemake
snakemake -s GeneScreener.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 60
