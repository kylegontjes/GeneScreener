# GeneScreener
Snakemake to quickly screen gene presence in genes. 

The current version (8/23) is set to run blastn using a set database on a collection of one or more isolates. However, this pipeline can altered to run any flavor of blast.

# Installation
git clone https://github.com/kylegontjes/GeneScreener.git

## Easy run
### Singularity and snakemake are necessary 
module load singularity

module load snakemake

### Do dry run to check ability to run
snakemake -s GeneScreener.smk --dryrun -p

### Full run
snakemake -s GeneScreener.smk --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 30

# To run many isolates at the same time (and possibly close the computer, try running this command)
srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash
