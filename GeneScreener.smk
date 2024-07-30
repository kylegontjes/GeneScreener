# Author: Kyle Gontjes

configfile: "config/config.yaml"

import pandas as pd
import os
import re
PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id']) 

if not os.path.exists("results/" + PREFIX):
    try:
        os.makedirs("results/" + PREFIX)
    except OSError as e:
        print(f"Error creating directory: {e}")

rule all:
    input:
        database_nbd = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(f".ndb")),
        blast_out =  expand("results/{prefix}/{sample}/{sample}_blast.out",prefix=PREFIX,sample=SAMPLE) 
        
# Step 1: Create dictionary for blast
rule create_db:
    input:
        database_fasta = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(f".fasta"))
    output:
        database_ndb = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(".ndb")) 
    params:   
        database_directory=config["database_directory"],
        database_name=config["database_name"],
        dbtype=config["dbtype"],
        dbversion=config["dbversion"]
    log:
        create_db_log = "create_blast_db.log"
    singularity:
        "docker://staphb/blast:2.15.0"
    shell:
        """
        cd {params.database_directory}
        makeblastdb -in {input.database_fasta} -out {params.database_name}  -dbtype {params.dbtype} -blastdb_version {params.dbversion} &> {log.create_db_log}
        """

rule blast:
    input: 
        genome_fasta = lambda wildcards: expand(str(config["queries"])  + "/" + f"{wildcards.sample}.fasta"),
        blast_database = lambda wildcards: expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(f".fasta")),
        blast_ndb = lambda wildcards: expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(f".ndb"))
    output:
        blast_out = f"results/{{prefix}}/{{sample}}/{{sample}}_blast.out"
    params:     
        outfmt=config["outfmt"], 
        evalue=config["evalue"], 
        max_target_seqs=config["max_target_seqs"]
    log:
        blast_log = "logs/{prefix}/{sample}/{sample}_blast.log"
    singularity:
        "docker://staphb/blast:2.15.0"
    shell: 
        "blastn -query {input.genome_fasta} -db {input.blast_database} -out {output.blast_out} -outfmt {params.outfmt} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} &> {log.blast_log}"    