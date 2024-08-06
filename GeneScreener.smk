# Author: Kyle Gontjes

configfile: "config/config.yaml"

import pandas as pd
import os
import re 

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id']) 

if not os.path.exists("results/"):
    try:
        os.makedirs("results/")
    except OSError as e:
        print(f"Error creating directory: {e}")

if not os.path.exists("logs"):
    try:
        os.makedirs("logs")
    except OSError as e:
        print(f"Error creating directory: {e}")

rule all:
    input:
        database_nbd = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".ndb"),
        blast_out = expand("results/{sample}_blast.out",sample=SAMPLE),
        presence_mat = expand("results/{sample}_blast_mat.tsv",sample=SAMPLE)
        
# Step 1: Create dictionary for blast
rule create_db:
    input:
        database_fasta = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".fasta")
    output:
        database_ndb = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + str(".ndb")) 
    params:   
        database_directory=config["database_directory"],
        database_name=config["database_name"],
        dbtype=config["dbtype"],
        dbversion=config["dbversion"] 
    singularity:
        "docker://staphb/blast:2.15.0"
    shell:
        """
        cd {params.database_directory}
        makeblastdb -in {input.database_fasta} -out {params.database_name}  -dbtype {params.dbtype} -blastdb_version {params.dbversion}
        """

rule blast:
    input: 
        genome_fasta = str(config['queries']) + "/" + "{sample}",
        blast_ndb = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".ndb")
    output:
        blast_out = f"results/{{sample}}_blast.out" 
    params:     
        outfmt=config["outfmt"], 
        evalue=config["evalue"], 
        max_target_seqs=config["max_target_seqs"],
        blast_database=str(config["database_directory"]) + "/" + str(config["database_name"])
    log:
        blast_log = "logs/{sample}_blast.log"
    singularity:
        "docker://staphb/blast:2.15.0"
    shell: 
        "blastn -query {input.genome_fasta} -db {params.blast_database} -out {output.blast_out} -outfmt {params.outfmt} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} &> {log.blast_log}"    
    
rule clean_blast:
    input: 
        blast_out = str("results/" + "{sample}" + "_blast.out"), 
        blast_db = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".fasta")
    output: 
        presence_mat = f"results/{{sample}}_blast_mat.tsv",
        blast_annots = f"results/{{sample}}_blast_annots.tsv"
    conda:
        "envs/R.yaml"
    script:
        "scripts/curate_blast_results.R"