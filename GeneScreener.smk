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

rule all:
    input:
        database_nbd = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".ndb"),
        blast_out = expand("results/{sample}_blast.out",sample=SAMPLE),
        presence_mat = expand("results/{sample}_blast_clean.tsv",sample=SAMPLE)
        
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
        genome_fasta = str(config['queries']) + "/" + "{sample}" + ".fasta",
        blast_ndb = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".ndb")
    output:
        blast_out = f"results/{{sample}}_blast.out" 
    params:     
        outfmt=str(config["outfmt"]), 
        evalue=config["evalue"],  
        culling_limit=config["culling_limit"],
        blast_database=str(config["database_directory"]) + "/" + str(config["database_name"]) 
    singularity:
        "docker://staphb/blast:2.15.0"
    shell: 
        "blastn -query {input.genome_fasta} -db {params.blast_database} -out {output.blast_out} -outfmt {params.outfmt} -evalue {params.evalue} -culling_limit {params.culling_limit}"    
    
rule clean_blast:
    input: 
        blast_out = str("results/" + "{sample}" + "_blast.out"),
        database_fasta = expand(str(config["database_directory"]) + "/" + str(config["database_name"]) + ".fasta")
    output: 
        blast_annots = f"results/{{sample}}_blast_clean.tsv" 
    params:
        directory=config["directory"],
        outfmt=str(config["outfmt"])
    singularity:
        "docker://rocker/tidyverse"
    script:
        "scripts/curate_blast_results.R"
