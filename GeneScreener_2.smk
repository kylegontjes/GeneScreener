# Author: Kyle Gontjes

configfile: "config/config_2.yaml"

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
        database_nbd = expand(str(config["database_directory"]) + "/" + "{sample}" + ".pdb",sample=SAMPLE),
        blast_out = expand("results/{sample}_blast.out",sample=SAMPLE),
        presence_mat = expand("results/{sample}_blast_clean.tsv",sample=SAMPLE)
        
# Step 1: Create dictionary for blast
rule create_db:
    input:
        database_fasta = str(config['sample_file_path']) + "/" + "{sample}"
    output:
        database_ndb = expand(str(config["database_directory"]) + "/" + str("{{sample}}") + str(".pdb")) 
    params:   
        database_directory=config["database_directory"], 
        dbtype=config["dbtype"],
        dbversion=config["dbversion"],
        dbname=str("{sample}")
    singularity:
        "docker://staphb/blast:2.15.0"
    shell:
        """
        cd {params.database_directory}
        makeblastdb -in {input.database_fasta} -out {params.dbname}  -dbtype {params.dbtype} -blastdb_version {params.dbversion}
        """

rule blast:
    input: 
        query = str(config["query"]),
        blast_ndb = str(config["database_directory"]) + "/" + str("{sample}") + ".pdb"
    output:
        blast_out = f"results/{{sample}}_blast.out" 
    params:     
        outfmt=str(config["outfmt"]), 
        evalue=config["evalue"],  
        culling_limit=config["culling_limit"],  
        blastdb=str(config["database_directory"]) + "/" + str("{sample}") 
    singularity:
        "docker://staphb/blast:2.15.0"
    shell: 
        "blastx -query {input.query} -db {params.blastdb} -out {output.blast_out} -outfmt {params.outfmt} -evalue {params.evalue} -culling_limit {params.culling_limit} -subject_besthit"    
    
rule clean_blast:
    input: 
        blast_out = str("results/" + "{sample}" + "_blast.out"),
        query = str(config["query"])
    output: 
        blast_annots = f"results/{{sample}}_blast_clean.tsv" 
    params:
        directory=config["directory"],
        outfmt=str(config["outfmt"]) 
    singularity:
        "docker://rocker/tidyverse"
    script:
        "scripts/curate_blast_results.R"