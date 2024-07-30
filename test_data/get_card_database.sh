#!/bin/sh

# Create directories
## Download CARD Dataset 
cd ../database
wget https://card.mcmaster.ca/download/0/broadstreet-v3.2.6.tar.bz2
tar -xvf broadstreet-v3.2.6.tar.bz2
 
grep 'KPC beta-lactamase' aro_index.tsv | cut -f1 | cut -d: -f2 > kpc_aro_accessions
## Fasta files
touch kpc_nuc_db.fasta

for accession in `cat kpc_aro_accessions`;
do
  sed -n "/$accession/,/>/{p;}" nucleotide_fasta_protein_homolog_model.fasta | sed '$d' >> kpc_nuc_db.fasta
done