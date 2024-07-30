#!/bin/sh

# Create master contig list
cd ../queries
touch crkp_contigs.fasta

for isolate in `ls /nfs/esnitkin/Project_Penn_KPC/Sequence_data/assembly/2021-02-15_assembly/PCMP_H10*fasta`
  do
  echo $isolate
  cat $isolate >> crkp_contigs_1.fasta
done

for isolate in `ls /nfs/esnitkin/Project_Penn_KPC/Sequence_data/assembly/2021-02-15_assembly/PCMP_H20*fasta`
  do
  echo $isolate
  cat $isolate >> crkp_contigs_2.fasta
done