#!/bin/bash
set -euo pipefail
ACC="GCA_019448385.1"
cd data/raw
datasets download genome accession "$ACC" --include genome --filename dlutheri.zip
unzip -o -q dlutheri.zip
cp ncbi_dataset/data/"$ACC"/*.fna dlutheri_genome.fa
cd ../..
seqkit stats -a -T data/raw/dlutheri_genome.fa > results/genome_stats.tsv
seqkit fx2tab -n -l -g data/raw/dlutheri_genome.fa | sort -t$'\t' -k2,2nr > results/contig_lengths.tsv
echo "== Whole genome =="; cat results/genome_stats.tsv
echo "== 5 longest contigs =="; head -n 5 results/contig_lengths.tsv
