#!/usr/bin/env bash
# Ingests a csv file as the first argument, formats into a Primer3 compatible file,
# designs primers and outputs them as a fasta file, then analyzes them for dimers, secondary structure,
# and off-targeting.
# Dependencies: Primer3, MFEPrimer3.1 with human reference genome indexed (i.e., hg38.fa)
# Author: blake hummer
# Date: 2022/11/18
# Example: ./PrimerDesignerandQC.sh ./input.txt

input=$1

# Reformat to Primer3 compatible file
sed '/^$/d' $input | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
"PRIMER_PICK_LEFT_PRIMER=1" "\n" \
"PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
"PRIMER_OPT_SIZE=18" "\n" \
"PRIMER_MIN_SIZE=15" "\n" \
"PRIMER_MAX_SIZE=28" "\n" \
"PRIMER_MAX_NS_ACCEPTED=1" "\n" \
"PRIMER_PRODUCT_SIZE_RANGE=115-150" "\n" \
"P3_FILE_FLAG=1" "\n" \
"PRIMER_EXPLAIN_FLAG=1" "\n" \
"PRIMER_MAX_END_GC=3" "\n" \
"PRIMER_OPT_TM=65" "\n" \
"PRIMER_MIN_TM=63" "\n" \
"PRIMER_MAX_TM=66" "\n" \
"="}' | sed '$ s/$//' > Inputfile

# Create fasta with all sequences
sed '/^$/d' $input | sed 's/^/>/g' | sed 's/,/\n/g' > sequence.fa

# Design primers, grab the highest scoring primers for each given sequence and reformat into fasta file
primer3_core Inputfile | grep 'PRIMER_LEFT_0_SEQUENCE' | nl -w 2 -s _  | sed 's/^/>/g' | sed 's/\>[[:space:]]/>/g' | sed 's/=/\n/g' > LeftPrimer.fa
primer3_core Inputfile | grep 'PRIMER_RIGHT_0_SEQUENCE' | nl -w 2 -s _ | sed 's/^/>/g' | sed 's/\>[[:space:]]/>/g' | sed 's/=/\n/g' > RightPrimer.fa

cat LeftPrimer.fa RightPrimer.fa > result.fa

# Ingest fasta file and do analyses
mfeprimer dimer -i result.fa -d hg38.fa > results.txt
