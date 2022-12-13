#!/usr/bin/env zsh
# This script ingests a csv or vcf file, formats into a primer3 compatible file,
# designs nested primers and outputs them as a fasta file, then analyzes them for dimers and secondary structure, and finally blasts for potential amplicons.
# Additionally, you can designate if you want your primer designs to be for NGS or PCR based
# Dependencies: Primer3, MFEPrimer3.1, BLASTn2.2.18, bedtools
# See README file for further information on setup and dependencies
# Author: blake hummer
# Date: 2022/12/12

# Help documentation
help() {
    echo "Usage: $(basename $0) [-dsvglh]"
    echo '   -d   Designates the design schema: 1 designs NGS primers that flank your ROI, 2 designs PCR primers that overlap with your ROI, 3 designs taqman based primers with an internal probe that interescts your ROI'
    echo '   -s   Requires a comma separated file with a sequence and associated name, for example: seq1,atcgatcg'
    echo '   -v   Requires a variant call file'
    echo '   -g   Requires a reference genome as a fasta file (e.g., hg38.fa)'
    echo '   -l   Requires a chromosome length file (e.g., hg38.chrom.sizes)'
    echo '   -h   Shows the help documentation. Example of command: ./PrimerDesigner.sh -d 2 -v file.vcf -g hg38.fa -l hg38.chrom.sizes'
    exit 1
}

# Establish options
while getopts ":d:s:v:g:l:h" opt
   do
     case $opt in
        d ) design=$OPTARG;;
        s ) seq=$OPTARG;;
        v ) vcf=$OPTARG;;
        g ) genome=$OPTARG;;
        l ) lengths=$OPTARG;;
        h ) help;;
        ? ) help;;
     esac
done

# Handle all potential inpt files
INPUT=`awk '{OFS=","}{print $1,$2}' <(bedtools getfasta -fi $genome -bed <(bedtools slop -i <(tail -n +24 $vcf | awk -F"\t" '{OFS="\t"}{print $1,$2,$2}') -g $lengths -b 75) -tab -name)`
echo $INPUT > inputfile.txt
: | cat $seq >> inputfile.txt

# Create fasta with all sequences
sed '/^$/d' inputfile.txt | sed 's/^/>/g' | sed 's/,/\n/g' > sequence.fa

# Designate design schema
if [[ $design -eq '1' ]]
then
    echo "you chose NGS"

# Parse input csv and formats into Primer3 compatible format for inner primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "SEQUENCE_TARGET=70,10" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=20,50,80,50" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=50-140" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=3" "\n" \
    "PRIMER_OPT_TM=55" "\n" \
    "PRIMER_MIN_TM=55" "\n" \
    "PRIMER_MAX_TM=60" "\n" \
    "="}' | sed '$ s/$//' > InnerInputfile

# Parse input csv and formats into Primer3 compatible format for outer primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "SEQUENCE_TARGET=70,10" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,35,115,35" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=80-150" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=3" "\n" \
    "PRIMER_OPT_TM=55" "\n" \
    "PRIMER_MIN_TM=50" "\n" \
    "PRIMER_MAX_TM=60" "\n" \
    "="}' | sed '$ s/$//' > OuterInputfile

elif [[ $design -eq '2' ]]
then
    echo "you chose PCR"

# Parse input csv and formats into Primer3 compatible format for inner primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "SEQUENCE_TARGET=90,10" "\n" \
    "SEQUENCE_OVERLAP_JUNCTION_LIST=75" "\n" \
    "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=3" "\n" \
    "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=3" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=35,45,80,45" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=50-140" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=4" "\n" \
    "PRIMER_OPT_TM=66" "\n" \
    "PRIMER_MIN_TM=64" "\n" \
    "PRIMER_MAX_TM=68" "\n" \
    "="}' | sed '$ s/$//' > InnerInputfile

# Parse input csv and formats into Primer3 compatible format for outer primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "SEQUENCE_TARGET=95,10" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,35,115,50" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=80-150" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=3" "\n" \
    "PRIMER_OPT_TM=66" "\n" \
    "PRIMER_MIN_TM=57" "\n" \
    "PRIMER_MAX_TM=68" "\n" \
    "="}' | sed '$ s/$//' > OuterInputfile

elif [[ $design -eq '3' ]]
then
    echo "you chose internal probe PCR"

# Parse input csv and formats into Primer3 compatible format for inner primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "PRIMER_PICK_INTERNAL_OLIGO=1" "\n" \
    "SEQUENCE_TARGET=70,10" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=20,50,80,50" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=50-140" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=3" "\n" \
    "PRIMER_OPT_TM=55" "\n" \
    "PRIMER_MIN_TM=55" "\n" \
    "PRIMER_MAX_TM=60" "\n" \
    "="}' | sed '$ s/$//' > InnerInputfile

# Parse input csv and formats into Primer3 compatible format for outer primers
    sed '/^$/d' inputfile.txt | awk -F, '{print "SEQUENCE_ID=" $1, "\n" "SEQUENCE_TEMPLATE=" $2}' | awk '1; NR%2 == 0 {print "PRIMER_TASK=generic" "\n" \
    "PRIMER_PICK_LEFT_PRIMER=1" "\n" \
    "PRIMER_PICK_RIGHT_PRIMER=1" "\n" \
    "SEQUENCE_TARGET=70,10" "\n" \
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,35,115,35" "\n" \
    "PRIMER_OPT_SIZE=18" "\n" \
    "PRIMER_MIN_SIZE=15" "\n" \
    "PRIMER_MAX_SIZE=28" "\n" \
    "PRIMER_MAX_NS_ACCEPTED=1" "\n" \
    "PRIMER_PRODUCT_SIZE_RANGE=80-150" "\n" \
    "P3_FILE_FLAG=1" "\n" \
    "PRIMER_EXPLAIN_FLAG=1" "\n" \
    "PRIMER_MAX_END_GC=3" "\n" \
    "PRIMER_OPT_TM=55" "\n" \
    "PRIMER_MIN_TM=50" "\n" \
    "PRIMER_MAX_TM=60" "\n" \
    "="}' | sed '$ s/$//' > OuterInputfile

fi

# Design inner primers, grab the highest scoring primers for each given sequence and reformat into fasta file
primer3_core InnerInputfile

for i in *.for
do
    echo `basename "$i"` | sed 's/\ .*//g' | sed 's/^/>for_/g'
    if [[ $(wc -l < $i) -gt 5 ]]
    then
        sed -n '4p' $i | awk '{print $2}'
    else
        echo 'nnnn'
    fi
done > InnerLeftPrimer.fa

for i in *.rev
do
    echo `basename "$i"` | sed 's/\ .*//g' | sed 's/^/>rev_/g'
    if [[ $(wc -l < $i) -gt 5 ]]
    then
        sed -n '4p' $i | awk '{print $2}'
    else
        echo 'nnnn'
    fi
done > InnerRightPrimer.fa

# Merge the forward and reverse primer files into one
cat InnerLeftPrimer.fa InnerRightPrimer.fa > InnerPrimers.fa

mkdir I_Primers

if [[ $design -eq '3' ]]
then
    for i in *.int
    do
        echo `basename "$i"` | sed 's/\ .*//g' | sed 's/^/>int_/g'
        if [[ $(wc -l < $i) -gt 5 ]]
        then
           sed -n '4p' $i | awk '{print $2}'
        else
            echo 'nnnn'
        fi
    done >> InnerPrimers.fa
    mv *.int ./I_Primers
else
    echo ""
fi

# Then hide the source files containing the remaining primer options
mv *.for ./I_Primers
mv *.rev ./I_Primers

# Repeat the above logic but with outer primers
primer3_core OuterInputfile

for i in *.for
do
    echo `basename "$i"` | sed 's/\ .*//g' | sed 's/^/>for_/g'
    if [[ $(wc -l < $i) -gt 5 ]]
    then
        sed -n '4p' $i | awk '{print $2}'
    else
        echo 'nnnn'
    fi
done > OuterLeftPrimer.fa

for i in *.rev
do
    echo `basename "$i"` | sed 's/\ .*//g' | sed 's/^/>rev_/g'
    if [[ $(wc -l < $i) -gt 5 ]]
    then
        sed -n '4p' $i | awk '{print $2}'
    else
        echo 'nnnn'
    fi
done > OuterRightPrimer.fa

cat OuterLeftPrimer.fa OuterRightPrimer.fa > OuterPrimers.fa

mkdir O_Primers
mv *.for ./O_Primers
mv *.rev ./O_Primers

# Move all the input files into their own directory
mkdir InputFiles
mv sequence.fa ./InputFiles

# Ingest fasta file and do dimer and hairpin analysis for inner primers
mfeprimer dimer -i InnerPrimers.fa > InnerResults.txt
mfeprimer hairpin -i InnerPrimers.fa >> InnerResults.txt

# Ingest fasta file and do dimer and hairpin analysis for outer primers
mfeprimer dimer -i OuterPrimers.fa > OuterResults.txt
mfeprimer hairpin -i OuterPrimers.fa >> OuterResults.txt

# Blast inner and outer primers WITHOUT MBs for potential amplicons
blastn -task blastn-short -db $genome -query InnerPrimers.fa > InnerPrimerBlast.txt
blastn -task blastn-short -db $genome -query OuterPrimers.fa > OuterPrimerBlast.txt

# Move all result files into thier own directory
mkdir ResultFiles
mv *Result*.* ./ResultFiles
mv *Blast.txt ./ResultFiles

mkdir Primers
mv *Primer*.* ./Primers