#!/bin/bash
bam=$1
genomeCoverageBed=$2
chrmSizes=$3
wigToBigWig=$4
species=$5

$genomeCoverageBed -split -bg -ibam $1 -g ~/mount/publicdata/$species/chrmSizes.$species > $1.bedgraph
$genomeCoverageBed -split -strand + -bg -ibam $1 -g ~/mount/publicdata/$species/chrmSizes.$species > $1.plus.bedgraph
$genomeCoverageBed -split -strand - -bg -ibam $1 -g ~/mount/publicdata/$species/chrmSizes.$species > $1.minus.bedgraph
NUM_FRAGS=`samtools view -c $1`
cat $1.bedgraph | awk -v NUM_FRAGS=$NUM_FRAGS 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/(NUM_FRAGS/1000000)}' > $1.normalised.bedgraph
cat $1.plus.bedgraph | awk -v NUM_FRAGS=$NUM_FRAGS 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/(NUM_FRAGS/1000000)}' > $1.plus.normalised.bedgraph
cat $1.minus.bedgraph | awk -v NUM_FRAGS=$NUM_FRAGS 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/(NUM_FRAGS/1000000)}' > $1.minus.normalised.bedgraph
rm -f $1.bedgraph
rm -f $1.plus.bedgraph
rm -f $1.minus.bedgraph
$wigToBigWig $1.normalised.bedgraph ~/mount/publicdata/$species/chrmSizes.$species $1.bigWig
$wigToBigWig $1.plus.normalised.bedgraph ~/mount/publicdata/$species/chrmSizes.$species $1.plus.bigWig
$wigToBigWig $1.minus.normalised.bedgraph ~/mount/publicdata/$species/chrmSizes.$species $1.minus.bigWig
rm -f $1.normalised.bedgraph
rm -f $1.plus.normalised.bedgraph
rm -f $1.minus.normalised.bedgraph
