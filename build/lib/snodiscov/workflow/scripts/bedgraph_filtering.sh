#!/usr/bin/env bash

bedgraphFile=$1  # coco bedgraph file
bedFile=$2 # bedgraph file
output=$3 # output intesect results

awk '{OFS="\t"; print ("chr"$1,$2,$3,$4,$5,$6)}' "$bedFile" \
    | bedtools intersect \
        -wb \
        -a "$bedgraphFile" \
        -b stdin \
         > "$output"
