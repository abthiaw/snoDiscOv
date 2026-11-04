#!/usr/bin/env bash

gtfFile=$1
outGtf=$2

sed 's/""\([^"]*\)""/"\1"/g' "$gtfFile" \
   | sed 's/"\(.*\)"$/\1/;s/^"\(.*\)"/\1/' \
   > "$outGtf"
