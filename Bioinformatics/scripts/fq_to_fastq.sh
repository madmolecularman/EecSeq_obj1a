#!/bin/bash
set -eu 
 
# Find all .fq file and rename them to .fastq
find . -name "*.fq" -exec rename .fq .fastq {} \;