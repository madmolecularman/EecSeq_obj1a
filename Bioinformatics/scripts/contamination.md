Looking for potential contamination in Capture1_N1 samples

Retaining ambiguous barcode samples and dropped reads from sequencing run 1

```bash
declare -a StringArray=("Capture1_N1" "Capture2_N1" "Capture3_N1" "Capture4_N1")
for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}.R2.fastq.gz -o $WORKING_DIR/01_process_reads/contaminanttest/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r -D
done
```

Retaining ambiguous barcode samples and dropped reads from sequencing run 2 
```bash
declare -a StringArray=("Capture1_N1" "Capture2_N1" "Capture3_N1" "Capture4_N1")
for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}.reseq.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}.reseq.R2.fastq.gz -o $WORKING_DIR/01_process_reads/contaminanttest/${i}reseq/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r -D
done
```

Some helpful lines of code from Jon

```bash
zcat sequence_file_here | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h

zcat sequence_file_here | head -100000 | mawk 'NR%4==1' | cut -f2 -d " " | cut -f4 -d ":" | sort | uniq -c | mawk '$1 > 20' | sort -h
```

Count number of line in each discard file in the first run that contains the "GTCGAC" SAII cut site

```bash
awk 'NR%4==2' Capture1_N1.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 7319840
awk 'NR%4==2' Capture1_N1.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 7288151
awk 'NR%4==2' Capture2_N1.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 10751251
awk 'NR%4==2' Capture2_N1.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 10656175
awk 'NR%4==2' Capture3_N1.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 3716168
awk 'NR%4==2' Capture3_N1.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 3681856
awk 'NR%4==2' Capture4_N1.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 1560735
awk 'NR%4==2' Capture4_N1.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 1545056
```

Count number of line in each discard file in the second run that contains the "GTCGAC" SAII cut site

```bash
awk 'NR%4==2' Capture1_N1.reseq.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 251903
awk 'NR%4==2' Capture1_N1.reseq.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 250024
awk 'NR%4==2' Capture2_N1.reseq.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 1984140
awk 'NR%4==2' Capture2_N1.reseq.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 1962958
awk 'NR%4==2' Capture3_N1.reseq.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 12324214
awk 'NR%4==2' Capture3_N1.reseq.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 12225601
awk 'NR%4==2' Capture4_N1.reseq.R1.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 6644491
awk 'NR%4==2' Capture4_N1.reseq.R2.fastq.gz.discards | cut -c 8-14 | grep -c "GTCGAC"
> 6591971
```