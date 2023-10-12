# Probe Analysis

## Justification for the presence of probes

- Mixed files were a molecular protocol error G3 and N2
- Overabundant reads in N1 files run1 and run2
    - I think we have probes in here
    - Every CaptureX+N1 read in the discard file has the same FASTQ header:
        - ATTACTCG+AGGCTATA which corresponds to the 501/701 primers inline barcode
        - These primers were used on the 150 bp cDNA libraries for amplification
        - Why would Capture3 and Capture4 have these same barcodes?
    
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
    
- What we are seeing here is the discarded read files from process_shortreads program. These reads were binned into the N1 files but are not N1 gDNA reads. They are pooled cDNA probes. I was able to pull out the exact cut site for the SAII enzyme from all of them.
    - Using a random BLAST of these sequences they are oyster reads.
- Why did our probes get sequenced?
    - Our probes were not fully treated with the restriction enzyme which means they still contained the adapter for a sequenced run
        - Improper or ineffective cutting of this region during the molecular protocol could be the issue
    - Other options to discuss
        - Our probes were not biotinylated fully meaning they couldn’t be captured by the streptavidin beads
            - One thing to note here is that Biointylated probes
        - We did not add enough Streptavidin beads or the reaction was efficient (beads not being viable) which lead to probes not being captured and still retained in the effleunt of the Hybridization section.
    - Where do we go from here
        - Initially it was dissapointing to know that our gDNA reads were compromised by the presence of probes. We would have gotten a lot more data about our gDNA if these were not there.
        - The presence of these probe reads means that we could possibly use them, but to what end?
            - We are concerned with how our different size probes are able to enrich our gDNA regions. Knowing what probes are in this reaction could potentially beenfit this process
            - We could clean and map these probes to the Oyster genome
                - If they overlap with the regions that were enrich we know that these reads were potentially hybridized
                - If they do not overlap we know then that these reads did not enrich these regions
        - Options
            - Move forward with subsample reads analysis w/o discarded sequenced probes
            - Increase RE enzyme from 4 to 24 hours
            - Switch out RE enzyme
            - Map probes and compare with our gDNA reads
            - Take a look at Capture and UCE literature for probe, bead, and gDNA ratio
    - What does this mean for Obj1b (probe and gDNA diversity)?
        - Something in our protocol isnt working right. I am not confident that I haven’t made the same mistake again and our currently synthesized probes could:
            
            1) Contain the cut site after RE treatment
            
            2) Biotinylation may not have been effective
            
            3) Our current hybridization protocol may need some work
            
            4) There are reagents in this pipeline that are problematic
            
        - We could run a higher percentage agarose gel of pre and psot RE treatment libraries to try and separate them out
        - We could also review the tapestations below to try and tease this out but the resolution and large modal peaks, with biointylated probes may make this difficult
        - We could run a small URI sequencing run pre and post RE to figure this out (expensive)
    - What direction do I move in now?

## Bioinformatic analysis of probes

In each contaminanttest directory contains probes sequences within the discard files. We need to pull out only those fastq files that have the adapter cutsite in both the forward and reverse reads?

```bash
cat Capture1_N1.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture1_probes.F.fq
cat Capture1_N1.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture1_probes.R.fq
cat Capture1_N1.reseq.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture1_reseq_probes.F.fq
cat Capture1_N1.reseq.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture1_reseq_probes.R.fq
cat Capture2_N1.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture2_probes.F.fq
cat Capture2_N1.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture2_probes.R.fq
cat Capture2_N1.reseq.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture2_reseq_probes.F.fq
cat Capture2_N1.reseq.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture2_reseq_probes.R.fq
cat Capture3_N1.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture3_probes.F.fq
cat Capture3_N1.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture3_probes.R.fq
cat Capture3_N1.reseq.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture3_reseq_probes.F.fq
cat Capture3_N1.reseq.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture3_reseq_probes.R.fq
cat Capture4_N1.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture4_probes.F.fq
cat Capture4_N1.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture4_probes.R.fq
cat Capture4_N1.reseq.R1.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture4_reseq_probes.F.fq
cat Capture4_N1.reseq.R2.fastq.gz.discards | grep -E -B 1 -A 2 '^.{7}GTCGAC' > Capture4_reseq_probes.R.fq
```