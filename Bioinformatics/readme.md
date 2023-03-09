# This code repeats the analyses done in EAGER Obj1a 2022

This code assumes that there are two subdirectories in this working directory ./DNA which contain the raw sequencing files from the EAGER 2022 Obj1a sequenicng project and  ./Genome directory which contains files from the eastern oyster genome.

## Setup conda environment
```bash
conda update conda
mamba create eagerobj1a ddocent
conda activate eagerobj1a
```
## Copying files to new directory and renaming to fit with ddocent naming convention
```bash
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_1.fq.gz Capture1_B3.R1.fq.gz
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_2.fq.gz Capture1_B3.R2.fq.gz
```

## Enmasse copying files and renaming
```bash
# Set string array with the prefix name of each read file.
declare -a StringArray=("Capture1_B3" "Capture1_G5" "Capture1_M3" "Capture1_N2" "Capture2_G3" "Capture2_K4" "Capture2_N1" "Capture3_B4" "Capture3_K3" "Capture3_M4" "Capture4_B3" "Capture4_G5" "Capture4_M3" "Capture4_N2" "Capture1_B4" "Capture1_K3" "Capture1_M4" "Capture2_B3" "Capture2_G5" "Capture2_M3" "Capture2_N2" "Capture3_G3" "Capture3_K4" "Capture3_N1" "Capture4_B4" "Capture4_K3" "Capture4_M4" "Capture1_G3" "Capture1_K4" "Capture1_N1" "Capture2_B4" "Capture2_K3" "Capture2_M4" "Capture3_B3" "Capture3_G5" "Capture3_M3" "Capture3_N2" "Capture4_G3" "Capture4_K4" "Capture4_N1")

## Check naming convention
for i in "${StringArray[@]}"
do
echo $i
done

## Copy Forward 1 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L1_1.fq.gz ${i}.R1.fastq.gz
done

## Copy Forward 2 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L2_1.fq.gz ${i}.R1.fastq.gz
done

## Copy Forward 3 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L3_1.fq.gz ${i}.reseq.R1.fastq.gz
done

## Copy Forward 4 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L4_1.fq.gz ${i}.reseq.R1.fastq.gz
done

## Copy Reverse 1 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L1_2.fq.gz ${i}.R2.fastq.gz
done

## Copy Reverse 2 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L2_2.fq.gz ${i}.R2.fastq.gz
done

## Copy Reverse 3 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L3_2.fq.gz ${i}reseq.R2.fastq.gz
done

## Copy Reverse 4 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L4_2.fq.gz ${i}reseq.R2.fastq.gz
done
```

## Run fastqc and multiqc to visualize files
```bash
fastqc -f fastq -t 20 *.fastq.gz
multiqc .
```

## Results
Multiqc plots show uneven distribution of sequence counts specifically in samples Capture1_N1.R1, Capture1_N1.R2, Capture2_N1.R1, Capture2_N1.R2, Capture3_N1.reseq.R1, Capture3_N1.reseq.R2, Capture4_N1.reseq.R1, Capture4_N1.reseq.R2. See [multiqc report](https://github.com/madmolecularman/EecSeq_obj1a/blob/main/Bioinformatics/01_process_reads/multiqc_report_20230301_run1.html)

## Renaming reads

I need to be careful with some changes to the reads. 
1) Capture 2 = Capture 1
2) Capture 1 = Capture 2
3) Capture 3 = Capture 3
4) Capture 4 = Capture 4

Move capture 1 and capture 2 files into different directories
```bash
mkdir capture1
mkdir capture2
mv Capture2* capture1/
mv Capture1* capture2/
  ```

Rename Capture1 preffix to Capture2 preffix
```bash
rename Capture1 Capture2 Capture1*
rename Capture2 Capture1 Capture2*
```

Redo fastqc and multiqc
```bash
multiqc .
```

## Fixing mixed files

Capture1_G3.R1, Capture1_G3.R2, Capture1_G3.reseq.R1, Capture1_G3.reseq.R2 are actually capture2

Capture2_G3.R1, Capture2_G3.R2, Capture2_G3.reseq.R1, Capture2_G3.reseq.R2 are actually capture1

Capture2_N2.R1, Capture2_N2.R2, Capture2_N2.reseq.R1, Capture2_N2.reseq.R1 are actually capture1

Capture1_N2.R1, Capture2_N2.R2, Capture2_N2.reseq.R1, Capture2_N2.reseq.R1 are actually capture2

```bash
rename Capture2 Capture1 Capture2*
rename Capture1 Capture2 Capture1*
```

Now we have all the files organized which you can see in this [multiqc](https://github.com/madmolecularman/EecSeq_obj1a/blob/main/Bioinformatics/01_process_reads/multiqc_report_20230301_run3.html) report. Most of the files are demultiplexed. Yet Capture[1-4]_N1 and paired reseq files have very messy 

## Processing short reads on files that are not demultiplexed

## Demultiplexing capture files

Use process_shortreads to demultiplex files from raw/ to the demux/samples/ folder. One of the major problems is that we are not sure if novogene was able to separate out the different capture pools specifically for Capture[1-4]_N1 samples and Capture[1-4]_N1 resequenced. We need to create a demultiplexing run that will filter the index and the inline barcode into the proper file. The question here is can I demultiplex just the index and then the adapter inline barcode?

Its important to understand what ![process_shortreads](https://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php) is doing and how to create the barcode files ![Stacks manual](https://catchenlab.life.illinois.edu/stacks/manual/index.php#clean)

First I have to create a barcode file that is tab delimited

```bash
nano combined_barcodes.txt
```

```text
ATCGCG	CGTGAT	Capture1
CGATGT	ACATCG	Capture2
TTAGGC	GCCTAA	Capture3
TGGCCA	TGGCCA	Capture4
```

Be careful to create this files in nane and don't edit it in the VScode text editor. VScode will change the endline with a "M^" character that is incompatible with process_shortreads program. To view if there are there use the following

```bash
cat -A <barcodes.tsv>
```

Make directories for all of the files

```bash
mkdir samples
cd samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Renaming samples directory to combined barcode samples directory

```bash
mv samples/ combined_barcode_samples
```

Try process_shortreads on one pair of files (ran from 01_process_reads/raw/ directory)

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/combined_barcode_samples/Capture1_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture2_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture3_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture4_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
```

Make smaller process_shortread log files for easier viewing

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
echo ${i}
done

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

See 01_process_reads/demux/combined_barcode_samples with each directories process_shortreads.short.log file. What we can see are large number of ambiguous reads being discarded amny of which are in the sequences not recorded. 

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 648            | 0           | 29671776           | 0             | 0                         | 29672424 |


| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 29671776 |
| Low Quality          | 0        |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 648      |

| Barcode       | Filename | Total | Retained |
|---------------|----------|-------|----------|
| ATCGCG-CGTGAT | Capture1 | 34    | 34       |
| CGATGT-ACATCG | Capture2 | 310   | 310      |
| TTAGGC-GCCTAA | Capture3 | 122   | 122      |
| TGGCCA-TGGCCA | Capture4 | 182   | 182      |

**Sequences not recorded**

| Barcode       | Total    |
|---------------|----------|
| CGATGT-CGATGT | 14620908 |
| CGATGT-GGGGGG | 32992    |
| CGATGT-AGATGT | 22422    |

This is the case for Capture[1-3]_N1 and their reseq files. Capture4_N1 does not have this issue as its adapter barcode is a palindrome. It seems like the forward and reverse reads contain the barcode not its compliment. Going to test a few things

## Using index_null process short reads to investigate issue with ambiguous reads
New directory for single inline barcodes. I am wanting to test out how we retain reads with single barcodes.

Make new barcode file

``bash
nano single_barcode.txt
```

```text
ATCGCG	Capture1
CGATGT	Capture2
TTAGGC	Capture3
TGGCCA	Capture4
```

```bash
mkdir single_barcode_samples
cd single_barcode_samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Try process_shortreads on one pair of files (ran from 01_process_reads/raw/ directory)

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/single_barcode_samples/Capture1_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture2_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture3_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture4_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
```

Make smaller process_shortread log files for easier viewing (run from single_barcode_sample)

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
echo ${i}
done

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

See 01_process_reads/demux/single_barcode_samples with each Capture1_N1 process_shortreads.short.log file. We retain more reads! Yay! But we have 14 million ambiguous barcodes? Is it not reading on of the paired files? We can also see the the ambiguous barcodes are not recorded.

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 15129309       | 32341       | 14510774           | 0             | 0                         | 29672424 |


| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 14510774 |
| Low Quality          | 32341    |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 15129309 |

| Barcode | Filename | Total    | Retained |
|---------|----------|----------|----------|
| ATCGCG  | Capture1 | 60438    | 60311    |
| CGATGT  | Capture2 | 14938002 | 14906109 |
| TTAGGC  | Capture3 | 92194    | 92011    |
| TGGCCA  | Capture4 | 71016    | 70878    |


Sequences not recorded
| Barcode | Total |
|---------|-------|
| GTTTGG  | 15032 |
| GTGTTC  | 14364 |


Now lets place exact copies of the barcodes doubling them in the file and try to pull out both the forward and reverse files with that.

## Using duplicated barcodes and inline_inline to invetigate ambigous barcode drop for R2 files

New directory for double inline barcodes. I am wanting to test out how we retain reads with double barcodes.

Make new barcode file

``bash
nano double_barcode.txt
```

```text
ATCGCG  ATCGCG	Capture1
CGATGT  CGATGT	Capture2
TTAGGC  TTAGGC	Capture3
TGGCCA  TGGCCA	Capture4
```

```bash
mkdir double_barcode_samples
cd double_barcode_samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
```

Make smaller process_shortread log files for easier viewing (run from single_barcode_sample)

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
echo ${i}
done

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

I am now more confused. So here we see a drop in about half of our reads in the ambiguous reads. We have lost a lot of the retain barcodes from the other capture compaired to the single approach.

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 14658692       | 31070       | 14982662           | 0             | 0                         | 29672424 |

| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 14982662 |
| Low Quality          | 31070    |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 14658692 |

| Barcode       | Filename | Total    | Retained |
|---------------|----------|----------|----------|
| ATCGCG-ATCGCG | Capture1 | 34       | 34       |
| CGATGT-CGATAT | Capture2 | 14687394 | 14656325 |
| TTAGGC-TTAGGC | Capture3 | 2152     | 2151     |
| TGGCCA-TGGCCA | Capture4 | 182      | 182      |

Sequences not recorded
| Barcode       | Total |
|---------------|-------|
| CGATGT-GGGGGG | 32992 |
| CGATGT-AGATGT | 22422 |

Looking at files barcodes

```bash
Capture1_N1_R1.fastq.gz
zcat Capture1_N1_R1.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1_R2.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1.reseq.R1.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1.reseq.R2.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
```

After meeting with Jon on 3/2/2023 we have a few ideas of what the issue are regarding these samples.

1. Probe contamination. The elongated barcodes that are visualized in the multiqc report show a GT overthang compared to the other "clean" samples.
2. Adapter ligation. The adapters could have been ligated weirdly. But since we found clean samples in those pools this was discarded.
3. Sample contamination. Since we shared a lane with other samples our could've gotten contaminated. Novogene sent us back the samples demultipelxed and may have mistakenly sent back these reads or there was an issue where they were localized to on the sequencer. 

To move forward in processing the reads we are going to do the following:

1. Concatenate all the demultiplexed samples into one read. 
2. Process_shortread on each of the contaminated reads.
  * Then concatenate them into one file.
3. Explore other contamination programs
4. Normalization

I will keep 3 and 4 on the back burner in order to begin ddocent mapping, snp calling, and coverage.

## Decontaminating our Capture[1-4]_N1 reads

Use double barcode file located here: 

> /home/jgreen/EAGER_OBJ1a/01_process_reads/barcodes/double_barcode.txt

Code for process_shortreads on the Capture[1-4]_N1 and reseq files

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture2_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture3_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture4_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
```

## Fastqc and multiqc demux files

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
fastqc -f fastq -t 20 ${i}/*.fq.gz
done

for i in "${StringArray[@]}"
do
multiqc ${i}/. -o ${i}
done
```

## Cat reseq reads to original reads

```bash
declare -a StringArray=("Capture1_B3" "Capture1_G5" "Capture1_M3" "Capture1_N2" "Capture2_G3" "Capture2_K4" "Capture2_N1" "Capture3_B4" "Capture3_K3" "Capture3_M4" "Capture4_B3" "Capture4_G5" "Capture4_M3" "Capture4_N2" "Capture1_B4" "Capture1_K3" "Capture1_M4" "Capture2_B3" "Capture2_G5" "Capture2_M3" "Capture2_N2" "Capture3_G3" "Capture3_K4" "Capture3_N1" "Capture4_B4" "Capture4_K3" "Capture4_M4" "Capture1_G3" "Capture1_K4" "Capture1_N1" "Capture2_B4" "Capture2_K3" "Capture2_M4" "Capture3_B3" "Capture3_G5" "Capture3_M3" "Capture3_N2" "Capture4_G3" "Capture4_K4" "Capture4_N1")

for i in "${StringArray[@]}"
do
cp ${i}.R1.fastq.gz ../clean/
cp ${i}.R2.fastq.gz ../clean/
cat ${i}.reseq.R1.fastq.gz >> ../clean/${i}.R1.fastq.gz
cat ${i}.reseq.R2.fastq.gz >> ../clean/${i}.R2.fastq.gz
done
```

Need to pull Capture[1-4]_N1 original and reseq reads back into the clean directory

Names for the files: Capture1_N1.R1.fastq.gz Capture1_N1.R2.fastq.gz Capture2_N1.R1.fastq.gz Capture2_N1.R2.fastq.gz Capture3_N1.R1.fastq.gz Capture3_N1.R2.fastq.gz Capture4_N1.R1.fastq.gz Capture4_N1.R2.fastq.gz

```bash
#Capture1_N1
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.R1.fastq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R2.fastq.gz
#Capture1_N1_reseq 
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.R1.fastq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R2.fastq.gz
#Capture2_N1
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.R1.fastq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R2.fastq.gz
#Capture2_N1_reseq
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.R1.fastq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R2.fastq.gz
#Capture3_N1
cat Capture3.1.fq.gz >> ../../../clean/Capture3_N1.R1.fastq.gz
cat Capture3.2.fq.gz >> ../../../clean/Capture3_N1.R2.fastq.gz
#Capture3_N1_reseq
cat Capture3.1.fq.gz >> ../../../clean/Capture3_N1.R1.fastq.gz
cat Capture3.2.fq.gz >> ../../../clean/Capture3_N1.R2.fastq.gz
#Capture4_N1
cat Capture4.1.fq.gz >> ../../../clean/Capture4_N1.R1.fastq.gz
cat Capture4.2.fq.gz >> ../../../clean/Capture4_N1.R2.fastq.gz
#Capture4_N1_reseq
cat Capture4.1.fq.gz >> ../../../clean/Capture4_N1.R1.fastq.gz
cat Capture4.2.fq.gz >> ../../../clean/Capture4_N1.R2.fastq.gz
```

## Analyze clean files

## Run fastqc and multiqc to visualize files
```bash
fastqc -f fastq -t 20 *.fastq.gz
multiqc .
```

## Rename files for ddocent analysis

```bash
rename R1 F *R1.fastq.gz
rename R2 R *R2.fastq.gz
rename fastq fq *.fastq.gz
```
## Run ddocent on all individual files

Make config.file

```bash
nano config.file

Number of Processors
24
Maximum Memory
0
Trimming
yes
Assembly?
no
Type_of_Assembly
PE
Clustering_Similarity%
0.90
Minimum within individual coverage level to include a read for assembly (K1)
2
Minimum number of individuals a read must be present in to include for assembly (K2)
2
Mapping_Reads?
yes
Mapping_Match_Value
1
Mapping_MisMatch_Value
3
Mapping_GapOpen_Penalty
5
Calling_SNPs?
yes
Email
gree9242@uri.edu
```

## Download new dDocent script in scripts directory
```bash
wget "https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/scripts/dDocent_ngs"
```

I receive an error when downloading this as most of the file is in an html format.

```bash 
nano dDocent_ngs.sh
```

Copy and paste the whole script from the previous location

```bash
chmod +x dDocent_ngs.sh
```

## Run the dDocent script in the directory containing our newly named files

```bash
../../scripts/dDocent_ngs.sh config.file
```

## Masking old bed files with haplotig bed file

Copy old bed files

```bash
cp /home/Genomic_Resources/C_virginica/bed_and_GFF/*bed .
```

Pull haplotig bed from github repo https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/Haplotig_Masking/Haplotig_Masking.md

```bash
wget https://raw.githubusercontent.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/main/Haplotig_Masking/Output/haplotigs.bed
```

Use bedtools subtract to remove haplotigs present inthe bed file from the other bed files

```bash
declare -a StringArray=("sorted.ref3.0.CDS.sc" "sorted.ref3.0.exon.sc" "sorted.ref3.0.gene.sc" "sorted.ref3.0.UTR.sc" )
for i in "${StringArray[@]}"
do
bedtools subtract -A -a ${i}.bed -b haplotigs.bed > ${i}.hmask.bed
done
```

```text
sorted.ref3.0.CDS.sc.hmask.bed  sorted.ref3.0.exon.sc.hmask.bed  sorted.ref3.0.gene.sc.hmask.bed  sorted.ref3.0.UTR.sc.hmask.bed
```
# Create reference.file from reference.fasta

```bash
mawk -v OFS='\t' {'print $1,$2'} reference.fasta.fai > reference.file
```

# Get total coverage counts per exon for each capture

sc = single copy gene bed file
hmask = haplotype masked bed file

We will be using the new haplotype masked bed file that Jon has made after the analysis finished up on the C. virginica genome

Capture 1
```bash

declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

#CDS
#Test one one capture
bedtools coverage -hist -b $WORKING_DIR/03_mapping/Capture1_B3.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.CDS.sc.hmask.bed -g $WORKING_DIR/Genome/reference.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/Capture1_B3.hist.AllCDS.all.split.txt
#Looping for all samples
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.CDS.sc.hmask.bed -g $WORKING_DIR/Genome/reference.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/${i}.hist.AllCDS.all.split.txt
done
#Exon
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/reference.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/${i}.hist.AllExon.all.split.txt
done
#Gene
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.gene.sc.hmask.bed -g $WORKING_DIR/Genome/reference.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/${i}.hist.AllGene.all.split.txt
done
#UTR
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.UTR.sc.hmask.bed -g $WORKING_DIR/Genome/reference.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/${i}.hist.AllUTR.all.split.txt
done
#Intergenic
#bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllIntergenic.all.split.txt
```

R plot for one file

```R
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

print(files <- list.files(pattern="Capture1_B3.hist"))
(files <- c("Capture1_B3.hist.AllCDS.all.split.txt", "Capture1_B3.hist.AllExon.all.split.txt", "Capture1_B3.hist.AllGene.all.split.txt", "Capture1_B3.hist.AllUTR.all.split.txt"))
print(files)

labs <- c("CDS","Exon","Gene","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p1 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(1,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2_Capture1_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p1
dev.off()

# Capture 2 B3
print(files <- list.files(pattern="Capture2_B3.hist"))
(files <- c("Capture2_B3.hist.AllCDS.all.split.txt", "Capture2_B3.hist.AllExon.all.split.txt", "Capture2_B3.hist.AllGene.all.split.txt", "Capture2_B3.hist.AllUTR.all.split.txt"))
print(files)

labs <- c("CDS","Exon","Gene","UTR")

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p2 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(1,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2_Capture2_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p2
dev.off()

# Capture 3 B3
print(files <- list.files(pattern="Capture3_B3.hist"))
(files <- c("Capture3_B3.hist.AllCDS.all.split.txt", "Capture3_B3.hist.AllExon.all.split.txt", "Capture3_B3.hist.AllGene.all.split.txt", "Capture3_B3.hist.AllUTR.all.split.txt"))
print(files)

labs <- c("CDS","Exon","Gene","UTR")

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p3 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(1,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2_Capture3_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p3
dev.off()

# Capture 4 B3
print(files <- list.files(pattern="Capture4_B3.hist"))
(files <- c("Capture4_B3.hist.AllCDS.all.split.txt", "Capture4_B3.hist.AllExon.all.split.txt", "Capture4_B3.hist.AllGene.all.split.txt", "Capture4_B3.hist.AllUTR.all.split.txt"))
print(files)

labs <- c("CDS","Exon","Gene","UTR")

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p4 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(1,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2_Capture4_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p4
dev.off()

library(gridExtra)
library(grid)
p5 <- grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
p5
```

These values are a lot lower than when I first did the analysis on just the single copy bed files. Additionally the 0 depth the % bases > 1 is 1. I don't know if this makes sense/ I am also missing a bam file for Capture4_n1.

# Generate data for Table with mapped read stats

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

for i in "${StringArray[@]}"
do 
nom=$(samtools view -@32 $WORKING_DIR/03_mapping/${i}.F.bam -c -L $WORKING_DIR/Genome/mtDNA.bed); denom=$(samtools view -@32 $WORKING_DIR/03_mapping/${i}.F.bam -c); dup=$(mawk '/Unknown/' $WORKING_DIR/03_mapping/${i}_dup_metrics.txt | cut -f9); paste <(echo $i) <(echo $(( `zcat $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz | wc -l` /2 ))) <(echo $(( `zcat $WORKING_DIR/01_process_reads/clean/${i}.R1a.fq.gz | wc -l` /2 ))) <(samtools view -@ 32 $WORKING_DIR/03_mapping/${i}-RGmd.bam -c) <(python -c "print(round("$dup" * 100,2))") <(echo $denom) <(python -c "print(round("$nom"/"$denom" *100,2))") 
done > data.table2
echo -e "Pool\tRaw_Reads\tFiltered_Reads\tMapped_Reads\t%_Duplicate\tFiltered_Mapped_Reads\t%_mapping_to_mitochondrial_genome" > header
cat header data.table2 > table2.txt
```

## Generate data for Figure 3

Calculate Exon percentiles

Get total coverage counts per exon

Set string array
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")
```

Get each capture name with .F.bam using string array
```bash
for i in "${StringArray[@]}"
do
echo $i.F.bam
done
```

Merge capture pools into one file
```bash
#Capture1
samtools merge -@20 capture1_filter.merged.bam Capture1_B3.F.bam Capture1_B4.F.bam Capture1_G3.F.bam Capture1_G5.F.bam Capture1_K3.F.bam Capture1_K4.F.bam Capture1_M3.F.bam Capture1_M4.F.bam Capture1_N2.F.bam

#Capture2
samtools merge -@20 capture2_filter.merged.bam Capture2_B3.F.bam Capture2_B4.F.bam Capture2_G3.F.bam Capture2_G5.F.bam Capture2_K3.F.bam Capture2_K4.F.bam Capture2_M3.F.bam Capture2_M4.F.bam Capture2_N1.F.bam Capture2_N2.F.bam

#Capture3
samtools merge -@20 capture3_filter.merged.bam Capture3_B3.F.bam Capture3_B4.F.bam Capture3_G3.F.bam Capture3_G5.F.bam Capture3_K3.F.bam Capture3_K4.F.bam Capture3_M3.F.bam Capture3_M4.F.bam Capture3_N1.F.bam Capture3_N2.F.bam

#Capture4
samtools merge -@20 capture3_filter.merged.bam Capture4_B3.F.bam Capture4_B4.F.bam Capture4_G3.F.bam Capture4_G5.F.bam Capture4_K3.F.bam Capture4_K4.F.bam Capture4_M3.F.bam Capture4_M4.F.bam Capture4_N1.F.bam Capture4_N2.F.bam
```

Get total coverage counts per exon per sample
```bash
for i in "${StringArray[@]}"
do
bedtools coverage -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hamsk.bed -g $WORKING_DIR/Genome/reference.file -mean -split > ${i}cov.mean.filtered.exon.stats
done
```

Get total coverage counts per exon
```bash
declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
bedtools coverage -b _${i}filter.merged.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.bed -sorted -g $WORKING_DIR/Genome/reference.file -mean > ${i}.cov.mean.filtered.merged.exon.stats
done
```

Goal is to compare capture1 against all others

Remove mtDNA from stat files

```bash
declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
mawk '!/NC_007175.2/' ${i}.cov.mean.filtered.merged.exon.stats > ${i}.DNA.mean.stats
```

Calculate lower 10th percentile of exon sizes
```bash
declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}cov.mean.filtered.merged.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
done
```

Calculate upper 10th percentile of exon sizes
```bash
declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}cov.mean.filtered.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
done
```

Mark exons into size classes based on size distribution and create data table
```bash
declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
mawk '{if ( $3 -$2 > 517 ) print $0 "\tUpper"; else if ( $3 - $2 < 59 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 518) print $0 "\tMiddle" }' ${i}.cov.mean.filtered.C1.exon.stats > ${i}.mean.cov.stats.class
done

echo -e "Chrom\tStart\tEnd\tDNA_Coverage\tExon_Size_Class" > header

declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
cat header ${i}.mean.cov.stats.class > ${i}.ExonMeanCoverage.txt
```

Next step is to find exons with minimum thresholds of gDNA.  These will be our "target" sets along with confidence intervals.  Based on overal DNA coverage, we chose 35X as our "target" set and choose 15X boundaries around that number for confidence intervals.  We will create three `bed` files from our RNA exon coverage stats.

```bash

declare -a StringArray=("capture1" "capture2" "capture3" "capture4")
for i in  "${StringArray[@]}"
do
mawk 'BEGIN { FS = "\t" } ; $11 > 19' ${i}.cov.mean.filtered.exon.stats > .EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 49' ${i}.cov.mean.filtered.exon.stats > .EiRc50.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 34' ${i}.cov.mean.filtered.exon.stats > .EiRc35.bed
done
```

### Calculating data for table 3

We will use a BASH function to automate this for us:

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $1.F.bam -a sorted.ref3.0.exon.sc.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc20.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc50.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc35.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat sorted.ref3.0.exon.sc.bed | wc -l )
X35X=$(cat C1.EiRc20.bed | wc -l)
X20X=$(cat C1.EiRc50.bed | wc -l)
X50X=$(cat C1.EiRc35.bed | wc -l)

#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target C1_Sample1) <(counts_per_target C2_Sample1) <(counts_per_target C3_Sample1) <(counts_per_target C4_Sample1) > Table3.txt
```

## Graveyard

# Used for orginally renaming files for demultiplexing

Rename fq suffix to fastq
```bash
for file in *.fq.gz; do mv -- "$file" "${file%fq.gz}fastq.gz"; done
```

Rename R1 files to R2 to reflect reverse read naming convention
```bash
for file in *R1.fastq.gz; do mv -- "$file" "${file%R1.fastq.gz}R2.fastq.gz"; done
for file in *R1_fastqc.zip; do mv -- "$file" "${file%R1_fastqc.zip}R2_fastqc.zip"; done
for file in *R1_fastqc.html; do mv -- "$file" "${file%R1_fastqc.html}R2_fastqc.html"; done
```

Rename F1 and F2 files to R1 to reflect forward read naming convention

```bash
# For F1 file -> R1 file naming
for file in *F1.fastq.gz; do mv -- "$file" "${file%F1.fastq.gz}R1.fastq.gz"; done
for file in *F1_fastqc.zip; do mv -- "$file" "${file%F1_fastqc.zip}R1_fastqc.zip"; done
for file in *F1_fastqc.html; do mv -- "$file" "${file%F1_fastqc.html}R1_fastqc.html"; done

# For F2 file -> R1 file naming
for file in *F2.fastq.gz; do mv -- "$file" "${file%F2.fastq.gz}R1.fastq.gz"; done
for file in *F2_fastqc.zip; do mv -- "$file" "${file%F2_fastqc.zip}R1_fastqc.zip"; done
for file in *F2_fastqc.html; do mv -- "$file" "${file%F2_fastqc.html}R1_fastqc.html"; done

```


Making directories for all original samples
```bash
#declare -a StringArray=("Capture1_B3" "Capture1_G5" "Capture1_M3" "Capture1_N2" "Capture2_G3" "Capture2_K4" "Capture2_N1" "Capture3_B4" "Capture3_K3" "Capture3_M4" "Capture4_B3" "Capture4_G5" "Capture4_M3" "Capture4_N2" "Capture1_B4" "Capture1_K3" "Capture1_M4" "Capture2_B3" "Capture2_G5" "Capture2_M3" "Capture2_N2" "Capture3_G3" "Capture3_K4" "Capture3_N1" "Capture4_B4" "Capture4_K3" "Capture4_M4" "Capture1_G3" "Capture1_K4" "Capture1_N1" "Capture2_B4" "Capture2_K3" "Capture2_M4" "Capture3_B3" "Capture3_G5" "Capture3_M3" "Capture3_N2" "Capture4_G3" "Capture4_K4" "Capture4_N1")

#for i in "${StringArray[@]}"; do; mkdir ${i}; done
```

## Continuing with global capture files

For right now I am going to keep the capture files as one files. Meaning all individuals will be group into one category based on their size selection of the probe and gDNA fragment length. The Multiqc shows great stats across the board and I will be able to do an initial pass for coverage differences across all the treatment irrespective of individual differences or comparisons. 

## Going to concatenate all the files from each capture into one file
```bash
for i in 1:4
do
cat ../raw/Capture${i}_*.F*.fq.gz >> capture${i}_F.fq.gz
done

for i in 1:4
do
cat ../raw/Capture${i}_*.R*.fq.gz >> capture${i}_R.fq.gz
done

## Run Fastqq and multiqc
fastqc -f fastq -t 20 *.fq.gz
multiqc .
```

## Rename files

```bash
mv capture1_F.fq.gz C1_Sample1.F.fq.gz
mv capture1_R.fq.gz C1_Sample1.R.fq.gz
mv capture2_F.fq.gz C2_Sample1.F.fq.gz
mv capture2_R.fq.gz C2_Sample1.R.fq.gz
mv capture3_F.fq.gz C3_Sample1.F.fq.gz
mv capture3_R.fq.gz C3_Sample1.R.fq.gz
mv capture4_F.fq.gz C4_Sample1.F.fq.gz
mv capture4_R.fq.gz C4_Sample1.R.fq.gz
```

```bash
#CDS
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.CDS.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.exon.sc.hmask.bed -g reference.file -sorted -split | grep ^all > C1.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.gene.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.UTR.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllUTR.all.split.txt
#Intergenic
#bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllIntergenic.all.split.txt
```

Capture 2
```bash
#CDS
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C2.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllIntergenic.all.split.txt
```

Capture 3
```bash
#CDS
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C3.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllIntergenic.all.split.txt
```

Capture 4
```bash
#CDS
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C4.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllIntergenic.all.split.txt
```

## Generate data for Figure 3

Calculate Exon percentiles

Get total coverage counts per exon

```bash
#Capture 1
bedtools coverage -b C1_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean -split > cov.mean.filtered.C1.exon.stats
#Capture 2
bedtools coverage -b C2_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean --split> cov.mean.filtered.C2.exon.stats
#Capture 3
bedtools coverage -b C3_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean > cov.mean.filtered.C3.exon.stats
#Capture 4
bedtools coverage -b C4_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean > cov.mean.filtered.C4.exon.stats
```

```bash
samtools merge -@20 filter.merged.bam C1_Sample1.F.bam C2_Sample1.F.bam C3_Sample1.F.bam C4_Sample1.F.bam
```

Get total coverage counts per exon
```bash
bedtools coverage -b filter.merged.bam -a sorted.ref3.0.exon.sc.bed -sorted -g reference.file -mean > capture.cov.mean.filtered.merged.exon.stats
```

Not sure what to do with these stats files since we dont have RNA coverage. Think this through. Maybe use percentiles as x - axis

Remove mtDNA from stat files

```bash
mawk '!/NC_007175.2/' capture.cov.mean.filtered.merged.exon.stats > capture.DNA.mean.stats
```

Calculate lower 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' cov.mean.filtered.C1.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
```
Resut: `59`

Calculate upper 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' cov.mean.filtered.C1.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
```
Result: `517`

Mark exons into size classes based on size distribution and create data table
```bash
mawk '{if ( $3 -$2 > 517 ) print $0 "\tUpper"; else if ( $3 - $2 < 59 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 518) print $0 "\tMiddle" }' cov.mean.filtered.C1.exon.stats > C1.rnd.mean.cov.stats.class
echo -e "Chrom\tStart\tEnd\tDNA_Coverage\tExon_Size_Class" > header
cat header C1.rnd.mean.cov.stats.class > C1.ExonMeanCoverage.txt
```

Next step is to find exons with minimum thresholds of gDNA.  These will be our "target" sets along with confidence intervals.  Based on overal DNA coverage, we chose 35X as our "target" set and choose 15X boundaries around that number for confidence intervals.  We will create three `bed` files from our RNA exon coverage stats.

```bash
mawk 'BEGIN { FS = "\t" } ; $11 > 19' cov.mean.filtered.C1.exon.stats > C1.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 49' cov.mean.filtered.C1.exon.stats > C1.EiRc50.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 34' cov.mean.filtered.C1.exon.stats > C1.EiRc35.bed
```
### Calculating data for table 3

We will use a BASH function to automate this for us:

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $1.F.bam -a sorted.ref3.0.exon.sc.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc20.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc50.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc35.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat sorted.ref3.0.exon.sc.bed | wc -l )
X35X=$(cat C1.EiRc20.bed | wc -l)
X20X=$(cat C1.EiRc50.bed | wc -l)
X50X=$(cat C1.EiRc35.bed | wc -l)

#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target C1_Sample1) <(counts_per_target C2_Sample1) <(counts_per_target C3_Sample1) <(counts_per_target C4_Sample1) > Table3.txt
```