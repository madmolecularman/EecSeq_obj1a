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
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_1.fq.gz Capture1_B3.F.fq.gz
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_2.fq.gz Capture1_B3.R.fq.gz
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
Multiqc plots show uneven distribution of sequence counts specifically in samples Capture1_N1.R1, Capture1_N1.R2, Capture2_N1.R1, Capture2_N1.R2, Capture3_N1.reseq.R1, Capture3_N1.reseq.R2, Capture4_N1.reseq.R1, Capture4_N1.reseq.R2. See [multiqc report]()

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

New directory for single inline barcodes. I am wanting to test out how we retain reads with single barcodes.

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
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r -c -q
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

## Run the dDoccent script
```bash
## Make a config.file

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

## Download new dDocent script
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

## Run the dDocent script in the directory containing our newly names files

```bash
./dDocent_ngs.sh congfig.file
```

# Get total coverage counts per exon for each capture

sc = single copy

Capture 1
```bash
#CDS
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C1.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllIntergenic.all.split.txt
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