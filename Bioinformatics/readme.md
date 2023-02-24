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
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L4_2.fq.gz ${i}reseq.R3.fastq.gz
done
```

## Run fastqc and multiqc to visualize files
```bash
fastqc -f fastq -t 20 *.fastq.gz
multiqc .
```

## Results
Multiqc plots show uneven distribution of sequence counts specifically in samples Capture1_N1.F1, Capture1_N1.R1, Capture2_N1.F2, Capture2_N1.R2, Capture3_N1.F3, Capture3_N1.R3, Capture4_N1.F4, Capture4_N1.R4


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

## Renaming files for demultiplexing

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

## Demultiplexing capture files

Use process_shortreads to demultiplex files. In the DNA/demux folder. One of the major problems is that we are not sure if novogene was able to separate out the different capture pools. We need to create a demultiplexing run that will filter the index and the inline barcode into the proper file. The question here is can I demultiplex just the index and then the adapter inline barcode?

Its important to understand what ![process_shortreads](https://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php) is doing and how to create the barcode files ![Stacks manual](https://catchenlab.life.illinois.edu/stacks/manual/index.php#clean)

First I have to create a barcode file that is tab delimited

```bash
nano barcodes.tsv
```

```text
ATCACG	ATCACG  Capture1
CGATGT	CGATGT  Capture2
TTAGGC	TTAGGC  Capture3
TGGCCA	TGGCCA  Capture4
```

Make directories for all of the files

```bash
mkdir samples
cd samples

declare -a StringArray=("Capture1_B3" "Capture1_G5" "Capture1_M3" "Capture1_N2" "Capture2_G3" "Capture2_K4" "Capture2_N1" "Capture3_B4" "Capture3_K3" "Capture3_M4" "Capture4_B3" "Capture4_G5" "Capture4_M3" "Capture4_N2" "Capture1_B4" "Capture1_K3" "Capture1_M4" "Capture2_B3" "Capture2_G5" "Capture2_M3" "Capture2_N2" "Capture3_G3" "Capture3_K4" "Capture3_N1" "Capture4_B4" "Capture4_K3" "Capture4_M4" "Capture1_G3" "Capture1_K4" "Capture1_N1" "Capture2_B4" "Capture2_K3" "Capture2_M4" "Capture3_B3" "Capture3_G5" "Capture3_M3" "Capture3_N2" "Capture4_G3" "Capture4_K4" "Capture4_N1")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Try it on one file (ran from DNA directory)

```bash
process_shortreads -P -1 ./raw/Capture1_N1.R1.fastq.gz -2 ./raw/Capture1_N1.R2.fastq.gz -o ./samples/Capture1_B3/ -b ./barcodes/barcodes.tsv --inline_inline -c -q -r --barcode_dist_2 3
```

This result is in ./samples/Capture1_N1/process_shortreads.log. There are some weird results including a mixture of the different barcodes.

<center>

|Sequence Category|Sequence Number|
|-----------------|---------------|
|Total Sequences|	35392266|
|Ambiguous Barcodes|	35309540|
|Low Quality|	273|
|Trimmed Reads|	0|
|Orphaned Paired-ends|	0|
|Retained Reads|	82453|


|Barcode|	Total|	Retained|
|-------|--------|----------|
|ATCACG-CGTGAT|	8128|	8112|
|CGATGT-ACATCG|	52322|	52119|
|TTAGGC-GCCTAA|	17622|	17576|
|TGGCCA-TGGCCA|	4654|	4646|

</center>

### Sequences not recorded (subset)

<center>

|Barcode|	Total|
|-------|--------|
|CGATGT-CGATGT|	13471732|
|CGATGT-GGGGGG|	34484|
|CGATGT-CGAGGT|	10924|

</center>


```bash
# Using process shortreads on a whole directory containing these files
# Cannot use yet
process_shortreads -P -p ./raw/ -o ./samples/ -b ./barcodes/barcodes.tsv --inline_inline -c -q -r --barcode_dist 3
```

## Continuing with global capture files

For right now I am going to keep the capture files as one files. Meaning all individuals will be group into one category based on their size selection of the probe and gDNA fragment length. The Multiqc shows great stats across the board and I will be able to do an initial pass for coverage differences across all the treatment irrespective of individual differences or comparisons. 

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