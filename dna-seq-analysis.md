---
description: script for 202104 DNA-seq data analysis
---

# DNA-seq analysis

## 1 Check raw data

* md5sum

```bash
// md5sum -c md5.md5
```

## 2 Qc

```bash
// fastqc /home/jiangb/rna_seq_202103/raw/N2101557_LL_80-600237629_eukRNASEQ/combined/*.gz -t 6 -o qc
// multiqc .
```

## 3 Trim #5' 10bp

Based on qc report, we find that sequencing quality of 5' 10bp is not stable, then we cut it down

```bash
// #!/bin/bash

echo -e "\033[1;36mTrim finished at $(date)\033[0m"

for filename in *.gz;
        do
        # Get export file name
        exportFilename=`echo "${filename}" | sed -E 's/-LDD[0-9]{4}_combined//g' | sed -E 's/(R[0-9]).fastq/\1_trimmed10.fastq/g'`
        # Specify path to export
        exportPath="/home/jiangb/dna_seq_202104/trimmed/$exportFilename"
        echo -e "Reading \033[1;36m${filename}\033[0m; Export to \033[1;31m${exportPath}\033[0m"
        # Run commands
        time zcat $filename | fastx_trimmer -f 11 -z > "$exportPath" && echo "**trim done**"
done

echo -e "\033[1;36mTrim finished at $(date)\033[0m"
```

## 4 Clean reads

filter read<100 qc<30 and remove adaptors

```bash
// #!/bin/bash

# get all filename in specified path
path="/home/jiangb/dna_seq_202104/raw/N2104223_LL_80-637653917_DNASEQ/210416-Nova"
files=`ls ${path}/*.gz` 
for filename in $files;
do
	basename "${filename}" | sed -E 's/-LDD[0-9]{4}_combined_(R[0-9]).fastq.gz//g' >> filename.txt
done

#remove duplicate filename
sort -n filename.txt | uniq >> uniq_filename.txt
rm filename.txt


#filter bad quality reads and remove adaptors

#mkdir clean

echo -e "\033[1;36mTrim_galore started at $(date)\033[0m"

cat uniq_filename.txt | while read i;
do
	echo -e "\033[1;36mNow for $i\033[0m"
	echo $i | time trim_galore -output_dir clean --paired --fastqc --length 100 --quality 30 --gzip --stringency 5 \
	 trimmed/${i}_R1_trimmed10.fastq.gz trimmed/${i}_R2_trimmed10.fastq.gz && echo "**clean done**"
done

echo -e "\033[1;36mTrim_galore finished at $(date)\033[0m"

#after clean, qc again
multiqc . -o qc
```

## 5 Mapping

```bash
// 
#reference index
bwa index reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
#.fasta.fai
samtools faidx reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
#dict for reference
picard CreateSequenceDictionary R=reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa  O=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.dict

#mapping
#!/bin/bash
echo -e "\033[1;36mbwa star at $(date)\033[0m"

mkdir mapping/sorted
cat uniq_filename.txt | while read i;
do
        tag="@RG\tID:${i}\tSM:${i}\tPL:Illumina\tLB:library1"
        echo $i |time bwa mem -t 20 -M  -R "${tag}" reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa clean/${i}_R1_trimmed10_val_1.fq.gz clean/${i}_R2_trimmed10_val_2.fq.gz |samtools view -bS - > mapping/${i}.bam && echo "**bwa mapping done **"
        time samtools sort -@6 -o mapping/sorted/${i}.bam mapping/${i}.bam && echo "**samtools sorting done **"
        rm mapping/${i}.bam && echo "**raw bam removing done **"
done

echo -e "\033[1;36mbwa mapping finish at $(date)\033[0m"
```

## 6 Mark duplicate and sort

```bash
//#!/bin/bash
  
echo -e "\033[1;36m started at $(date)\033[0m"

mkdir 2_variant_calling/markdup_sort
mkdir 2_variant_calling/markdup_sort/markdup_matrixs

cat uniq_filename.txt | while read i;
do
        echo -e "\033[1;36mNow for $i\033[0m"
        echo $i |time gatk-4.2.0.0/gatk MarkDuplicates I=1_bwa_mapping/mapping/sorted/${i}.bam O=2_variant_calling/markdup_sort/${i}.sorted.markdup.bam M=2_variant_calling/markdup_sort/markdup_matrixs/${i}.sorted.markdup.txt REMOVE_DUPLICATES=false && echo "**markduplicates done**"
        #make index for bam
        echo $i |time samtools index 2_variant_calling/markdup_sort/${i}.sorted.markdup.bam && echo "**index done**"

done

echo -e "\033[1;36m finished at $(date)\033[0m"
```

## 7 Gatk variantion calling

```bash
/#!/bin/bash

echo -e "\033[1;36mvariant calling started at $(date)\033[0m"

cat name0.txt | while read i;
do
       echo -e "\033[1;36mNow for $i\033[0m"
       echo $i| time gatk-4.2.0.0/gatk  HaplotypeCaller \
                -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
                -ERC GVCF \
                -ploidy 1 \
                -I mapping/${i}.sorted.markdup.bam \
                -O variant_calling/gvcf/${i}.g.vcf      && echo "**get gvcf done**"

done



find  variant_calling/gvcf/ -name "*.g.vcf" >input.list

#Consolidate GVCFs `CombineGVCFs`
time gatk-4.2.0.0/gatk CombineGVCFs \
        -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -V input.list \
        -O variant_calling/combined.g.vcf && echo "**merge gvcf done**"


#variant_calling by joint genotyping
time gatk-4.2.0.0/gatk GenotypeGVCFs \
        -ploidy 1 \
        -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -V variant_calling/combined.g.vcf \
        -O variant_calling/variant.vcf && echo "**variant calling done**"

#select snp
time gatk-4.2.0.0/gatk SelectVariants\
        -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -V variant_calling/variant.vcf \
        -select-type SNP \
        -O variant_calling/snp.vcf && echo "**snp select done**"

#select indel
time gatk-4.2.0.0/gatk SelectVariants\
        -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -V variant_calling/variant.vcf \
        -select-type INDEL \
        -O variant_calling/indel.vcf && echo "**indel select done**"

#hard filter snp
time gatk-4.2.0.0/gatk VariantFiltration \
        -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -V variant_calling/snp.vcf \
        --cluster-window-size 10  \
        --cluster-size 3 \
        --missing-values-evaluate-as-failing \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
         -O variant_calling/snp.hardfiltered.vcf && echo "**snp filter done**"

#select PASS snp
#pass snp including those without two info MQRankSum and ReadPosRankSum
awk '/^#/||$7=="PASS"||$7=="MQRankSum-12.5;ReadPosRankSum-8"||$7=="MQRankSum-12.5"||$7=="ReadPosRankSum-8"' variant_calling/snp.hardfiltered.vcf  > variant_calling/snp.hardfiltered.PASS.vcf


echo -e "\033[1;36mgvcf calling finished at $(date)\033[0m"


#filter indel
#!/bin/bash
  
#hard filter indel
time gatk-4.2.0.0/gatk VariantFiltration \
                -R reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
                -V variant_calling/indel.vcf \
                -filter "QD < 2.0" --filter-name "QD2" \
                -filter "FS > 200.0" --filter-name "FS200" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 10.0" --filter-name "SOR10" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                -O variant_calling/indel.hardfiltered.vcf && echo "**indel filter done**"

#select PASS indel
#pass indel including those without two info MQRankSum and ReadPosRankSum
awk '/^#/||$7=="PASS"||$7=="MQRankSum-12.5;ReadPosRankSum-8"||$7=="MQRankSum-12.5"||$7=="ReadPosRankSum-8"' variant_calling/indel.hardfiltered.vcf  > variant_calling/indel.hardfiltered.PASS.vcf


echo -e "\033[1;36mindel filter finished at $(date)\033[0m"
```

## 8 snpEff annotate vcf

```bash
#gft2gff3
gffread reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.101.gtf -o-> reference/ensembl/gff3fromgtfbygffread/Saccharomyces_cerevisiae.R64-1-1.101.gff3

#snpEff annotate vcf
 unzip snpEff_latest_core.zip 
 cd snpEff/
 mkdir data
 cd data
 mkdir genomes ;copy fa to this file
 mkdir sc ;copy gff3 to this file
 cd snpEff/
 echo "sc.genome:sc" >> snpEff.config 
 cd ..
 java -jar snpEff/snpEff.jar build -gff3 sc


#annotate pass_snp and pass_indel
java -jar snpEff/snpEff.jar -v sc 2_variant_calling/snp.hardfiltered.PASS.vcf > 2_variant_calling/annotate_by_snpEff/snp.hardfiltered.PASS.anno.vcf 

java -jar snpEff/snpEff.jar -v sc 2_variant_calling/indel.hardfiltered.PASS.vcf >2_variant_calling/annotate_by_snpEff/indel.hardfiltered.PASS.anno.vcf
```

## 9 Get readcount table

```bash
##get readcount table of sna-seq samples for checking of markers genes
featureCounts -T 10 -p -t exon -g gene_id  -a reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.101.gtf -o readcount/dna_readcount.txt \
mapping/S0d-ADE1-1.sorted.bam \
mapping/S0d-ADE1-2.sorted.bam \
mapping/S0d-ADE1-6.sorted.bam \
mapping/S0d-HAP4-1.sorted.bam \
mapping/S0d-HAP4-2.sorted.bam \
mapping/S0d-HAP4-3.sorted.bam \
mapping/S0d-WT.sorted.bam \
mapping/14d-ADE1-1.sorted.bam \
mapping/14d-ADE1-2.sorted.bam \
mapping/14d-ADE1-6.sorted.bam \
mapping/14d-HAP4-1.sorted.bam \
mapping/14d-HAP4-2.sorted.bam \
mapping/14d-HAP4-3.sorted.bam \
mapping/14d-WT1.sorted.bam \
mapping/14d-WT2.sorted.bam \
mapping/14d-WT3.sorted.bam \
mapping/28d-ADE1-1.sorted.bam \
mapping/28d-ADE1-2.sorted.bam \
mapping/28d-ADE1-6.sorted.bam \
mapping/28d-HAP4-1.sorted.bam \
mapping/28d-HAP4-2.sorted.bam \
mapping/28d-HAP4-3.sorted.bam \
mapping/28d-WT1.sorted.bam \
mapping/28d-WT2.sorted.bam \
mapping/28d-WT3.sorted.bam >featurecount_dna.log 2>&1 &
```

## 10 CNVkit

* install

```bash
git clone https://github.com/etal/cnvkit
cd cnvkit/
pip install -e .

pip install matplotlib -i https://pypi.douban.com/simple

conda activate cnv
```

* usage

```bash
python3 cnvkit/cnvkit.py batch \
-m wgs \
--haploid-x-reference \
-p 20 \
--output-dir 3_cnv \
-n 1_bwa_mapping/mapping/sorted/0d-WT.bam \
-f  reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
--annotate reference/ensembl/Saccharomyces_cerevisiae.R64-1-1.101.gtf \
1_bwa_mapping/mapping/sorted/0d-ADE1-1.bam \
1_bwa_mapping/mapping/sorted/0d-ADE1-2.bam \
1_bwa_mapping/mapping/sorted/0d-ADE1-6.bam \
1_bwa_mapping/mapping/sorted/0d-HAP4-1.bam \
1_bwa_mapping/mapping/sorted/0d-HAP4-2.bam \
1_bwa_mapping/mapping/sorted/0d-HAP4-3.bam \
1_bwa_mapping/mapping/sorted/14d-ADE1-1.bam \
1_bwa_mapping/mapping/sorted/14d-ADE1-2.bam \
1_bwa_mapping/mapping/sorted/14d-ADE1-6.bam \
1_bwa_mapping/mapping/sorted/14d-HAP4-1.bam \
1_bwa_mapping/mapping/sorted/14d-HAP4-2.bam \
1_bwa_mapping/mapping/sorted/14d-HAP4-3.bam \
1_bwa_mapping/mapping/sorted/14d-WT1.bam \
1_bwa_mapping/mapping/sorted/14d-WT2.bam \
1_bwa_mapping/mapping/sorted/14d-WT3.bam \
1_bwa_mapping/mapping/sorted/28d-ADE1-1.bam \
1_bwa_mapping/mapping/sorted/28d-ADE1-2.bam \
1_bwa_mapping/mapping/sorted/28d-ADE1-6.bam \
1_bwa_mapping/mapping/sorted/28d-HAP4-1.bam \
1_bwa_mapping/mapping/sorted/28d-HAP4-2.bam \
1_bwa_mapping/mapping/sorted/28d-HAP4-3.bam \
1_bwa_mapping/mapping/sorted/28d-WT1.bam \
1_bwa_mapping/mapping/sorted/28d-WT2.bam \
1_bwa_mapping/mapping/sorted/28d-WT3.bam>log/cnv.log 2>&1 &

 #-y Use or assume a male reference --male-reference, --haploid-x-reference
```

* plot

```
cnvkit.py scatter -h
cnvkit.py diagram -h
cnvkit.py heatmap -h
```

```bash
###overview
#cnv overview by bin-level log2 coverages or copy ratios (.cnn, .cnr)
#it would take a much long time 
python3 cnvkit/cnvkit.py heatmap 3_cnv/*.cnr \
--desaturate \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cnr_d.pdf

#cnv overview by segements which is by far more clear
python3 cnvkit/cnvkit.py heatmap 3_cnv/cns/*.cns \
--desaturate \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cns_d.pdf
```

```bash
###cnv specific chr overview
#chr IV
python3 cnvkit/cnvkit.py heatmap 3_cnv/cns/*.cns \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cns_IV_d.pdf \
-c IV \
--desaturate 
#A samller window
python3 cnvkit/cnvkit.py heatmap 3_cnv/cns/*.cns \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cns_IV_500000-580000.pdf \
-c IV:500000-580000 \
--desaturate 

#chr XII
python3 cnvkit/cnvkit.py heatmap 3_cnv/cns/*.cns \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cns_XII_d.pdf \
-c XII \
--desaturate 
#A samller window
python3 cnvkit/cnvkit.py heatmap 3_cnv/cns/*.cns \
--haploid-x-reference \
-o 3_cnv/plot/pheatmap_cns_XII_460000-500000.pdf \
-c XII:460000-500000 \
--desaturate 


#scatter chr XII
#plotted with CNVkit’s “scatter” command, 
#the size of the plotted datapoints is proportional to each bin’s weight 
#a relatively small point indicates a less reliable bin.
python3 cnvkit/cnvkit.py scatter \
-s 3_cnv/cns/0d-ADE1-6.cns \
3_cnv/0d-ADE1-6.cnr \
--title 0d-ADE1-6_XII \
-c XII \
-t \
-o 3_cnv/plot/0d-ADE1-6_XII_scatter.pdf

python3 cnvkit/cnvkit.py scatter \
-s 3_cnv/cns/14d-ADE1-6.cns \
3_cnv/14d-ADE1-6.cnr \
--title 14d-ADE1-6_XII \
-c XII \
-t \
-o 3_cnv/plot/14d-ADE1-6_XII_scatter.pdf

python3 cnvkit/cnvkit.py scatter \
-s 3_cnv/cns/28d-ADE1-6.cns \
3_cnv/28d-ADE1-6.cnr \
--title 28d-ADE1-6_XII \
-c XII \
-t \
-o 3_cnv/plot/28d-ADE1-6_XII_scatter.pdf



##plot raw log2 depth 
python3 cnvkit/cnvkit.py scatter \
3_cnv/0d-WT.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/0d-WT_log2_depth.pdf \
--title depth_0d-WT:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 

#0d-ADE1-6
python3 cnvkit/cnvkit.py scatter \
3_cnv/0d-ADE1-6.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/0d-ADE1-6_log2_depth.pdf \
--title depth_0d-ADE1-6:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 

#14d-ADE1-6
python3 cnvkit/cnvkit.py scatter \
3_cnv/14d-ADE1-6.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/14d-ADE1-6_log2_depth.pdf \
--title depth_14d-ADE1-6:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 
#28d-ADE1-6
python3 cnvkit/cnvkit.py scatter \
3_cnv/28d-ADE1-6.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/28d-ADE1-6_log2_depth.pdf \
--title depth_28d-ADE1-6:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 

#0d-ADE1-2
python3 cnvkit/cnvkit.py scatter \
3_cnv/0d-ADE1-2.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/0d-ADE1-2_log2_depth.pdf \
--title depth_0d-ADE1-2:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 
#14d-HAP4-2
python3 cnvkit/cnvkit.py scatter \
3_cnv/14d-HAP4-2.targetcoverage.cnn \
-c XII:460000-500000 \
-o 3_cnv/plot/14d-HAP4-2_log2_depth.pdf \
--title depth_14d-HAP4-2:XII_460000-500000 \
-g YLR155C,YLR156W,RDN5-4,YLR159W \
-t 
```
