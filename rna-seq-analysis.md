# RNA-seq analysis

## 1 Qc

```
echo "fastqc started at $(date)"
fastqc *.gz -t <thred_num> -o <out.dir> -f <input_format> <input_file_1> <input_file_2> ...
fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN
```

对二代测序数据质量快速检验的工具，可以输入fastq（fastq.gz）、sam或者bam文件。查看输出结果解释。可以联合multiqc使用，查看多个qc的报告。

主要是包括前面的各种选项和最后面的可以加入N个文件

\-o --outdir FastQC生成的报告文件的储存路径，生成的报告的文件名是根据输入来定的

\--extract 生成的报告默认会打包成1个压缩文件，使用这个参数是让程序不打包

\-t --threads 选择程序运行的线程数，每个线程会占用250MB内存，越多越快咯

\-c --contaminants 污染物选项，输入的是一个文件，格式是Name \[Tab] Sequence，里面是可能的污染序列，如果有这个选项，FastQC会在计算时候评估污染的情况，并在统计的时候进行分析，一般用不到

\-a --adapters 也是输入一个文件，文件的格式Name \[Tab] Sequence，储存的是测序的adpater序列信息，如果不输入，目前版本的FastQC就按照通用引物来评估序列时候有adapter的残留

\-q --quiet 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况。

```
#merge qc reports
multiqc . #对当前文件夹
multqc -o <out.path> *.fastqc.zip --ignore *.html #指定文件夹
echo "fastqc finished at $(date)"
```

## 2 Read a README

## 3 Trim 5' \~15bp

```
mkdir trimmed #建文件夹
#!/bin/bash
for filename in *.gz; do
    # Get export file name
    exportFilename=`echo "${filename}" | sed -E 's/S-([^-]+).+combined_(R[0-9])/\1_\2_trimmed15/'`
    
    # Specify path to export
    exportPath="/home/jiangb/trimmed/$exportFilename"
    echo -e "Reading \033[1;36m${filename}\033[0m; Export to \033[1;31m${exportPath}\033[0m"

    # Run commands
    zcat $filename | fastx_trimmer -f 16 -z > "$exportPath"
done
```

## 4 Clean

\#使用trim\_galore前需要确认需先安装fastqc和cutadapt是否已经安装

```
#remove bad quality reads and remove adaptors
mkdir clean

#!/bin/bash
echo " trim_galore cut adapters started at $(date)"
for i in B1 B2 B3 0H2 0H3 0H4 14H2 14H3 14H4 ; do
	trim_galore -output_dir clean --paired --fastqc --length 100 --quality 30 --gzip --stringency 5 \
	trimmed/${i}_R1_trimmed15.fastq.gz trimmed/${i}_R2_trimmed15.fastq.gz
done
echo "trim_galore cut adapters finished at $(date)"
```

{% embed url="https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md" %}

\--quality:设定Phred quality score阈值，默认为20

\--phred33：：选择-phred33或者-phred64，表示测序平台使用的Phred quality score。

\--adapter：输入adapter序列。也可以不输入，Trim Galore!会自动寻找可能性最高的平台对应的adapter。自动搜选的平台三个，也直接显式输入这三种平台，即--illumina、--nextera和--small\_rna。

\--stringency：设定可以忍受的前后adapter重叠的碱基数，默认为1（非常苛刻）。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到。

\--length：设定输出reads长度阈值，小于设定值会被抛弃。

\--paired：对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个会被同样抛弃，而不管是否达到标准。

\--retain\_unpaired：对于双端测序结果，一对reads中，如果一个read达到标准，但是对应的另一个要被抛弃，达到标准的read会被单独保存为一个文件。

\--gzip和--dont\_gzip：清洗后的数据zip打包或者不打包。

\--output\_dir：输入目录。需要提前建立目录，否则运行会报错。

\--trim-n : 移除read一端的reads

\--fastqc :完成质控后进行fastqc

\--e ：Maximum allowed error rate (no. of errors divided by the length of the matching region) Default: 0.1

```
#再用fastqc和multiqc观察以下过滤的结果
fastqc *.gz #这一步可以省略因为trim_galore可以直接进行fastqc
multiqc .
```

## 5 Download reference.fa & reference.gtf

## 6 STAR 2-pass mode mapping

```
#!/bin/bash
####mapping-STAR_2_pass_mode

mkdir mapping
mkdir 2-pass_mapping
mkdir 2-pass_index

#mapping，output bam
echo -e "\033[1;36mSTAR mapping started at $(date)\033[0m"

mkdir index_S288c_gtf

#index 1st：
time STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir index_S288c_gtf \
--genomeFastaFiles reference/latest/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
--sjdbGTFfile reference/latest/Saccharomyces_cerevisiae.R64-1-1.104.gtf \
--sjdbOverhang 134 \
--genomeSAindexNbases 10 && echo "**STAR index 1st done**"

#mapping 1st:
cat uniq_filename1.txt | while read i;
do
       echo -e "\033[1;36mNow for $i\033[0m"    
        time STAR \
        --genomeDir index_S288c_gtf \
        --runThreadN 12 \
        --readFilesIn clean/${i}_R1_trimmed15_val_1.fq.gz clean/${i}_R2_trimmed15_val_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix mapping/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 12 \
        --limitBAMsortRAM 1244948663 && echo "**STAR mapping 1st done**"


#index 2nd： 
        mkdir 2-pass_index/${i}_S288c_gtf_2-pass
        echo -e "\033[1;36m2_pass mapping for $i\033[0m"
        time STAR \
        --runThreadN 12 \
        --runMode genomeGenerate \
        --genomeDir 2-pass_index/${i}_S288c_gtf_2-pass \
        --genomeFastaFiles reference/latest/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
        --sjdbGTFfile reference/latest/Saccharomyces_cerevisiae.R64-1-1.104.gtf \
        --sjdbOverhang 134 \
        --sjdbFileChrStartEnd mapping/${i}_SJ.out.tab \
        --genomeSAindexNbases 10 && echo "**STAR index 2nd done**"


#2_pass mapping based on the 2nd index
        time STAR \
        --genomeDir 2-pass_index/${i}_S288c_gtf_2-pass \
        --runThreadN 12 \
        --readFilesIn clean/${i}_R1_trimmed15_val_1.fq.gz clean/${i}_R2_trimmed15_val_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix 2-pass_mapping/${i}_2-pass_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 12 \
        --limitBAMsortRAM 1244948663 && echo "**STAR mapping 2nd done**"
done

echo -e "\033[1;36mSTAR mapping finished at $(date)\033[0m"
```

```
#用samtools查看bam结果文件
samtools view ${i}_Aligned.sortedByCoord.out.bam |head -n 5


#用samtools查看比对率：
samtools flagstat BY_Aligned.sortedByCoord.out.bam
```

## 7 Get readcount table

```
```
