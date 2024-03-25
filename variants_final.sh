
conda init bash
conda activate base
conda install -c bioconda trimmomatic
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda varscan
conda install -c bioconda bedtools
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt
tar -zxvf patient7.tar.gz
mkdir index
wget -P ./index http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz

normal_tumor=("TCRBOA7-N-WEX-chr16_r1F.fastq.gz" "TCRBOA7-T-WEX-chr16_r1F.fastq.gz")

for F in "${normal_tumor[@]}"; do

    base_name=$(basename "$F" "_r1F.fastq.gz")
    #eliminate low quality sequences at read extermities
    #If the final read is shorter than 50, it is deleted
    trimmomatic PE "./patient7.exome/${F}" "./patient7.exome/${base_name}_r2F.fastq.gz" -baseout "${base_name}" \
    LEADING:20 TRAILING:20 MINLEN:50

    #create the BWA index 
    bwa index -a bwtsw ./index/chr16.fa.gz


    #mapping with bwa
    bwa mem -M -t 2 -A 2 -E 1 ./index/chr16.fa.gz "./${base_name}_1P" "./${base_name}_2P" > "./${base_name}_output.sam"
#     Sam 2 Bam


    samtools view -S -b "./${base_name}_output.sam"  > "./${base_name}_output.bam"

    # flagstats


    samtools flagstat "./${base_name}_output.bam" 

    #Sort Bam

   
    samtools sort "./${base_name}_output.bam" > "./${base_name}_output_sorted.bam" 

    #Index bam file
 

    samtools index "./${base_name}_output_sorted.bam"
    gunzip ./index/chr16.fa.gz
#     Convert to Mpileup


    samtools mpileup -B -A -f ./index/chr16.fa  "./${base_name}_output_sorted.bam" > \
    "./${base_name}_mpileupfile"
done

#we use varscan to create variant lists from Normal and Tumor Pileup files:
varscan somatic ./TCRBOA7-N-WEX-chr16_mpileupfile ./TCRBOA7-T-WEX-chr16_mpileupfile ./output_varscan

grep -i 'somatic' ./output_varscan.indel > ./output_varscan_filtered.indel
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
./output_varscan_filtered.indel > ./indel.bed
grep -i 'somatic' ./output_varscan.snp > ./output_varscan_filtered.snp
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
   ./output_varscan_filtered.snp > ./snp.bed

# Download annotation file with wget from this URL and uncompress it:
wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
gunzip gencode.v24lift37.basic.annotation.gtf.gz 
conda install -n base -c bioconda bedtools

bedtools intersect -a ./gencode.v24lift37.basic.annotation.gtf  -b ./indel.bed > ./indel.intersect
grep '\sgene\s' ./indel.intersect | awk '{print " " $1 " " $4 " " $5 " " $16}'


bedtools intersect -a ./gencode.v24lift37.basic.annotation.gtf  -b ./snp.bed > ./snp.intersect
grep '\sgene\s' ./snp.intersect | awk '{print " " $1 " " $4 " " $5 " " $16}'