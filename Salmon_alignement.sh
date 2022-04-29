###############################
###############################

# Salmon indexes


# # https://github.com/COMBINE-lab/salmon/issues/603
# # similar to here https://github.com/COMBINE-lab/salmon/issues/214
# 
# mkdir /opt/genomes/mouse_mm10/Salmon_homemade/
# cd /opt/genomes/mouse_mm10/Salmon_homemade/
# 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
# 
# gunzip gencode.vM23.transcripts.fa.gz
# gunzip GRCm38.primary_assembly.genome.fa.gz
# 
# grep "^>" GRCm38.primary_assembly.genome.fa | cut -d " " -f 1 > GRCm38.decoys.txt
# sed -i 's/>//g' GRCm38.decoys.txt
# cat gencode.vM23.transcripts.fa GRCm38.primary_assembly.genome.fa | gzip > GRCm38.gentrome.fa.gz
# salmon index -t GRCm38.gentrome.fa.gz -d GRCm38.decoys.txt -p 12 -i salmon_index --gencode

# 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
# 
# gunzip gencode.v38.transcripts.fa.gz
# gunzip GRCh38.primary_assembly.genome.fa.gz
# 
# grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > GRCh38.decoys.txt
# sed -i 's/>//g' GRCh38.decoys.txt
# cat gencode.v38.transcripts.fa GRCh38.primary_assembly.genome.fa | gzip > GRCh38.gentrome.fa.gz
# salmon index -t GRCh38.gentrome.fa.gz -d GRCh38.decoys.txt -p 12 -i salmon_index --gencode


###############################
###############################

# Data downloaded from SRA for each dataset using the accession list files:
  
parallel --eta -j 20 "
fastq-dump --gzip --split-3 {}
" < SRR_Acc_List.txt

# each placed in its own folder


###############################
###############################

#GSE148224

cd ./GSE148224/

# Running QC FastQC
mkdir ./Output
mkdir ./Output/QC_Raw
find *.fastq.gz | parallel --eta "fastqc -o ./Output/QC_Raw {}"

# Running Trim-Galore
mkdir ./Output/trim_fastq
find *.fastq.gz | sed 's/_..fastq.gz//' | sort -n | uniq | parallel --eta "trim_galore --paired {}_1.fastq.gz {}_2.fastq.gz --fastqc --phred33 --gzip -o ./Output/trim_fastq "

# Alignement
mkdir ./Salmon_homemade_hg38_bootstraped/

while read p; 
do
echo "$p"
salmon quant -i /opt/genomes/human_hg38/Salmon_homemade/salmon_index -l A --validateMappings --numBootstraps 100 -1 <(zcat ./Output/trim_fastq/"$p"_1_val_1.fq.gz) -2 <(zcat ./Output/trim_fastq/"$p"_2_val_2.fq.gz) -o ./Salmon_homemade_hg38_bootstraped/"$p" 
done < ./SRR_Acc_List.txt


# QC
multiqc ./Salmon_homemade_hg38_bootstraped/ -o ./Salmon_homemade_hg38_bootstraped/ -d -f -s -v -p

cd ../

###############################
###############################


#GSE114853

cd ./GSE114853/

# Running QC FastQC
mkdir ./Output
mkdir ./Output/QC_Raw
find *.fastq.gz | parallel "fastqc -o ./Output/QC_Raw {}"


# Running Trim-Galore
mkdir ./Output/trim_fastq
find *.fastq.gz | parallel "trim_galore --fastqc --phred33 --gzip -o ./Output/trim_fastq {}"

# Alignement
mkdir ./Salmon_homemade_hg38_bootstraped/

ls ./Output/trim_fastq/*_trimmed.fq.gz >  ./sample_list_trimmedfastq.txt
while read p;
do
echo "$p"
salmon quant -i /opt/genomes/human_hg38/Salmon_homemade/salmon_index -l A --validateMappings --numBootstraps 100 -r <(zcat "$p") -o ./Salmon_homemade_hg38_bootstraped/"$p"
done < ./sample_list_trimmedfastq.txt


# QC
multiqc ./Salmon_homemade_hg38_bootstraped/ -o ./Salmon_homemade_hg38_bootstraped/ -d -f -s -v -p

cd ../

###############################
###############################


#GSE114913

cd ./GSE114913/


# Running QC FastQC
mkdir ./Output
mkdir ./Output/QC_Raw
find *.fastq.gz | parallel "fastqc -o ./Output/QC_Raw {}"


# Running Trim-Galore
mkdir ./Output/trim_fastq
find *.fastq.gz | parallel "trim_galore --fastqc --phred33 --gzip -o ./Output/trim_fastq {}"

# Alignement
mkdir ./Salmon_homemade_mm10_bootstraped/
  
ls ./Output/trim_fastq/*_trimmed.fq.gz >  ./sample_list_trimmedfastq.txt
while read p;
do
echo "$p"
salmon quant -i /opt/genomes/mouse_mm10/Salmon_homemade/salmon_index -l A --validateMappings --numBootstraps 100 -r <(zcat "$p") -o ./Salmon_homemade_mm10_bootstraped/"$p"
done < ./sample_list_trimmedfastq.txt

# QC
multiqc ./Salmon_homemade_mm10_bootstraped/ -o ./Salmon_homemade_mm10_bootstraped/ -d -f -s -v -p

cd ../

###############################
###############################


#GSE193866 

cd ./GSE193866/

# Running QC FastQC
mkdir ./Output
mkdir ./Output/QC_Raw
find *.fastq.gz | parallel "fastqc -o ./Output/QC_Raw {}"


# Running Trim-Galore
mkdir ./Output/trim_fastq
find *.fastq.gz | parallel "trim_galore --fastqc --phred33 --gzip -o ./Output/trim_fastq {}"

# Alignement
mkdir ./Salmon_homemade

ls ./Output/trim_fastq/*_trimmed.fq.gz >  ./sample_list_trimmedfastq.txt

while read p;
do
echo "$p"
salmon quant -i /opt/genomes/mouse_mm10/Salmon_homemade/salmon_index  -l SF --noLengthCorrection --validateMappings --numBootstraps 100  -r <(zcat "$p") -o ./Salmon_homemade/"$p"
done < ./sample_list_trimmedfastq.txt


# QC
multiqc ./Salmon_homemade/ -o ./Salmon_homemade/ -d -f -s -v -p

  
  

