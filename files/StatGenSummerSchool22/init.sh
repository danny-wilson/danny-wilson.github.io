#!/usr/bin/env bash
#
# Statistical Genomics Summer School
# initialization script
# Staphylococcus aureus antimicrobial resistance example with 992 genomes

mkdir -p ~/Wednesday/amr
cd ~/Wednesday/amr

# Download the pre-computed kmer_pipeline results
echo 'Downloading pre-computed data for 992 Staphylococcus aureus genomes ...'
DATA=https://www.well.ox.ac.uk/bioinformatics/training/whg_training_resources/data/genome_wide_association_studies/pathogen_gwas
wget --quiet $DATA/saur.tar && tar -xf saur.tar && rm saur.tar && chmod -R 755 take1/nucleotidekmer31_patternbatches && rm -rf take1/nucleotidekmer31_patternbatches

# Download the phenotypes
echo 'Downloading antimicrobial resistance phenotypes for 992 S. aureus genomes ...'
wget --quiet $DATA/pheno.tar -O pheno.tar && tar -xf pheno.tar && rm pheno.tar

# Customize the prompt
echo 'PS1="\[\e]0;\u@Singularity: \w\a\]${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@Singularity\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "' >> ~/.bashrc

# Custom directory for scripts, to fix a bug in plotManhattanbowtie.Rscript
mkdir -p ~/Wednesday/amr/scripts
ln -s /usr/local/bin/* ~/Wednesday/amr/scripts/
rm ~/Wednesday/amr/scripts/plotManhattanbowtie.Rscript
sed 's/as.integer(args\[13\]))/as.integer(args\[13\])/g' /usr/local/bin/plotManhattanbowtie.Rscript > ~/Wednesday/amr/scripts/plotManhattanbowtie.Rscript
sed -i 's/minor_allele_threshold = 0, ,/minor_allele_threshold = 0,/g' ~/Wednesday/amr/scripts/plotManhattanbowtie.Rscript
chmod +x ~/Wednesday/amr/scripts/plotManhattanbowtie.Rscript

# Install Entrez Direct
cd ~
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

# Download the S. aureus reference genome and genbank annotation
cd ~/Wednesday/amr
wget https://raw.githubusercontent.com/danny-wilson/danny-wilson.github.io/main/files/StatGenSummerSchool22/BX571856.1.esearch.txt
#~/edirect/esearch -db nucleotide -query "BX571856.1" | ~/edirect/efetch -format fasta > ~/Wednesday/amr/R00000022.fa
#~/edirect/esearch -db nucleotide -query "BX571856.1" | ~/edirect/efetch -format gb > ~/Wednesday/amr/R00000022.gbk
cat BX571856.1.esearch.txt | ~/edirect/efetch -format fasta > ~/Wednesday/amr/R00000022.fa
cat BX571856.1.esearch.txt | ~/edirect/efetch -format gb > ~/Wednesday/amr/R00000022.gbk

# Rename the FASTA file to that used in the precomputed steps
sed -i 's/BX571856.1/R00000022 BX571856.1/g' ~/Wednesday/amr/R00000022.fa

# Update the pipeline software file
sed 's$scriptpath\t/usr/local/bin$scriptpath\t/home/jovyan/Wednesday/amr/scripts$g' /usr/share/kmer_pipeline/example/pipeline_software_location.txt > ~/Wednesday/amr/pipeline_software_location.txt

# Signal finished
echo 'Done !'
