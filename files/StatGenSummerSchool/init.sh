#!/usr/bin/env bash
#
# Statistical Genomics Summer School
# initialization script
# Staphylococcus aureus antimicrobial resistance example with 992 genomes

mkdir -p ~/BacterialGWAS/amr
cd ~/BacterialGWAS/amr

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
mkdir -p ~/BacterialGWAS/amr/scripts
ln -s /usr/local/bin/* ~/BacterialGWAS/amr/scripts/
rm ~/BacterialGWAS/amr/scripts/plotManhattanbowtie.Rscript
sed 's/as.integer(args\[13\]))/as.integer(args\[13\])/g' /usr/local/bin/plotManhattanbowtie.Rscript > ~/BacterialGWAS/amr/scripts/plotManhattanbowtie.Rscript
sed -i 's/minor_allele_threshold = 0, ,/minor_allele_threshold = 0,/g' ~/BacterialGWAS/amr/scripts/plotManhattanbowtie.Rscript
chmod +x ~/BacterialGWAS/amr/scripts/plotManhattanbowtie.Rscript

# Download the S. aureus reference genome and genbank annotation
cd ~/BacterialGWAS/amr
wget -O ~/BacterialGWAS/amr/R00000022.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz && gunzip ~/BacterialGWAS/amr/R00000022.fa.gz
wget -O ~/BacterialGWAS/amr/R00000022.gbk.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gbff.gz && gunzip ~/BacterialGWAS/amr/R00000022.gbk.gz

# Rename the FASTA file to that used in the precomputed steps
sed -i 's/NC_007795.1/R00000022 NC_007795.1/g' ~/BacterialGWAS/amr/R00000022.fa

# Update the pipeline software file
sed 's$scriptpath\t/usr/local/bin$scriptpath\t/home/jovyan/BacterialGWAS/amr/scripts$g' /usr/share/kmer_pipeline/example/pipeline_software_location.txt > ~/BacterialGWAS/amr/pipeline_software_location.txt

# Signal finished
echo 'Done !'
