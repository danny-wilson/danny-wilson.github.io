#!/usr/bin/env bash
#
# Statistical Genomics Summer School
# initialization script
# Staphylococcus aureus antimicrobial resistance example with 992 genomes

mkdir -p ~/Wednesday/amr
cd ~/Wednesday/amr

# Download the pre-computed kmer_pipeline results
echo 'Downloading pre-computed data for 992 Staphylococcus aureus genomes ...'
#wget --quiet --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1uQ08isFiCKLQD7Ma7ABpcy6mZc7U7m-5' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1uQ08isFiCKLQD7Ma7ABpcy6mZc7U7m-5" -O saur.tar && rm -rf /tmp/cookies.txt && 
wget --quiet https://git.ecdf.ed.ac.uk/hbecher/2206-osgss/-/raw/main/saur.tar && tar -xf saur.tar && rm saur.tar && chmod -R 755 take1/nucleotidekmer31_patternbatches && rm -rf take1/nucleotidekmer31_patternbatches

# Download the phenotypes
echo 'Downloading antimicrobial resistance phenotypes for 992 S. aureus genomes ...'
wget --quiet --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ki-ebZgbONlKRJjnBP1G1RRg_IJGJsSD' -O pheno.tar && tar -xf pheno.tar && rm pheno.tar

# Customize the prompt
echo 'PS1="\[\e]0;\u@Singularity: \w\a\]${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@Singularity\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "' >> ~/.bashrc
echo 'Done !'

