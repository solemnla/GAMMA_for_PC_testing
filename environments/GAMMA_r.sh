#!/bin/bash

# Create base R environment
echo "Creating base R environment..."
conda create -n GAMMA_r r-base=4.3.1 -c conda-forge -y

# Activate environment
echo "Activating GAMMA_r environment..."
conda activate GAMMA_r

# Install base packages batch 1
echo "Installing base packages (batch 1)..."
conda install -c conda-forge mkl r-essentials r-biocmanager -y

# Install data manipulation packages batch 2
echo "Installing data manipulation packages (batch 2)..."
conda install -c conda-forge r-data.table r-stringr r-dplyr r-r.utils -y

# Install development and tidyverse packages batch 3
echo "Installing tidyverse packages (batch 3)..."
conda install -c conda-forge r-devtools r-tidyverse -y

# Install statistical packages batch 4
echo "Installing statistical packages (batch 4)..."
conda install -c conda-forge r-glmnet r-rcpparmadillo r-rcppgsl -y

# Install Bioconductor packages batch 5
echo "Installing Bioconductor packages (batch 5)..."
conda install -c conda-forge jq -y
conda install -c bioconda bioconductor-genomeinfodb bioconductor-iranges bioconductor-s4vectors -y

# Install reporting packages batch 6
echo "Installing reporting packages (batch 6)..."
conda install -c conda-forge r-knitr r-rmarkdown -y

# Install more Bioconductor packages batch 7
echo "Installing more Bioconductor packages (batch 7)..."
conda install -c bioconda bioconductor-genomicranges bioconductor-regioner bioconductor-repitools bioconductor-qvalue -y

# Install analysis packages batch 8
echo "Installing analysis packages (batch 8)..."
conda install -c conda-forge r-optparse r-coloc r-mgcv r-survey r-ncmisc -y

# Install Mendelian Randomization packages batch 9
echo "Installing MR packages (batch 9)..."
conda install -c conda-forge r-mendelianRandomization r-mrmix r-mrpresso r-penalized r-mr.raps -y

echo "Environment setup complete!"
