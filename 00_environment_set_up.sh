#!/bin/bash
# You should have conda installed before running this script.

conda env create -f ./environments/GAMMA_depict.yml && \
    conda env create -f ./environments/GAMMA_pops.yml && \
    conda env create -f ./environments/GAMMA_ldsc.yml && \
    conda env create -f ./environments/GAMMA_NetworkX.yml && \
    conda env create -f ./environments/GAMMA_ML.yml

# conda env create -f ./environments/GAMMA_r.yml && \
#     /bin/bash -c "source activate GAMMA_r && R -e 'devtools::install_github(\"ZikunY/CARMA\")' && \
#     R -e 'remotes::install_github(\"gabraham/plink2R\", subdir=\"plink2R\", ref=\"master\")' && \
#     R -e 'remotes::install_github(\"cotsapaslab/jlim\", subdir=\"jlimR\", ref=\"master\")'"

bash ./environments/GAMMA_r.sh && \
    /bin/bash -c "source activate GAMMA_r && \
    R -e 'install.packages(c(\"htmltools\", \"remotes\"), repos=\"https://cran.r-project.org\")' && \
    R -e 'devtools::install_github(\"ZikunY/CARMA\")' && \
    R -e 'remotes::install_github(\"gabraham/plink2R\", subdir=\"plink2R\", ref=\"master\")' && \
    R -e 'remotes::install_github(\"cotsapaslab/jlim\", subdir=\"jlimR\", ref=\"master\")'"

# Download yq release
wget https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 -O ~/.local/bin/yq
chmod +x ~/.local/bin/yq
yq --version
