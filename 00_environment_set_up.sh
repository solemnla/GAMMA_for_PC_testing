#!/bin/bash
# You should have conda installed before running this script.

conda env create -f ./environments/depict.yml && \
    conda env create -f ./environments/pops.yml && \
    conda env create -f ./environments/ldsc.yml && \
    conda env create -f ./environments/NetworkX.yml && \
    conda env create -f ./environments/py39.yml

conda env create -f ./environments/r.yml && \
    /bin/bash -c "source activate GAMMA_r && R -e 'devtools::install_github(\"ZikunY/CARMA\")' && \
    R -e 'remotes::install_github(\"gabraham/plink2R\", subdir=\"plink2R\", ref=\"master\")' && \
    R -e 'remotes::install_github(\"cotsapaslab/jlim\", subdir=\"jlimR\", ref=\"master\")'"