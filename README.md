# GAMMA Demo

## Overview

This is a **DEMO version of GAMMA** designed for testing purposes. The demo has been modified to run only on **chromosome 11** with a limited number of xQTLs to facilitate quick testing and evaluation.

## Prerequisites

- Ensure you have the necessary system requirements installed
- Sufficient disk space for test data and results

## Demo Data

Download the test data from: 
https://yanglab.westlake.edu.cn/data/gamma-download/GAMMA_test_data.tar.gz

## Tutorial

Follow these steps to run the GAMMA demo:

### 0. Code initialization

```bash
git clone https://github.com/solemnla/GAMMA_for_PC_testing.git
cd GAMMA_for_PC_testing
chmod 754 ./00_environment_set_up.sh ./01_run_pipeline.sh
```

### 1. Download demo data

Download the test data from the link above and extract the data directory under the `GAMMA_for_PC_testing` directory.

```bash
# download demo data (T2D chr11 GWAS and other functional annotation data used in GAMMA analysis)
wget https://yanglab.westlake.edu.cn/data/gamma-download/GAMMA_test_data.tar.gz
tar -xvzf GAMMA_test_data.tar.gz
```

### 2. Install Required Environment

Set up the required environment using the provided script:

```bash
bash ./00_environment_set_up.sh
```

### 3. Configure GAMMA Home Path

Adjust the `deploy/GAMMA.yaml` file:
- **Line 2**: Update the `GAMMA_HOME` path to the absolute path where your `GAMMA_for_PC_testing` directory is located

```yaml

GAMMA_HOME=change/this/to/path/to/your/GAMMA_for_PC_testing
CONFIG=${GAMMA_HOME}/deploy/GAMMA.yaml

yq -i ".input.GAMMA_HOME = \"$GAMMA_HOME\"" "$CONFIG"
```

### 4. Run the Pipeline

Execute the pipeline script:

```bash
bash ./01_run_pipeline.sh ${CONFIG}
```

A `results` directory will be generated under `GAMMA_for_PC_testing`.

### 5. View Results

After successful execution, you can find:

- **GAMMA summary data**: `results/GAMMA/feature/T2D_chr11_GAMMA.feature`
- **Machine learning scores**: `results/GAMMA/score/ML_score.csv`

## Directory Structure

```
GAMMA_for_PC_testing/
├── data/                          # Extracted test data
├── deploy/
│   └── GAMMA.yaml                 # Configuration file
├── environments/
├── results/                       # Generated after running pipeline
│   └── GAMMA/
│       ├── feature/
│           └── T2D_chr11_GAMMA.feature               # GAMMA summary data
│       └── score/
│           └── ML_score.csv      # Machine learning scores
├── scripts/
├── softwares/
├── 00_environment_set_up.sh      # Environment setup script
└── 01_run_pipeline.sh            # Main pipeline script
```

## Notes

- This demo is limited to chromosome 11 and a subset of xQTLs for demonstration purposes
- For full functionality, please refer to the complete GAMMA implementation

## Support

For issues or questions, please open an issue on the GitHub repository.
