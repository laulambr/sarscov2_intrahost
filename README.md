# sarscov2_intrahost
 SARS-CoV-2 intrahost pipeline
 
 Scripts used to perform analysis of (intrahost) SARS-CoV-2 NGS data for both Illumina and Nanopore data.

**Table of contents**
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Example](#quick-start)

# System Requirements

## Hardware requirements
The pipeline can be run on either a standard computer or a HPC server. We tested the pipeline on a standard desktop computer with the following specifications:

RAM: 16+ GB
CPU: 4+ cores, @1.90 GHz

## Software requirements
### OS requirements
The pipeline has been tested on several Linux operating systems including the following systems:

Linux: Ubuntu 16.04, Ubuntu 18.04, Ubuntu 20.04

### Dependencies

A conda package manager like [Miniconda3](https://docs.conda.io/en/latest/miniconda.html). Instructions on how to install:
1.  Download the latest miniconda installation script
   ```
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```
2. Make the miniconda installation script executable
  ```
  chmod +x Miniconda3-latest-Linux-x86_64.sh
  ```
3. Run miniconda installation script
  ```
  bash ./Miniconda3-latest-Linux-x86_64.sh
  ```

# Installation
## Oxford Nanopore 
Install the [ARTIC pipeline](https://github.com/artic-network/artic-ncov2019). 

1.  Download the initial installation file 
```
git clone https://github.com/artic-network/artic-ncov2019.git
```
2. Install the ARTIC pipeline  
```cd artic-ncov2019
conda env remove -n artic-ncov2019
conda env create -f environment.yml
```


## Illumina
Time to install the pipeline took less than 15 minutes on standard desktop computer.

1.  Download the initial environment installation file 
   ```
   wget https://raw.githubusercontent.com/laulambr/virus_assembly/main/scripts/install_env.sh
   ```
2. Run the script in the terminal 
  ```
   bash ./install_env.sh
   ```
3. Check if installation worked
  ```
   conda activate virus_assembly
  ```
# Usage
   ```
 Pipeline: NGS pipeline for viral assembly.
usage: virus_assembly [-h -v -p -q] (-i dir -m value -t value )
(-s string) 
with:
    -h  Show help text
    -v  Version of the pipeline
    -n  Name of RUN.
    -i  Input directory
    -s  Viral species [HIV, RSV, RRV, HMPV]
    -c  Perform clipping of primers
    -q  Perform quality check using fastQC
    -m  Memory
    -t  Number of threads
   ```
## Quick start
 1. Activate environment.
 
   ```
   conda activate virus_assembly
  ```
 2. Head to the directory where you will perform the analysis.
 3. Place the raw fastq.gz files in a directory called source.  
 4. Create a list holding the sample names from you sequencing files called IDs.list and place it in the main directory.

