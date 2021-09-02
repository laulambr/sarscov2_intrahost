# sarscov2_intrahost
 SARS-CoV-2 intrahost pipeline
 
 Scripts used to perform analysis of (intrahost) SARS-CoV-2 NGS data for both Illumina and Nanopore data.

**Table of contents**
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Usage](#usage)

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


Following tools should also be installed

1. [VarScan.v2.4.3](https://github.com/dkoboldt/varscan)
2. [ngmlr.v0.2.7](https://github.com/philres/ngmlr)
3. [Sniffles](https://github.com/fritzsedlazeck/Sniffles)

# Installation of work environments
## Oxford Nanopore 
Install the [ARTIC pipeline](https://github.com/artic-network/artic-ncov2019). For more information on the pipeline look at this [website](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

1.  Download the initial installation file 
```
git clone https://github.com/artic-network/artic-ncov2019.git
```
2. Install the ARTIC pipeline  
```
cd artic-ncov2019
conda env remove -n artic-ncov2019
conda env create -f environment.yml
```


## Illumina
Time to install the Illumina pipeline took less than 15 minutes on standard desktop computer.

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
  ## Github repo
Install primer schemes

1. Activate environment
 ```
   conda activate virus_assembly
  ```

2. Clone git repository into current directory
  ```
   git clone  https://github.com/laulambr/sarscov2_intrahost.git
  ```
3. The primer schemes can be found within the cloned directory sarscov2_intrahost/primer_schemes
  
# Data availability
Raw data can be downloaded from the Sequencing Read Archive under the following Bioproject PRJNA724859.
The [ARTIC pipeline](https://github.com/artic-network/artic-ncov2019) has included stimulated reads data on which can be used to proof run the scripts.
  
# Usage
## Oxford Nanopore 

To analyse SARS-CoV-2 data generated with ONT. The following steps can be followed:

1.  Activate work environment

  ```
  conda activate artic-ncov2019
  ```
  
2. Adapt variables in the script (Based on used primer scheme select *_Midnight.sh or *_Eden.sh)

  ```
  ...
   wrk=/path/to/output/location
   prj=name_of_output_folder
# location of demultiplexed fastq ONT data
	fastq_loc=/path/to/fastq
# location of fast5 ONT data
	fast5_loc=/path/to/fast5
# location of sequencing sumary ONT data
	seq_sum_loc=/path/to/summary
 
# location of primer schemes (installed via cloning of the sarscov2_intrahost repository)
	prim=path/to/primer_schemes
 
# install location of varscan
	instal_varscan=/path/to/varscan/
# install location of ngmlr
	instal_ngmlr=/path/to/ngmlr-0.2.7/
# install location of sniffles
	instal_sniffles=/path/to/Sniffles-maste
   ...
  ```
3. Run the script with adapted variables
  For Midnight primers:
  ```
   bash SARSCOV2_intra_ONT_MIDNIGHT.sh
  ```
  For Eden primers:
  ```
   bash SARSCOV2_intra_ONT_EDEN.sh
  ```
## Illumina
To analyse SARS-CoV-2 data generated with Illumina. The following script were followed used:
1. Activate environment
 ```
   conda activate virus_assembly
  ```
2. Adapt variables in the script
   ```
   bash *.sh
   ```
