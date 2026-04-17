# Welcome to the "Cool Name for Pipeline Here" Github

This Snakemake pipeline is intended to analyze ChIP-seq and ATAC-seq data for the species Plasmodium falciparum.





## Clone the repository
We recommend cloning this to your machine. Follow the instructions below. 

Before cloning, configure your Git identity:
```bash
git config --global user.name 'example'
```
and
```bash
git config --global user.email 'example@example.com'
```

Then clone the repository using this command:
```bash
git clone https://github.com/hcallachor/ChIP-Seq_Project
```


---

 Before using the pipeline, please ensure you have the following tools downloaded. This document includes instructions on how to download the tools.

## Tools used in this pipeline:
| Tool | Version |
|------|--------|
| Snakemake | 7.32.4 |
| PyYAML | python module |
| fasterq-dump | 3.0.3 |
| Trimmomatic | 0.39 |
| BWA | 0.7.17-r1188 |
| Samtools | 1.19.2 |
| Picard (MarkDuplicates) | 3.4.0 |
| MACS3 | 3.0.3 |
| Bedtools | 2.31.1 |
| bedToBigBed | 2.10 |
| Conda | 26.1.1 |
| Python | 3.9.25 |
| Java | 21.0.10 |

---

## Directions to Install Softwares

### Install Miniconda (Conda)

install conda using this command:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
```

set conda path:

```bash
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
conda init bash
```

---

### Install OpenJDK from conda-forge
This is a version of java you can obtain through conda

```bash
conda install -c conda-forge openjdk

# Verify installation
java -version
```

---

### Trimmomatic

Requires Java.

```bash
java -jar trimmomatic-0.39.jar
```

Please see the link below for further details and to download trimmomatic-0.39.jar:
https://github.com/usadellab/Trimmomatic

---

### Picard MarkDuplicates

Requires Java.

```bash
java -jar picard.jar
```

Please see the link below for further details and to download:
https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

---

### bedToBigBed

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
```
After a bigBed file is created, it needs to be hosted on a web server in order to be compatible with the UCSC Genome Browser. A convenient option is to upload the bigBed file to GitHub, then copy the raw file URL and paste it into the custom tracks function found under the "My Data" tab of the genome browser website.

---

### Create environment and intsall macs3

```bash
conda create -n macs3_env python=3.9
conda activate macs3_env
conda install -c bioconda -c conda-forge macs3
```

---


## To Begin Using the Pipeline

### 1. Configure your samples
   
Edit the file:

```
YourSamplesHere.yaml
```

This yaml connects to the sample download code to begin downloading the data you would like to use. 

You should not need to make any edits to the SampleDownload.py code itself, except for adding the reference genome if you plan to work with something other than Plasmodium falciparum.

Add the **SRA IDs** you want to analyze.  
An example format is already included in the file.

---

### Sample test dataset

For quick testing of the pipeline, use:

**SRR34020775**

This dataset is very small and should run in ~2 minutes.

---

### 2. Install Python dependency

To install the YAML functionality, use the following command:
```bash
pip install PyYAML
```

---

### 3. Download sample data
Run the download script:

```bash
python ./SampleDownload.py
```

### What it does:
- Creates a directory called **initial_data/**
- Reads SRA IDs from `YourSamplesHere.yaml`
- Downloads sequencing data using SRA tools
- Generates the config file for Snakemake

**Implementation notes:**
- Uses Python `subprocess` to run command-line tools
- Checks and creates directories using the operating system
- SRA list can be modified inside the script (line 19)
- Reference genome settings are defined in the output dictionary

---

## Running the pipeline

### Activate Conda environment
(if conda has not yet been installed, see directions at the end)

```bash
conda activate macs3_env
```

### Run Snakemake (background mode)

```bash
nohup snakemake -s Snakefile -c 4 --configfile CompProjectconfig.yaml > snakemake.log 2>&1 &
```

---

### Run cleanup workflow
```bash
snakemake cleanup -c 1 --configfile CompProjectconfig.yaml
```

---

### List available rules

```bash
snakemake --list
```

---


## Notes

- Java is required for Trimmomatic and Picard
- Always activate the conda environment before running Snakemake
- Ensure `YourSamplesHere.yaml` is properly formatted
- Modify reference genome settings in `SampleDownload.py` if needed

