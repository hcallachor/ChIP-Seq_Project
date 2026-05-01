# Welcome to the PlasmoPeak Pipeline Github

This Snakemake pipeline is intended to analyze ChIP-seq and ATAC-seq data for the species *Plasmodium falciparum*.






## Clone the repository
We recommend cloning this to your machine. Follow the instructions below. 

Before cloning, configure your Git identity:
```bash
git config --global user.name 'example'
git config --global user.email 'example@example.com'
```

Then clone the repository using this command:
```bash
git clone https://github.com/hcallachor/ChIP-Seq_Project
```

Then, change your working directory to the cloned repository:
```bash
cd ChIP-Seq_Project
```

---

Before running the pipeline, ensure Conda is installed. All required tools will be installed automatically via the provided environment below.

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
| bedGraphToBigWig | 2.10 |
| Conda | 26.1.1 |
| Python | 3.9.25 |
| Java | 21.0.10 |

---

## Directions to Install Softwares

### Install Miniconda (Conda)
This pipeline is designed to run with a conda environment that contains the softwares and tools needed. 

install conda using this command:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
```

-Note: if you have already previously installed miniconda, you only need to initialize the conda:
```bash
conda init bash
```

If this is your first time using miniconda, set conda path:

```bash
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
Initialize Conda:
```bash
conda init bash
```
Activate Conda:
```bash
conda activate environment
```
---
# Create the environment from the provided YAML file
```bash
conda env create -f environment.yaml -n environment
```
### Note: If you add new tools to the YAML later, update the environment using:
```bash
env update -n environment -f environment.yaml
```



### Install PyYAML

```bash
conda install -c conda-forge pyyaml
```

---

## To Begin Using the Pipeline

### 1. Configure your samples
   
Edit the file:

```
YourSamplesHere.yaml
```

This yaml connects to the sample download code to begin downloading the data you would like to use. 

You should not need to make any edits to the SampleDownload.py code itself, except for adding the reference genome if you plan to work with something other than *Plasmodium falciparum*. (ref genome as of May 2026)

Add the **SRA IDs** you want to analyze.  
An example format is already included in the file.

---

### Sample test dataset

For quick testing of the pipeline, use:

**SRR34020775**

This dataset is very small and should run in ~2 minutes if using 4 cores.

---


### 2. Download sample data
Run the download script:

```bash
python ./SampleDownload.py
```

### What it does:
- Creates a directory called **initial_data/**
- Reads SRA IDs from `YourSamplesHere.yaml`
- Downloads sequencing data using SRA tools
- Generates the config file for Snakemake and the PlasmoPeak file

**Implementation notes:**
- Uses Python `subprocess` to run command-line tools
- Checks and creates directories using the operating system
- SRA list can be modified inside the script (line 19)
- Reference genome settings are defined in the output dictionary

---

## Running the pipeline

### Run Snakemake

```bash
snakemake -s PlasmoPeak -c 4 --configfile config.yaml --use-conda
```

### Run Snakemake (background mode)

```bash
nohup snakemake -s PlasmoPeak -c 4 --configfile config.yaml --use-conda &
```

---

### Run cleanup workflow
```bash
snakemake cleanup -c 1 --configfile config.yaml
```

---

### List available rules

```bash
snakemake --list
```

---

## Viewing your results
- bigWig files can be viewed on the UCSC Genome Browser: https://genome.ucsc.edu/
- bigWig files must be hosted on a web-accessible URL (such as on GitHub)


### Pushing to GitHub
- First, create a repository on your GitHub account, and clone the repository following the instructions listed at the top of this README.
- Move the bigWig output files to your own repository, and push them to GitHub using git add, commit, and push.

### Using the Genome Browser
- Open your GitHub repository to the bigWig file, and right-click the "raw" button in the top right, in order to copy the raw data URL. (this can also be accomplished by changing "blob" to "raw" in the URL)
- Under the "Genomes" tab on the browser website, search for your reference genome (in this case GCF_000002765.6), and select "go".
- Next, paste the URL to your bigWig file in the "Custom Tracks" section found under "My Data", submit and select "Go To First Annotation".
- Viewing settings must be manually configured by selecting the gear icon on the left side of your sample track.
- Simply change the "dense" setting to "full" and submit changes, then the tracks will show peak data.

---

## Notes
- Ensure `YourSamplesHere.yaml` is properly formatted
- Modify reference genome settings in `SampleDownload.py` if needed
- We chose to set this up in a conda environment, but all of these tools could be downloaded and their calls changed within the pipeline to reflect the versions individually downloaded. 
- You can rename the main pipeline file if you make changes, but will need to change the name within the snakemake calls

```bash
snakemake -s NewName -c 4 -config --configfile config.yaml --use-conda
```

