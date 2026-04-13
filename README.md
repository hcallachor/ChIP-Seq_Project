Hello thank you for using our repo!
This snakemake pipeline is intended to analyze ChIP-seq and ATAC-seq data for the species Plasmodium falciparum.

To clone this repo, first log into git using:
```bash
git config --global user.name 'example'
```
and
```bash
git config --global user.email 'example@example.com'
```

then clone this repo using the command:
```bash
git clone https://github.com/hcallachor/ChIP-Seq_Project
```

To Begin:
You need to add your samples to the "YourSamplesHere.yaml"
First: open the YourSamplesHere.yaml
Next: add the SRAs you plan to use, there is a formatting example inside the YourSamplesHere.yaml
This yaml connects to the sample download code to begin downloading the data you would like to use. 
You should not need to make any edits to the SampleDownloadV2.py code itself, except for adding the reference genome if
you plan to work with something other than Plasmodium falciparum

For quick testing of the pipeline, we recommend using "SRR34020775" in "YourSamplesHere.yaml", as it is a very small sample file capable of running in ~2 minutes.

To install the yaml functionality, use the following command:
```bash
pip install PyYAML
```

SampleDownload.py:
This Script is made to download the initial data to run through the pipeline and create the config file needed for the pipeline.
It will create a directory to store samples called 
"initial_data"
The function will use the operating system to check for the directory. 
Subproccess allows for the running of commandline tool from within Python. 
To adapt the Sample Download code to other samples, the list at line 19 can be changed to other SRR/SRA numbers
Please check that SampleDownloadV2.py is in your current working directory before using ./
Otherwise list full path
Additionally, if you need to change the reference genome you wish to use, this is under the output dictionary 

Example to run:
```bash
python ./SampleDownloadV2.py
```
To run the Snakefile in background with a log file outputted
Must run with the conda Macs3 environment activated
```bash
conda activate macs3_env
```
```bash
nohup snakemake -s Snakefile -c 4 --configfile CompProjectconfig.yaml > snakemake.log 2>&1 &
```

To run the Snakefile cleanup
```bash
snakemake cleanup -c 1 --configfile CompProjectconfig.yaml
```

To See What Rules your Snakefile Can See:
```bash
snakemake --list
```

Trimmomatic Notes:
Your environment needs Java to run Trimmomatic, an example that downloads Trimmomatic is as follows:

java -jar trimmomatic-0.39.jar

Please See https://github.com/usadellab/Trimmomatic/blob/main/README.md For Further Details and to Download trimmomatic-0.39.jar

trimmomatic was downloaded using this zip file:
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

Picard MarkDupilcates Notes:
Your environment needs Java to run Picard MarkDuplicates, an example that downloads Picard MarkDuplicates is as follows:

java -jar picard.jar 

Picard MarkDupilcates was downloaded using this file:
https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

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

Install Macs3:
```bash
conda install -c bioconda macs3
```

Install bedToBigBed:
```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
```
After a bigBed file is created, it needs to be hosted on a web server in order to be compatible with the UCSC Genome Browser. A convenient option is to upload the bigBed file to GitHub, then copy the raw file URL and paste it into the custom tracks function found under the "My Data" tab of the genome browser website.

Create and activate conda environment:
```bash
conda create -n macs3_env python=3.9
conda activate macs3_env
conda install -c bioconda -c conda-forge macs3
```

