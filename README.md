Hello thank you for using our repo!
This snakemake pipeline is intended to analyze ChIP-seq and ATAC-seq data for the species Plasmodium falciparum.


SampleDownloadPFal.py:
This Script is made to download the initial data to run through the pipeline. 
It will create a directory to store samples called 
"initial_data"
The function will use the operating system to check for the directory. 
Subproccess allows for the running of commandline tool from within Python. 
To adapt the Sample Download code to other samples, the list at line 19 can be changed to other SRR/SRA numbers
Please check that SampleDownloadPFal.py is in your current working directory before using ./
Otherwise list full path


Example to run:
```bash
python ./SampleDownloadPFal.py
```
To run the Snakefile in background with a log file outputted
```bash
nohup snakemake -s snakefile -c 4 --configfile CompProjectconfig.yaml > snakemake.log 2>&1 &
```

To run the Snakefile cleanup
```bash
snakemake -c 1 cleanup
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
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

