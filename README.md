Hello thank you for using our repo!
This snakemake pipeline is intended to analyze ChIP-seq and ATAC-seq data for the species Plasmodium falciparum.

To run the Snakefile in background with a log file outputted
```bash
nohup snakemake -s snakefile -c 4  > snakemake.log 2>&1 &
```

To run the Snakefile clenup
```bash
snakemake -c 1 cleanup
```
